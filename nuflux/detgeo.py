import numpy as np
import matplotlib.pyplot as plt
import DarkNews as dn
from scipy.stats import uniform
import copy

import importlib

import sys
sys.path.append("/home/luc/Research/BIN_MC/nuflux/detector_geometries")
import helpers
import useful_data
import time
from numba import njit

def D3distance(point1, point2):
    return np.sqrt((point1[:,0] - point2[:,0])**2 + (point1[:,1] - point2[:,1])**2 + (point1[:,2] - point2[:,2])**2)

def get_cs(E, part):
    sigmanue,sigmanuebar,sigmanumu,sigmanumubar= useful_data.cs_interp()
    if part == "nue":
        return sigmanue(E)
    elif part == "numu":
        return sigmanumubar(E)
    
def get_quantities(sim):
    return helpers.cc(sim)

class SimulateDetector():

    def __init__(self, coord_object, geom, particle):
        self.time = np.zeros(6)
        self.time[0] = time.time()
        self.cc = coord_object
        self.dec_pos = self.cc.p
        self.w = self.cc.weights.reshape((self.cc.weights.shape[0],1)) # already normalized
        self.Geometry = geom
        self.Nmu = self.cc.Nmu
        self.sample_size = self.dec_pos.shape[0]
        self.particle = particle
        if self.particle=='numu':
            self.momenta = self.cc.mnumu
            self.E = self.cc.Enumu
        elif self.particle=='nue':
            self.momenta = self.cc.mnue
            self.E = self.cc.Enue

        geom = importlib.import_module(self.Geometry)
        self.objects = geom.OBJECTS
        self.object_ids = np.array([obj.id for obj in self.objects])
        self.zbeginning = geom.zbeginning
        self.rmax = geom.rmax
        self.rbp = geom.rbp
        self.zending = geom.zending
        self.initials = geom.INITIALIZERS #should not contain minus_one
        self.decayer = geom.DECAYER #is the minus_one
        self.outside = geom.OUTSIDE # parent class of outside-going faces
        self.face_dict = geom.facedict
        self.outside_ids = np.array([obj.id for obj in self.outside])
        self.iterations = geom.iterations #size of the interactions_points array, one more than distances and densities
        

        self.initialize_quantities()


    def find_info(self): 
        #first interaction points
        #self.time[1] = time.time()
        count = 1
        delta_z = self.zbeginning - self.dec_pos[:,2]
        to = delta_z / self.momenta[:,2]
        self.t = to.reshape((self.sample_size, 1))
        ip = self.dec_pos + self.t * self.momenta #interaction points
        self.r_values = np.sqrt(ip[:,0]**2+ip[:,1]**2)
        
        cond_1 = ~(((self.intersection_points[:,0,2] < self.zending) & (self.intersection_points[:,0,2] > self.zbeginning) & (self.intersection_points[:,0,1] > -1 * self.rmax) & (self.intersection_points[:,0,1] < self.rmax)))
        cond_2 = ((self.r_values > self.rmax) | (to < 0))
        self.mask_not_accepted = cond_2 & cond_1
        
        not_accepted_indices = np.where(self.mask_not_accepted)[0] #indices of this new array 
        self.location[not_accepted_indices, 1] = 0
        self.intersection_points[not_accepted_indices, 1, :] = self.intersection_points[not_accepted_indices, 0, :]
        self.time[1] = time.time()
        self.mask_accepted = ~self.mask_not_accepted
        for obj in self.initials:
            new_indices = obj.check_in(self.r_values,self.t,  self.mask_accepted)
            self.location[new_indices, 1] = obj.id
            self.intersection_points[new_indices, 1, :] = ip[new_indices,:]
        self.time[2] = time.time()
        
        #at this point, there should be -1s (decaying in detector), 0s (did not reach detector), and ids of initials
        mask_decay = (self.location[:,1] == -1) & ~((self.dec_pos[:,2] > self.zbeginning) & (self.dec_pos[:,2] < self.zending) & (np.sqrt(self.dec_pos[:,0]**2 + self.dec_pos[:,1]**2) > self.rbp) & (np.sqrt(self.dec_pos[:,0]**2 + self.dec_pos[:,1]**2) < self.rmax))
        self.update_intersections(self.decayer[0], count, mask_decay)

        #treating decays within the detector
        weird_decay_mask = (self.location[:,1] == -1)
        self.update_intersections(self.decayer[1], count, weird_decay_mask)
        self.location[weird_decay_mask,0] = -2
        #self.time[2] = time.time()
        #there should be no more -1s, just 0 - x
        count = 2
        while count < self.iterations:
            
            locs = np.unique(self.location[:, count - 1])
            objects = [self.objects[i] for i in locs]
            for obj in objects:
                
                part_mask = (self.location[:, count -1] == obj.id)
                particles = np.where(part_mask)[0]#these are indices of the particles at that face

                if obj.id in self.outside_ids:
                    self.intersection_points[particles, count, :] = self.intersection_points[particles, count-1, :]
                    self.location[particles, count] = self.location[particles, count - 1]
                    continue
                
                if np.any(particles):
                    self.update_intersections(obj, count, part_mask)
            count+=1

        self.time[3] = time.time()

    def get_probs(self):

        for i in range(self.iterations - 1):
            self.distances[:,i] = D3distance(self.intersection_points[:,i,:],self.intersection_points[:,i+1,:]) #these input arrays will be two dimensional, i.e., (sample_size, 3); distances will be i.e., (sample_size, 14)
        self.part_line_integrals = self.distances * self.densities # huge array of line_ints at each component for each particle
        self.line_integrals = np.sum(self.part_line_integrals, axis = 1)
        self.cs = get_cs(self.E, self.particle)
        self.probs = 1 - np.exp(-1*self.cs * self.line_integrals)

        self.factors = np.full((self.sample_size,1), self.Nmu) * self.w

        self.counts = self.factors * self.probs.T.reshape((self.sample_size,1)) # (sample_size, 1)
        self.total_count = np.sum(self.counts)
        self.time[4] = time.time()

    def get_event_positions(self):
        '''weights have already been given'''

        self.mask = (self.counts>0)
        self.mask = self.mask[:, 0]
        scales = self.distances[self.mask]
        temporary = np.empty((1,np.sum(self.mask), self.iterations - 1))
        temporary[:,:, :] = [uniform.rvs(loc = 0, scale = scales)]
        self.events_position, self.part_face_counts = get_probs_njit(self.dec_pos, self.momenta, self.mask, self.distances, self.sample_size, self.iterations, self.counts, self.cs, self.part_line_integrals, temporary)
        self.time[5] = time.time()
        return


    def initialize_quantities(self):

        self.intersection_points = np.full((self.sample_size, self.iterations, 3), 1e4) #1e4 is arbitrary
        self.densities = np.zeros((self.sample_size, self.iterations - 1))
        self.location = np.full((self.sample_size, self.iterations), -1) #starting at initials
        self.intersection_points[:,0,:] = self.dec_pos
        self.distances = np.zeros((self.sample_size, self.iterations - 1))


    def update_params(self):

        if self.particle =='numu':
            self.momenta = self.cc.mnumu

            self.E = self.cc.Enumu
        elif self.particle == 'nue':
            self.momenta = self.cc.mnue
            self.E = self.cc.Enue

    
    def update_intersections(self, obj, count, mask):
        for neighbor in obj.next_ids:
            neigh = self.objects[neighbor]
            new_indices, ips = neigh.check_intersection(self.intersection_points[:,count - 1,:], self.momenta, mask)
                        
            mask_1 = (ips[:,2] < self.intersection_points[new_indices, count, 2])& (ips[:,2] > self.intersection_points[new_indices, 0, 2])
            accepted_ix = new_indices[mask_1]
            self.intersection_points[accepted_ix, count, :] = ips[mask_1]
            self.location[accepted_ix, count] = neigh.id
            self.densities[accepted_ix, count-1] = obj.density

    def run(self):
        if (self.particle == "nue") | (self.particle == "numu"):
            self.find_info()
            self.get_probs()
            self.get_event_positions()
            print("{:.3g} {} events".format(self.total_count, self.particle))
            #print('time: {:.3g} for initialization;\n{:.3g} for initial objects;\n{:.3g} for other objects;\n{:.3g} for probs;\n{:.3g} for event positions;\n{:.3g} for total time.'.format(
                #self.time[1] - self.time[0],self.time[2] - self.time[1], 
                #self.time[3] - self.time[2], self.time[4] - self.time[3], 
                #self.time[5] - self.time[4],self.time[5] - self.time[0]))
            return self
        
        elif self.particle=='both':
            self.particle="numu"
            self.update_params()
            self.run()
            sim2 = copy.deepcopy(self)
            sim2.particle="nue"
            sim2.update_params()
            sim2.initialize_quantities()
            sim2.run()
            return self, sim2
        
        else:
            raise ValueError("No particle of that name incldued in the detector!")

    def get_face_counts(self, arg, sim2):
        self.facecounts = copy.deepcopy(self.face_dict)
        locs = self.location[:, :-1]
        for key in self.face_dict.keys():
            total = 0
            for face in self.face_dict[key]:
                mask = (locs == face)
                total += np.sum(self.part_face_counts[mask])
            self.facecounts[key] = total
            if arg != 'both':
                print('{:.3g} decays in {}'.format(total, key))

        if arg =='both':
            facecounts = copy.deepcopy(sim2.face_dict)
            locs = sim2.location[:, :-1]
            for key in facecounts.keys():
                total = 0
                for face in facecounts[key]:
                    mask = (locs == face)
                    total += np.sum(sim2.part_face_counts[mask])
                self.facecounts[key] += total
                print('{:.3g} decays in {}'.format(self.facecounts[key], key))
            
            print('{:.3g} total decays'.format(np.sum(self.part_face_counts) + np.sum(sim2.part_face_counts)))
        
        else:
            print('{:.3g} total decays'.format(np.sum(self.part_face_counts)))
        

def plot_sim(geom):

    if geom == 'approximate_muon_detector_1':
        T1 = [[-231, 231,231,28,-28,-231,-231], [150,150,24,3,3,24,150]]

        ECAL1=[[-231, -231, 231, 231, -231],[150, 170, 170, 150, 150]]
        ECAL2 = [[231, 231, 251, 251, 231],[24, 170, 170, 26, 24]]
        ECAL3 = [[-1 * i for i in ECAL2[0]], ECAL2[1]]

        HCAL1 = [[-251, -251, 251, 251, -251],[170, 348, 348, 170, 170]]
        HCAL2= [[251, 251, 418, 418,251, 251,251],[170, 348, 348, 43,26, 170,170]]
        HCAL3 = [[-1 * i for i in HCAL2[0]], HCAL2[1]]

        SOLENOID =[[-418, -418, 418, 418, -418],[348, 446, 446, 348, 348]]

        MD1 = [[-418, -418, 418, 418, -418],[446, 645, 645, 446, 446]]
        MD2 = [[418, 418, 564, 564, 418],[43, 645, 645, 58, 43]]
        MD3 = [[-1 * i for i in MD2[0]], MD2[1]]

        CONE1 = [[28,564,564,28],[3,58,3,3]]
        CONE2 = [[-1 * i for i in CONE1[0]], CONE1[1]]

        BL = [[-564, -564, 564, 564, -564], [-3, 3, 3, -3, -3]]

        dets = [T1, ECAL1, ECAL2, ECAL3, HCAL1, HCAL2, HCAL3,SOLENOID, MD1, MD2, MD3, CONE1, CONE2, BL]
        cols = ['lightgrey']*1 + ['dimgrey']*3 + ['grey']*3 + ['darkgrey'] + 3*['grey'] + ['black']*2 + ['white']
        plt.figure(figsize = (20,12))
        for i, det in enumerate(dets):
            plt.plot(det[0],det[1], color = cols[i])
            plt.fill_between(det[0], det[1], color = cols[i], alpha=0.7)
            new_y =[-1*k for k in det[1]]
            plt.plot(det[0], new_y, color = cols[i])
            plt.fill_between(det[0], new_y, color = cols[i], alpha=0.7)

        plt.xlabel("z-coordinate (cm)")
        plt.ylabel("r-coordinate (cm)")

    
    elif (geom == 'approximate_muon_detector_2') | (geom == 'approximate_muon_detector_3'):
        ECAL1=[[-221, -221, 221, 221, -221],[150, 170.2, 170.2, 150, 150]]
        ECAL2 = [[230.7, 230.7, 250.9, 250.9, 230.7],[31, 170, 170, 33.9, 31]]
        ECAL3 = [[-1 * i for i in ECAL2[0]], ECAL2[1]]

        HCAL1 = [[-221, -221, 221, 221, -221],[174, 333, 333, 174, 174]]
        HCAL2= [[235.4, 235.4, 412.9, 412.9,250.9, 250.9,235.4],[170, 324.6, 324.6, 56.8,33.9, 170,170]]
        HCAL3 = [[-1 * i for i in HCAL2[0]], HCAL2[1]]

        SOLENOID =[[-412.9, -412.9, 412.9, 412.9, -412.9],[348.3, 352.3, 352.3, 348.3, 348.3]]
        SOLENOID_2 =[[-412.9, -412.9, 412.9, 412.9, -412.9],[364.9, 399.3, 399.3, 364.9, 364.9]]
        SOLENOID_3 =[[-412.9, -412.9, 412.9, 412.9, -412.9],[425, 429, 429, 425, 425]]

        MD1 = [[-563.8, -563.8, 563.8, 563.8, 417.9, 417.9, -417.9, -417.9, -563.8],[78.2, 645, 645, 78.2, 57.5, 446.1, 446.1, 57.5, 78.2]]


        CONE1 = [[6.5,230.7, 250.9, 412.9, 417.9, 563.8,563.8, 6.5],[2.2,31,33.9, 56.8, 57.5 ,78.2, 2.2, 2.2]]
        CONE2 = [[-1 * i for i in CONE1[0]], CONE1[1]]

        BL = [[-563.8, -563.8, 563.8, 563.8, -563.8], [-2.2, 2.2, 2.2, -2.2, -2.2]]

        dets = [ECAL1, ECAL2, ECAL3, HCAL1, HCAL2, HCAL3,SOLENOID, SOLENOID_2, SOLENOID_3, MD1, CONE1, CONE2, BL]
        cols = ['dimgrey']*3 + ['darkgrey']*3 + ['grey'] + ['lightgrey'] + ['gray'] + ['grey'] + ['black']*2 + ['white']
        #plt.figure(figsize = (20,12))
        for i, det in enumerate(dets):
            plt.plot(det[0],det[1], color = cols[i])
            plt.fill_between(det[0], det[1], color = cols[i], alpha=0.7)
            new_y =[-1*k for k in det[1]]
            plt.plot(det[0], new_y, color = cols[i])
            plt.fill_between(det[0], new_y, color = cols[i], alpha=0.7)

        plt.xlabel("z-coordinate (cm)")
        plt.ylabel("y-coordinate (cm)")


    else:
        print("this geometry has not been implemented yet!")

#finally, no use of njit; numpy does this efficiently already
def get_probs_njit(dec_pos, momenta, mask, distances, sample_size, iterations,counts, cs, part_line_integrals, temporary):
    part_face_counts =  np.zeros(distances.shape) #counts at each face for each particle
    int_parameters = np.zeros(distances.shape) #t values starting at the last interaction point for each simulated interaction
    events_position = np.empty((sample_size, iterations - 1, 3))

    pweights = np.zeros(distances.shape) #these are probabilities, NOT WEIGHTS
    temp = np.zeros((sample_size, iterations - 1)) #probs of interaction in each segment
    temp2 = np.ones((sample_size, iterations - 1)) #probs of non-interaction in each segment
        
    middle = cs[mask]
    middle2 = part_line_integrals[mask]
    temp2[mask, :] = np.exp(-1 *middle[:, np.newaxis] * middle2)
    temp[mask,:] = 1 - temp2[mask,:] # probs of interaction on each segment

    mi = np.prod(temp2, axis = 1)
    mid = mi[:, np.newaxis] / temp2
    pweights[mask,:] = temp[mask,:] * mid[mask,:]
    max_index = np.argmax(pweights[mask], axis=1)

    pweights[mask, :] = pweights[mask] / pweights[mask, max_index][:, np.newaxis]

    mid  = np.sum(pweights, axis = 1)

    pweights[mask,:] = pweights[mask] / mid[mask, np.newaxis] #good multinomial probabilities
    part_face_counts[:,:] = counts * pweights # these are the face counts (expected)
        

    #getting the distances to the event points for each ray
    int_parameters[mask, :] = temporary
    cumulative_distances = np.cumsum(distances, axis = 1)
    int_parameters[:,1:] += cumulative_distances[:,:-1] 
        
    #time to get the event points for each ray
    normed_m = momenta / np.linalg.norm(momenta, axis = 1)[:, np.newaxis]
    events_position[:,:, :] = dec_pos[:,np.newaxis,:] + int_parameters[:,:,np.newaxis] * normed_m[:, np.newaxis, :]
    return events_position, part_face_counts
    
