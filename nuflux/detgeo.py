import numpy as np
import matplotlib.pyplot as plt
import DarkNews as dn
from scipy.stats import uniform
import copy

import importlib

import sys
import os
cfp = os.path.dirname(os.path.abspath(__file__))
td = os.path.join(cfp,'detector_geometries')
sys.path.append(td)
cfp = os.path.dirname(os.path.abspath(__file__))
sys.path.append(cfp)
import helpers
import useful_data
import time
from numba import njit, jit
from prettytable import PrettyTable
from memory_profiler import profile
import data
import gc

def D3distance(point1, point2):
    return np.sqrt((point1[:,0] - point2[:,0])**2 + (point1[:,1] - point2[:,1])**2 + (point1[:,2] - point2[:,2])**2)

def get_cs(E, part):
    sigmanue,sigmanuebar,sigmanumu,sigmanumubar= useful_data.cs_interp()
    if part == "nue":
        return sigmanue(E)
    elif part == "numu":
        return sigmanumubar(E)

##@profile
def get_quantities(param, N_evals):
    dt = list(data.get_particles(param, N_evals))
    sim =  helpers.cc(R = dt[0], w =  dt[1], sample_size = dt[2], Enumu = dt[3], Enue = dt[4], N_mu = dt[5], pnumu_ar = dt[6], pnue_ar = dt[7], pos_at = dt[8])

    del dt
    return sim

class SimulateDetector():

    ##@profile
    def __init__(self, coord_object, f, geom, particle):
        self.time = np.zeros(7)
        self.f = f
        self.Geometry = geom
        geom = importlib.import_module(self.Geometry)
        self.objects = geom.OBJECTS
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
        del geom 
        
        self.cc = copy.deepcopy(coord_object)
        self.cc.straight_segment_at_detector(self.f, self.rmax)
        
        self.w = self.cc.weights.reshape((self.cc.weights.shape[0],1)) # already normalized (does not include all decays!!!)
        self.Nmu = self.cc.Nmu
        self.sample_size = self.cc.sample_size
        self.particle = particle
        if self.particle=='numu':
            self.momenta = self.cc.pnumu[:,1:]
            self.E = self.cc.pnumu[:,0]
        elif self.particle=='nue':
            self.momenta = self.cc.pnue[:,1:]
            self.E = self.cc.pnue[:,0]

        

        self.initialize_quantities()

    #@profile
    def find_info(self): 
        
        #first interaction points
        count = 1

        delta_z = self.zbeginning - self.intersection_points[:,0,2]
        to = delta_z / self.momenta[:,2]
        ip = self.intersection_points[:,0,:] + to[:, np.newaxis] * self.momenta #intersection points
        r_values = np.sqrt(ip[:,0]**2+ip[:,1]**2)
        
        cond_1 = ~(((self.intersection_points[:,0,2] < self.zending) & (self.intersection_points[:,0,2] > self.zbeginning) & (self.intersection_points[:,0,1] > -1 * self.rmax) & (self.intersection_points[:,0,1] < self.rmax)))
        cond_2 = ((r_values > self.rmax) | (to < 0))
        mask_not_accepted = cond_2 & cond_1
        
        not_accepted_indices = np.where(mask_not_accepted)[0] #indices of this new array 
        self.location[not_accepted_indices, 1] = 0
        self.intersection_points[not_accepted_indices, 1, :] = self.intersection_points[not_accepted_indices, 0, :]
        self.time[1] = time.time()

        #resolved! efficient now
        mask_accepted = (~mask_not_accepted) & (to > 0)
        for obj in self.initials:
            mask_accepted = mask_accepted & (self.location[:,1]==-1)
            new_indices = obj.check_in(r_values, mask_accepted) #gives indices (numbers) of those that go in one of the initials

            self.location[new_indices, 1] = obj.id
            self.intersection_points[new_indices, 1, :] = ip[new_indices,:]
        
        
        del self.initials
        
        #at this point, there should be -1s (decaying in detector), 0s (did not reach detector), and ids of initials
        mask_decay = (self.location[:,1] == -1) & ~(np.sqrt(self.intersection_points[:,0,0]**2 + self.intersection_points[:,0,1]**2) > self.rbp)
        
        #this one takes a bunch of time! will replace by common cylinder solving. Need to get new ips; new face id; own density

        a,b,c = helpers.barrel_get_pols(self.rbp, self.intersection_points[mask_decay,0,:], self.momenta[mask_decay])

        roots = helpers.get_roots(a,b,c)

        new_mask_1 = (np.round(roots[:,0],decimals =12) > 0)
        new_mask_2 = (np.round(roots[:,1], decimals = 12) > 0)
        t = roots.copy() #need t to get the correct root
    
        doubles = new_mask_1 & new_mask_2
    
        t[new_mask_2,0] = roots[new_mask_2,1]

        if np.any(doubles):
            t[doubles,0] = np.min(roots[doubles])
        
        
        
        ip2 = self.intersection_points[mask_decay,0,:] + t[:,0][:, np.newaxis]*self.momenta[mask_decay]
        for neighb in self.decayer[0].next_ids:
            neighbor = self.objects[neighb]
            loc_mask  = (self.location[mask_decay, 1]==-1)
            new_indices, new_mask = self.decayer[0].check_in(neighbor, ip2[:,2], mask_decay, loc_mask) #gives indices (numbers)
            self.location[new_indices, 1] = neighb
            ips = ip2[new_mask,:]
            self.intersection_points[new_indices, 1, :] = ips


        

        #treating decays within the detector
        weird_decay_mask = (self.location[:,1] == -1)
        self.update_intersections(self.decayer[1], count, weird_decay_mask)
        self.location[weird_decay_mask,0] = -2
        
        del self.decayer
        
        #there should be no more -1s, just 0 - x
        count = 2

        self.time[2] = time.time()
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
        del self.objects
        del self.outside
        del self.outside_ids
    
    #@profile
    def get_probs(self):

        for i in range(self.iterations - 1):
            self.distances[:,i] = D3distance(self.intersection_points[:,i,:],self.intersection_points[:,i+1,:]) #these input arrays will be two dimensional, i.e., (sample_size, 3); distances will be i.e., (sample_size, 14)
        
        self.dec_pos = self.intersection_points[:,0,:]
        del self.intersection_points
        self.part_line_integrals = self.distances * self.densities # huge array of line_ints at 
        line_integrals = np.sum(self.part_line_integrals, axis = 1)
        self.cs = get_cs(self.E, self.particle)
        probs = 1 - np.exp(-1*self.cs * line_integrals)

        factors = np.full((self.sample_size,1), self.Nmu) * self.w
        del self.w
        self.counts = factors * probs.T.reshape((self.sample_size,1)) # (sample_size, 1)
        self.total_count = np.sum(self.counts)
        self.time[4] = time.time()
        del self.densities 
    
    #@profile
    def get_event_positions(self):
        '''weights have already been given'''

        self.mask = (self.counts>0)
        self.mask = self.mask[:, 0]
        
        self.momenta = self.momenta[self.mask]
        self.dec_pos = self.dec_pos[self.mask]
        self.distances = self.distances[self.mask]
        self.part_line_integrals = self.part_line_integrals[self.mask]
        self.counts = self.counts[self.mask]
        self.location = self.location[self.mask]
        self.E = self.E[self.mask]
        self.cs = self.cs[self.mask]
        
        temporary = np.empty((1,np.sum(self.mask), self.iterations - 1)) #same size as other big arrays now
        temporary[:,:, :] = [uniform.rvs(loc = 0, scale = self.distances)]
        self.time[5] = time.time()
        
        
        
        self.events_position, self.part_face_counts = get_probs_njit(self.dec_pos, self.momenta, self.distances, self.iterations, self.counts, self.cs, self.part_line_integrals, temporary)
        self.time[6] = time.time()
        
        return

    ##@profile
    def initialize_quantities(self):

        self.intersection_points = np.full((self.sample_size, self.iterations, 3), 1e4) #1e4 is arbitrary
        self.densities = np.zeros((self.sample_size, self.iterations - 1))
        self.location = np.full((self.sample_size, self.iterations), -1) #starting at initials
        self.intersection_points[:,0,:] = self.cc.p
        self.distances = np.zeros((self.sample_size, self.iterations - 1))


    def update_params(self):

        if self.particle =='numu':
            self.momenta = self.cc.pnumu[:,1:]
            self.E = self.cc.pnumu[:,0]

        elif self.particle == 'nue':
            self.momenta = self.cc.pnue[:,1:]
            self.E = self.cc.pnue[:,0]

    
    def update_intersections(self, obj, count, mask):
 
        for neighbor in obj.next_ids:
            neigh = self.objects[neighbor]
            new_indices, ips = neigh.check_intersection(self.intersection_points[:,count - 1,:], self.momenta, mask)
                        
            mask_1 = (ips[:,2] < self.intersection_points[new_indices, count, 2])& (ips[:,2] > self.intersection_points[new_indices, count - 1, 2])
            accepted_ix = new_indices[mask_1]
            self.intersection_points[accepted_ix, count, :] = ips[mask_1]
            self.location[accepted_ix, count] = neigh.id
            self.densities[accepted_ix, count-1] = obj.density

    #@profile
    def run(self):
        if (self.particle == "nue") | (self.particle == "numu"):
            self.time[0] = time.time()
            self.find_info()
            self.get_probs()
            self.get_event_positions()
            
            print(f'sim time: {(self.time[6] - self.time[0]):.3g} ({(self.time[1] - self.time[0]):.3g} init; {(self.time[2] - self.time[1]):.3g} init obj; {(self.time[3] - self.time[2]):.3g} obj; {(self.time[4] - self.time[3]):.3g} get_probs; {(self.time[5] - self.time[4]):.3g} uniform rvs; {(self.time[6] - self.time[5]):.3g} event pos)')
            
            
            if self.particle =='numu':
                return self, None
            else:
                return None, self
        
        elif self.particle=='both':
            self.particle="numu"
            self.update_params()
            self.run()
            _, sim2 = SimulateDetector(self.cc, 0, self.Geometry, 'nue').run()
            
            
            
            self.get_face_counts('both', sim2)
            
            self.clear_mem()
            sim2.clear_mem()
            
            return self, sim2
        
        else:
            raise ValueError("No particle of that name incldued in the detector!")
    
    ##@profile
    def get_face_counts(self, arg, sim2):

        # Example usage within your class
        self.facecounts = calculate_facecounts(self.face_dict, self.location, self.part_face_counts)

        if arg == 'both':
            self.facecounts2 = calculate_facecounts(sim2.face_dict, sim2.location, sim2.part_face_counts)

        names = list(self.face_dict.keys())
        names.append('TOTAL')
        
        data = [[self.facecounts2[key], self.facecounts[key], self.facecounts[key] + self.facecounts2[key]] for key in self.face_dict.keys()]
        data.append([sum(list(self.facecounts2.values())), sum(list(self.facecounts.values())), sum(list(self.facecounts2.values())) + sum(list(self.facecounts.values()))])
        table = PrettyTable()
        table.field_names = ['Detector Parts', 'ν_e events', 'ν_μ events ', 'Total events']
        
        for name, row in zip(names, data):
            formatted_row = [f"{x:.3e}" for x in row]
            table.add_row([name] + formatted_row)

        print(table)
    
    ##@profile
    def clear_mem(self):
        deletables=['cc','momenta','location','distances','counts','part_line_integrals','time','f','Geometry','zbeginning','rmax','rbp','zending','iterations','Nmu','sample_size','particle', 'face_dict', 'facecounts','facecounts2']
        for att in deletables:
            if hasattr(self, att):
                delattr(self, att)
               


@jit(nopython = False, forceobj=True)
def calculate_facecounts(face_dict, location, part_face_counts):
    facecounts = {}
    locs = location[:, :-1]
    
    for key, faces in face_dict.items():
        mask = np.isin(locs, faces)
        facecounts[key] = np.sum(part_face_counts[mask])

    return facecounts

#can't njit this since a bunch of np operations along axes
#@profile
def get_probs_njit(dec_pos, momenta, distances, iterations,counts, cs, part_line_integrals, temporary):
    #initialize
    #pweights these are probabilities, NOT WEIGHTS
    
    temporary = temporary.reshape((temporary.shape[1], temporary.shape[2]))
    temp2 = np.ones((distances.shape[0], iterations - 1)) #probs of non-interaction in each segment
    

    temp2 = np.exp(-1 *cs[:, np.newaxis] * part_line_integrals)
    mi = np.prod(temp2, axis = 1)
    
    #to njit
    pweights = get_events_njit1(temp2, mi[:, np.newaxis])
    max_index = np.argmax(pweights[:], axis=1)
 
    mask = np.ones_like(pweights, dtype=bool)[:,0]
    # Fetch maximum values from pweights using max_index
    max_p = pweights[mask, max_index]

    # Ensure max_p is reshaped correctly for broadcasting
    max_p = max_p[:, np.newaxis]  # Convert to column vector

    pweights = pweights / max_p
    mid  = np.sum(pweights, axis = 1)

    

    cumulative_distances = np.cumsum(distances, axis = 1)
    normed_m = momenta / np.linalg.norm(momenta, axis = 1)[:, np.newaxis]
    dec_pos = dec_pos[:,np.newaxis,:]
    normed_m = normed_m[:,np.newaxis, :]

    
    info =  get_events_njit2(distances.shape, pweights, mid[:,np.newaxis], counts, dec_pos, normed_m, cumulative_distances, temporary)

    return info[0], info[1]


@njit
def get_events_njit1(temp2, mi):
    
    res = (1 - temp2) * mi / temp2
    
    return res


@njit
def get_events_njit2(shape, pweights, mid, counts, dec_pos, normed_m, cumulative_distances, temporary):
    # shape is (:, 25); pweights is (:, 25); mid is (:,1); mask is (:,); counts is (:,1); dec_pos is (:, 1, 3); normed_m is [:, 1, 3]; cumulative_distances is [:,25]; temporary is [:, 25]
    
    pweights = pweights / mid #good multinomial probs

    temporary[:,1:] += cumulative_distances[:,:-1]
    
    part_face_counts = counts * pweights # these are the face counts (expected)
    
    events_position = dec_pos + temporary[:,:,np.newaxis] * normed_m
    
    return events_position, part_face_counts


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

@profile
def che():
    print('mem')