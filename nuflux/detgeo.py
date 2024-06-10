import numpy as np
import matplotlib.pyplot as plt
import DarkNews as dn
from scipy import interpolate
from scipy.stats import binom, multinomial, uniform
import copy

import importlib

import sys
sys.path.append("/home/luc/Research/BIN_MC/nuflux/detector_geometries")
import helpers
import useful_data
from DarkNews import const
from DarkNews import Cfourvec as Cfv

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

    def __init__(self, coord_object, parameters):
        self.cc = coord_object
        self.dec_pos = self.cc.p
        self.w = self.cc.weights.reshape((self.cc.weights.shape[0],1)) / np.sum(self.cc.weights) #normalized
        self.Geometry = parameters["geom"]
        self.Nmu = parameters["Nmu"]
        self.sample_size = self.dec_pos.shape[0]
        self.iterations = parameters["iterations"] #size of the interactions_points array, one more than distances and densities
        self.particle = parameters["particle"]
        if self.particle=='numu':
            self.momenta = self.cc.mnumu
            self.E = self.cc.Enumu
        elif self.particle=='nue':
            self.momenta = self.cc.mnue
            self.E = self.cc.Enue

        geom = importlib.import_module(self.Geometry)
        self.objects = geom.OBJECTS
        self.zbeginning = geom.zbeginning
        self.rmax = geom.rmax
        self.zending = geom.zending
        self.initials = geom.INITIALIZERS #should not contain minus_one
        self.decayer = geom.DECAYER #is the minus_one
        self.outside = geom.OUTSIDE # parent class of outside-going faces
        

        self.initialize_quantities()


    def find_info(self): 
        #first interaction points
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

        self.mask_accepted = ~self.mask_not_accepted
        for obj in self.initials:
            new_indices = obj.check_in(self.r_values,self.t,  self.mask_accepted)
            self.location[new_indices, 1] = obj.id
            self.intersection_points[new_indices, 1, :] = ip[new_indices,:]

        
        #at this point, there should be -1s (decaying in detector), 0s (did not reach detector), and ids of initials
        for neighbor in self.decayer[0].next_ids: #decayer is id -1; could be more added, just needa loop over the other decayers
            mask_decay = (self.location[:,1] == -1) & ~((self.dec_pos[:,2] > self.zbeginning) & (self.dec_pos[:,2] < 564) & (np.sqrt(self.dec_pos[:,0]**2 + self.dec_pos[:,1]**2) > 3) & (np.sqrt(self.dec_pos[:,0]**2 + self.dec_pos[:,1]**2) < 6e2))
            neigh = self.objects[neighbor]
            new_indices, ips = neigh.check_intersection(self.intersection_points[:,0,:], self.momenta, mask_decay)

            int_mask1 = (ips[:,2] < self.intersection_points[new_indices, count, 2])
            int_mask2 = (ips[:,2] > self.intersection_points[new_indices, 0, 2])

            mask_1 = int_mask1 & int_mask2
            accepted_ix = new_indices[mask_1]

            self.intersection_points[accepted_ix, count, :] = ips[mask_1]
            self.location[accepted_ix, 1] = neigh.id

        weird_decay_mask = (self.location[:,1] == -1)
        for neighbor in [18,17,19,20, 27,36]: #decayer is id -1; could be more added, just needa loop over the other decayers
            
            neigh = self.objects[neighbor]
            
            
            new_indices, ips = neigh.check_intersection(self.intersection_points[:,0,:], self.momenta, weird_decay_mask)

            mask_1 = (ips[:,2] < self.intersection_points[new_indices, count, 2]) & (ips[:,2] > self.intersection_points[new_indices, 0, 2])
            accepted_ix = new_indices[mask_1]

            self.intersection_points[accepted_ix, count, :] = ips[mask_1]
            self.location[accepted_ix, 1] = neigh.id

        #there should be no more -1s, just 0 - x
        count = 2
        while count < self.iterations:
            
            for obj in self.objects:
                
                part_mask = (self.location[:, count -1] == obj.id)
                particles = np.where(part_mask)[0]#these are indices of the particles at that face

                if obj in self.outside:
                    self.intersection_points[particles, count, :] = self.intersection_points[particles, count-1, :]
                    self.location[particles, count] = self.location[particles, count - 1]
                    continue
                
                if np.any(particles):
                    for neighbor in obj.next_ids:
                        neigh = self.objects[neighbor]
                        new_indices, ips = neigh.check_intersection(self.intersection_points[:,count - 1,:], self.momenta, part_mask)
                        
                        mask_1 = (ips[:,2] < self.intersection_points[new_indices, count, 2])& (ips[:,2] > self.intersection_points[new_indices, 0, 2])
                        accepted_ix = new_indices[mask_1]
                        self.intersection_points[accepted_ix, count, :] = ips[mask_1]
                        self.location[accepted_ix, count] = neigh.id
                        self.densities[accepted_ix, count-1] = obj.parent.density
            count+=1


    def get_probs(self):
        for i in range(self.iterations - 1):
            self.distances[:,i] = D3distance(self.intersection_points[:,i,:],self.intersection_points[:,i+1,:]) #these input arrays will be two dimensional, i.e., (sample_size, 3); distances will be i.e., (sample_size, 14)
        self.part_line_integrals = self.distances * self.densities # huge array of line_ints at each component for each particle
        self.line_integrals = np.sum(self.part_line_integrals, axis = 1)
        self.cs = get_cs(self.E, self.particle)
        self.probs = 1 - np.exp(-1*self.cs * self.line_integrals)

        self.factors = np.full(( self.sample_size,1), self.Nmu) * self.w

        self.counts = self.factors * self.probs.T.reshape((self.sample_size,1)) # (sample_size, 1)
        self.total_count = np.sum(self.counts)


    def get_event_positions(self):
        '''weights have already been given'''
        self.part_face_counts =  np.zeros(self.distances.shape) #counts at each face for each particle
        self.int_parameters = np.zeros(self.distances.shape) #t values starting at the last interaction point for each simulated interaction
        self.events_position = np.empty((self.sample_size, self.iterations - 1, 3))

        self.pweights = np.zeros(self.distances.shape) #these are probabilities, NOT WEIGHTS
        self.mask = (self.counts > 0)
        self.mask = self.mask.reshape((self.mask.shape[0],))
        temp = np.zeros((self.sample_size, self.iterations - 1)) #probs of interaction in each segment
        temp2 = np.ones((self.sample_size, self.iterations - 1)) #probs of non-interaction in each segment
        
        middle = self.cs[self.mask]
        middle2 = self.part_line_integrals[self.mask]
        temp2[self.mask, :] = np.exp(-1 *middle.reshape(middle.shape[0],1) * middle2)
        temp[self.mask,:] = 1 - temp2[self.mask,:] # probs of interaction on each segment

        mi = np.prod(temp2, axis = 1)
        mid = mi.reshape((mi.shape[0], 1)) / temp2
        self.pweights[self.mask,:] = temp[self.mask,:] * mid[self.mask,:]
        self.pweights[self.mask, :] = self.pweights[self.mask]
        max_index = np.argmax(self.pweights[self.mask], axis=1)

        self.pweights[self.mask, :] = self.pweights[self.mask] / self.pweights[self.mask, max_index].reshape((self.pweights[self.mask].shape[0], 1))

        mid  = np.sum(self.pweights, axis = 1)

        self.pweights[self.mask,:] = self.pweights[self.mask] / mid.reshape((mid.shape[0], 1))[self.mask] #good multinomial probabilities
        self.part_face_counts[:,:] = self.counts * self.pweights # these are the face counts (expected)
        

        #getting the distances to the event points for each ray
        temporary = np.empty((1,np.sum(self.mask), self.iterations - 1))
        self.scales = self.distances[self.mask]
        temporary[:,:, :] = [uniform.rvs(loc = 0, scale = self.scales)]
        self.int_parameters[self.mask, :] = temporary
        self.cumulative_distances = np.cumsum(self.distances, axis = 1)
        self.int_parameters[:,1:] += self.cumulative_distances[:,:-1] 
        
        #time to get the event points for each ray
        self.normed_m = self.momenta / np.linalg.norm(self.momenta, axis = 1)[:, np.newaxis]
        self.events_position[:,:, :] = self.dec_pos[:,np.newaxis,:] + self.int_parameters[:,:,np.newaxis] * self.normed_m[:, np.newaxis, :]
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
    def run(self):
        if (self.particle == "nue") | (self.particle == "numu"):
            self.find_info()
            self.get_probs()
            self.get_event_positions()
            print("{:.3g} {} events".format(self.total_count, self.particle))
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







