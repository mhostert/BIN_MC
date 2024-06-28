import importlib
import numpy as np
from DarkNews import Cfourvec as Cfv
import copy
from numba import njit
import sys
AVOGADRO = 6.02214e23
class face:
    def __init__(self, density):
        self.density = density #density in which they're going

class cap:

    def __init__(self,parent,id, next_ids, zpos, rbeg, rend): #what to do for the end
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.zpos = zpos
        self.rbeg = rbeg
        self.rend = rend

    def check_intersection(self, position,momenta, mask): #position and momentum will be in (sample_size, 3) shape, 0,0,0 is center of detector
        return cap_check_i(self.zpos, self.rbeg, self.rend, position, momenta, mask)


class barrel:

    def __init__(self, parent, id, next_ids,rpos, zbeg, zend):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.rpos = rpos
        self.zbeg = zbeg
        self.zend = zend

    def check_intersection(self, position,momenta, mask):
        indices = np.where(mask)[0] #array of indices of particles that we will consider
        a,b,c = barrel_get_pols(self.rpos, position, momenta, indices)

        coeffs = np.vstack((a,b,c)).T #( indices size), 3)
        roots=np.empty((a.shape[0], 2))
        for i,poly in enumerate(coeffs):
            root = np.roots(poly)
            root1 = np.zeros(2)
            if isinstance(root[0],complex):
                 root1[0] = -1
            else:
                root1[0] = root[0]  
            if isinstance(root[1], complex):
                root1[1] = -1
            else:
                root1[1] = root[1]    
            roots[i,:] = root1# (indices size, 2) THIS MIGHT NOT ALWAYS have size two
        info =  barrel_check_i(self.zbeg, self.zend, position, momenta, indices, roots)
        return info[0], info[1]
                    

class conic:

    def __init__(self, parent, id, next_ids, tan_theta, zbeg, zend, rsmall, direction):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.tan_theta = tan_theta # positive
        self.rsmall  = rsmall
        self.zbeg = zbeg # big opening
        self.zend = zend # smallest radius point
        self.direction = direction
        if self.direction=='towards':
            self.zcenter = zend + rsmall/tan_theta
        else:
            self.zcenter = zbeg - rsmall/tan_theta

    def check_intersection(self, position,momenta, mask):
        info =  conic_check_i(self.tan_theta, self.zcenter, self.zbeg, self.zend, position, momenta, mask)
        return info[0], info[1]
        

class initializer:

    def __init__(self,parent, id, next_ids, rsmall, rbig):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.rsmall = rsmall
        self.rbig = rbig
    
    def check_in(self, r,t, mask):
        return init_check_in(self.rbig, self.rsmall, r, t.reshape((t.shape[0],)), mask)

class decayer:

    def __init__(self, parent, id, next_ids):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids

@njit
def cap_check_i(zpos, rbeg, rend, position, momenta, mask):
    indices = np.where(mask)[0] #array of indices of particles that we will consider
    delta_z = zpos - position[indices][:,2]
    t = delta_z / momenta[indices][:,2]
    ip = position[indices] + t[:, np.newaxis] * momenta[indices] #interaction points
    r_values = np.sqrt(ip[:,0]**2 + ip[:,1]**2) 
    mask_new = (r_values < rend) & (r_values > rbeg) & (t > 0) #particles that have a correct ip
    kept_indices = indices[np.where(mask_new)[0]] #indexing the indices to get the number of the particles
    return kept_indices, ip[mask_new]

@njit
def init_check_in(rbig, rsmall, r, to, mask):
    indices = np.where(mask)[0]
    new_mask = (r[indices] < rbig) & (r[indices]> rsmall) & (to[indices] > 0)
    return indices[np.where(new_mask)[0]]

@njit
def conic_check_i(tan_theta, zcenter, zbeg, zend, position, momenta, mask):
    indices = (np.where(mask)[0]) #indices of particles that we will consider
    delta_z = zcenter - position[indices,2]
    a = -1 * momenta[indices,2]**2 * tan_theta**2 + momenta[indices,1]**2 + momenta[indices,0]**2
    b = 2*position[indices,0] * momenta[indices, 0] + 2*position[indices,1]*momenta[indices,1] + tan_theta**2 * 2*zcenter*momenta[indices,2] - tan_theta**2 * 2*position[indices,2]*momenta[indices,2]
    c = -1*tan_theta**2 * delta_z**2 + position[indices,0]**2 + position[indices,1]**2
        
    coeffs = np.vstack((a,b,c)).T # (sample_size, 3)
    
    roots=np.empty((coeffs.shape[0], 2))
        
    for i,poly in enumerate(coeffs):
        root = np.roots(poly)
        if isinstance(root[0],complex):
            root[0] = -1
        if isinstance(root[1], complex):
            root[1] = -1
        roots[i,:] = root# (indices size, 2) THIS MIGHT NOT ALWAYS have size two

    ip_1 = position[indices] + roots[:,0][:, np.newaxis] * momenta[indices]
    ip_2 = position[indices] + roots[:,1][:, np.newaxis] * momenta[indices]
    #conditions: assert r between rsmall and (zcenter - zbeg)*tan_theta, z between zend and zbeg (first), root is positive
    new_mask_1 = (ip_1[:,2] > zbeg) & (ip_1[:,2] < zend) & (np.round(roots[:,0], decimals = 12) > 0)
    new_mask_2 = (ip_2[:,2] > zbeg) & (ip_2[:,2] < zend) & (np.round(roots[:,1], decimals = 12) > 0)
            
    t = roots.copy() #need t to get the correct root
    any = new_mask_1 | new_mask_2
    doubles = new_mask_1 & new_mask_2
            

    t[new_mask_2,0] = roots[new_mask_2,1]
    if np.any(doubles):
        t[doubles,0] = np.min(roots[doubles])
            
    ip = position[indices][any] + t[:,0][any, np.newaxis] * momenta[indices][any]
    kept_indices = indices[np.where(any)[0]] # indexing the indices to get the correct particle numbers (ids)
        
    return kept_indices, ip

@njit
def barrel_get_pols(rpos, position, momenta, indices):
    a = (momenta[indices,0])**2 + (momenta[indices,1])**2
    b = 2 * (position[indices,0]*momenta[indices,0] + position[indices,1]*momenta[indices,1])
    c = position[indices,1]**2 + position[indices,0]**2 - rpos**2
    return a,b,c

@njit
def barrel_check_i(zbeg, zend, position,momenta, indices, roots):

    ip_1 = position[indices,2] + roots[:,0]*momenta[indices,2] # only z
    ip_2 = position[indices,2] + roots[:,1]*momenta[indices,2] # only z

    new_mask_1 = (ip_1 > zbeg) & (ip_1 < zend) & (np.round(roots[:,0],decimals =12) > 0)
    new_mask_2 = (ip_2 > zbeg) & (ip_2 < zend) & (np.round(roots[:,1], decimals = 12) > 0)
    t = roots.copy() #need t to get the correct root
    any = new_mask_1 | new_mask_2

    doubles = new_mask_1 & new_mask_2
    
    t[new_mask_2,0] = roots[new_mask_2,1]

    if np.any(doubles):
        t[doubles,0] = np.min(roots[doubles])
            
    ip = position[indices][any] + t[any,0][:, np.newaxis]*momenta[indices][any]
                    
    kept_indices = indices[np.where(any)[0]] # indexing the indices to get the correct particle numbers (ids)
            
    return kept_indices, ip


class cc:

    def __init__(self, sim):
        
        self.sim = sim
        self.weights = self.sim.w / np.sum(self.sim.w)
        self.Enumu = self.sim.Enumu
        self.Enue =self.sim.Enue
        self.Nmu = self.sim.N_mu

        old_m = self.sim.pnumu_ar[:,1:4] # after initial rotation within program, no need to translate
        old_m2 = self.sim.pnue_ar[:,1:4]
        po = self.sim.pos_at.T #need to translate
        old_p = np.copy(po)
        old_p[:,2] = po[:,2] - self.sim.Racc

        self.mnumu = np.empty(old_m.shape)
        self.mnue = np.empty(old_m2.shape)
        self.p = np.empty(old_p.shape)

        old = [old_m, old_m2, old_p]
        new = [self.mnumu, self.mnue, self.p]
        for i, coord in enumerate(new):
            coord[:,0] = -1 * old[i][:,0]
            coord[:,1] = old[i][:,2]
            coord[:,2] = old[i][:, 1]
        
        self.pnumu = np.empty((self.mnumu.shape[0], 4))
        self.pnue = np.empty((self.mnue.shape[0], 4))

        self.pnumu[:,0] = self.Enumu
        self.pnue[:,0] = self.Enue
        self.pnumu[:,1:] = self.mnumu
        self.pnue[:,1:] = self.mnue
    
    def completely_circular(self):
        return self

    def straight_segment_at_detector(self, ssl):
        #ssl is straight segment length, ENTIRE, not half
        self.L = ssl/2
        self.d = np.sqrt(self.sim.Racc**2 - self.L**2)

        #mask - change everything that is on the straight segment
        self.sample_size = self.p.shape[0]
        self.mask = (self.p[:, 2] > -1*self.L) & (self.p[:,2] < self.L) & (self.p[:,1] > - self.sim.Racc)


        self.new_p = np.copy(self.p)
        self.new_mnumu = np.copy(self.mnumu)
        self.new_mnue = np.copy(self.mnue)

        #lower dimension quantities
        self.dphi = np.arcsin(self.p[:,2] / self.sim.Racc)
        self.lengths = self.dphi * self.sim.Racc
        self.new_p[self.mask,2] = self.lengths[self.mask]
        self.tantheta = self.L * np.sqrt(1 / (self.sim.Racc**2 - self.L**2))
    
        self.new_p[self.mask,0] = self.p[self.mask, 0]
        self.new_p[self.mask,1] = np.zeros((np.sum(self.mask), )) #approximation; these don't actually have y= 0, they are a little off the beam

        #these are on lower dimension already
        self.maskleft = (self.dphi  < 0) & (self.mask)
        self.maskright = (self.dphi > 0) & (self.mask)


        
        self.new_mnumu[self.maskleft, :] = Cfv.rotationx(self.pnumu[self.maskleft], -1*(self.dphi[self.maskleft]))[:,1:] #lower dimension
        self.new_mnumu[self.maskright, :] = Cfv.rotationx(self.pnumu[self.maskright], (self.dphi[self.maskright]))[:,1:] #lower dimension
        self.new_mnue[self.maskleft, :] = Cfv.rotationx(self.pnue[self.maskleft], -1 * (self.dphi[self.maskleft]))[:, 1:] #lower dimension
        self.new_mnue[self.maskright, :] = Cfv.rotationx(self.pnue[self.maskright], (self.dphi[self.maskright]))[:, 1:] #lower dimension

        self.notmask = ~self.mask
        self.new_p[self.notmask, 1] = self.p[self.notmask, 1] + self.sim.Racc - self.L/self.tantheta

        self.p = np.copy(self.new_p)
        self.mnumu = np.copy(self.new_mnumu)
        self.mnue = np.copy(self.new_mnue)

        return self
    
    def completely_linear(self):
        pass


class material:
    def __init__(self):
        self.N = 0

class subs(material):
    def __init__(self, density, am, A):
        super().__init__()
        self.density = density
        self.am = am
        self.A = A
        self.N = AVOGADRO * self.density / self.am * self.A

class comp(material):
    def __init__(self, table):
        super().__init__()
        for row in table:
            self.N += row[0].N * row[1] #row[0] is N of an element; row[1] is the percentage of it that occupies the total material