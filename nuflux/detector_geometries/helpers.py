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
        a,b,c = barrel_get_pols(self.rpos, position[indices], momenta[indices])

        '''coeffs = np.vstack((a,b,c)).T #( indices size), 3)
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
            roots[i,:] = root1# (indices size, 2) THIS MIGHT NOT ALWAYS have size two'''
        
        info = get_sols(a,b,c, self.zbeg, self.zend, position[mask], momenta[mask], mask)
            
        #info =  barrel_check_i(self.zbeg, self.zend, position, momenta, indices, roots)
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
        
        a,b,c = conic_get_pols(self.zcenter, self.tan_theta, position[mask], momenta[mask])
        info = get_sols(a,b,c, self.zbeg, self.zend, position[mask], momenta[mask], mask)
        
        #info =  conic_check_i(self.tan_theta, self.zcenter, self.zbeg, self.zend, position, momenta, mask)
        return info[0], info[1]
        

class initializer:

    def __init__(self,parent, id, next_ids, rsmall, rbig):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.rsmall = rsmall
        self.rbig = rbig
    
    def check_in(self, r, mask):
        indices = np.where(mask)[0]
        new_mask =  (r[indices] < self.rbig)
        return indices[np.where(new_mask)[0]]

class decayer:

    def __init__(self, parent, id, next_ids, last):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.last = last
    
    def check_in(self, neighbor, zs, mask_decay, loc_mask):
        indices = np.where(mask_decay)[0]
        if neighbor.id == self.last:
            return indices[np.where(loc_mask)[0]], loc_mask
        local_mask = (zs > 0) & ((zs < neighbor.zend) & (loc_mask))
        return indices[np.where(local_mask)[0]], local_mask
        

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
def barrel_get_pols(rpos, position, momenta):
    a = (momenta[:,0])**2 + (momenta[:,1])**2
    b = 2 * (position[:,0]*momenta[:,0] + position[:,1]*momenta[:,1])
    c = position[:,1]**2 + position[:,0]**2 - rpos**2
    return a,b,c

@njit
def conic_get_pols(zcenter, tan_theta, position, momenta):
    delta_z = zcenter - position[:,2]
    a = -1 * momenta[:,2]**2 * tan_theta**2 + momenta[:,1]**2 + momenta[:,0]**2
    b = 2*position[:,0] * momenta[:, 0] + 2*position[:,1]*momenta[:,1] + tan_theta**2 * 2*zcenter*momenta[:,2] - tan_theta**2 *2*position[:,2]*momenta[:,2]
    c = -1*tan_theta**2 * delta_z**2 + position[:,0]**2 + position[:,1]**2
    return a,b,c

@njit
def get_roots(a,b,c):
    disc = b**2 - 4*a*c
    roots = np.empty((a.shape[0], 2))
    
    mask_1 = (disc > 0)
    mask_2 = (disc == 0)
    mask_3 = (disc < 0)
    
    div = 2*a
    r = -1*b/div
    
    r2 = np.sqrt(disc[mask_1])/div[mask_1]
    roots[mask_1,0] = r[mask_1] + r2
    roots[mask_1,1] = r[mask_1] - r2
    
    roots[mask_2,0] = r[mask_2]
    roots[mask_2,1] = r[mask_2]
    
    roots[mask_3, :] = -1
    return roots

@njit
def get_sols(a,b,c,zbeg, zend, position, momenta, mask):
    indices = np.where(mask)[0]
    #input position and momenta already sliced to indices
    #need to get roots in efficient way
    roots = get_roots(a,b,c)
    
    ip_1 = position[:,2] + roots[:,0]*momenta[:,2] # only z
    ip_2 = position[:,2] + roots[:,1]*momenta[:,2] # only z

    new_mask_1 = (ip_1 > zbeg) & (ip_1 < zend) & (np.round(roots[:,0],decimals =12) > 0)
    new_mask_2 = (ip_2 > zbeg) & (ip_2 < zend) & (np.round(roots[:,1], decimals = 12) > 0)
    t = roots.copy() #need t to get the correct root
    any = new_mask_1 | new_mask_2

    doubles = new_mask_1 & new_mask_2
    
    t[new_mask_2,0] = roots[new_mask_2,1]

    if np.any(doubles):
        t[doubles,0] = np.min(roots[doubles])
            
    ip = position[any] + t[any,0][:, np.newaxis]*momenta[any]
                    
    kept_indices = indices[np.where(any)[0]] # indexing the indices to get the correct particle numbers (ids)
    
    return kept_indices, ip


class cc:

    def __init__(self, R, w, sample_size, Enumu, Enue, N_mu, pnumu_ar, pnue_ar, pos_at):
        
        self.Racc = R
        nw = np.copy(w)
        del w
        self.weights = nw / np.sum(nw)
        self.sample_size = sample_size
        self.Nmu = N_mu

        old_m = pnumu_ar[:,1:4] # after initial rotation within program, no need to translate
        old_m2 = pnue_ar[:,1:4]
        po = np.copy(pos_at).T #need to translate
        
        old_p = np.copy(po)
        old_p[:,2] = po[:,2] - self.Racc

        mnumu = np.empty(old_m.shape)
        mnue = np.empty(old_m2.shape)
        self.p = np.empty(old_p.shape)

        old = [old_m, old_m2, old_p]
        new = [mnumu, mnue, self.p]
        for i, coord in enumerate(new):
            coord[:,0] = -1 * old[i][:,0]
            coord[:,1] = old[i][:,2]
            coord[:,2] = old[i][:, 1]
        
        self.pnumu = np.empty((self.sample_size, 4))
        self.pnue = np.empty((self.sample_size, 4))
        
        self.pnumu[:,0] = np.copy(Enumu)
        self.pnue[:,0] = np.copy(Enue)
        self.pnumu[:,1:] =mnumu
        self.pnue[:,1:] = mnue

    def straight_segment_at_detector(self, f, h):
        #ssl is straight segment length, ENTIRE, not half
        dphi = np.arcsin(self.p[:,2] / self.Racc)
<<<<<<< HEAD
        
        if f==0:
            mask_not_accepted = (dphi > np.pi/100) | (self.p[:,1] < -1*self.Racc) | (dphi < -1*np.pi/8)
            mask_acc = ~mask_not_accepted
            
        elif f==-1:
            return
=======
        
        if f==0:
            mask_not_accepted = (dphi > np.pi/100) | (self.p[:,1] < -1*self.Racc) | (dphi < -1*np.pi/8)
            mask_acc = ~mask_not_accepted
        
>>>>>>> 502e55835199c64fe1919866c5dab403db34d53d
        else:
            self.Lc = 2*  np.sqrt((-1*h**2 + np.sqrt(h**4 + 4*h**2 * self.Racc**2))/2)
            L = (f * self.Lc)/2
            d = np.sqrt(self.Racc**2 - L**2)

            #mask - change everything that is on the straight segment
    
            mask = (self.p[:, 2] > -1*L) & (self.p[:,2] < L) & (self.p[:,1] > - self.Racc)


            new_p = np.copy(self.p)
            new_mnumu = np.copy(self.pnumu[:,1:])
            new_mnue = np.copy(self.pnue[:,1:])

        #lower dimension quantities
            lengths = dphi * self.Racc
            new_p[mask,2] = lengths[mask]
            tantheta = L * np.sqrt(1 / (self.Racc**2 - L**2))
        
            new_p[mask,0] = self.p[mask, 0]
            new_p[mask,1] = np.zeros((np.sum(mask), )) #approximation; these don't actually have y= 0, they are a little off the beam

        #these are on lower dimension already
            maskleft = (dphi  < 0) & (mask)
            maskright = (dphi > 0) & (mask)


        
            new_mnumu[maskleft, :] = Cfv.rotationx(self.pnumu[maskleft], -1*(dphi[maskleft]))[:,1:] #lower dimension
            new_mnumu[maskright, :] = Cfv.rotationx(self.pnumu[maskright], (dphi[maskright]))[:,1:] #lower dimension
            new_mnue[maskleft, :] = Cfv.rotationx(self.pnue[maskleft], -1 * (dphi[maskleft]))[:, 1:] #lower dimension
            new_mnue[maskright, :] = Cfv.rotationx(self.pnue[maskright], (dphi[maskright]))[:, 1:] #lower dimension

            notmask = ~mask
            new_p[notmask, 1] = self.p[notmask, 1] + self.Racc - L/tantheta

            self.p = np.copy(new_p)
            self.pnumu[:,1:] = np.copy(new_mnumu)
            self.pnue[:,1:] = np.copy(new_mnue)
        
            #to remove all decays that are too far ~15/16
<<<<<<< HEAD
            co = f*self.Lc/self.Racc
            if (co >= 1/2) and (f>=2): #arbitrary; should be very careful
                if co>=2:
                    raise ValueError('This length for straight segment is way too big. In fact, it is over twice the radius.')
                mask_not_accepted = (dphi > np.pi/100) | (self.p[:,1] < -1*L/tantheta)| (dphi < -1*np.pi/4*co)
            else:
                mask_not_accepted = (dphi > np.pi/100) | (self.p[:,1] < -1*L/tantheta)| (dphi < -1*np.pi/8)
            
=======
            mask_not_accepted = (dphi > np.pi/100) | (self.p[:,1] < -1*L/tantheta) | (dphi < -1*np.pi/8)
>>>>>>> 502e55835199c64fe1919866c5dab403db34d53d
            mask_acc = ~mask_not_accepted
        
        self.p = self.p[mask_acc]
        self.pnumu = self.pnumu[mask_acc]
        self.pnue = self.pnue[mask_acc]
        self.weights = self.weights[mask_acc]
        self.sample_size = np.sum(mask_acc)
    
    def completely_linear(self):
        pass

    def clear_mem(self):
        deletables = ['Racc', 'p', 'pnumu', 'pnue', 'weights', 'sample_size', 'Nmu']
        for att in deletables:
            if hasattr(self, att):
                delattr(self, att)

class material:
    def __init__(self):
        self.N = 0
        self.density = 0

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
            self.density += row[0].density * row[1]
            self.N += row[0].N * row[1] #row[0] is N of an element; row[1] is the percentage of it that occupies the total material