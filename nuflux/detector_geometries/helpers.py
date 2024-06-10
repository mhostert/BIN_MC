import importlib
import numpy as np

class face:
    def __init__(self, density):
        self.density = density #density in which they're going

    class cap:

        def __init__(self,parent,id, next_ids, zpos, rbeg, rend): #what to do for the end
            self.parent = parent
            self.id = id
            self.next_ids = next_ids
            self.zpos = zpos
            self.rbeg = rbeg
            self.rend = rend

        def check_intersection(self, position,momenta, mask): #position and momentum will be in (sample_size, 3) shape, 0,0,0 is center of detector
            indices = np.where(mask)[0] #array of indices of particles that we will consider

            delta_z = self.zpos - position[indices][:,2]
            t = delta_z / momenta[indices][:,2]
            ip = position[indices] + t.reshape((t.shape[0],1)) * momenta[indices] #interaction points
            r_values = np.sqrt(ip[:,0]**2 + ip[:,1]**2)
            
            mask_new = (r_values < self.rend) & (r_values > self.rbeg) & (t > 0) #particles that have a correct ip
            kept_indices = indices[np.where(mask_new)[0]] #indexing the indices to get the number of the particles
            
            return kept_indices, ip[mask_new]


    class barrel:

        def __init__(self, parent, id, next_ids,rpos, zbeg, zend):
            self.parent = parent
            self.id = id
            self.next_ids = next_ids
            self.rpos = rpos
            self.zbeg = zbeg
            self.zend = zend

        def check_intersection(self, position,momenta, mask):
            indices = np.where(mask)[0] #array of indices of particles that we will consider

            a = (momenta[indices,0])**2 + (momenta[indices,1])**2
            b = 2 * (position[indices,0]*momenta[indices,0] + position[indices,1]*momenta[indices,1])
            c = position[indices,1]**2 + position[indices,0]**2 - self.rpos**2
            
            coeffs = np.vstack((a,b,c)).T #( indices size), 3)
            roots=np.empty((a.shape[0], 2))
            for i,poly in enumerate(coeffs):
                 root = np.roots(poly)
                 if isinstance(root[0],complex):
                     root[0] = -1
                 if isinstance(root[1], complex):
                    root[1] = -1
                 roots[i,:] = root# (indices size, 2) THIS MIGHT NOT ALWAYS have size two

            ip_1 = position[indices,2] + roots[:,0]*momenta[indices,2] # only z
            ip_2 = position[indices,2] + roots[:,1]*momenta[indices,2] # only z

            new_mask_1 = (ip_1 > self.zbeg) & (ip_1 < self.zend) & (np.round(roots[:,0],decimals =12) > 0)
            new_mask_2 = (ip_2 > self.zbeg) & (ip_2 < self.zend) & (np.round(roots[:,1], decimals = 12) > 0)
            t = roots.copy() #need t to get the correct root
            any = new_mask_1 | new_mask_2

            doubles = new_mask_1 & new_mask_2
    
            t[new_mask_2,0] = roots[new_mask_2,1]

            if np.any(doubles):
                t[doubles,0] = np.min(roots[doubles], axis=1)
            
            ip = position[indices][any] + t[any,0][:, np.newaxis]*momenta[indices][any]
                    
            kept_indices = indices[np.where(any)[0]] # indexing the indices to get the correct particle numbers (ids)
            
            return kept_indices, ip
                    

    class conic:

        def __init__(self, parent, id, next_ids, tan_theta, zbeg, zend, rsmall, direction):
            self.parent = parent
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
            indices = np.where(mask)[0] #indices of particles that we will consider

            delta_z = self.zcenter - position[indices,2]
            a = -1 * momenta[indices,2]**2 * self.tan_theta**2 + momenta[indices,1]**2 + momenta[indices,0]**2
            b = 2*position[indices,0] * momenta[indices, 0] + 2*position[indices,1]*momenta[indices,1] + self.tan_theta**2 * 2*self.zcenter*momenta[indices,2] - self.tan_theta**2 * 2*position[indices,2]*momenta[indices,2]
            c = -1*self.tan_theta**2 * delta_z**2 + position[indices,0]**2 + position[indices,1]**2
            
            coeffs = np.vstack((a,b,c)).T # (sample_size, 3)
            roots=np.empty((a.shape[0], 2))
            for i,poly in enumerate(coeffs):
                 root = np.roots(poly)
                 if isinstance(root[0],complex):
                     root[0] = -1
                 if isinstance(root[1], complex):
                    root[1] = -1
                 roots[i,:] = root# (indices size, 2) THIS MIGHT NOT ALWAYS have size two

            

            ip_1 = position[indices] + roots[:,0].reshape((roots.shape[0],1)) * momenta[indices]
            ip_2 = position[indices] + roots[:,1].reshape((roots.shape[0],1)) * momenta[indices]

            #conditions: assert r between rsmall and (zcenter - zbeg)*tan_theta, z between zend and zbeg (first), root is positive
            new_mask_1 = (ip_1[:,2] > self.zbeg) & (ip_1[:,2] < self.zend) & (np.round(roots[:,0], decimals = 12) > 0)
            new_mask_2 = (ip_2[:,2] > self.zbeg) & (ip_2[:,2] < self.zend) & (np.round(roots[:,1], decimals = 12) > 0)
            
            t = roots.copy() #need t to get the correct root
            any = new_mask_1 | new_mask_2
            doubles = new_mask_1 & new_mask_2
            

            t[new_mask_2,0] = roots[new_mask_2,1]
            if np.any(doubles):
                t[doubles,0] = np.min(roots[doubles])
            
            ip = position[indices][any] + t[:,0].reshape((t.shape[0],1))[any] * momenta[indices][any]
            kept_indices = indices[np.where(any)[0]] # indexing the indices to get the correct particle numbers (ids)
        
            return kept_indices, ip

    class initializer:

        def __init__(self,parent, id, next_ids, rsmall, rbig):
            self.parent = parent
            self.id = id
            self.next_ids = next_ids
            self.rsmall = rsmall
            self.rbig = rbig
        
        def check_in(self, r,t, mask):
            indices = np.where(mask)[0]
            to = t.reshape((t.shape[0],))
            new_mask = (r[indices] < self.rbig) & (r[indices]> self.rsmall) & (to[indices] > 0)
            return indices[np.where(new_mask)[0]]


class cc:

    def __init__(self, sim):

        self.weights = sim.w
        self.Enumu = sim.Enumu
        self.Enue =sim.Enue
        old_m = sim.pnumu_ar[:,1:4] # after initial rotation within program, no need to translate
        old_m2 = sim.pnue_ar[:,1:4]
        po = sim.pos_at.T #need to translate
        old_p = np.copy(po)
        old_p[:,2] = po[:,2] - sim.Racc

        self.mnumu = np.empty(old_m.shape)
        self.mnue = np.empty(old_m2.shape)
        self.p = np.empty(old_p.shape)

        old = [old_m, old_m2, old_p]
        new = [self.mnumu, self.mnue, self.p]
        for i, coord in enumerate(new):
            coord[:,0] = -1 * old[i][:,0]
            coord[:,1] = old[i][:,2]
            coord[:,2] = old[i][:, 1]
    
    def completely_circular(self):
        return self

    def straight_segment_at_detector(self):
        pass

    def completely_linear(self):
        pass