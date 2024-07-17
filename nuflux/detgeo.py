import sys, os
cfp = os.path.dirname(os.path.abspath(__file__))
td = os.path.join(cfp,'detector_geometries')
sys.path.append(td)
cfp = os.path.dirname(os.path.abspath(__file__))
sys.path.append(cfp)

import copy, gc, importlib, time
import numpy as np
import matplotlib.pyplot as plt
from numba import njit, jit
from prettytable import PrettyTable
from memory_profiler import profile
from scipy.stats import uniform

import DarkNews as dn
import data
import helpers
import useful_data as ud

@profile
def check_mem():
    '''For memory-consuming-checking processes.'''
    z = 1+2
    return

def D3distance(point1, point2):
    '''3D Euclidian distance between two points.'''
    
    return np.sqrt((point1[:,0] - point2[:,0])**2 + (point1[:,1] - point2[:,1])**2 + (point1[:,2] - point2[:,2])**2)

def get_cs(E, part):
    '''Wrapper for detector_geometries.useful_data's cross section interpolation function.'''
    sigmanue,sigmanuebar,sigmanumu,sigmanumubar= ud.cs_interp()
    
    if part == "nue":
        
        return sigmanue(E)
    
    elif part == "nuebar":
        
        return sigmanuebar(E)
    
    elif part == "numu":
        
        return sigmanumu(E)
    
    elif part == "numubar":
        
        return sigmanumubar(E)

def SimulateDecays(param = 'mutristan_small', N_evals = 1e5, alr_loaded=False, dt = None):
    '''Wrapper for data's Monte Carlo generation of muon decays. Alr_loaded if dt pre-generated, in which case you need to specify dt.'''
    
    t0 = time.time()
    
    text = '(pre-generated dataset)'
    
    if not alr_loaded:
        dt = list(data.get_particles(param, N_evals))
        text = ''
        
    if dt:
        sim =  helpers.cc(C = dt[0], w =  dt[1], sample_size = dt[2], N_mu = dt[3], pnumu = dt[4], pnue = dt[5], pos = dt[6], name = dt[7], Emu = dt[8])
    
    else:
        raise ValueError('No dt provided.')
    
    t1 = time.time() - t0
    
    print(f'Simulation: {param} parameter set with {N_evals:.3e} evaluations {text}')
    print(f'{dt[2]:.3e} MC generations; took {t1:.3} s')
        
    return sim


muonic_neutrinos = ['numu', 'numubar']
electronic_neutrinos = ['nue', 'nuebar']
part_names = {'nue': 'ν_e', 'nuebar': 'anti ν_e', 'numu': 'ν_μ', 'numubar': 'anti ν_μ'}

class SimulateDetector():
    '''Detector Simulation.'''
    def __init__(self, coord_object, geom, particle = 'nue', Lss = 0):
        self.time = np.zeros(6)
        
        #Detector-related quantities
        self.Geometry = geom
        geom = importlib.import_module(self.Geometry)
        self.objects = geom.OBJECTS
        self.zbeginning = geom.zbeginning
        self.rmax = geom.rmax
        self.rbp = geom.rbp
        self.zending = geom.zending
        self.initials = geom.INITIALIZERS #should not contain minus_one,minus_two
        self.decayer = geom.DECAYER #is the minus_one,minus_two
        self.outside = geom.OUTSIDE # parent class of outside-going faces
        self.face_dict = geom.facedict
        self.outside_ids = np.array([obj.id for obj in self.outside])
        self.iterations = geom.iterations #size of the interactions_points array, one more than distances and densities
        self.detname = geom.name
        del geom 
        
        self.particle = particle
        
        #Decay-related quantities
        self.cc = copy.deepcopy(coord_object)
        self.cc.straight_segment_at_detector(self.zending, Lss = Lss)
        self.mutimes = self.cc.times
        self.w = self.cc.weights.reshape((self.cc.weights.shape[0],1)) # already normalized (does not include all decays!!!)
        self.Nmu = self.cc.Nmu
        self.sample_size = self.cc.sample_size
        self.paramname = self.cc.name
        
        if self.cc.L < 100:
            self.Lval = 0
            
        else:
            self.Lval = self.cc.L
        
        if self.particle in muonic_neutrinos:
            self.momenta = self.cc.pnumu[:,1:]
            self.E = self.cc.pnumu[:,0]
            
        elif self.particle in electronic_neutrinos:
            self.momenta = self.cc.pnue[:,1:]
            self.E = self.cc.pnue[:,0]

        #Detector-related quantities
        self.initialize_quantities()

    def find_info(self): 
        '''Getting all intersection points of neutrinos with detector components.'''
        
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

        mask_accepted = (~mask_not_accepted) & (to > 0)
        
        for obj in self.initials:
            mask_accepted = mask_accepted & (self.location[:,1]==-1)
            new_indices = obj.check_in(r_values, mask_accepted) #gives indices (numbers) of those that go in one of the initials

            self.location[new_indices, 1] = obj.id
            self.intersection_points[new_indices, 1, :] = ip[new_indices,:]
        
        del self.initials
        
        #at this point, there should be -1s (decaying in detector), 0s (did not reach detector), and ids of initials
        #treating decays in beampipe
        mask_decay = (self.location[:,1] == -1) & ~(np.sqrt(self.intersection_points[:,0,0]**2 + self.intersection_points[:,0,1]**2) > self.rbp)    
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
            self.densities[new_indices, 0] = self.decayer[0].density


        #treating decays within the detector components themselves
        weird_decay_mask = (self.location[:,1] == -1)
        self.update_intersections(self.decayer[1], count, weird_decay_mask)
        self.location[weird_decay_mask,0] = -2
        
        del self.decayer
        
        #there should be no more -1s, just 0 - x.
        #now iteratively finding intersection points
        self.time[2] = time.time()
        
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
        
        del self.objects
        del self.outside
        del self.outside_ids
    
    def get_probs(self):
        '''Determine the weights of each interacting neutrino.'''
        
        for i in range(self.iterations - 1):
            self.distances[:,i] = D3distance(self.intersection_points[:,i,:],self.intersection_points[:,i+1,:]) #these input arrays will be two dimensional, i.e., (sample_size, 3); distances will be i.e., (sample_size, 14
        self.dec_pos = self.intersection_points[:,0,:]
        
        #might have been easier to just save the t parameters throughout...
        del self.intersection_points
        
        self.part_line_integrals = self.distances * self.densities # huge array of distances in eahc component per ray
        for_mask = np.sum(self.part_line_integrals[:,1:], axis = 1)
        line_integrals = for_mask + self.part_line_integrals[:,0]
        self.cs = get_cs(self.E, self.particle)
        probs = 1 - np.exp(-1*self.cs * line_integrals)

        factors = np.full((self.sample_size,1), self.Nmu) * self.w
        
        del self.w #can be retrieved by np.sum(sim.part_face_counts, axis=1)
        
        self.counts = factors *probs[:, np.newaxis] # (sample_size, 1)
        self.mask = (self.counts>0)
        self.mask = self.mask[:, 0]
        self.counts = self.counts[self.mask]
        
        self.time[4] = time.time()
        
        del self.densities 
    
    #@profile
    def get_event_positions(self):
        '''Monte Carlo generation of event positions'''
        
        #freeing memory
        self.momenta = self.momenta[self.mask]
        self.dec_pos = self.dec_pos[self.mask]
        self.distances = self.distances[self.mask]
        self.part_line_integrals = self.part_line_integrals[self.mask]
        self.location = self.location[self.mask]
        self.E = self.E[self.mask]
        self.cs = self.cs[self.mask]
        self.mutimes = self.mutimes[self.mask, np.newaxis]
        
        self.times = self.distances / helpers.LIGHT_SPEED - self.mutimes
        
        temporary = np.empty((1,np.sum(self.mask), self.iterations - 1)) #same size as other big arrays now
        temporary[:,:, :] = [uniform.rvs(loc = 0, scale = self.distances)]        
        
        self.events_position, self.part_face_counts = get_probs_njit(self.dec_pos, self.momenta, self.distances, self.iterations, self.counts, self.cs, self.part_line_integrals, temporary)
        
        self.time[5] = time.time()
        
        return

    ##@profile
    def initialize_quantities(self):
        '''Initializes detector-related quantities.'''
        self.intersection_points = np.full((self.sample_size, self.iterations, 3), 1e4) #1e4 is arbitrary
        self.densities = np.zeros((self.sample_size, self.iterations - 1))
        self.densities[:,0] = ud.EARTH_DENSITY
        self.location = np.full((self.sample_size, self.iterations), -1) #starting at initials
        self.intersection_points[:,0,:] = self.cc.p
        self.distances = np.zeros((self.sample_size, self.iterations - 1))


    def update_params(self):
        '''In case one needs to run another simulation but cannot re-initialize.'''
        
        if self.particle in muonic_neutrinos:
            self.momenta = self.cc.pnumu[:,1:]
            self.E = self.cc.pnumu[:,0]

        elif self.particle in electronic_neutrinos:
            self.momenta = self.cc.pnue[:,1:]
            self.E = self.cc.pnue[:,0]

    
    def update_intersections(self, obj, count, mask):
        '''For a single object, with particles on it, finds their next component, iteratively on each neighbor, at a specific step of the sim.'''
        
        for neighbor in obj.next_ids:
            neigh = self.objects[neighbor]
            new_indices, ips = neigh.check_intersection(self.intersection_points[:,count - 1,:], self.momenta, mask)
            mask_1 = (ips[:,2] < self.intersection_points[new_indices, count, 2])& (ips[:,2] > self.intersection_points[new_indices, count - 1, 2])
            accepted_ix = new_indices[mask_1]
            
            self.intersection_points[accepted_ix, count, :] = ips[mask_1]
            self.location[accepted_ix, count] = neigh.id
            self.densities[accepted_ix, count-1] = obj.density

    def run(self, show_components = 1, show_time = 1, collision = 'mu+mu+'):
        '''Runs the whole simulation. show_components for distribution within different detector components. show_time for time summary. collision is the type of collision.'''        
        
        acc_colls_types = ['mu+e-', 'mu+mu+', 'mu+mu-']
        acc_colls_dict = {'mu+e-': 'μ+e-', 'mu+mu+': 'μ+μ+', 'mu+mu-': 'μ+μ-'}
        colls_types_to_part = {'mu+e-': ['numubar', 'nue'], 'mu+mu+': ['numubar', 'nue'], 'mu+mu-':['numubar', 'nue', 'numu', 'nuebar']}
        
        if (self.particle in muonic_neutrinos) | (self.particle in electronic_neutrinos):
            
            self.time[0] = time.time()
            
            self.find_info()
            self.get_probs()
            self.get_event_positions()
       
            return self
        
        elif collision in acc_colls_types:
            
            t0 = time.time()
            
            #have to run both (or four) sims.
            sim2 = None
            sim3 = None
            sim4 = None
            sims = [sim2, sim3, sim4]
            parts = colls_types_to_part[collision]
            nsims = len(parts)
            
            #first is actual instance of the simulation.
            self.particle = parts.pop(0)
            self.update_params()
            self.run()
            self.part_face_counts *= 2/nsims
            self.facecounts = calculate_facecounts(self.face_dict,self.location, self.part_face_counts)
            self.clear_mem(arg = 'first')
            self.collision = collision
            
            #others are distinct simulations.
            for i, part in enumerate(parts):
                sims[i] = SimulateDetector(self.cc, self.Geometry, part, Lss=-1).run()
                sims[i].part_face_counts *= 2/nsims
                sims[i].facecounts = calculate_facecounts(sims[i].face_dict, sims[i].location, sims[i].part_face_counts)
                sims[i].clear_mem()
            
            self.clear_mem(arg = 'last')
            
            if nsims == 2:
                sims = [self, sims[0]]
                total_count = np.sum(sims[0].w) + np.sum(sims[1].w)
        
            elif nsims == 4:
                sims = [self, sims[0], sims[1], sims[2]]
                total_count = np.sum(sims[0].w) + np.sum(sims[1].w) + np.sum(sims[2].w) + np.sum(sims[3].w)
            
            else:
                raise ValueError("Something went wrong with the number of simulatioms.")
            
            t1 = time.time() - t0
            
            print(f'Simulation: {self.paramname} ({acc_colls_dict[collision]}) at L = {self.Lval/100:.2f} m with {self.detname} as a detector')
            print(f'Total Count: {total_count:.2e} events; took {t1:.3} s')
            
            if show_components:
                get_face_counts(sims)
            
            if show_time:
                get_timetable(sims)
            
            return sims
        
        else:
            raise ValueError("No particle of that name included in this simulation!")
    
    def clear_mem(self, arg = 'another'):
        '''Freeing up memory.'''
        
        if arg == 'last':
            deletables = ['cc']
        
        else:
            mask = (self.location[:,0] == -1)
            self.part_face_counts[mask,0] = 0
            self.w  = self.part_face_counts.flatten()
            
            del self.part_face_counts
        
            new_mask = (self.w > 0)
            self.arrx = self.events_position[:,:,0].flatten()[new_mask]
            self.arry = self.events_position[:,:,1].flatten()[new_mask]
            self.arrz = self.events_position[:,:,2].flatten()[new_mask]
            
            del self.events_position
        
            self.times = self.times.flatten()[new_mask]
            self.w    = self.w[new_mask]
            
            del new_mask, mask
        
            deletables=['mutimes','momenta','location','distances','counts','part_line_integrals','f','zbeginning','rbp','iterations','Nmu','sample_size', 'cs', 'mask', 'dec_pos']
        
            if arg != 'first':
                deletables.append('cc')
        
        for att in deletables:
            
            if hasattr(self, att):
                delattr(self, att)

                
def plot(sims, nbins = 200, cmin = 1, orientation = 'z-y', give_data = False, savefig = None, fs = (20,12), cmap = 'viridis', ax = None, title = True, xl = True, yl = True):
    '''Plotting the detector event distribution as a hist2d instance, with the detector geometry behind.'''
    
    if not ax:
        fig, ax = plt.subplots(figsize = fs)
        
    bs = np.linspace(-1* sims[0].zending, sims[0].zending, nbins)
    bs2 = np.linspace(-1*sims[0].rmax, sims[0].rmax, nbins)
        
    if sims[0].collision == 'mu+mu-':
        x = np.concatenate((sims[0].arrx, sims[1].arrx, -1*sims[2].arrx, -1*sims[3].arrx))
        y = np.concatenate((sims[0].arry, sims[1].arry, sims[2].arry, sims[3].arry))
        z = np.concatenate((sims[0].arrz, sims[1].arrz, -1*sims[2].arrz, -1*sims[3].arrz))
        w = np.concatenate((sims[0].w, sims[1].w, sims[2].w, sims[3].w))
        
    elif sims[0].collision == 'mu+e-':
        x = np.concatenate((sims[0].arrx, sims[1].arrx))
        y = np.concatenate((sims[0].arry, sims[1].arry))
        z = np.concatenate((sims[0].arrz, sims[1].arrz))
        w = np.concatenate((sims[0].w, sims[1].w))
    
    elif sims[0].collision == 'mu+mu+':
        x = np.concatenate((sims[0].arrx, sims[1].arrx, -1*sims[0].arrx, -1*sims[1].arrx))
        y = np.concatenate((sims[0].arry, sims[1].arry, sims[0].arry, sims[1].arry))
        z = np.concatenate((sims[0].arrz, sims[1].arrz, -1*sims[0].arrz, -1*sims[1].arrz))
        w = np.concatenate((sims[0].w, sims[1].w, sims[0].w/2, sims[1].w/2))
    
    else:
        raise ValueError('This collision plotting has not been implemented yet!')
    
    lbl4 =
    lbl3 = "Collision: {}" 
    lbl = r"$L_{ss} = $" + f"{sims[0].Lval/100:.0f} m"
    lbl2 = r"$N_{events} = $" + f"{np.sum(w):.3e} events/yr"
    
    if orientation == 'z-y':
        ax.hist2d(z, y, alpha = 1, zorder = 30, bins = (bs, bs2), weights = w, cmin = cmin, cmap = cmap)
        plot_det(sims[0].Geometry, ax, xl  = xl, yl = yl)
        
    elif orientation == 'z-x':
        print(lbl)
        ax.hist2d(z, x, alpha = 1, zorder = 30, bins = (bs, bs2), weights = w, cmin = cmin, cmap = cmap)
        plot_det(sims[0].Geometry, ax, xl = xl, yl = yl)
            
    elif orientation == 'x-y':
        ax.hist2d(x, y, alpha = 1, zorder = 30, bins = (bs, bs2), weights = w, cmin = cmin, cmap = cmap)
        plot_det(sims[0].Geometry, ax, orientation = 'x-y', xl = xl, yl = yl)
        ax.set_xlim(-1* sims[0].rmax*10/12, sims[0].rmax *10/12)
        ax.set_ylim(0, sims[0].rmax)
        ax.legend([lbl], loc='lower right').set_zorder(50)

    else:
        raise ValueError('Only orientations are z-y, x-y, and z-x!')  
        
    ax.legend([lbl], loc='lower right').set_zorder(50)
    
    if title:
        ax.set_title(f'Event Distribution: {orientation}')
    
    if savefig:
        plt.savefig(savefig, bbox_inches = 'tight', dpi = 300)
            
    if give_data:
        return x,y,z,w

def event_timing(sims, fs = (20,12), histtype = 'barstacked', nbins = 100, give_data = False, savefig = None, label = '', legend = False):
    '''Wrapper to plot a hist of the neutrino interaction times.'''
    
    if fs:
        plt.figure(figsize  = fs)
    
    if len(sims) == 2:
        times = np.concatenate((sims[0].times, sims[1].times))
        w = np.concatenate((sims[0].w, sims[1].w))
    
    elif len(sims) == 4:
        times = np.concatenate((sims[0].times, sims[1].times, sims[2].times, sims[3].times))
        w = np.concatenate((sims[0].w, sims[1].w, sims[2].w, sims[3].w))
        
    else:
        raise ValueError('Wrong number of sims. Something went wrong.')
    
    plt.xlabel('Time before collision (s)')
    plt.ylabel(r'$N_{events}$')
    plt.title('Event Timing (with respect to collision time)')
    plt.hist(times, weights = w, histtype = histtype, bins = nbins, label = label)
    
    if legend:
        plt.legend(loc='best')
    
    if savefig:
        plt.savefig(savefig, bbox_inches = 'tight', dpi = 300)
            
    if give_data:
        return times, w
    
def phi_distribution(sims, fs = (20,12), histtype = 'step', nbins = 100, give_data = False, savefig = None, ylog = True, label='', legend = False):
    '''Wrapper to plot the phi distribution of neutrino events.'''
    
    if fs:
        plt.figure(figsize = fs)
    
    if len(sims) == 2:
        x = np.concatenate((sims[0].arrx, sims[1].arrx))
        y = np.concatenate((sims[0].arry, sims[1].arry))
        w = np.concatenate((sims[0].w, sims[1].w))
    
    elif len(sims) == 4:
        x = np.concatenate((sims[0].arrx, sims[1].arrx, -1*sims[2].arrx, -1*sims[3].arrx))
        y = np.concatenate((sims[0].arry, sims[1].arry, sims[2].arry, sims[3].arry))
        w = np.concatenate((sims[0].w, sims[1].w, sims[2].w, sims[3].w))
    
    phi = np.arctan(x/y)
    
    plt.hist(phi, weights = w, histtype = histtype, bins = nbins, label = label)
    plt.ylabel(r'$N_{events}/yr$')
    plt.xlabel(r'$\phi$ Distribution (rad)')
    
    if legend:
        plt.legend(loc='best')
        
    if ylog:
        plt.yscale('log')
    
    if savefig:
        plt.savefig(savefig, bbox_inches = 'tight', dpi = 300)
            
    if give_data:
        return phi, w

def get_timetable(sims):
    '''Prints the table of detailed time durations for the simulation.'''
    
    if len(sims) == 2:
        cols = ['ν_e time', 'anti ν_μ time']
        data = [[sims[1].time[i+1] - sims[1].time[i], sims[0].time[i+1] - sims[0].time[i], sims[0].time[i+1] - sims[0].time[i] + sims[1].time[i+1] - sims[1].time[i]] for i in range(0, len(sims[0].time)-1)]
        data.append([sims[1].time[-1] - sims[1].time[0], sims[0].time[-1] - sims[0].time[0], sims[1].time[-1] - sims[0].time[0]])

    elif len(sims) == 4:
        cols = ['ν_e time', 'anti ν_e time', 'ν_μ time', 'anti ν_μ time']
        data = [[sims[1].time[i+1] - sims[1].time[i], sims[3].time[i+1] - sims[3].time[i], sims[2].time[i+1] - sims[2].time[i], sims[0].time[i+1] - sims[0].time[i], sims[0].time[i+1] - sims[0].time[i] + sims[1].time[i+1] - sims[1].time[i] + sims[2].time[i+1] - sims[2].time[i] + sims[3].time[i+1] - sims[3].time[i]] for i in range(0, len(sims[0].time)-1)]
        data.append([sims[1].time[-1] - sims[1].time[0], sims[3].time[-1] - sims[3].time[0],sims[2].time[-1] - sims[2].time[0], sims[0].time[-1] - sims[0].time[0], sims[3].time[-1] - sims[0].time[0]])
    
    table = PrettyTable()
    names = ['init','init obj','other obj','probs','events', 'totals']
    labels = ['Simulation Time']
    labels.extend(cols)
    labels.append('TOTAL')
    table.field_names = labels
    
    for name, row in zip(names, data):
        formatted_row = [f'{x:.3e}' for i,x in enumerate(row)]
        table.add_row([name] + formatted_row)
        
    print(table)

def get_face_counts(sims):
    '''Prints the table of detailed distribution of events in detector components.'''
    names = list(sims[0].face_dict.keys())
    names.append('TOTAL')
    total_count = 0
    totals = []
    
    for sim in sims:
        totals.append(sum(list(sim.facecounts.values())))
        
    total_count = sum(totals)
    
    if len(sims) == 4:
        data = [[sims[1].facecounts[key], sims[3].facecounts[key], sims[2].facecounts[key], sims[0].facecounts[key], sims[0].facecounts[key] + sims[1].facecounts[key] + sims[2].facecounts[key] + sims[3].facecounts[key]] for key in sims[0].face_dict.keys()]
        data.append([totals[1], totals[3], totals[2], totals[0], total_count])
        parts = ['ν_e events','anti ν_e events', 'ν_μ events','anti ν_μ events']
        
    elif len(sims) == 2:
        data = [[sims[1].facecounts[key], sims[0].facecounts[key], sims[0].facecounts[key] + sims[1].facecounts[key]] for key in sims[0].face_dict.keys()]
        data.append([totals[1], totals[0], total_count])
        parts = ['ν_e events', 'ν_μ events']
        
    else:
        raise ValueError("Number of simulations not accepted; something went wrong.")
        
    table = PrettyTable()
    labels = ['Detector Parts']
    labels.extend(parts)
    labels.append('Total Events')
    table.field_names = labels
        
    for name, row in zip(names, data):
        formatted_row = [f"{x:.3e}" for x in row]
        table.add_row([name] + formatted_row)

    print(table)
    

@jit(nopython = False, forceobj=True)
def calculate_facecounts(face_dict, location, part_face_counts):
    '''Gets the facecounts for each detector component in the simulation.'''
    facecounts = {}
    locs = location[:, :-1]
    
    for key, faces in face_dict.items():
        mask = np.isin(locs, faces)
        facecounts[key] = np.sum(part_face_counts[mask])

    return facecounts


def get_probs_njit(dec_pos, momenta, distances, iterations,counts, cs, part_line_integrals, temporary):
    '''Gets the number of events (weights) of each interacting neutrino.'''
    #pweights these are probabilities, NOT WEIGHTS
    temporary = temporary.reshape((temporary.shape[1], temporary.shape[2]))
    temp2 = np.ones((distances.shape[0], iterations - 1)) #probs of non-interaction in each segment
    temp2 = np.exp(-1 *cs[:, np.newaxis] * part_line_integrals)
    mi = np.prod(temp2, axis = 1)
    
    pweights = get_events_njit1(temp2, mi[:, np.newaxis])
    
    max_index = np.argmax(pweights[:], axis=1)
    mask = np.ones_like(pweights, dtype=bool)[:,0]
    max_p = pweights[mask, max_index]
    max_p = max_p[:, np.newaxis]

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
    '''Gets weight of each event along the segment of a single ray.'''
    res = (1 - temp2) * mi / temp2
    
    return res


@njit
def get_events_njit2(shape, pweights, mid, counts, dec_pos, normed_m, cumulative_distances, temporary):
    '''Gets weight of each event generated by all interacting neutrinos.'''
    # shape is (:, 25); pweights is (:, 25); mid is (:,1); mask is (:,); counts is (:,1); dec_pos is (:, 1, 3); normed_m is [:, 1, 3]; cumulative_distances is [:,25]; temporary is [:, 25]
    pweights = pweights / mid #good multinomial probs
    temporary[:,1:] += cumulative_distances[:,:-1]
    part_face_counts = counts * pweights # these are the face counts (expected)
    events_position = dec_pos + temporary[:,:,np.newaxis] * normed_m
    
    return events_position, part_face_counts


def plot_det(geom, ax, orientation ='z-y', xl = True, yl = True):
    '''Plots the detector geometry behind a sim plot.'''

    if geom == 'det_v1':
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
        
        for i, det in enumerate(dets):
            ax.plot(det[0],det[1], color = cols[i])
            ax.fill_between(det[0], det[1], color = cols[i], alpha=0.7)
            new_y =[-1*k for k in det[1]]
            ax.plot(det[0], new_y, color = cols[i])
            ax.fill_between(det[0], new_y, color = cols[i], alpha=0.7)

        ax.set_xlabel("z-coordinate (cm)")
        ax.set_ylabel("r-coordinate (cm)")

    
    elif (geom == 'det_v2') | (geom == 'zero_density_test'):
        
        if orientation=='x-y':
            MDET = 645
            SPS1 = 446.1
            SOL_1= 429
            SPS2 = 425
            SOL_2= 399.3
            SPS3 = 364.9
            SOL_3= 352.3
            SPS4 = 348.3
            HCAL = 333
            SPS5 = 174
            ECAL = 170.2
            SPS6 = 150
            CONE = 31
            BL   = 2.2
            dets = [MDET, SPS1, SOL_1, SPS2, SOL_2, SPS3, SOL_3, SPS4, HCAL, SPS5, ECAL, SPS6, CONE, BL]
            cols = ['grey', 'white', 'gray', 'white', 'lightgrey', 'white','gray','white','darkgrey','white','dimgrey','white', 'black', 'white']
            
            for i, det in enumerate(dets):
                circle = plt.Circle((0,0), det, zorder = i, alpha = 0.7, edgecolor = cols[i], facecolor=cols[i])
                ax.add_artist(circle)
                
            ax.scatter(MDET,MDET, color = cols[3])
            
            if xl:
                ax.set_xlabel('x-coordinate (cm)')
            
            if yl:
                ax.set_ylabel('y-coordinate (cm)')
                
        else:
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
            cols = ['dimgrey']*3 + ['darkgrey']*3 + ['gray'] + ['lightgrey'] + ['gray'] + ['grey'] + ['black']*2 + ['white']
            
            for i, det in enumerate(dets):
                ax.plot(det[0],det[1], color = cols[i])
                ax.fill_between(det[0], det[1], color = cols[i], alpha=0.7)
                new_y =[-1*k for k in det[1]]
                ax.plot(det[0], new_y, color = cols[i])
                ax.fill_between(det[0], new_y, color = cols[i], alpha=0.7)

            if xl:
                ax.set_xlabel("z-coordinate (cm)")
            
            if yl:
                ax.set_ylabel("y-coordinate (cm)")

    else:
        print("this geometry has not been implemented yet!")
