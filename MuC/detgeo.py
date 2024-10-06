import os

import copy
import importlib
import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from numba import njit
from prettytable import PrettyTable
from scipy.stats import uniform
from matplotlib.colors import LogNorm

import MuC
from MuC import mudecay_tools as mudec
from MuC import detector_tools as det
from MuC import collider_tools as col
from MuC import xsecs

# from memory_profiler import profile

# @profile
# def check_mem():
#     """For memory-consuming-checking processes."""
#     z = 1 + 2
#     return


def D3distance(point1, point2):
    """3D Euclidian distance between two points."""

    return np.sqrt(
        (point1[:, 0] - point2[:, 0]) ** 2
        + (point1[:, 1] - point2[:, 1]) ** 2
        + (point1[:, 2] - point2[:, 2]) ** 2
    )


def get_total_xsec(E, part):
    """Wrapper for detector_geometries.useful_data's cross section interpolation function."""
    if part == "nue":

        return xsecs.total_sigmanue(E)

    elif part == "nuebar":

        return xsecs.total_sigmanuebar(E)

    elif part == "numu":

        return xsecs.total_sigmanumu(E)

    elif part == "numubar":

        return xsecs.total_sigmanumubar(E)


def SimulateDecays(design, N_evals=1e5, alr_loaded=False, dt=None):
    """Wrapper for data's Monte Carlo generation of muon decays. Alr_loaded if dt pre-generated, in which case you need to specify dt."""

    t0 = time.time()

    text = "(pre-generated dataset)"
    if not alr_loaded:
        attempt = 0
        max_attempts = 10

        while attempt < max_attempts:

            try:
                dt = list(mudec.get_particles(design, N_evals))
                max_attempts = 0

            except Exception as e:
                attempt += 1
                print(
                    f"Failed event generation. Attempt number {attempt}. Trying again. Exception: {e}",
                )

        if not dt:
            raise RuntimeError(
                "Max attempts reached. Could not generate dataset for decays."
            )

        text = ""

    if dt:
        sim = det.cc(
            C=dt[0],
            w=dt[1],
            sample_size=dt[2],
            N_mu=dt[3],
            pnumu=dt[4],
            pnue=dt[5],
            pos=dt[6],
            name=dt[7],
            Emu=dt[8],
        )

    else:
        raise ValueError("No dt provided.")

    t1 = time.time() - t0

    print(f"{dt[7]} parameter set with {N_evals:.3e} evaluations {text}.")
    print(f"{dt[2]:.3e} MC generations; took {t1:.3} s.")

    return sim


class SimNeutrinos:
    """Detector Simulation."""

    # @profile
    def __init__(self, coord_object, geom, particle=None, direc="left"):
        self.time = np.zeros(6)

        self.d = direc
        # Detector-related quantities
        self.objects = geom.OBJECTS
        self.zbeginning = geom.zbeginning
        self.rmax = geom.rmax
        self.rbp = geom.rbp
        self.zending = geom.zending
        self.initials = geom.INITIALIZERS  # should not contain minus_one,minus_two
        self.decayer = geom.DECAYER  # is the minus_one,minus_two
        self.outside = geom.OUTSIDE  # parent class of outside-going faces
        self.face_dict = geom.facedict
        self.outside_ids = np.array([obj.id for obj in self.outside])
        self.iterations = (
            geom.iterations
        )  # size of the interactions_points array, one more than distances and densities
        self.detname = geom.name
        del geom

        self.particle = particle

        # Decay-related quantities
        self.mutimes = coord_object.times
        self.sample_size = coord_object.sample_size
        self.Nmu = coord_object.Nmu
        self.weights = coord_object.weights.reshape(
            (coord_object.weights.shape[0], 1)
        )  # already normalized (does not include all decays!!!)
        self.weights = np.full((self.sample_size, 1), self.Nmu) * self.weights
        self.paramname = coord_object.name

        if coord_object.L < 100:
            self.Lval = 0

        else:
            self.Lval = coord_object.L

        if self.particle in MuC.anti_neutrinos:
            self.momenta = coord_object.pnumu[:, 1:]
            self.E = coord_object.pnumu[:, 0]

        elif self.particle in MuC.neutrinos:
            self.momenta = coord_object.pnue[:, 1:]
            self.E = coord_object.pnue[:, 0]

        # Detector-related quantities
        self.initialize_quantities(coord_object.p)

    # @profile
    def find_info(self):
        """Getting all intersection points of neutrinos with detector components."""

        # first interaction points
        count = 1
        delta_z = self.zbeginning - self.intersection_points[:, 0, 2]
        to = delta_z / self.momenta[:, 2]
        ip = (
            self.intersection_points[:, 0, :] + to[:, np.newaxis] * self.momenta
        )  # intersection points
        r_values = np.sqrt(ip[:, 0] ** 2 + ip[:, 1] ** 2)
        cond_1 = ~(
            (
                (self.intersection_points[:, 0, 2] < self.zending)
                & (self.intersection_points[:, 0, 2] > self.zbeginning)
                & (self.intersection_points[:, 0, 1] > -1 * self.rmax)
                & (self.intersection_points[:, 0, 1] < self.rmax)
            )
        )
        cond_2 = (r_values > self.rmax) | (to < 0)
        mask_not_accepted = cond_2 & cond_1
        not_accepted_indices = np.where(mask_not_accepted)[
            0
        ]  # indices of this new array

        self.location[not_accepted_indices, 1] = 0
        self.intersection_points[not_accepted_indices, 1, :] = self.intersection_points[
            not_accepted_indices, 0, :
        ]

        self.time[1] = time.time()

        mask_accepted = (~mask_not_accepted) & (to > 0)

        for obj in self.initials:
            mask_accepted = mask_accepted & (self.location[:, 1] == -1)
            new_indices = obj.check_in(
                r_values, mask_accepted
            )  # gives indices (numbers) of those that go in one of the initials
            self.location[new_indices, 1] = obj.id
            self.intersection_points[new_indices, 1, :] = ip[new_indices, :]
            self.new_indices = new_indices

        del self.initials

        # at this point, there should be -1s (decaying in detector), 0s (did not reach detector), and ids of initials
        # treating decays in beampipe
        mask_decay = (self.location[:, 1] == -1) & ~(
            np.sqrt(
                self.intersection_points[:, 0, 0] ** 2
                + self.intersection_points[:, 0, 1] ** 2
            )
            > self.rbp
        )
        a, b, c = det.barrel_get_pols(
            self.rbp,
            self.intersection_points[mask_decay, 0, :],
            self.momenta[mask_decay],
        )
        roots = det.get_roots(a, b, c)

        new_mask_1 = np.round(roots[:, 0], decimals=12) > 0
        new_mask_2 = np.round(roots[:, 1], decimals=12) > 0
        t = roots.copy()  # need t to get the correct root
        doubles = new_mask_1 & new_mask_2
        t[new_mask_2, 0] = roots[new_mask_2, 1]

        if np.any(doubles):
            t[doubles, 0] = np.min(roots[doubles])

        ip2 = (
            self.intersection_points[mask_decay, 0, :]
            + t[:, 0][:, np.newaxis] * self.momenta[mask_decay]
        )

        if self.decayer:

            for neighb in self.decayer[0].next_ids:
                neighbor = self.objects[neighb]
                loc_mask = self.location[mask_decay, 1] == -1
                new_indices, new_mask = self.decayer[0].check_in(
                    neighbor, ip2[:, 2], mask_decay, loc_mask
                )  # gives indices (numbers)
                self.location[new_indices, 1] = neighb
                ips = ip2[new_mask, :]
                self.intersection_points[new_indices, 1, :] = ips
                self.densities[new_indices, 0] = self.decayer[0].density

            # treating decays within the detector components themselves
            weird_decay_mask = self.location[:, 1] == -1
            self.update_intersections(self.decayer[1], count, weird_decay_mask)
            self.location[weird_decay_mask, 0] = -2

        del self.decayer

        # there should be no more -1s, just 0 - x.
        # now iteratively finding intersection points
        self.time[2] = time.time()

        count = 2
        while count < self.iterations:

            locs = np.unique(self.location[:, count - 1])
            objects = [self.objects[i] for i in locs]

            for obj in objects:
                part_mask = self.location[:, count - 1] == obj.id
                particles = np.where(part_mask)[
                    0
                ]  # these are indices of the particles at that face

                if obj.id in self.outside_ids:
                    self.intersection_points[particles, count, :] = (
                        self.intersection_points[particles, count - 1, :]
                    )
                    self.location[particles, count] = self.location[
                        particles, count - 1
                    ]
                    continue

                if np.any(particles):
                    self.update_intersections(obj, count, part_mask)

            count += 1

        self.time[3] = time.time()

        del self.objects
        del self.outside
        del self.outside_ids

    # @profile
    def get_probs(self):
        """Determine the weights of each interacting neutrino."""

        for i in range(self.iterations - 1):
            self.distances[:, i] = D3distance(
                self.intersection_points[:, i, :], self.intersection_points[:, i + 1, :]
            )  # these input arrays will be two dimensional, i.e., (sample_size, 3); distances will be i.e., (sample_size, 14
        self.dec_pos = self.intersection_points[:, 0, :]

        # might have been easier to just save the t parameters throughout...
        # del self.intersection_points

        self.part_line_integrals = (
            self.distances * self.densities
        )  # huge array of distances in eahc component per ray
        for_mask = np.sum(self.part_line_integrals[:, 1:], axis=1)
        line_integrals = for_mask + self.part_line_integrals[:, 0]

        # Total event rate (total nu cross section)
        self.cs = get_total_xsec(self.E, self.particle)
        probs = 1 - np.exp(-1 * self.cs * line_integrals)

        self.counts = self.weights * probs[:, np.newaxis]  # (sample_size, 1)
        self.mask = self.counts > 0
        self.mask = self.mask[:, 0]
        self.counts = self.counts[self.mask]

        self.time[4] = time.time()

        del self.densities

    # @profile
    def get_probs_njit(self):
        """Gets the number of events (weights) of each interacting neutrino."""
        # pweights these are probabilities, NOT WEIGHTS
        self.t_values = self.t_values.reshape(
            (self.t_values.shape[1], self.t_values.shape[2])
        )
        temp2 = np.ones(
            (self.distances.shape[0], self.iterations - 1)
        )  # probs of non-interaction in each segment
        temp2 = np.exp(-1 * self.cs[:, np.newaxis] * self.part_line_integrals)
        mi = np.prod(temp2, axis=1)

        pweights = get_events_njit1(temp2, mi[:, np.newaxis])

        max_index = np.argmax(pweights[:], axis=1)
        mask = np.ones_like(pweights, dtype=bool)[:, 0]
        max_p = pweights[mask, max_index]
        max_p = max_p[:, np.newaxis]

        pweights = pweights / max_p
        mid = np.sum(pweights, axis=1)
        cumulative_distances = np.cumsum(self.distances, axis=1)
        normed_m = self.momenta / np.linalg.norm(self.momenta, axis=1)[:, np.newaxis]
        self.dec_pos = self.dec_pos[:, np.newaxis, :]
        self.costheta = normed_m[:, 2]
        normed_m = normed_m[:, np.newaxis, :]

        self.events_position, self.part_face_counts = get_events_njit2(
            pweights,
            mid[:, np.newaxis],
            self.counts,
            self.dec_pos,
            normed_m,
            cumulative_distances,
            self.t_values,
        )

    # @profile
    def get_event_positions(self):
        """Monte Carlo generation of event positions"""
        # freeing memory
        self.momenta = self.momenta[self.mask]
        self.dec_pos = self.dec_pos[self.mask]
        self.distances = self.distances[self.mask]
        self.part_line_integrals = self.part_line_integrals[self.mask]
        self.intersection_points = self.intersection_points[self.mask]
        self.location = self.location[self.mask]
        self.E = self.E[self.mask]
        self.cs = self.cs[self.mask]
        self.mutimes = self.mutimes[self.mask, np.newaxis]
        self.weights = self.weights[self.mask]

        self.t_values = np.empty(
            (1, np.sum(self.mask), self.iterations - 1)
        )  # same size as other big arrays now
        self.t_values[:, :, :] = [uniform.rvs(loc=0, scale=self.distances)]

        self.get_probs_njit()

        self.time[5] = time.time()

        return

    # @profile
    def initialize_quantities(self, p):
        """Initializes detector-related quantities."""
        self.intersection_points = np.full(
            (self.sample_size, self.iterations, 3), 1e4
        )  # 1e4 is arbitrary
        self.densities = np.zeros((self.sample_size, self.iterations - 1))
        self.densities[:, 0] = det.EARTH_DENSITY
        self.location = np.full(
            (self.sample_size, self.iterations), -1
        )  # starting at initials
        self.intersection_points[:, 0, :] = p
        self.distances = np.zeros((self.sample_size, self.iterations - 1))

    # @profile
    def update_intersections(self, obj, count, mask):
        """For a single object, with particles on it, finds their next component, iteratively on each neighbor, at a specific step of the sim."""

        for neighbor in obj.next_ids:
            neigh = self.objects[neighbor]
            new_indices, ips = neigh.check_intersection(
                self.intersection_points[:, count - 1, :], self.momenta, mask
            )
            mask_1 = (ips[:, 2] < self.intersection_points[new_indices, count, 2]) & (
                ips[:, 2] > self.intersection_points[new_indices, count - 1, 2]
            )
            accepted_ix = new_indices[mask_1]

            self.intersection_points[accepted_ix, count, :] = ips[mask_1]
            self.location[accepted_ix, count] = neigh.id
            self.densities[accepted_ix, count - 1] = obj.density

    # @profile
    def calculate_facecounts(self):
        """Gets the facecounts for each detector component in the simulation."""
        self.facecounts = {}

        for key in self.face_dict.keys():
            face_mask = self.get_face_masks(key)
            self.facecounts[key] = np.sum(self.w[face_mask])

    def get_face_masks(self, sec):
        """Generate a masks for a specific detector component."""

        if sec == "all":
            mask = np.ones_like(self.w, dtype=bool)

        else:
            mask = np.isin(self.location, self.face_dict[sec])

        return mask

    # @profile
    def get_lum_q(self, param):
        """Getting the luminosity parameters for a beam collision."""
        self.intersection_points = self.intersection_points[self.new_indices, 1, :]
        self.weights = self.weights[self.new_indices]
        self.Nnu = np.sum(self.weights) / col.parameters[param]["syr"]
        wax = np.average(self.intersection_points[:, 0], weights=self.weights[:, 0])
        way = np.average(self.intersection_points[:, 1], weights=self.weights[:, 0])
        self.nusdx = np.sqrt(
            np.average(
                (self.intersection_points[:, 0] - wax) ** 2, weights=self.weights[:, 0]
            )
        )
        print(f"{self.nusdx*10:.3} mm for sigma x")
        self.nusdy = np.sqrt(
            np.average(
                (self.intersection_points[:, 1] - way) ** 2, weights=self.weights[:, 0]
            )
        )
        print(f"{self.nusdy*10:.3} mm for sigma y")

    def run(self):
        """Runs the simulation for a single neutrino species."""

        if (self.particle in MuC.anti_neutrinos) | (self.particle in MuC.neutrinos):

            self.time[0] = time.time()

            self.find_info()
            self.get_probs()
            self.get_event_positions()
            # print('ran it')
            return self

        else:
            raise ValueError("No particle of that name included in this simulation!")

    # @profile
    def clear_mem(self, arg=None):
        """Freeing up memory."""

        if arg:

            if hasattr(self, arg):
                delattr(self, arg)
                return

        # print('good')
        mask = self.location[:, 0] == -1
        self.part_face_counts[mask, 0] = 0
        self.w = self.part_face_counts.flatten()
        self.times = self.t_values / det.LIGHT_SPEED - self.mutimes
        self.t_values = self.t_values.flatten()
        # print('good')
        # print(self.E.shape, self.weights.shape, self.part_face_counts.shape)
        self.E = self.E[:, np.newaxis] * np.ones(self.part_face_counts.shape)
        self.weights = self.weights * np.ones(self.part_face_counts.shape)
        self.costheta = self.costheta[:, np.newaxis] * np.ones(
            self.part_face_counts.shape
        )
        # print('good')
        del self.part_face_counts
        # print('good')
        new_mask = self.w > 0
        # print('mask ok')
        # check_mem()
        # print(self.E.shape, self.weights.shape)
        self.E = self.E.flatten()[new_mask]
        self.weights = self.weights.flatten()[new_mask]
        # print('good')
        self.costheta = self.costheta.flatten()[new_mask]
        self.arrx = self.events_position[:, :, 0].flatten()[new_mask]
        self.arry = self.events_position[:, :, 1].flatten()[new_mask]
        self.arrz = self.events_position[:, :, 2].flatten()[new_mask]
        self.t_values = self.t_values[new_mask]
        self.distances = self.distances.flatten()[new_mask]
        # print('pli?')
        nl = self.part_line_integrals.flatten()[new_mask]
        self.wnl = nl * self.weights
        # print('w r fine')
        ones = np.ones_like(self.times)
        self.mutimes = ones * self.mutimes
        self.mutimes = self.mutimes.flatten()[new_mask]
        self.location = self.location[:, :-1].flatten()[new_mask]

        del self.events_position

        self.times = self.times.flatten()[new_mask]
        self.w = self.w[new_mask]

        del new_mask, mask

        if self.d == "right":
            # print('got one right')
            self.arrx = -1 * self.arrx
            self.arrz = -1 * self.arrz

        deletables = [
            "momenta",
            "counts",
            "part_line_integrals",
            "f",
            "zbeginning",
            "rbp",
            "iterations",
            "Nmu",
            "sample_size",
            "cs",
            "mask",
            "mutimes",
            "distances",
            "t_values",
            "weights",
            "dec_pos",
        ]

        for att in deletables:

            if hasattr(self, att):
                delattr(self, att)


class SimulateDetector:
    """Simulating the whole detector interactions.
    Hierarchy is as follows: Initializing SimulateDetector calls SimulateDecays; running simulates many SimNeutrinos based on the collision type,
      which is a single-neutrino-species MC generation of events within a detector, which are all saved in a list in the .sims attribute.
    """

    def __init__(
        self,
        design,
        N_evals=1e5,
        alr_loaded=False,
        dt=None,
        geom="det_v2",
        save_mem=True,
    ):
        """Initializes from muon decay sim."""

        self.save_mem = save_mem
        self.design = design
        self.N_evals = N_evals
        self.cco = SimulateDecays(
            design=design, N_evals=N_evals, alr_loaded=alr_loaded, dt=dt
        )

        # Detector geometry
        self.geom = getattr(MuC, geom)
        self.comps = list(self.geom.facedict.keys())
        self.zending = self.geom.zending
        self.rmax = self.geom.rmax
        cc1, cc2 = self.cco.straight_segment_at_detector(
            self.geom.zending, Lss=self.design["Lss"], two=True
        )

        self.ccs = [cc1, cc1, cc2, cc2]
        if self.save_mem:
            del cc1
            del cc2

    # @profile
    def run(
        self,
        show_components=0,
        show_time=0,
    ):
        """Runs the whole simulation, based on a storage ring geometry, detector geometry, and collision.

        Args:
            show_components (bool): for distribution within different detector components.
            show_time (bool): for time summary.
            geom (str): the detector version. Latest is det_v2; uniform block is block_test; zero_density_test is exactly that; and det_v1 is simpler.
            Lss (float): Length of the straight segment upon which the IR detector is centered.
        """

        t0 = time.time()

        # Attempts to perform a muon decay simulation
        self.ntimes = 6  # should be constant, except during debugging/improvements.

        # Neutrino flux types -- determines how many muon decay simulations we need
        self.parts = [
            part[0] for part in col.colls_types_to_part[self.design["collision_type"]]
        ]
        self.nsims = len(self.parts)

        # Cotainer for all simulations
        self.sims = [[None]] * self.nsims

        for i, part in enumerate(self.parts):

            print(f"Simulating muon decays for {part}")
            self.sims[i] = SimNeutrinos(
                self.ccs[i], self.geom, part, MuC.directions[i]
            ).run()

            self.sims[i].clear_mem()
            # NOTE: Where is this factor of 2 coming from? Two beams? I dont get it.
            self.sims[i].w *= 2 / self.nsims
            self.sims[i].calculate_facecounts()

        self.total_count = np.sum([np.sum(self.sims[i].w) for i in range(self.nsims)])
        t1 = time.time() - t0
        extra = ""

        if t1 > 1.75 * self.nsims / 2 * self.N_evals / 1e5:
            extra = " (numba pre-compilation needed)"

        self.name = self.design["name"]
        print(f"Successfully simulated neutrino event rates within {self.geom.name}:")
        print(
            f"{self.name} ({col.acc_colls_dict[self.design['collision_type']]}) at L = {self.design['Lss']:.2f} m."
        )
        print(f"Total count: {self.total_count:.2e} events; took {t1:.3} s{extra}.\n")

        if show_components:
            self.get_face_counts()

        if show_time:
            self.get_timetable()

        sim = copy.deepcopy(self)
        sim.sims = self.sims

        if self.save_mem:
            del self.sims
            del sim.cco
            del sim.geom

        return sim

    def NuNuLuminosity(self, particle1="nue", particle2="numu"):
        """Computing the luminosity of a neutrino collision."""

        cc1, cc2 = self.cco.straight_segment_at_detector(
            0, Lss=self.design["Lss"], two=False
        )
        self.collision = "mu+mu-"
        self.ntimes = 6  # should be constant, except during debugging/improvements.
        sim1 = None
        sim2 = None
        sim3 = None
        sim4 = None
        sims = [sim1, sim2, sim3, sim4]
        self.parts = [part[0] for part in col.colls_types_to_part[self.collision]]
        nsims = len(self.parts)
        geom = importlib.import_module("MuC.detector_geometries.nunulum")
        self.zending = geom.zending
        self.rmax = geom.rmax
        mask = np.array(self.parts) == "lol"
        mask[np.array(self.parts) == particle1] = True
        mask[np.array(self.parts) == particle2] = True
        indices = np.where(mask)[0]
        cc1, cc2 = self.cco.straight_segment_at_detector(
            geom.zending, Lss=self.design["Lss"], two=True
        )
        ccs = [cc1, cc1, cc2, cc2]
        for i, part in enumerate(self.parts):

            if i not in indices:
                continue

            sims[i] = SimNeutrinos(ccs[i], geom, part, MuC.directions[i])
            sims[i].find_info()
            sims[i].weights *= 2 / nsims
            sims[i].get_lum_q(self.design)

        sims = [sims[i] for i in range(nsims)]
        p1 = indices[0]
        p2 = indices[1]

        return sims, sims[p1].Nnu * sims[p2].Nnu / 4 / np.pi / np.sqrt(
            sims[p1].nusdx ** 2 + sims[p2].nusdx ** 2
        ) / np.sqrt(sims[p1].nusdy ** 2 + sims[p2].nusdy ** 2)

    def get_timetable(self):
        """Prints the table of detailed time durations for the simulation."""

        cols = [
            MuC.part_names[part[0]] + " time (" + part[1] + ")"
            for part in col.colls_types_to_part[self.collision]
        ]
        data = []

        for i in range(0, self.ntimes - 1):
            row = [sim.time[i + 1] - sim.time[i] for sim in self.sims]
            row.append(sum(row))
            data.append(row)

        trow = [sim.time[-1] - sim.time[0] for sim in self.sims]
        trow.append(self.sims[-1].time[-1] - self.sims[0].time[0])
        data.append(trow)

        table = PrettyTable()
        names = ["init", "init obj", "other obj", "probs", "events", "totals"]
        labels = ["Simulation Time"]
        labels.extend(cols)
        labels.append("TOTAL")
        table.field_names = labels

        for name, row in zip(names, data):
            formatted_row = [f"{x:.3e}" for i, x in enumerate(row)]
            table.add_row([name] + formatted_row)

        print("Time Distribution:\n", table)

    def get_face_counts(self, percentage=0):
        """Prints the table of detailed distribution of events in detector components."""

        f = 1

        if percentage:
            f = 100 / self.total_count

        names = copy.deepcopy(self.comps)
        names.append("TOTAL")
        total_count = 0
        totals = []

        for sim in self.sims:
            totals.append(sum(list(sim.facecounts.values())) * f)

        total_count = sum(totals)

        data = []

        for key in self.comps:
            row = [sim.facecounts[key] * f for sim in self.sims]
            row.append(np.sum(row))
            data.append(row)

        trow = [totals[i] for i in range(len(self.sims))]
        trow.append(total_count)
        data.append(trow)

        parts = [
            MuC.part_names[part[0]] + " events (" + part[1] + ")"
            for part in col.colls_types_to_part[self.design['collision_type']]
        ]

        table = PrettyTable()
        labels = ["Detector Parts"]
        labels.extend(parts)
        labels.append("Total Events")
        table.field_names = labels

        for name, row in zip(names, data):

            if percentage:
                formatted_row = [f"{x:.1f}" for x in row]

            else:
                formatted_row = [f"{x:.3e}" for x in row]

            table.add_row([name] + formatted_row)

        print("Event Distribution:\n", table)
        return np.array(data)
    
    def get_data(self, sec="all", part="all", genie=0):
        """Retrieving data from the sims object.

        Args:
            sec (str): the detector section (component) that one wants to single out. This usually fall in two categories: endcaps (ec) and barrels.
                Set to: all, muon_detector_ec, muon_detector_barrel, ecal_ec, ecal_barrel, hcal_ec, hcal_barrel, solenoid_borders, solenoid_mid, or nozzles.
            part (str): a particle one would want to single out. Can be either nue, nuebar, numu, or numubar.
            genie (bool): to get the n*l factor in too."""

        if (sec != "all") & (sec not in self.comps):
            raise ValueError(
                "This is not an implemented detector component. Choices are: "
                + str(self.comps)
            )

        if (part != "all") & (part not in self.parts):
            raise ValueError(
                "This is not a valid particle! Choices are: " + str(self.parts)
            )

        elif part == "all":
            maskp = np.ones_like(self.sims, dtype=bool)

        else:
            maskp = [(sim.particle == part) for sim in self.sims]

        sims = [item for item, keep in zip(self.sims, maskp) if keep]

        x = np.concatenate([sim.arrx for sim in sims])
        y = np.concatenate([sim.arry for sim in sims])
        z = np.concatenate([sim.arrz for sim in sims])
        w = np.concatenate([sim.w for sim in sims])
        times = np.concatenate([sim.times for sim in sims])
        E = np.concatenate([sim.E for sim in sims])
        costheta = np.concatenate([sim.costheta for sim in sims])
        mask = np.concatenate([sim.get_face_masks(sec) for sim in sims])

        if genie:
            w = np.concatenate([sim.wnl for sim in sims])

        return x[mask], y[mask], z[mask], w[mask], times[mask], E[mask], costheta[mask]

    def plot(
        self,
        nbins=100,
        cmin=1,
        orientation="z-y",
        savefig=None,
        fs=(20, 12),
        cmap="viridis",
        ax=None,
        title=True,
        xl=True,
        yl=True,
        vmin=None,
        vmax=None,
        h=False,
        sec="all",
        colorbar=True,
        part="all",
    ):
        """Plotting the detector event distribution as a hist2d instance, with the detector geometry behind.

        Args:
            orientation (str): coordinate representation for the plot. Can be z-y (default), z-x, or x-y.
            savefig (str): name of file to save plot to.
            fs (tuple): figsize of plot. Can be None to display on the same plot that is being worked on.
            ax (plt.axes obj): if one wants to plot on a specific ax (i.e., subplot), specify.
            title (bool): if one wants to display the pre-generated title.
            xl/yl (bool): to display the x-label and y-label.
            vmin/vmax (float): to change color displays.
            h (bool): to return the plot obj.
            sec (str): the component of the detector one wants to single out. Options are the sames as those written in get_data() method description.
            part (str): a particle one would want to single out. Can be either nue, nuebar, numu, or numubar.
        """

        if not ax:
            fig, ax = plt.subplots(figsize=fs)

        bs = np.linspace(-1 * self.zending, self.zending, nbins)
        bs2 = np.linspace(-1 * self.rmax, self.rmax, nbins)

        x, y, z, w, _, _, _ = self.get_data(sec=sec, part=part)

        # lbl4 = f"Experiment: {self.design}"
        # lbl3 = f"Collision: {col.acc_colls_dict[self.collision]}"
        # lbl = r"$L_{ss} = $" + f"{self.L:.0f} m"
        # lbl2 = r"$N_{events} = $" + f"{np.sum(w):.3e} events"
        # ax.plot(
        #     [
        #         -1 * self.zending,
        #         -1 * self.zending,
        #         self.zending,
        #         self.zending,
        #         -1 * self.zending,
        #     ],
        #     [-1 * self.rmax, self.rmax, self.rmax, -1 * self.rmax, -1 * self.rmax],
        #     color="black",
        #     zorder=50,
        # )
        c_dict = {"z-y": [z, y], "z-x": [z, x], "x-y": [x, y]}
        ha = ax.hist2d(
            c_dict[orientation][0],
            c_dict[orientation][1],
            alpha=1,
            zorder=2.5,
            bins=(bs, bs2),
            weights=w,
            cmin=cmin,
            cmap=cmap,
            norm=LogNorm(),
        )
        ax.clear()
        vmin, vmax = ha[-1].get_clim()
        print(f"vmin: {vmin}, vmax: {vmax}")
        ha = ax.hist2d(
            c_dict[orientation][0],
            c_dict[orientation][1],
            alpha=1,
            zorder=2.5,
            bins=(bs, bs2),
            weights=w,
            cmin=cmin,
            cmap=cmap,
            norm=LogNorm(vmin=vmax / 1e5, vmax=vmax),
        )

        plot_det(self.Geometry, ax, orientation=orientation, xl=xl, yl=yl)

        # ax.legend([lbl4, lbl3, lbl, lbl2], loc='lower right').set_zorder(50)

        if title:
            ax.set_title(f"Event Distribution: {orientation}")

        if savefig:
            plt.savefig(savefig, bbox_inches="tight", dpi=300)

        if colorbar:
            fig.colorbar(ha[3], fraction=0.025, pad=0.04, aspect=30)

        if h:
            return ha

    def event_timing(
        self,
        fs=(20, 12),
        histtype="barstacked",
        nbins=100,
        savefig=None,
        legend=False,
        title=True,
        sec="all",
        part="all",
    ):
        """Wrapper to plot neutrino interaction times.

        Args:
            savefig (str): name of file to save plot to.
            fs (tuple): figsize of plot. Can be None to display on the same plot that is being worked on.
            title (bool): if one wants to display the pre-generated title.
            legend (bool): to display the pre-generated legend.
            sec (str): the component of the detector one wants to single out. Options are the sames as those written in get_data() method description.
            part (str): a particle one would want to single out. Can be either nue, nuebar, numu, or numubar.
        """

        if fs:
            plt.figure(figsize=fs)

        _, _, _, w, times, _, _ = self.get_data(sec=sec, part=part)

        label = f"{sec}; " + r"$N_{events}$" + f": {np.sum(w):.3e}"

        if sec == "all":
            label = r"$L_{ss} = $" + f"{self.design['Lss']:.0f} m"

        times *= 1e9

        plt.xlabel("Time before collision (ns)")
        plt.ylabel(r"$N_{events}$")

        if title:
            plt.title(
                f"Event Timing (wrt bunch crossing); ({self.name} at L = {self.design['Lss']:.2f})"
            )

        plt.hist(times, weights=w, histtype=histtype, bins=nbins, label=label)

        if legend:
            plt.legend(loc="best")

        plt.yscale("log")

        if savefig:
            plt.savefig(savefig, bbox_inches="tight", dpi=300)

    def phi_distribution(
        self,
        fs=(20, 12),
        histtype="step",
        nbins=100,
        savefig=None,
        ylog=True,
        label="",
        legend=False,
        sec="all",
        part="all",
    ):
        """Wrapper to plot the phi distribution of neutrino events.

        Args:
            savefig (str): name of file to save plot to.
            histtype (str): plt histtype options.
            fs (tuple): figsize of plot. Can be None to display on the same plot that is being worked on.
            legend (bool): to display the legend.
            ylog (bool): to set the y-scale as log.
            label (str): Label of the plot.
            sec (str): the component of the detector one wants to single out. Options are the sames as those written in get_data() method description.
            part (str): a particle one would want to single out. Can be either nue, nuebar, numu, or numubar.
        """

        if fs:
            plt.figure(figsize=fs)

        x, y, _, w, _, _, _ = self.get_data(sec=sec, part=part)

        phi = np.arctan(x / y)

        plt.hist(phi, weights=w, histtype=histtype, bins=nbins, label=label)
        plt.ylabel(r"$N_{events}/yr$")
        plt.xlabel(r"$\phi$ Distribution (rad)")

        if legend:
            plt.legend(loc="best")

        if ylog:
            plt.yscale("log")

        if savefig:
            plt.savefig(savefig, bbox_inches="tight", dpi=300)

    def energies(
        self,
        fs=(20, 12),
        histtype="step",
        nbins=100,
        savefig=None,
        label="",
        legend=False,
        linestyle="-",
        sec="all",
        part="all",
    ):
        """Wrapper to plot energy flux of particles.

        Args:
            savefig (str): name of file to save plot to.
            histtype (str): plt histtype options.
            fs (tuple): figsize of plot. Can be None to display on the same plot that is being worked on.
            legend (bool): to display the legend.
            linestyle (str): plt linestyle options.
            sec (str): the component of the detector one wants to single out. Options are the sames as those written in get_data() method description.
            part (str): a particle one would want to single out. Can be either nue, nuebar, numu, or numubar.
        """

        if fs:
            plt.figure(figsize=fs)

        _, _, _, w, _, E, _ = self.get_data(sec=sec, part=part)

        plt.hist(
            E,
            weights=w,
            histtype=histtype,
            bins=nbins,
            label=label,
            linestyle=linestyle,
        )
        plt.ylabel(r"$N_{events}/yr$")
        plt.xlabel(r"$E_{\nu}$ (GeV)")

        if legend:
            plt.legend(loc="best")

        if savefig:
            plt.savefig(savefig, bbox_inches="tight", dpi=300)

    def get_GENIE_flux(self, sec, part, nbins=100, filename=None):
        """Creates a flux .data file for GENIE simulation of events.

        Args:
            filename (str): name of file to be saved in the fluxes/ folder (opt).
            nbins (int): number of bins for the histogram.
            sec (str): detector component to consider.
            part (str): a particle one would want to single out. Can be either nue, nuebar, numu, or numubar.
        """

        _, _, _, w, _, E, _ = self.get_data(sec=sec, part=part, genie=1)
        h = np.histogram(E, weights=w, bins=100)
        flux = h[0]
        eg = h[1]
        egs = [(eg[i] + eg[1 + i]) / 2 for i in range(len(eg) - 1)]

        assert len(egs) == len(flux), "Lists must have the same length"

        if filename:
            fn = f"{filename}"

        else:
            fn = f"fluxes/{self.design['short_name']}_{MuC.compsto2[sec]}_{part}.data"

        with open(fn, "w") as file:
            for item1, item2 in zip(egs, flux):
                file.write(f"{item1}\t{item2}\n")

        print(f"Data has been written to {fn}")


@njit
def get_events_njit1(temp2, mi):
    """Gets weight of each event along the segment of a single ray."""
    res = (1 - temp2) * mi / temp2

    return res


@njit
def get_events_njit2(
    pweights, mid, counts, dec_pos, normed_m, cumulative_distances, t_values
):
    """Gets weight of each event generated by all interacting neutrinos."""
    # pweights is (:, 25); mid is (:,1); mask is (:,); counts is (:,1); dec_pos is (:, 1, 3); normed_m is [:, 1, 3]; cumulative_distances is [:,25]; t_values is [:, 25]
    pweights = pweights / mid  # good multinomial probs
    t_values[:, 1:] += cumulative_distances[:, :-1]
    part_face_counts = counts * pweights  # these are the face counts (expected)
    events_position = dec_pos + t_values[:, :, np.newaxis] * normed_m
    return events_position, part_face_counts


def plot_det(geom, ax, orientation="z-y", xl=True, yl=True):
    """Plots the detector geometry behind a sim plot."""

    if geom == "det_v1":
        T1 = [[-231, 231, 231, 28, -28, -231, -231], [150, 150, 24, 3, 3, 24, 150]]

        ECAL1 = [[-231, -231, 231, 231, -231], [150, 170, 170, 150, 150]]
        ECAL2 = [[231, 231, 251, 251, 231], [24, 170, 170, 26, 24]]
        ECAL3 = [[-1 * i for i in ECAL2[0]], ECAL2[1]]

        HCAL1 = [[-251, -251, 251, 251, -251], [170, 348, 348, 170, 170]]
        HCAL2 = [[251, 251, 418, 418, 251, 251, 251], [170, 348, 348, 43, 26, 170, 170]]
        HCAL3 = [[-1 * i for i in HCAL2[0]], HCAL2[1]]

        SOLENOID = [[-418, -418, 418, 418, -418], [348, 446, 446, 348, 348]]

        MD1 = [[-418, -418, 418, 418, -418], [446, 645, 645, 446, 446]]
        MD2 = [[418, 418, 564, 564, 418], [43, 645, 645, 58, 43]]
        MD3 = [[-1 * i for i in MD2[0]], MD2[1]]

        CONE1 = [[28, 564, 564, 28], [3, 58, 3, 3]]
        CONE2 = [[-1 * i for i in CONE1[0]], CONE1[1]]

        BL = [[-564, -564, 564, 564, -564], [-3, 3, 3, -3, -3]]

        components = [
            T1,
            ECAL1,
            ECAL2,
            ECAL3,
            HCAL1,
            HCAL2,
            HCAL3,
            SOLENOID,
            MD1,
            MD2,
            MD3,
            CONE1,
            CONE2,
            BL,
        ]
        cols = (
            ["lightgrey"] * 1
            + ["dimgrey"] * 3
            + ["grey"] * 3
            + ["darkgrey"]
            + 3 * ["grey"]
            + ["black"] * 2
            + ["white"]
        )

        for i, component in enumerate(components):
            ax.plot(component[0], component[1], color=cols[i])
            ax.fill_between(component[0], component[1], color=cols[i], alpha=0.7)
            new_y = [-1 * k for k in component[1]]
            ax.plot(component[0], new_y, color=cols[i])
            ax.fill_between(component[0], new_y, color=cols[i], alpha=0.7)

        ax.set_xlabel("z-coordinate (cm)")
        ax.set_ylabel("r-coordinate (cm)")

    elif (geom == "det_v2") | (geom == "zero_density_test"):

        if orientation == "x-y":
            MDET = 645
            SPS1 = 446.1
            SOL_1 = 429
            SPS2 = 425
            SOL_2 = 399.3
            SPS3 = 364.9
            SOL_3 = 352.3
            SPS4 = 348.3
            HCAL = 333
            SPS5 = 174
            ECAL = 170.2
            SPS6 = 150
            CONE = 31
            BL = 2.2
            components = [
                MDET,
                SPS1,
                SOL_1,
                SPS2,
                SOL_2,
                SPS3,
                SOL_3,
                SPS4,
                HCAL,
                SPS5,
                ECAL,
                SPS6,
                CONE,
                BL,
            ]
            cols = [
                "grey",
                "white",
                "gray",
                "white",
                "lightgrey",
                "white",
                "gray",
                "white",
                "darkgrey",
                "white",
                "dimgrey",
                "white",
                "black",
                "white",
            ]

            for i, component in enumerate(components):
                circle = plt.Circle(
                    (0, 0),
                    det,
                    zorder=i / len(components),
                    alpha=1,
                    edgecolor=cols[i],
                    facecolor=cols[i],
                )
                ax.add_artist(circle)

            ax.scatter(MDET, MDET, color=cols[3])
            ax.scatter(MDET, MDET, color=cols[3])
            ax.scatter(MDET, MDET, color=cols[3])
            ax.scatter(MDET, MDET, color=cols[3])

            ax.set_xlim(-1 * 645 * 10 / 12, 645 * 10 / 12)
            ax.set_ylim(0, 645)

            if xl:
                ax.set_xlabel("x-coordinate (cm)")

            if yl:
                ax.set_ylabel("y-coordinate (cm)")

        else:
            ECAL1 = [[-221, -221, 221, 221, -221], [150, 170.2, 170.2, 150, 150]]
            ECAL2 = [[230.7, 230.7, 250.9, 250.9, 230.7], [31, 170, 170, 33.9, 31]]
            ECAL3 = [[-1 * i for i in ECAL2[0]], ECAL2[1]]

            HCAL1 = [[-221, -221, 221, 221, -221], [174, 333, 333, 174, 174]]
            HCAL2 = [
                [235.4, 235.4, 412.9, 412.9, 250.9, 250.9, 235.4],
                [170, 324.6, 324.6, 56.8, 33.9, 170, 170],
            ]
            HCAL3 = [[-1 * i for i in HCAL2[0]], HCAL2[1]]

            SOLENOID = [
                [-412.9, -412.9, 412.9, 412.9, -412.9],
                [348.3, 352.3, 352.3, 348.3, 348.3],
            ]
            SOLENOID_2 = [
                [-412.9, -412.9, 412.9, 412.9, -412.9],
                [364.9, 399.3, 399.3, 364.9, 364.9],
            ]
            SOLENOID_3 = [
                [-412.9, -412.9, 412.9, 412.9, -412.9],
                [425, 429, 429, 425, 425],
            ]

            MD1 = [
                [-563.8, -563.8, 563.8, 563.8, 417.9, 417.9, -417.9, -417.9, -563.8],
                [78.2, 645, 645, 78.2, 57.5, 446.1, 446.1, 57.5, 78.2],
            ]

            CONE1 = [
                [6.5, 230.7, 250.9, 412.9, 417.9, 563.8, 563.8, 6.5],
                [2.2, 31, 33.9, 56.8, 57.5, 78.2, 2.2, 2.2],
            ]
            CONE2 = [[-1 * i for i in CONE1[0]], CONE1[1]]

            BL = [[-563.8, -563.8, 563.8, 563.8, -563.8], [-2.2, 2.2, 2.2, -2.2, -2.2]]

            components = [
                ECAL1,
                ECAL2,
                ECAL3,
                HCAL1,
                HCAL2,
                HCAL3,
                SOLENOID,
                SOLENOID_2,
                SOLENOID_3,
                MD1,
                CONE1,
                CONE2,
                BL,
            ]
            cols = (
                ["dimgrey"] * 3
                + ["darkgrey"] * 3
                + ["gray"]
                + ["lightgrey"]
                + ["gray"]
                + ["grey"]
                + ["black"] * 2
                + ["white"]
            )

            for i, component in enumerate(components):
                ax.plot(component[0], component[1], color=cols[i])
                ax.fill_between(component[0], component[1], color=cols[i], alpha=1)
                new_y = [-1 * k for k in component[1]]
                ax.plot(component[0], new_y, color=cols[i])
                ax.fill_between(component[0], new_y, color=cols[i], alpha=1)

            if xl:
                ax.set_xlabel("z-coordinate (cm)")

            if yl:
                ax.set_ylabel("y-coordinate (cm)")

            # ax.set_xlim(-564, 564)
            # ax.set_ylim(-645, 645)

    else:
        print("this geometry has not been implemented yet!")


def load_data(fil, direc=None, n_events=1e5, getQ=False):
    """Loads a GENIE analysis file from $GENANA. Adds weights. Note that SB is solenoid border while SM is solenoid middle (they differ in material).

    Args:
    n_events: number of events that the GENIE file had. Matheus, use 1e5."""

    if direc:
        os.system("source ~/.bashrc")
        directory = os.getenv(direc)
    else:
        directory = "luc_analysis"
    #     raise ValueError(f"Environment variable {direc} is not set. Ensure it's defined in your .bashrc and loaded correctly.")

    filename = os.path.join(directory, fil)

    with open(filename, "r") as file:

        for i, line in enumerate(file):

            if i == 0:
                continue

            elif i == 1:
                exp = (line.split(":")[1])[:-1]

            elif i == 2:
                particles = ((line.split(":")[1])[:-1]).split(",")

            elif i == 3:
                comps = (line.split(":")[1])[:-1]

            else:
                break

    expname = (col.weights[exp])["Name"]
    parts = [MuC.partn_names[part] for part in particles]

    particlenames = ", ".join(parts)
    t = comps.replace(",", ", ")

    print(f"Loading generated data for a {expname} experiment;")
    print(
        f"It includes interactions from {particlenames} within the {t} of the muon detector."
    )

    data = pd.read_csv(filename, sep="\s+", skiprows=5)

    print("Adding weights...")

    data["w"] = data.apply(
        lambda row: get_weight(row["DComp"], row["IncL"], exp, n_events=n_events),
        axis=1,
    )

    if getQ:
        print("Computing Q2...")
        data["Q2"] = data.apply(
            lambda row: get_Q2(row["nu_E"], row["E"], row["pz"]), axis=1
        )

    print("Done!")

    return data


def get_Q2(nu_E, E, pz):
    """Getting Q squared from the generated events."""
    return -1 * ((nu_E - E) ** 2 - (nu_E - pz) ** 2)


def get_weight(comp, p, exp, n_events):
    """Getting the weight of a genie particle from its detector component."""
    try:
        return (
            (col.weights[exp])["total_count"]
            * ((col.weights[exp])[MuC.pdg2names[str(p)]])[comp]
            / 100
            / n_events
        )  # the last factor depends on how many generated events there are in the files. It only supports same n files across detectors.
    except KeyError:
        return 0
