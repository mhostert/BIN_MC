import os
import pickle
import time

import numpy as np
import pandas as pd
from numba import njit

from scipy.stats import uniform
from scipy.optimize import bisect

from prettytable import PrettyTable

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import MuC
from MuC import const
from MuC import mudecay_tools as mudec
from MuC import detector_tools as det
from MuC import collider_tools as col
from MuC import xsecs

from DarkNews import Cfourvec as Cfv

import MuC.plot_tools

N_LIFETIMES = 1  # number of beam lifetimes to consider


def D3distance(point1, point2):
    """3D Euclidian distance between two points."""

    return np.sqrt(
        (point1[:, 0] - point2[:, 0]) ** 2
        + (point1[:, 1] - point2[:, 1]) ** 2
        + (point1[:, 2] - point2[:, 2]) ** 2
    )


def to_opt(x, K):
    """Optimizing provides root for transformed radius, R', of storage ring for nonzero Lss."""
    return np.sin(x) - K * (np.pi - x)


# Precompute func_radius for all theta values
def func_radius(theta):
    """This is the regular parametrization; straight segment comes after. DO NOT CHANGE unless altering the actual curve parametrization."""
    radius = 1630e2
    matching_angular_distance = 25 / 180 * np.pi
    tan_theta = np.abs(np.tan(theta))
    condition = tan_theta > np.tan(matching_angular_distance)
    arc_tan_value = np.where(condition, np.arctan(tan_theta), matching_angular_distance)
    result = (
        (arc_tan_value - matching_angular_distance)
        / (np.pi / 2 - matching_angular_distance)
        * (np.pi / 2)
    )
    return radius + (-radius * 1.2 / 10) * np.sin(result) ** 2


class DetectorSimulator:
    """
    Detector class for MuC

    """

    def __init__(
        self,
        mustorage,
        det_geom,
        save_mem=True,
    ):

        self.mustorage = mustorage
        self.geom = det_geom

        # Detector-related quantities
        self.detname = det_geom.name
        self.iterations = det_geom.iterations
        self.objects = det_geom.OBJECTS
        self.zbeginning = det_geom.zbeginning
        self.rmax = det_geom.rmax
        self.beam_pipe_radius = det_geom.beam_pipe_radius
        self.zending = det_geom.zending
        self.initials = det_geom.INITIALIZERS
        self.decayer = det_geom.DECAYER  # is the minus_one,minus_two
        self.outside = det_geom.OUTSIDE  # parent class of outside-going faces
        self.face_dict = det_geom.facedict
        self.outside_ids = np.array([obj.id for obj in self.outside])

        """ Processing the mustorage simulation """

        self.sample_size = mustorage.sample_size
        self.pos = mustorage.pos
        self.pnu = mustorage.pnu
        self.nuflavor = mustorage.nuflavor
        self.direction = mustorage.direction
        self.x_CM = mustorage.x_CM
        self.costheta_CM = mustorage.costheta_CM

        self.intersection_points = np.full(
            (self.sample_size, self.iterations, 3), 1e4
        )  # 1e4 is arbitrary

        self.x_CM = np.reshape(
            np.repeat(self.x_CM, self.iterations), (self.sample_size, self.iterations)
        )
        self.costheta_CM = np.reshape(
            np.repeat(self.costheta_CM, self.iterations),
            (self.sample_size, self.iterations),
        )

        self.densities = np.zeros((self.sample_size, self.iterations - 1))
        # self.densities[:, 0] = det.EARTH_DENSITY # NOTE: this is probably old?
        self.component_id = np.full(
            (self.sample_size, self.iterations), -1
        )  # starting at initials
        self.intersection_points[:, 0, :] = self.mustorage.pos
        self.distances = np.zeros((self.sample_size, self.iterations - 1))

        self.save_mem = save_mem

    # @profile
    def find_intersections(self):
        """
        Getting all intersection points of neutrinos with detector components.

        Arrays are extended to dimension, such that, e.g.,

            self.intersection_points: (sample_size, iterations, 3)

        with `iterations` the number of iterations of the algorithm

        """

        """ first iteration points """
        count = 1

        # Distance to face of detector
        delta_z = self.zbeginning - self.intersection_points[:, 0, 2]

        # auxiliary parameter for intersection point
        to = delta_z / self.pnu[:, 3]

        # intersection points on the plane of the face of the detector
        ip_at_detface = (
            self.intersection_points[:, 0, :] + to[:, np.newaxis] * self.pnu[:, 1:]
        )

        # ip_at_detface[:, 0] = delta_z * tan(theta) * cos(phi)
        # ip_at_detface[:, 1] = delta_z * tan(theta) * sin(phi)
        # ip_at_detface[:, 2] = delta_z

        # First treat decays inside the detector barrel
        started_outside_barrel = ~(
            (self.intersection_points[:, 0, 2] < self.zending)
            & (self.intersection_points[:, 0, 2] > self.zbeginning)
            & (self.intersection_points[:, 0, 1] > -1 * self.rmax)
            & (self.intersection_points[:, 0, 1] < self.rmax)
        )

        # X-Y of exit point should be inside barrel
        radius_at_detface = np.sqrt(ip_at_detface[:, 0] ** 2 + ip_at_detface[:, 1] ** 2)

        # Projected position in the plane of detector face outside of detector radius or neutrino pointing backwards (NOTE:??)
        hit_outside_detface = radius_at_detface > self.rmax
        pointing_backwards = self.pnu[:, 3] < 0  # NOTE: used to be (to < 0), but why?

        # Events where the neutrino missed altogether
        missed_detface = started_outside_barrel & (
            hit_outside_detface | pointing_backwards
        )

        # indices of this new array
        miss_indices = np.where(missed_detface)[0]

        # If particles do not reach detector, set (component_id = 0) and (intersection_point = parent muon) at 1st iteration
        self.component_id[miss_indices, 1] = 0
        self.intersection_points[miss_indices, 1, :] = self.intersection_points[
            miss_indices, 0, :
        ]

        # NOTE: Again, why to > 0? Is this a requirement on pnu or on deltaz?
        hit_detface = ~missed_detface  # & (to > 0)

        # Loop over initial faces of the detector
        # All these initial faces are on the same plane (z = zbeginning), so we can treat them on the same footing
        for obj in self.initials:

            # get particles in detector that are still left uninitialized
            hit_detface = hit_detface & (self.component_id[:, 1] == -1)

            # indices of those that cross this initial face/component
            in_initial_obj = obj.check_in(radius_at_detface, hit_detface)

            # updating location of the particle to the id of the component
            self.component_id[in_initial_obj, 1] = obj.id

            # position is given directly by the intersection point ip calculated before
            self.intersection_points[in_initial_obj, 1, :] = ip_at_detface[
                in_initial_obj, :
            ]

        if self.save_mem:
            del self.initials

        """ At this point, `component_id` should be one of the following:
                -1 (decaying in detector)
                0  (did not reach detector)
                id of initial face (e.g., id = 1, 2, or 3 for det_v2)
            
            Now we move on to treat decays in beampipe.
        
        """

        in_detector_decay_pipe = (self.component_id[:, 1] == -1) & (
            self.intersection_points[:, 0, 0] ** 2
            + self.intersection_points[:, 0, 1] ** 2
            < self.beam_pipe_radius**2
        )

        # solve for intersection points with beampipe
        a, b, c = det.barrel_get_polynomial(
            self.beam_pipe_radius,
            self.intersection_points[in_detector_decay_pipe, 0, :],
            self.pnu[in_detector_decay_pipe, 1:],
        )
        roots = det.get_roots(a, b, c)

        # NOTE: What is going on here? Why are we rounding to 12 decimals?
        new_mask_1 = np.round(roots[:, 0], decimals=12) > 0
        new_mask_2 = np.round(roots[:, 1], decimals=12) > 0
        t = roots.copy()  # need t to get the correct root
        doubles = new_mask_1 & new_mask_2
        t[new_mask_2, 0] = roots[new_mask_2, 1]

        if np.any(doubles):
            t[doubles, 0] = np.min(roots[doubles])

        # These are intersection points on the plane where nus cross the beampipe edge
        ip_at_beam_pipe_exit = (
            self.intersection_points[in_detector_decay_pipe, 0, :]
            + t[:, 0][:, np.newaxis] * self.pnu[in_detector_decay_pipe, 1:]
        )

        # If events can start inside beampipe or nozzles:
        if self.decayer:

            for neighb in self.decayer[0].next_ids:
                neighbor = self.objects[neighb]
                loc_mask = self.component_id[in_detector_decay_pipe, 1] == -1

                # check if particles hit a neighboring component (otherwise, down the pipe it goes)
                new_indices, in_neighbor_comp = self.decayer[0].check_in(
                    neighbor,
                    ip_at_beam_pipe_exit[:, 2],
                    in_detector_decay_pipe,
                    loc_mask,
                )

                # Update the component indices for the ones that do
                self.component_id[new_indices, 1] = neighb
                ips = ip_at_beam_pipe_exit[in_neighbor_comp, :]
                self.intersection_points[new_indices, 1, :] = ips
                self.densities[new_indices, 0] = self.decayer[0].density

            # NOTE: So, these decays should not happen, maybe we throw an error:
            # decays within the detector components themselves
            weird_decay_mask = self.component_id[:, 1] == -1
            self.update_intersections(self.decayer[1], count, weird_decay_mask)
            self.component_id[weird_decay_mask, 0] = -2
            if weird_decay_mask.sum() > 0:
                print(
                    "WARNING: Decays inside detector components found! Total of ",
                    weird_decay_mask.sum(),
                )

            if self.save_mem:
                del self.decayer

        # there should be no more -1s, just 0 - x.
        # now iteratively finding intersection points

        """ At this point, `component_id` should be one of the following:
                0  (did not reach detector)
                id of a component
                -2 (weird decays inside detector components)
            
            Now we move on to propagate the rays through the detector components.
        
        """

        # NOTE: CHECKED UP UNTIL HERE.
        count = 2
        while count < self.iterations:

            locs = np.unique(self.component_id[:, count - 1])
            objects = [self.objects[i] for i in locs]

            for obj in objects:
                part_mask = self.component_id[:, count - 1] == obj.id
                particles = np.where(part_mask)[
                    0
                ]  # these are indices of the particles at that face

                if obj.id in self.outside_ids:
                    self.intersection_points[particles, count, :] = (
                        self.intersection_points[particles, count - 1, :]
                    )
                    self.component_id[particles, count] = self.component_id[
                        particles, count - 1
                    ]
                    continue

                if np.any(particles):
                    self.update_intersections(obj, count, part_mask)

            count += 1

        if self.save_mem:
            del self.objects
            del self.outside
            del self.outside_ids

    # @profile
    def compute_interaction_rate(self):
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
        )  # huge array of distances in each component per ray
        for_mask = np.sum(self.part_line_integrals[:, 1:], axis=1)
        line_integrals = for_mask + self.part_line_integrals[:, 0]

        # Total event rate (total nu cross section)
        self.cs = xsecs.total_xsecs[self.nuflavor](self.pnu[:, 0])
        probs = 1 - np.exp(-1 * self.cs * line_integrals)

        self.counts = self.mustorage.weights * probs[:, np.newaxis]  # (sample_size, 1)
        self.in_acceptance = self.counts > 0
        if (self.counts < 0).sum() > 0:
            print("Negative counts found!")
        self.in_acceptance = self.in_acceptance[:, 0]
        self.counts = self.counts[self.in_acceptance]
        print(
            f"Efficiency of detector acceptance: {self.in_acceptance.sum() / self.sample_size:.2e}",
        )

        del self.densities

    def get_exclusive_rates(self):
        """
        Get the rate of exclusive interaction channels for this neutrino.

        Returns:
            dict: rate of rare events
        """
        self.exclusive_rates = {}
        Z = np.array([self.objects[comp].material.Z for comp in self.component_id])
        A = np.array([self.objects[comp].material.A for comp in self.component_id])
        for key in self.face_dict.keys():
            face_mask = self.get_face_masks(key)

            # NOTE: Fix-me to include other nuclei
            channels = xsecs.get_cross_sections(
                self.E[face_mask], self.nuflavor, Z[face_mask], A[face_mask]
            )

            for c, xsec in channels.items():
                self.exclusive_rates[(key, c)] = np.sum(
                    xsec / channels[f"{self.nuflavor}_total"] * self.w[face_mask]
                )
        return self.exclusive_rates

    def print_exclusive_counts(self, percentage=False):

        if not hasattr(self, "exclusive_rates"):
            self.get_exclusive_rates()

        sample_dict = self.exclusive_rates

        # Step 1: Identify unique rows and columns
        rows = sorted(set(key[0] for key in sample_dict.keys()))
        columns = sorted(set(key[1] for key in sample_dict.keys()))

        # Step 2: Initialize PrettyTable with columns
        table = PrettyTable([""] + columns)
        # Step 3: Populate the table with data
        for row in rows:
            row_data = [row]  # Start with the row label
            for column in columns:
                # Add the corresponding data from the dictionary, if it exists
                value = sample_dict[row, column]
                total = sample_dict[row, f"{self.nuflavor}_total"]
                if percentage and row != f"{self.nuflavor}_total":
                    row_data.append(f"{value/total*100:.1f}%")
                else:
                    row_data.append(f"{value:.1e}")

            table.add_row(row_data)

        return table

    # @profile
    def compute_interaction_rate_njit(self):
        """Gets the number of events (weights) of each interacting neutrino."""

        self.t_values = self.t_values.reshape(
            (self.t_values.shape[1], self.t_values.shape[2])
        )
        surv_prob_in_this_element = np.ones(
            (self.distances.shape[0], self.iterations - 1)
        )  # probs of non-interaction in each segment
        surv_prob_in_this_element = np.exp(
            -1 * self.cs[:, np.newaxis] * self.part_line_integrals
        )

        # product of survival probabilities along trajectory
        total_surv_prob = np.prod(surv_prob_in_this_element, axis=1)

        # Probability of interaction in a given segment
        pweights = (
            (1 - surv_prob_in_this_element)
            * total_surv_prob[:, np.newaxis]
            / surv_prob_in_this_element
        )

        # Where is the probability of interaction largest?
        max_index = np.argmax(pweights[:], axis=1)
        mask = np.ones_like(pweights, dtype=bool)[:, 0]
        max_p = pweights[mask, max_index]
        max_p = max_p[:, np.newaxis]

        # Normalize the probabilities
        pweights = pweights / max_p

        mid = np.sum(pweights, axis=1)
        cumulative_distances = np.cumsum(self.distances, axis=1)
        normed_m = (self.pnu[:, 1:].T / self.pnu[:, 0]).T
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
        self.pos = self.pos[self.in_acceptance]
        self.pnu = self.pnu[self.in_acceptance]
        self.dec_pos = self.dec_pos[self.in_acceptance]
        self.distances = self.distances[self.in_acceptance]
        self.part_line_integrals = self.part_line_integrals[self.in_acceptance]
        self.intersection_points = self.intersection_points[self.in_acceptance]
        self.component_id = self.component_id[self.in_acceptance]
        self.cs = self.cs[self.in_acceptance]
        self.mutimes = self.mustorage.mutimes[self.in_acceptance, np.newaxis]
        self.mutimes_to_bunchx = self.mustorage.mutimes_to_bunchx[
            self.in_acceptance, np.newaxis
        ]
        self.weights = self.mustorage.weights[self.in_acceptance]

        # decay properties
        self.x_CM = self.mustorage.x_CM[self.in_acceptance, np.newaxis]
        self.costheta_CM = self.mustorage.costheta_CM[self.in_acceptance, np.newaxis]

        self.t_values = np.empty(
            (1, np.sum(self.in_acceptance), self.iterations - 1)
        )  # same size as other big arrays now
        self.t_values[:, :, :] = [uniform.rvs(loc=0, scale=self.distances)]

        self.compute_interaction_rate_njit()

        return

    # @profile
    def update_intersections(self, obj, count, mask):
        """For a single object, with particles on it, finds their next component, iteratively on each neighbor, at a specific step of the sim."""

        for neighbor in obj.next_ids:
            neigh = self.objects[neighbor]
            new_indices, ips = neigh.check_intersection(
                self.intersection_points[:, count - 1, :], self.pnu[:, 1:], mask
            )
            mask_1 = (ips[:, 2] < self.intersection_points[new_indices, count, 2]) & (
                ips[:, 2] > self.intersection_points[new_indices, count - 1, 2]
            )
            accepted_ix = new_indices[mask_1]

            self.intersection_points[accepted_ix, count, :] = ips[mask_1]
            self.component_id[accepted_ix, count] = neigh.id
            self.densities[accepted_ix, count - 1] = obj.density

    # @profile
    def calculate_facecounts(self):
        """Gets the facecounts for each detector component in the simulation."""
        self.facecounts = {}

        for key in self.face_dict.keys():
            face_mask = self.get_face_masks(key)
            self.facecounts[key] = np.sum(self.w[face_mask])

    def reweigh_with_new_polarization(self, new_polarization):

        numerator = mudec.mudecay_matrix_element_sqr(
            self.x_CM,
            self.costheta_CM,
            new_polarization,
            self.mustorage.muon_charge,
            self.mustorage.nuflavor,
            self.mustorage.NLO,
        )
        denominator = mudec.mudecay_matrix_element_sqr(
            self.x_CM,
            self.costheta_CM,
            self.mustorage.muon_polarization,
            self.mustorage.muon_charge,
            self.mustorage.nuflavor,
            self.mustorage.NLO,
        )
        # NOTE: Changing polarization does not change the total integral, so this is a safe operation
        self.w_new_pol = self.w * numerator / denominator

        return self.w_new_pol

    def calculate_facecounts_new_pol(self, pol):
        """Gets the facecounts for each detector component in the simulation."""
        self.facecounts_new_pol = {}

        for key in self.face_dict.keys():
            face_mask = self.get_face_masks(key)
            self.facecounts_new_pol[key] = np.sum(
                self.reweigh_with_new_polarization(pol)[face_mask]
            )
        return self.facecounts_new_pol

    def get_face_masks(self, sec):
        """Generate a masks for a specific detector component."""

        if sec == "all":
            mask = np.ones_like(self.w, dtype=bool)

        else:
            mask = np.isin(self.component_id, self.face_dict[sec])

        return mask

    # # @profile
    # def get_lum_q(self, param):
    #     """Getting the luminosity parameters for a beam collision."""
    #     self.intersection_points = self.intersection_points[self.new_indices, 1, :]
    #     self.weights = self.weights[self.new_indices]
    #     self.Nnu = np.sum(self.weights) / col.parameters[param]["syr"]
    #     wax = np.average(self.intersection_points[:, 0], weights=self.weights[:, 0])
    #     way = np.average(self.intersection_points[:, 1], weights=self.weights[:, 0])
    #     self.nusdx = np.sqrt(
    #         np.average(
    #             (self.intersection_points[:, 0] - wax) ** 2, weights=self.weights[:, 0]
    #         )
    #     )
    #     print(f"{self.nusdx*10:.3} mm for sigma x")
    #     self.nusdy = np.sqrt(
    #         np.average(
    #             (self.intersection_points[:, 1] - way) ** 2, weights=self.weights[:, 0]
    #         )
    #     )
    #     print(f"{self.nusdy*10:.3} mm for sigma y")

    def run(self):
        """Runs the simulation for a single neutrino species."""

        if (self.nuflavor in MuC.anti_neutrinos) | (self.nuflavor in MuC.neutrinos):

            # !!!!!!!!!!
            # NOTE: let's try better variable names, more comments, less array dimension manipulations, etc.
            self.find_intersections()
            self.compute_interaction_rate()
            self.get_event_positions()
            # !!!!!!!!!!

            # NOTE: What is -1? What is mask doing?
            mask = self.component_id[:, 0] == -1
            self.part_face_counts[mask, 0] = 0
            self.w = self.part_face_counts.flatten()

            # Neutrino times (t = 0 is muon bunch crossing)
            self.times = self.t_values / const.c_LIGHT + self.mutimes_to_bunchx  # sec
            self.t_values = self.t_values.flatten()
            self.E = self.pnu[:, 0][:, np.newaxis] * np.ones(
                self.part_face_counts.shape
            )
            self.weights = self.weights * np.ones(self.part_face_counts.shape)
            self.costheta = self.costheta[:, np.newaxis] * np.ones(
                self.part_face_counts.shape
            )

            # NOTE: I was getting events far away from the detector. I am not sure why... maybe ray tracing was failing?
            # NOTE: I now enfore the two conditions below to avoid this.
            # NOTE: It could be that the first iteration is counted when it's just the decay of the muon...
            new_mask = (
                (self.w > 0)
                & (np.abs(self.events_position[:, :, 2].flatten()) < 10e2)
                & (np.abs(self.events_position[:, :, 1].flatten()) < 10e2)
            )

            # NOTE: We spent a lot of time tweaking array dimensions. Now we flatten it all.
            # NOTE: It would be much better to keep all arrays 1D from the start.
            self.E = self.E.flatten()[new_mask]
            self.weights = self.weights.flatten()[new_mask]

            self.costheta = self.costheta.flatten()[new_mask]
            self.arrx = self.events_position[:, :, 0].flatten()[new_mask]
            self.arry = self.events_position[:, :, 1].flatten()[new_mask]
            self.arrz = self.events_position[:, :, 2].flatten()[new_mask]
            # self.intersection_points[:, :, :] = self.intersection_points[new_mask, :, :]

            self.t_values = self.t_values[new_mask]
            self.distances = self.distances.flatten()[new_mask]
            nl = self.part_line_integrals.flatten()[new_mask]
            self.wnl = nl * self.weights

            ones = np.ones_like(self.times)
            self.mutimes = ones * self.mutimes
            self.mutimes = self.mutimes.flatten()[new_mask]

            self.mutimes_to_bunchx = ones * self.mutimes_to_bunchx
            self.mutimes_to_bunchx = self.mutimes_to_bunchx.flatten()[new_mask]

            self.component_id = self.component_id[:, :-1].flatten()[new_mask]

            self.times = self.times.flatten()[new_mask]
            self.w = self.w[new_mask]

            self.x_CM = self.x_CM * ones
            self.costheta_CM = self.costheta_CM * ones

            self.x_CM = self.x_CM.flatten()[new_mask]
            self.costheta_CM = self.costheta_CM.flatten()[new_mask]

            if self.direction == "right":
                self.arrx = -1 * self.arrx  # NOTE: Why are we flipping the x-axis?
                self.arrz = -1 * self.arrz

            return self

        else:
            raise ValueError("No particle of that name included in this simulation!")

    # @profile
    def clear_mem(self):
        """Freeing memory at the end of a sim, if necessary"""

        deletables = [
            "Racc",
            "weights",
            "sample_size",
            "counts",
            "part_line_integrals",
            "f",
            "zbeginning",
            "beam_pipe_radius",
            "iterations",
            "cs",
            "in_acceptance",
            "distances",
            "t_values",
            "dec_pos",
            "part_face_counts",
            "events_position",
        ]

        for att in deletables:
            if hasattr(self, att):
                delattr(self, att)


class MuStorageRingSimulator:
    """

    Monte Carlo generation of muon decays from a single beam

    """

    # @profile
    def __init__(
        self,
        design,
        nuflavor=None,
        direction="left",
        N_evals=1e5,
        preloaded_events=None,
        beam_lifetime=None,
        remove_ring_fraction=0,
        NLO=True,
    ):
        """
        Initialize the MuStorageRingSimulator class:

            This treats the decay of a single type of neutrino emitted by a circulating muon beam.

        Parameters:

        design (dict): Dictionary containing the design parameters.
        nuflavor (str, optional): Neutrino flavor, either 'numu' or 'nuebar'. Defaults to None.
        direction (str, optional): Direction of the beam, either 'left' or 'right'. Defaults to "left".
        N_evals (float, optional): Number of evaluations. Defaults to 1e5.
        preloaded_events (dict, optional): Dictionary of preloaded events. Defaults to None.
        beam_lifetime (float, optional): Lifetime of the beam -- different from muon lifetime if beam is dumped early or late. Defaults to muon lifetime in lab frames.
        remove_ring_fraction (float, optional): Fraction of the ring to be removed. Can be a tuple for beginning and end points to be removed or a float if it's the same. Defaults to 0 (entire ring).
        NLO (bool, optional): Next-to-leading order muon decay (radiative corrections). Defaults to True.

        Raises:

        ValueError: If muon_polarization is not between -1 and 1.
        """

        self.design = design
        self.N_evals = N_evals
        self.preloaded_events = preloaded_events
        self.remove_ring_fraction = remove_ring_fraction

        self.direction = direction
        self.nuflavor = nuflavor
        self.muon_charge = (
            -1 if (self.nuflavor == "numu" or self.nuflavor == "nuebar") else +1
        )
        self.muon_polarization = self.design["muon_polarization"]
        self.NLO = NLO

        self.muon_lifetime = const.tau0_mu * self.design["beam_p0"] / const.m_mu
        if beam_lifetime is not None:
            self.beam_lifetime = beam_lifetime
        else:
            self.beam_lifetime = self.muon_lifetime

        if abs(self.muon_polarization) > 1:
            raise ValueError(
                f"muon_polarization must be between -1 and 1, found design['muon_polarization'] = {self.muon_polarization}"
            )

        if self.preloaded_events:
            for var in self.preloaded_events.keys():
                setattr(self, var, preloaded_events[var])

    #######################################
    # Use vegas to simulate mu decays
    def simulate_decays(
        self,
        NLO=True,
        model="NLOmudecay_pol",
    ):
        """simulate_decays using get_MC_events() routine

        NOTE: This is for a given neutrino flavor only. Run this twice with different falvors to get both neutrinos.

        Correlations not implemented unless using a different `model`.

        Returns
        -------
        pd.DataFrame
            Dataframe containing all muon decay events.
        """

        std_tolerance = 4
        # muon helicities h #
        events = mudec.MC_events(
            model=model,
            Mparent=const.m_mu,
            Mdaughter=const.m_e,
            NLO=NLO,
            muon_polarization=self.muon_polarization,
            muon_charge=self.muon_charge,
            pmin=(1 - std_tolerance * self.design["beam_dpop"])
            * self.design["beam_p0"],  # minimum momentum
            pmax=(1 + std_tolerance * self.design["beam_dpop"])
            * self.design["beam_p0"],  # maximum momentum
            beam_p0=self.design["beam_p0"],  # average momentum (GeV)
            beam_dpop=self.design["beam_dpop"],  # little beam spread
            nuflavor=self.nuflavor,
            NINT=20,
            NINT_warmup=10,
            NEVAL=self.N_evals,
            NEVAL_warmup=self.N_evals / 10,
        )

        df_gen = events.get_MC_events()  # gets all events for +1 helicity

        # Muon 4-momenta
        self.pmu = df_gen["P_decay_mu"].to_numpy()  # all muon decaying momenta

        # Neutrino 4-momenta
        self.pnu = df_gen["P_decay_nu"].to_numpy()

        # rest frame variables in integration for ease of reweighting
        self.x_CM = df_gen["x_CM"].to_numpy()
        self.costheta_CM = df_gen["costheta_CM"].to_numpy()

        # Normalized weights. Adds up to 1 for a single muon decay.
        self.weights = df_gen["w_flux"].to_numpy() / np.sum(df_gen["w_flux"].to_numpy())
        self.weights *= self.design["Nmu_per_bunch"]
        self.weights = self.weights[:, np.newaxis]

        self.sample_size = np.size(self.pmu[:, 0])

        return df_gen

    def decay_muons(self):
        """
        Simulates the decay of muons and calculates their positions and velocities.

        If preloaded events are available, the function returns immediately.

        Otherwise, it simulates the decays using the specified NLO (Next-to-Leading Order) setting.

        Attributes:
            preloaded_events (bool): Flag indicating if events are preloaded.
            NLO (bool): Flag indicating if Next-to-Leading Order corrections are used.
            sample_size (int): Number of muon samples.
            pmu (np.ndarray): Array containing the momenta of the muons.
            pnu (np.ndarray): Array containing the momenta of the neutrinos.
            pos (np.ndarray): Array to store the positions of the muons.
            vmu (np.ndarray): Array to store the absolute velocities of the muons.
            momenta (np.ndarray): Array to store the momenta of the neutrinos.
            E (np.ndarray): Array to store the energies of the neutrinos.

        Returns:
            self: The instance of the class with updated attributes.
        """

        if self.preloaded_events:
            return self
        else:
            self.simulate_decays(NLO=self.NLO)

            # Muon decay position (initialized to 0)
            self.pos = np.zeros((self.sample_size, 3))

            self.vmu = const.c_LIGHT * np.sqrt(
                (1 - (const.m_mu / self.pmu[:, 0]) ** 2)
            )  # speed of muons

    def reweigh_with_new_polarization(self, new_polarization):

        numerator = mudec.mudecay_matrix_element_sqr(
            self.x_CM,
            self.costheta_CM,
            new_polarization,
            self.muon_charge,
            self.nuflavor,
            self.NLO,
        )
        denominator = mudec.mudecay_matrix_element_sqr(
            self.x_CM,
            self.costheta_CM,
            self.muon_polarization,
            self.muon_charge,
            self.nuflavor,
            self.NLO,
        )
        # NOTE: Changing polarization does not change the total integral, so this is a safe operation
        self.weights_new_pol = self.weights.flatten() * numerator / denominator

        return self.weights_new_pol

    def place_muons_on_simplified_ring(self, C, Lss, direction):
        """
        Changes the coordinate axis, position, and momenta of particles to fit a storage ring geometry.

        zlength is the length of the detector on one side.

        """

        if Lss == -1:
            return self

        elif Lss == 0:
            Lss = 0.01

        # meter to cm
        self.L = Lss

        K = Lss / (C - Lss)
        a = 1e-6
        b = np.pi - 1e-6

        # eta is the starting angle of the circular segment (end of straight section)
        eta = bisect(to_opt, a, b, args=(K,))

        # Racc is the radius of the circular segment
        self.Racc = Lss / 2 / np.sin(eta)  # cm

        # spread muons in time according to number of beam lifetimes desired
        self.mutimes = Cfv.random_generator(
            self.sample_size,
            0,
            N_LIFETIMES * self.beam_lifetime,
        )

        self.s_in_turn = (self.mutimes * self.vmu) % C

        # Now, if we want to increase our efficiency in the simulation, we better force particles to be close to the detector in some way.
        # Let's enforce this by checking what range of z_in_this_turn is most optimal:
        zacc_min = C * self.remove_ring_fraction[0]
        zacc_max = C * self.remove_ring_fraction[1]
        events_likely_within_acceptance = (self.s_in_turn <= zacc_min) | (
            self.s_in_turn > zacc_max
        )

        shifted_events = ~events_likely_within_acceptance
        deltaz_min = zacc_max - self.s_in_turn
        deltaz_max = C - self.s_in_turn
        shift_z = np.zeros(self.sample_size)
        shift_z[shifted_events] = (
            Cfv.random_generator(shifted_events.sum(), 0, 1)
            * (deltaz_max[shifted_events] - deltaz_min[shifted_events])
            + deltaz_min[shifted_events]
        )

        # Apply shift to s_in_this_turn and to the travel time of the muons
        self.s_in_turn = self.s_in_turn + shift_z
        self.mutimes += shift_z / self.vmu

        # Acceptance of simulated region
        self.weights[:, 0] = self.weights[:, 0] * (zacc_min + (C - zacc_max)) / C

        # Apply exponential suppression on total length travelled by muons
        self.weights[:, 0] = self.weights[:, 0] * (
            1 - np.exp(-self.mutimes / self.muon_lifetime)
        )

        # self.pos is currently entirely in z direction
        # turn_number = self.pos[:, 2] // C
        self.pos[:, 2] = self.s_in_turn  # self.mutimes * self.vmu

        if not Lss:
            on_straight_mask_left = [False] * self.sample_size
            on_straight_mask_right = [False] * self.sample_size

        else:
            on_straight_mask_left = self.s_in_turn > (C - Lss / 2)
            on_straight_mask_right = self.s_in_turn < Lss / 2

        # for straight segment decays
        on_circ = ~((on_straight_mask_right) | (on_straight_mask_left))

        # Origin of system is at the center of the straight segment
        self.pos[on_straight_mask_right, 2] = self.s_in_turn[on_straight_mask_right]
        self.pos[on_straight_mask_left, 2] = self.s_in_turn[on_straight_mask_left] - C

        # distance travelled by particles in the circular segment (not counting the straight segment)
        z_in_circ = self.s_in_turn[on_circ] - Lss / 2

        # angle travelled by particles in the circular segment
        phis = (z_in_circ) / self.Racc + eta

        ###################
        # Beam divergence
        theta_x = Cfv.random_normal(
            np.zeros(self.sample_size),
            self.design["beam_dtheta"] * np.ones(self.sample_size),
        )
        theta_y = Cfv.random_normal(
            np.zeros(self.sample_size),
            self.design["beam_dtheta"] * np.ones(self.sample_size),
        )

        # Rotate by beam divergence envolope
        self.pmu = Cfv.rotationx(self.pmu, -theta_x)
        self.pmu = Cfv.rotationy(self.pmu, -theta_y)

        self.pnu = Cfv.rotationx(self.pnu, -theta_x)
        self.pnu = Cfv.rotationy(self.pnu, -theta_y)

        # rotation of neutrino momenta
        self.pmu[on_circ, :] = Cfv.rotationx(self.pmu[on_circ], phis)
        self.pnu[on_circ, :] = Cfv.rotationx(self.pnu[on_circ], phis)

        # position of particles in the circular segment
        self.pos[on_circ, 2] = (self.Racc + self.pos[on_circ, 1]) * np.sin(phis)
        self.pos[on_circ, 1] = (self.Racc + self.pos[on_circ, 1]) * np.cos(
            phis
        ) - Lss / 2 / np.tan(eta)

        t_per_turn = C / self.vmu
        time_in_this_turn = self.mutimes % t_per_turn
        self.mutimes_to_bunchx = np.where(
            time_in_this_turn < t_per_turn / 2,  # in first half of turn
            time_in_this_turn,  # time wrt bunch xs is + (past crossing)
            time_in_this_turn - t_per_turn,  # time wrt bunch xs is - (future crossing)
        )

        # if right moving, then mirror trajectories through y axis
        if direction == "right":
            self.pos[:, 0] = -1 * self.pos[:, 0]
            self.pnu = Cfv.rotationx(self.pnu, np.full(self.sample_size, np.pi))
            self.pmu = Cfv.rotationx(self.pmu, np.full(self.sample_size, np.pi))

        return self

    def place_muons_on_lattice(self, direction, lattice):
        """Changes the coordinate axis, position, and momenta of particles to fit a storage ring geometry.


        For the lattice, x is horizontal and y is vertical directions in the transverse plane.

        lattice: dictionary with smooth function that describes the lattice (created from interpolation of .tfs files)

            each function returns the value of the lattice parameter as a function of a parameter `u`:
                `u` goes from 0 to 1, with 0 being the IP and 1 being the IP again after a full central orbit.

            lattice['x']: x position of the muons
            lattice['y']: y position of the muons
            lattice['s']: s displacement of the muons along lattice
            lattice['t']: time of the muons
            lattice['angle_of_central_p']: angle of the central momentum of the muons (tangent to the central orbit)
            lattice['beamsize_x']: beam size in x direction
            lattice['beamsize_y']: beam size in y direction
            lattice['beamdiv_x']: beam divergence in x direction
            lattice['beamdiv_y']: beam divergence in y direction
            lattice['dispersion_Dx']: dispersion in x direction
            lattice['dispersion_Dpx']: dispersion in angular x direction

            lattice['inv_s']: inverse function of s(u) -- it returns u(s) s.t. lattice['inv_s'](lattice['s'](u)) = u

        NOTE: the difference between coordinate systems in MuC and in the lattice

          Lattice reference:
         x -- horizontal (ring) plane
         y -- vertical plane
         z -- longitudinal along motion

         In MuC code:
         x -- normal to the ring (downwards when looking from the IP to the center of the ring)
         y -- horizontal plane (+y exits moves outwards from IP away from the ring)
         z -- along motion tangent to the ring at the IP (+z is a beam from left to right at the IP)

         So in our code, horizontal ~ y, vertical ~ x, and z ~ z

        """

        if not isinstance(lattice, dict):
            raise ValueError("Need a dictionary that describes the lattice.")

        # spread muons in time according to number of beam lifetimes desired
        self.mutimes = Cfv.random_generator(
            self.sample_size,
            0,
            N_LIFETIMES * self.beam_lifetime,
        )

        # get total length of the central orbit
        C = lattice["s"](1)  # cm

        self.s_in_turn = (self.mutimes * self.vmu) % C

        # Now, if we want to increase our efficiency in the simulation, we better force particles to be close to the detector in some way.
        # Let's enforce this by clipping the s_in_this_turn range:
        zacc_min = C * self.remove_ring_fraction[0]
        zacc_max = C * self.remove_ring_fraction[1]
        events_likely_within_acceptance = (self.s_in_turn <= zacc_min) | (
            self.s_in_turn > zacc_max
        )

        # Now we move muons along an extra distance to fall within acceptance, paying the price of a smaller survival probability
        shifted_events = ~events_likely_within_acceptance
        deltaz_min = zacc_max - self.s_in_turn
        deltaz_max = C - self.s_in_turn + zacc_min
        shift_z = np.zeros(self.sample_size)
        shift_z[shifted_events] = (
            Cfv.random_generator(shifted_events.sum(), 0, 1)
            * (deltaz_max[shifted_events] - deltaz_min[shifted_events])
            + deltaz_min[shifted_events]
        )

        # Apply shift to s_in_this_turn and to the travel time of the muons
        self.s_in_turn = self.s_in_turn + shift_z
        self.mutimes += shift_z / self.vmu
        self.s_in_turn = self.s_in_turn % C

        # Acceptance of simulated region
        self.weights[:, 0] = self.weights[:, 0] * (zacc_min + (C - zacc_max)) / C

        # Apply exponential suppression on total length travelled by muons
        self.weights[:, 0] = self.weights[:, 0] * (
            # 1 - np.exp(-self.mutimes / self.muon_lifetime)
            np.exp(-self.mutimes / self.muon_lifetime)
        )

        # Now deform locations to real space along the lattice

        # put everyone in the z axis
        self.pos[:, 2] = self.s_in_turn

        # parameter u that goes from 0 to 1 along the lattice
        u_parameter = lattice["inv_s"](self.s_in_turn)

        # Straight outta lattice parameterization
        self.pos[:, 2] = lattice["x"](u_parameter)
        self.pos[:, 1] = lattice["y"](u_parameter)
        self.pos[:, 0] = np.zeros(self.sample_size)

        # Rotation in 2D commutes, so as long as we only rotation in transverse plane
        theta_x = Cfv.random_normal(
            np.zeros(self.sample_size),
            lattice["beamdiv_x"](u_parameter),
        )
        theta_y = Cfv.random_normal(
            np.zeros(self.sample_size),
            lattice["beamdiv_y"](u_parameter),
        )

        # Rotate by beam divergence envolope
        self.pmu = Cfv.rotationx(self.pmu, -theta_x)
        self.pmu = Cfv.rotationy(self.pmu, -theta_y)

        self.pnu = Cfv.rotationx(self.pnu, -theta_x)
        self.pnu = Cfv.rotationy(self.pnu, -theta_y)

        # Rotate to central orbit
        theta_central_orbit = lattice["angle_of_central_p"](u_parameter)
        self.pnu = Cfv.rotationx(self.pnu, -theta_central_orbit)
        self.pmu = Cfv.rotationx(self.pmu, -theta_central_orbit)

        # Rotation in 2D commutes, so as long as we only rotation in transverse plane
        x_horizontal = Cfv.random_normal(
            np.zeros(self.sample_size),
            lattice["beamsize_x"](u_parameter),
        )
        x_vertical = Cfv.random_normal(
            np.zeros(self.sample_size),
            lattice["beamsize_y"](u_parameter),
        )

        # Now shift coordinated by location of beam envolope
        # vertical component is trivial
        self.pos[:, 0] = self.pos[:, 0] + x_vertical
        self.pos[:, 2] = self.pos[:, 2] + x_horizontal * np.sin(theta_central_orbit)
        self.pos[:, 1] = self.pos[:, 1] + x_horizontal * np.cos(theta_central_orbit)

        # Shift time so t = 0 is bunch crossing (NOTE: mutimes will always be negative.)
        t_per_turn = C / self.vmu
        time_in_this_turn = self.mutimes % t_per_turn
        self.mutimes_to_bunchx = np.where(
            time_in_this_turn < t_per_turn / 2,  # in first half of turn
            time_in_this_turn,  # time wrt bunch xs is + (past crossing)
            time_in_this_turn - t_per_turn,  # time wrt bunch xs is - (future crossing)
        )

        # if right moving, then mirror trajectories through y axis
        if direction == "right":
            self.pos[:, 0] = -1 * self.pos[:, 0]
            self.pnu = Cfv.rotationx(self.pnu, np.full(self.sample_size, np.pi))
            self.pmu = Cfv.rotationx(self.pmu, np.full(self.sample_size, np.pi))

        return self

    def get_flux_at_generic_location(
        self,
        det_location=[0, 0, 1e5],
        det_radius=1e2,
        ebins=100,
        acceptance=False,
        per_area=True,
        new_polarization=None,
        normalization=1,
    ):

        det_location = np.array(det_location)
        Ldet = np.sqrt(np.sum(det_location**2))

        # Determine normal vector and a point on the plane using the distance
        self.norm_det = det_location / np.linalg.norm(det_location)
        P0 = Ldet * self.norm_det

        # Calculate direction unit vectors
        d = self.pnu[:, 1:] / np.linalg.norm(self.pnu[:, 1:], axis=1, keepdims=True)

        # Calculate the parameter t for intersections
        numerator = np.dot((P0 - self.pos), self.norm_det)
        denominator = np.dot(d, self.norm_det)
        t = numerator / denominator  # Parameter for each neutrino trajectories

        # Calculate intersection points
        pos_on_det_plane = self.pos + t[:, np.newaxis] * d

        if isinstance(det_radius, tuple) or isinstance(det_radius, list):
            in_acceptance = (
                np.sqrt(np.sum((pos_on_det_plane - P0) ** 2, axis=1)) < det_radius[1]
            ) & (np.sqrt(np.sum((pos_on_det_plane - P0) ** 2, axis=1)) > det_radius[0])

            # area of detector
            area = np.pi * (det_radius[1] ** 2 - det_radius[0] ** 2)

        else:
            in_acceptance = (
                np.sqrt(np.sum((pos_on_det_plane - P0) ** 2, axis=1)) < det_radius
            )

            # area of detector
            area = np.pi * det_radius**2
        if new_polarization:
            self.reweigh_with_new_polarization(new_polarization)
            w = self.weights_new_pol
        else:
            w = self.weights[:, 0]
        nu_eff_ND = np.sum(w[in_acceptance]) / np.sum(w)

        if acceptance:
            print("Detector acceptance: {:.2e}".format(nu_eff_ND))
            return nu_eff_ND

        if nu_eff_ND > 0:
            Enu_ND, flux_nu_ND_p = get_flux(
                self.pnu[in_acceptance, 0], w[in_acceptance], ebins
            )

            if per_area:
                flux_nu_ND = normalization * flux_nu_ND_p / area / np.diff(Enu_ND)
            else:
                flux_nu_ND = normalization * flux_nu_ND_p / np.diff(Enu_ND)

            return Enu_ND, flux_nu_ND
        else:
            print("No flux through detector.")
            return 0, 0


class BINSimulator:
    """
      Main class of MuC library

      Simulates all beam-induced neutrino interactions in a given MuC detector.

      Hierarchy is as follows:
          Initializes by creating instances of DecaySimulation.

          self.run(), runs simulations of many SimNeutrinos based on the collision type,
    which is a single-neutrino-species MC generation of events within a detector, which are all saved in a list in the .sims attribute.
    """

    def __init__(
        self,
        design,
        N_evals=1e5,
        preloaded_events=None,
        det_geom="det_v2",
        save_mem=True,
        lattice=None,
        remove_ring_fraction=0,
    ):
        """Initializes a BIN simulation for a given detector geometry and collider lattice"""

        self.save_mem = save_mem
        self.design = design
        self.det_geom = getattr(MuC, det_geom)
        if isinstance(remove_ring_fraction, float) or isinstance(
            remove_ring_fraction, int
        ):
            self.remove_ring_fraction = [
                (1 - remove_ring_fraction) / 2,
                (1 + remove_ring_fraction) / 2,
            ]
        elif isinstance(remove_ring_fraction, list) or isinstance(
            remove_ring_fraction, tuple
        ):
            self.remove_ring_fraction = [
                (1 - remove_ring_fraction[0]) / 2,
                (1 + remove_ring_fraction[1]) / 2,
            ]

        else:
            print("No remove_ring_fraction specified. Using all ring.")
            self.remove_ring_fraction = [0, 0]

        self.lattice = lattice
        if isinstance(self.lattice, str):
            try:
                with open(lattice, "rb") as f:
                    self.lattice = pickle.load(f)
            except Exception as errormessage:
                print(errormessage)
                raise ValueError(f"Could not load lattice from file: {lattice}.")
        elif isinstance(self.lattice, dict):
            self.lattice = lattice
            assert (
                "x" in self.lattice.keys()
            ), f"Need x(u) in lattice dictionary. Found keys: {self.lattice.keys()}"
            assert (
                "inv_s" in self.lattice.keys()
            ), "Need inverse function of s(u) in lattice dictionary."
        else:
            print("No lattice specified. Using simplified ring geometry.")
            self.lattice = None

        # Total length of the ring in cm
        self.C = self.lattice["s"](1) if self.lattice is not None else self.design["C"]
        self.beam_lifetime = 1 / self.design["finj"]
        self.n_turns = self.beam_lifetime / (self.C / const.c_LIGHT)
        self.bunchx_in_a_year = (
            self.n_turns
            * self.design["duty_factor"]
            * self.design["finj"]
            * self.design["bunch_multiplicity"]
            * np.pi
            * 1e7  # seconds in a year
        )

        # Detector geometry
        self.comps = list(self.det_geom.facedict.keys())
        self.zending = self.det_geom.zending
        self.rmax = self.det_geom.rmax

        self.N_evals = N_evals

        # Container for all simulations
        self.mustorage_sims = []
        self.beam_cases = col.colls_types_to_beam_cases[self.design["collision_type"]]
        self.nuflavors = [part[0] for part in self.beam_cases]
        for nuflavor, direction in self.beam_cases:
            self.mustorage_sims.append(
                MuStorageRingSimulator(
                    design=design,
                    nuflavor=nuflavor,
                    direction=direction,
                    N_evals=N_evals,
                    beam_lifetime=self.beam_lifetime,
                    preloaded_events=preloaded_events,
                    remove_ring_fraction=self.remove_ring_fraction,
                )
            )

        # Number of simulations
        self.nsims = len(self.mustorage_sims)

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

        # Attempts to perform a muon decay simulation
        self.sims = []
        for mu_sim in self.mustorage_sims:

            # Decay muons
            mu_sim.decay_muons()

            # Place muon along the collider ring
            if self.lattice is not None:
                mu_sim.place_muons_on_lattice(
                    direction="left",  # mu_sim.direction,
                    lattice=self.lattice,
                )
            else:
                mu_sim.place_muons_on_simplified_ring(
                    C=self.C,
                    Lss=self.design["Lss"],
                    direction="left",  # mu_sim.direction,
                )
            # Now onto the detector simulations
            det_sim = DetectorSimulator(mu_sim, self.det_geom, save_mem=self.save_mem)
            det_sim.run()

            # NOTE: This weight is per bunch
            det_sim.w *= (
                self.design["finj"]
                * self.design["duty_factor"]
                * self.design["bunch_multiplicity"]
                * np.pi
                * 1e7  # s in a year
            )

            det_sim.calculate_facecounts()

            if self.save_mem:
                det_sim.clear_mem()

            self.sims.append(det_sim)

        self.total_count = np.sum([np.sum(self.sims[i].w) for i in range(self.nsims)])

        self.name = self.design["name"]
        # print(f"Successfully simulated neutrino event rates within {self.geom.name}:")
        # print(
        # f"{self.name} ({col.acc_colls_dict[self.design['collision_type']]}) at L = {self.design['Lss']:.2f} m."
        # )
        print(f"Total count: {self.total_count:.2e} events;\n")

        self.facecounts = self.get_face_counts()
        self.get_exclusive_rates()

        if show_components:
            self.print_face_counts()

        if show_time:
            self.print_timetable()

        return self

    def reweigh_with_new_polarization(self, new_polarization):

        for sim in self.sims:
            sim.calculate_facecounts_new_pol(new_polarization)
        return self.get_face_counts_new_pol()

    def get_exclusive_rates(self):

        # Append exclusive rates all neutrino species
        self.exclusive_rates = {}
        for s in self.sims:
            s.get_exclusive_rates()
            for (comp, channel), rate in s.get_exclusive_rates().items():
                key = s.nuflavor, comp, channel.replace(s.nuflavor + "_", "")
                if key in self.exclusive_rates.keys():
                    self.exclusive_rates[key] += rate
                else:
                    self.exclusive_rates[key] = rate

        # Exclusive rates from all neutrino species within ECAL and HCAL
        self.exclusive_rates_combined = {}
        for (flavor, comp, channel), rate in self.exclusive_rates.items():
            if comp == "hcal" or comp == "ecal":
                if (flavor, channel) in self.exclusive_rates_combined.keys():
                    self.exclusive_rates_combined[flavor, channel] += rate
                else:
                    self.exclusive_rates_combined[flavor, channel] = rate
        for test_flavor in ["nue", "numu", "nuebar", "numubar"]:
            for (flavor, comp, channel), rate in self.exclusive_rates.items():
                if (test_flavor, channel) in self.exclusive_rates_combined.keys():
                    continue
                else:
                    self.exclusive_rates_combined[test_flavor, channel] = 0

                # else:
                #     print(test_flavor, flavor, channel)
                #     self.exclusive_rates_combined[test_flavor, channel] = 0

    # def NuNuLuminosity(self, particle1="nue", particle2="numu"):
    #     """Computing the luminosity of a neutrino collision."""

    #     cc1, cc2 = self.cco.straight_segment_at_detector(
    #         0, Lss=self.design["Lss"], two=False
    #     )
    #     self.collision = "mu+mu-"
    #     self.ntimes = 6  # should be constant, except during debugging/improvements.
    #     sim1 = None
    #     sim2 = None
    #     sim3 = None
    #     sim4 = None
    #     sims = [sim1, sim2, sim3, sim4]
    #     self.parts = [part[0] for part in col.colls_types_to_part[self.collision]]
    #     nsims = len(self.parts)
    #     geom = importlib.import_module("MuC.detector_geometries.nunulum")
    #     self.zending = geom.zending
    #     self.rmax = geom.rmax
    #     mask = np.array(self.parts) == "lol"
    #     mask[np.array(self.parts) == particle1] = True
    #     mask[np.array(self.parts) == particle2] = True
    #     indices = np.where(mask)[0]
    #     cc1, cc2 = self.cco.straight_segment_at_detector(
    #         geom.zending, Lss=self.design["Lss"], two=True
    #     )
    #     ccs = [cc1, cc1, cc2, cc2]
    #     for i, part in enumerate(self.parts):

    #         if i not in indices:
    #             continue

    #         sims[i] = SimNeutrinos(ccs[i], geom, part, MuC.directions[i])
    #         sims[i].find_intersections()

    #         ########################################
    #         # NOTE: Let's check this.
    #         sims[i].weights *= 2 / nsims
    #         ########################################

    #         sims[i].get_lum_q(self.design)

    #     sims = [sims[i] for i in range(nsims)]
    #     p1 = indices[0]
    #     p2 = indices[1]

    #     return sims, sims[p1].Nnu * sims[p2].Nnu / 4 / np.pi / np.sqrt(
    #         sims[p1].nusdx ** 2 + sims[p2].nusdx ** 2
    #     ) / np.sqrt(sims[p1].nusdy ** 2 + sims[p2].nusdy ** 2)

    def get_face_counts(self):
        """Prints the table of detailed distribution of events in detector components."""
        rates = {}
        for sim in self.sims:
            for comp in self.comps:
                rates[MuC.compsto2[comp], f"{sim.nuflavor}_{sim.direction}"] = (
                    sim.facecounts[comp]
                )
            rates["Total", f"{sim.nuflavor}_{sim.direction}"] = np.sum(
                list(sim.facecounts.values())
            )
        rates["Total", "Total"] = 0
        for comp in self.comps:
            tot = 0
            for sim in self.sims:
                tot += rates[MuC.compsto2[comp], f"{sim.nuflavor}_{sim.direction}"]

            rates[MuC.compsto2[comp], "Total"] = tot
            rates["Total", "Total"] += tot

        return rates

    def get_face_counts_new_pol(self):
        """Prints the table of detailed distribution of events in detector components."""
        rates = {}
        for sim in self.sims:
            for comp in self.comps:
                rates[MuC.compsto2[comp], f"{sim.nuflavor}_{sim.direction}"] = (
                    sim.facecounts_new_pol[comp]
                )
            rates["Total", f"{sim.nuflavor}_{sim.direction}"] = np.sum(
                list(sim.facecounts_new_pol.values())
            )
        rates["Total", "Total"] = 0
        for comp in self.comps:
            tot = 0
            for sim in self.sims:
                tot += rates[MuC.compsto2[comp], f"{sim.nuflavor}_{sim.direction}"]

            rates[MuC.compsto2[comp], "Total"] = tot
            rates["Total", "Total"] += tot

        return rates

    def print_face_counts(self, percentage=0):
        """Prints the table of detailed distribution of events in detector components."""

        table = PrettyTable()
        sample_dict = self.facecounts

        # Step 1: Identify unique rows and columns
        rows = sorted(set(key[0] for key in sample_dict.keys())) + ["Total/bunchx"]
        columns = sorted(set(key[1] for key in sample_dict.keys()))

        # Step 2: Initialize PrettyTable with columns
        table = PrettyTable(["Det component"] + columns)
        # Step 3: Populate the table with data
        for row in rows:
            row_data = [
                (
                    "Total/bunchx"
                    if row == "Total/bunchx"
                    else MuC.comps_short_to_long[row]
                )
            ]
            for column in columns:
                if row == "Total/bunchx":
                    value = sample_dict["Total", column]
                    row_data.append(f"{value/self.bunchx_in_a_year:.2e}")
                else:
                    # Add the corresponding data from the dictionary, if it exists
                    value = sample_dict[row, column]
                    total = sample_dict["Total", column]
                    if percentage and row != "Total":
                        row_data.append(f"{value/total*100:.2g}%")
                    else:
                        row_data.append(f"{value:.2e}")
                        # if row == "Total" and column == "Total":
            table.add_row(row_data)

        return table

    def get_data(self, sec="all", nuflavor="all", genie=0):
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

        if (nuflavor != "all") & (nuflavor not in self.nuflavors):
            raise ValueError(
                f"Particle {nuflavor} is not a valid particle! Choices are: {str(self.nuflavors)}"
            )

        elif nuflavor == "all":
            maskp = np.ones_like(self.sims, dtype=bool)

        else:
            maskp = [(sim.nuflavor == nuflavor) for sim in self.sims]

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
        nuflavor="all",
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
            nuflavor (str): a particle one would want to single out. Can be either nue, nuebar, numu, or numubar.
        """

        if not ax:
            fig, ax = plt.subplots(figsize=fs)

        if orientation == "y-x":
            bs2 = np.linspace(-1 * self.rmax * 2, self.rmax * 2, nbins)
            bs = np.linspace(-1 * self.rmax * 2, self.rmax * 2, nbins)
        else:
            bs = np.linspace(-1 * self.zending * 2, self.zending * 2, nbins)
            bs2 = np.linspace(-1 * self.rmax * 2, self.rmax * 2, nbins)

        x, y, z, w, _, _, _ = self.get_data(sec=sec, nuflavor=nuflavor)

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
        c_dict = {"z-y": [z, y], "z-x": [z, x], "x-y": [x, y], "y-x": [y, x]}
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
            norm=LogNorm(vmin=vmax / 1e6, vmax=vmax),
        )

        MuC.plot_tools.plot_det("det_v2", ax, orientation=orientation, xl=xl, yl=yl)

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
        nuflavor="all",
    ):
        """Wrapper to plot neutrino interaction times.

        Args:
            savefig (str): name of file to save plot to.
            fs (tuple): figsize of plot. Can be None to display on the same plot that is being worked on.
            title (bool): if one wants to display the pre-generated title.
            legend (bool): to display the pre-generated legend.
            sec (str): the component of the detector one wants to single out. Options are the sames as those written in get_data() method description.
            nuflavor (str): a particle one would want to single out. Can be either nue, nuebar, numu, or numubar.
        """

        if fs:
            plt.figure(figsize=fs)

        _, _, _, w, times, _, _ = self.get_data(sec=sec, nuflavor=nuflavor)

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

    def get_GENIE_flux(self, sec, nuflavor, nbins=50):
        """
        Calculate the GENIE flux histogram for a given section and neutrino flavor.

        Parameters:
        sec (int): Detector component to consider.
        nuflavor (str): Neutrino flavor. Can be either nue, nuebar, numu, or numubar.
        nbins (int, optional): The number of bins for the histogram. Default is 50.

        Returns:
        tuple: A tuple containing the bin edges and the histogram values.
        """
        _, _, _, w, _, E, _ = self.get_data(sec=sec, nuflavor=nuflavor, genie=1)
        h = np.histogram(E, weights=w, bins=nbins)
        return h[1], h[0]

    def print_GENIE_flux_to_file(self, sec, nuflavor, nbins=50, filename=None):
        """
            Creates a flux .data file for GENIE simulation of events.

            sec (str): Detector component to consider.
            nuflavor (str): Neutrino flavor. Can be either nue, nuebar, numu, or numubar.
            nbins (int, optional): Number of bins for the histogram. Default is 50.
            filename (str, optional): Name of file to be saved in the fluxes/ folder. If not provided, a default name will be generated.

        Returns:
            None

        Saves:
            A .data file containing the flux information for GENIE simulation.
        """

        if filename:
            fn = f"{filename}"

        else:
            fn = f"fluxes/{self.design['short_name']}_{MuC.compsto2[sec]}_{nuflavor}.data"

        bins, flux = self.get_GENIE_flux(sec, nuflavor, nbins=nbins)
        bin_centers = bins[:-1] + np.diff(bins) / 2
        np.savetxt(fn, np.array([bin_centers, flux]).T)
        print(f"Flux file saved to {fn}")

    def plot_GENIE_flux(self, sec, nuflavor, nbins=50, ax=None):
        """
        Plots the GENIE flux for a given neutrino flavor and sector.

        Parameters:
            sec (str): The detector section for which the flux is to be plotted.
            nuflavor (str): Neutrino flavor. Can be either nue, nuebar, numu, or numubar.
            nbins (int), optional: The number of bins to use for the histogram (default is 50).
            ax (matplotlib.axes.Axes), optional: The axes on which to plot the histogram. If None, a new figure and axes are created (default is None).

        Returns:
            None
        """

        bins, flux = self.get_GENIE_flux(sec=sec, nuflavor=nuflavor, nbins=nbins)
        if ax is None:
            fig, ax = plt.subplots()
        _ = ax.hist(
            bins[:-1] + np.diff(bins) / 2,
            weights=flux,
            bins=bins,
            histtype="step",
            label=nuflavor + "_" + sec,
        )

    def get_GENIE_event_weight(self, comp, p, n_events):
        """Getting the weight of a genie particle from its detector component."""
        try:
            return (
                self.facecounts[comp, p + "_left"] / n_events
            )  # the last factor depends on how many generated events there are in the files. It only supports same n files across detectors.
        except KeyError:
            return self.facecounts[comp, p + "_right"] / n_events  #
        except KeyError:
            return 0

    def load_GENIE_file(self, filename, n_events=1e5):
        """Loads a single GENIE analysis file, adding respective weights based on this simulation's number of BIN interactions.

        Args:
            filename (str): name of the GENIE file to load.
            n_events: number of events that the GENIE file had.

        """

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

        expname = exp
        parts = [MuC.partn_names[part] for part in particles]

        particlenames = ", ".join(parts)
        t = comps.replace(",", ", ")

        print(f"Loading generated data for a {expname} experiment;")
        print(
            f"It includes interactions from {particlenames} within the {t} of the muon detector."
        )

        data = pd.read_csv(filename, sep="\s+", skiprows=5)

        print("Adding weights...")
        try:
            data["w"] = data.apply(
                lambda row: self.get_GENIE_event_weight(
                    row["DComp"], MuC.pdg2names[str(row["IncL"])], n_events=n_events
                ),
                axis=1,
            )
        except KeyError:
            data["w"] = data.apply(
                lambda row: self.get_GENIE_event_weight(
                    row["DComp"], str(row["Particle"]), n_events=n_events
                ),
                axis=1,
            )

        if "Q2" in data.columns:
            data["Q2"] = get_Q2(data["nu_E"], data["E"], data["pz"])

        print("Done!")

        return data

    def load_genie_events(self, filenames, n_events=1e6):
        """
        Load GENIE events from the specified filenames.
        Parameters:
        filenames (str or list of str): The filename(s) from which to load GENIE events.
                                        Can be a single filename or a list of filenames.
        n_events (int, optional): The number of events to load. Default is 1e6.

        Returns: None
            This method sets the following attributes:
                - genie_events: A DataFrame containing the loaded GENIE events.
                - genie_e: A boolean array indicating electron events.
                - genie_mu: A boolean array indicating muon events.
                - genie_tau: A boolean array indicating tau events.
                - genie_nue: A boolean array indicating electron neutrino events.
                - genie_numu: A boolean array indicating muon neutrino events.
                - genie_nutau: A boolean array indicating tau neutrino events.
                - genie_nuebar: A boolean array indicating electron antineutrino events.
                - genie_numubar: A boolean array indicating muon antineutrino events.
                - genie_nutaubar: A boolean array indicating tau antineutrino events.`
        """

        if isinstance(filenames, list):
            data_cases = []
            for filename in filenames:
                data_cases.append(
                    self.load_GENIE_file(f"{filename}", n_events=n_events)
                )
            self.genie_events = pd.concat(data_cases, axis=0)
        else:
            self.genie_events = self.load_GENIE_file(f"{filenames}", n_events=n_events)

        try:
            self.genie_e = np.abs(self.genie_events["OutL"]) == 11
            self.genie_mu = np.abs(self.genie_events["OutL"]) == 13
            self.genie_tau = np.abs(self.genie_events["OutL"]) == 15

            self.genie_nue = self.genie_events["IncL"] == 12
            self.genie_numu = self.genie_events["IncL"] == 14
            self.genie_nutau = self.genie_events["IncL"] == 16

            self.genie_nuebar = self.genie_events["IncL"] == -12
            self.genie_numubar = self.genie_events["IncL"] == -14
            self.genie_nutaubar = self.genie_events["IncL"] == -16
        except KeyError:
            self.genie_e = (self.genie_events["Name"] == "e-") | (
                self.genie_events["Name"] == "e+"
            )
            self.genie_mu = (self.genie_events["Name"] == "mu-") | (
                self.genie_events["Name"] == "mu+"
            )
            self.genie_tau = (self.genie_events["Name"] == "tau-") | (
                self.genie_events["Name"] == "tau+"
            )

            self.genie_nue = self.genie_events["Particle"] == "nue"
            self.genie_numu = self.genie_events["Particle"] == "numu"
            self.genie_nutau = self.genie_events["Particle"] == "nutau"

            self.genie_nuebar = self.genie_events["Particle"] == "nuebar"
            self.genie_numubar = self.genie_events["Particle"] == "numubar"
            self.genie_nutaubar = self.genie_events["Particle"] == "nutaubar"

        return self.genie_events


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


def get_Q2(nu_E, E, pz):
    """Getting Q squared from the generated events."""
    return -1 * ((nu_E - E) ** 2 - (nu_E - pz) ** 2)


def get_flux(x, w, nbins):
    hist1 = np.histogram(
        x, weights=w, bins=nbins, density=False, range=(np.min(x), np.max(x))
    )

    ans0 = hist1[1]
    ans1 = hist1[0]  # /(ans0[1]-ans0[0])
    return ans0, ans1
