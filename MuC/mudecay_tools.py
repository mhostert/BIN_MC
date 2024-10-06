# from memory_profiler import profile

import numpy as np
import random
import pandas as pd
import vegas as vg
import gc

from scipy.stats import truncexpon
from particle import literals as lp

# DarkNews modules
from DarkNews import const
from DarkNews import Cfourvec as Cfv
from DarkNews.MC import run_vegas
from DarkNews.MC import get_samples

# MuC modules
from MuC import integrands


def get_flux(x, w, nbins):
    hist1 = np.histogram(
        x, weights=w, bins=nbins, density=False, range=(np.min(x), np.max(x))
    )

    ans0 = hist1[1][:nbins]
    ans1 = hist1[0]  # /(ans0[1]-ans0[0])
    comb1 = ans1 / np.sum(ans1)
    return ans0, comb1


def get_particles(design, N_evals=1e5):
    """Monte Carlo generation of muon decays"""
    mdb = MuonDecay(N_mu=design["Nmu"])

    std_tolerance = 4
    mdb.simulate_decays(
        pmin=(1 - std_tolerance * design["beam_dpop"])
        * design["beam_p0"],  # minimum momentum
        pmax=(1 + std_tolerance * design["beam_dpop"])
        * design["beam_p0"],  # maximum momentum
        beam_p0=design["beam_p0"],  # average momentum (GeV)
        beam_dpop=design["beam_dpop"],  # little beam spread
        beam_theta0=0,
        beam_dtheta=design["beam_dtheta"],
        theta_max=design["beam_dtheta"] * std_tolerance,
        theta_min=0,
        Rpm=design[
            "muon_polarization"
        ],  # Fraction of total muons that are plus (or minus?)
        NINT=10,  # for MC integration
        NINT_warmup=10,
        NEVAL=N_evals,
        NEVAL_warmup=N_evals / 10,
    )

    mdb.propagate_to_detector(C=design["C"], circular=design["circular"], get_int=False)

    # summary quantities for det simulation
    C = mdb.C
    w = mdb.w
    sample_size = mdb.sample_size
    N_mu = mdb.N_mu
    pnumu = mdb.pnumu
    pnue = mdb.pnue
    pos = mdb.rvec_mu
    Emu = mdb.pmu[:, 0]

    # free memory
    del mdb
    gc.collect()

    return C, w, sample_size, N_mu, pnumu, pnue, pos, design["name"], Emu


class MuonDecay(object):
    def __init__(self, N_mu=1e18, save_mem=True):
        self.N_mu = N_mu  # number of muons to be decayed
        self.save_mem = save_mem

    #######################################
    # Use vegas to simulate mu decays
    def simulate_decays(
        self,
        Rpm=0.5,  # fraction that are plus helicity
        model="LOmudecay_unpol",
        pmin=0.010,  # GeV
        pmax=10.0,  # GeV
        theta_min=0.0,  # rad
        theta_max=1.0,  # rad
        beam_p0=3.8,  # GeV
        beam_dpop=0.001,  #
        beam_theta0=0.0,  # rad
        beam_dtheta=0.005,  # rad
        NINT=10,
        NINT_warmup=10,
        NEVAL=1e6,
        NEVAL_warmup=1e4,
        Nmu=2e19,
    ):
        """simulate_decays _summary_

        Parameters
        ----------
        Rpm : float, optional
            _description_, by default 0.5
        model : str, optional
            _description_, by default "LOmudecay_unpol"
        pmin : float, optional
            _description_, by default 0.010
        pmax : float, optional
            _description_, by default 10.0
        theta_min : float, optional
            _description_, by default 0.0
        NINT_warmup : int, optional
            _description_, by default 10
        NEVAL : _type_, optional
            _description_, by default 1e5
        NEVAL_warmup : _type_, optional
            _description_, by default 1e4

        Returns
        -------
        _type_
            _description_
        """

        Mparent = const.m_mu  # mass of the muon
        Mdaughter = const.m_e  # mass of the electron
        # muon helicities h #
        h_plus = MC_events(
            model=model,
            Mparent=Mparent,
            Mdaughter=Mdaughter,
            helicity=+1,  # For all h_plus muon events
            theta_min=theta_min,  # rad
            theta_max=theta_max,  # rad
            pmin=pmin,
            pmax=pmax,
            beam_p0=beam_p0,  # GeV
            beam_dpop=beam_dpop,
            beam_theta0=beam_theta0,
            beam_dtheta=beam_dtheta,
            NINT=NINT,
            NINT_warmup=NINT_warmup,
            NEVAL=NEVAL,
            NEVAL_warmup=NEVAL_warmup,
        )

        h_minus = MC_events(
            model=model,
            Mparent=Mparent,
            Mdaughter=Mdaughter,
            helicity=-1,
            theta_min=theta_min,  # rad
            theta_max=theta_max,  # rad
            pmin=pmin,
            pmax=pmax,
            beam_p0=beam_p0,  # GeV
            beam_dpop=beam_dpop,
            beam_theta0=beam_theta0,
            beam_dtheta=beam_dtheta,
            NINT=NINT,
            NINT_warmup=NINT_warmup,
            NEVAL=NEVAL,
            NEVAL_warmup=NEVAL_warmup,
        )

        df_plus = h_plus.get_MC_events()  # gets all events for +1 helicity
        df_minus = h_minus.get_MC_events()  # same but for all -1 helicity

        for w in ["w_flux", "w_decay_rate"]:
            df_plus[w] *= Rpm
            df_minus[w] *= 1 - Rpm

        df_gen = pd.concat([df_plus, df_minus], axis=0).reset_index(
            drop=True
        )  # adds all helicities together
        if self.save_mem:
            del df_plus
            del df_minus
            del h_plus
            del h_minus
            gc.collect()

        self.pmu = df_gen["P_decay_mu"].to_numpy()  # all muon decaying momenta
        self.pe = df_gen["P_decay_e"].to_numpy()  # all emitted electrons momenta
        self.pnue = df_gen[
            "P_decay_nu_e"
        ].to_numpy()  # all emitted electron neutrinos momenta
        self.pnumu = df_gen[
            "P_decay_nu_mu"
        ].to_numpy()  # all emitted muonic neutrinos momenta
        self.w = df_gen["w_flux"].to_numpy()  # flux weights

        if self.save_mem:
            del df_gen
        else:
            return df_gen

    #######################################
    # auxiliary func to decay particle along its trajectory
    def decay_position(self):
        Emu = self.pmu[:, 0]
        pmu = np.sqrt(Emu**2 - const.m_mu**2)  # gamma *m *v * c
        gammabetak = pmu / const.m_mu  # gamma * v/c
        lmudecay = gammabetak * const.c_LIGHT * (lp.mu_minus.lifetime * 1e-9)  # cm
        cmu1 = Cfv.get_cosTheta(
            self.pmu
        )  # costheta of the muon, not the decay particles, so always very close to 1

        # print(cmu1)

        if self.truncate_exp:
            self.rmu = truncexpon.rvs(
                (self.ZBEAMEND - self.ZBEAMEXIT) / lmudecay / cmu1,
                0.0,
                lmudecay,
                size=self.sample_size,
            )
        else:
            # Truncated exponential -- decay from Z=0 to Z=Z_ND
            self.rmu = np.random.exponential(scale=lmudecay, size=self.sample_size)
        # coordinates of decay of Klong in (cm,cm,cm)
        if self.include_beamdiv:
            self.rvec_mu = (Cfv.get_3direction(self.pmu).T * self.rmu).T
            # random x,y position within BEAM_HALFSIZE
            self.XPOS = np.linspace(
                -self.BEAM_HALFSIZE, self.BEAM_HALFSIZE, self.sample_size
            )
            self.rvec_mu[:, 0] += random.choices(
                self.XPOS, weights=Cfv.fit_x_position(self.XPOS), k=self.sample_size
            )

            self.YPOS = np.linspace(
                -self.BEAM_HALFSIZE, self.BEAM_HALFSIZE, self.sample_size
            )
            self.rvec_mu[:, 1] += random.choices(
                self.YPOS, weights=Cfv.fit_y_position(self.YPOS), k=self.sample_size
            )

            self.rvec_mu[:, 2] += self.ZBEAMEXIT
        else:
            self.rvec_mu = Cfv.put_in_z_axis(self.rmu)

        return self.rmu, self.rvec_mu

    # @profile
    def propagate_to_detector(
        self,
        ZBEAMEND=250e2,  # cm
        ZBEAMEXIT=0,  # cm
        BEAM_HALFSIZE=12,  # cm
        R_ND=[0, 0, 50e2],  # cm
        smearing=True,
        include_beamdiv=False,
        truncate_exp=False,
        circular=False,
        Racc=1e6,
        C=3e5,
        Ddetector=[3e2, 20],
        det_height=10e2,
        get_int=True,
    ):
        self.ZBEAMEND = ZBEAMEND  # cm
        self.ZBEAMEXIT = ZBEAMEXIT  # cm
        self.BEAM_HALFSIZE = BEAM_HALFSIZE  # cm
        self.R_ND = R_ND  # cm
        if circular:
            self.R_ND = [0, 0, Racc]
        self.Racc = Racc  # cm
        self.C = C  # cm
        self.Rdet = Ddetector[0]
        self.Rhole = Ddetector[1]
        self.det_height = det_height

        self.include_beamdiv = include_beamdiv
        self.truncate_exp = truncate_exp

        self.sample_size = np.size(self.pmu[:, 0])

        # Energies
        self.Emu = self.pmu[:, 0]
        self.Enumu = self.pnumu[:, 0]
        self.Enue = self.pnue[:, 0]
        self.Ee = self.pe[:, 0]

        # generate position of the mu decay along straight
        self.decay_position()
        # geometry - both for linear and circular

        # for circular
        if circular:
            self.delta = self.rvec_mu[:, 2] / self.Racc % 2 * np.pi

            # counterclockwise - new momentum
            self.pe[:, :] = Cfv.rotationx(self.pe, -1 * (self.delta + np.pi / 2))
            self.pnumu[:, :] = Cfv.rotationx(self.pnumu, -1 * (self.delta + np.pi / 2))
            self.pnue[:, :] = Cfv.rotationx(self.pnue, -1 * (self.delta + np.pi / 2))

            # translate coordinate axis - assign new coordinate positions based on deltay
            self.pos_at = np.zeros((3, self.sample_size))
            self.pos_at[0, :] = self.rvec_mu[:, 0]
            self.pos_at[1, :] = self.Racc * np.sin(self.delta)
            self.pos_at[2, :] = self.Racc * np.cos(self.delta)

            if get_int:
                self.momenta = [self.pe_ar, self.pnumu_ar, self.pnue_ar]

                """
                  This array will have all the coordinates for the intersection points, xyz:
                    first for y = 0,
                    then z= Racc + Rdet,
                    then y = det_height,
                    then x = Rdet,
                    then x = - Rdet;
                    then it will do it for all three particles

                """
                self.int_points = np.empty([self.sample_size, 5, 3, 3])
                # intersection point momentum - y=0 plane$ FOR ACCEPTANCE
                for i, p in enumerate(self.momenta):
                    self.int_points[:, 0, :, i] = self.pos_at.T - np.divide(
                        np.multiply(
                            self.pos_at[1, :].reshape((self.sample_size, 1)), p[:, 1:4]
                        ),
                        p[:, 2].reshape(self.sample_size, 1),
                    )
                    self.int_points[:, 1, :, i] = self.pos_at.T - np.divide(
                        np.multiply(
                            (
                                self.pos_at[2, :].reshape((self.sample_size, 1))
                                - np.full((self.sample_size, 1), self.Racc + self.Rdet)
                            ),
                            p[:, 1:4],
                        ),
                        p[:, 3].reshape(self.sample_size, 1),
                    )
                    self.int_points[:, 2, :, i] = self.pos_at.T - np.divide(
                        np.multiply(
                            (
                                self.pos_at[1, :].reshape((self.sample_size, 1))
                                - np.full((self.sample_size, 1), self.det_height)
                            ),
                            p[:, 1:4],
                        ),
                        p[:, 2].reshape(self.sample_size, 1),
                    )
                    self.int_points[:, 3, :, i] = self.pos_at.T - np.divide(
                        np.multiply(
                            (
                                self.pos_at[0, :].reshape((self.sample_size, 1))
                                - np.full((self.sample_size, 1), self.Rdet)
                            ),
                            p[:, 1:4],
                        ),
                        p[:, 1].reshape(self.sample_size, 1),
                    )
                    self.int_points[:, 4, :, i] = self.pos_at.T - np.divide(
                        np.multiply(
                            (
                                self.pos_at[0, :].reshape((self.sample_size, 1))
                                - np.full((self.sample_size, 1), -1 * self.Rdet)
                            ),
                            p[:, 1:4],
                        ),
                        p[:, 1].reshape(self.sample_size, 1),
                    )

            """planes for detector: Racc is the radius of accelerator, Rdet is the radius of the detector, Rhole is the radius of the hole.
                                    z_top = Racc + Rdet
                                    z_bottom = Racc - Rdet
                                    x_front = Rdet
                                    x_back = - Rdet
                                    y_front = 0
                                    y_back = det_height

                                    mask: accepted neutrinos already
                                    x cases:
                                    (1) = most common: y_front first, then z_top (if int. point with z_top has a 0 < y < det_height and - Rdet < x < Rdet)
                                    (2) = y_front, then y_back (if int. point with y_back has  - Rdet < x < Rdet)
                                    (3) = y_front, then x_front (if int. point with x_front has 0 < y < det_height)
                                    (4) = y_front, then x_back

                                    compute all lengths, then assert all lengths are below sqrt(8 * Rdet**2 + det_height**2)
            """

        else:

            # for linear
            # X,Y,Z coordinates of mu decay
            self.xmu = self.rvec_mu[:, 0]
            self.ymu = self.rvec_mu[:, 1]
            self.zmu = self.rvec_mu[:, 2]

            # distance to the detector
            self.X_ND = R_ND[0]
            self.Y_ND = R_ND[1]
            self.Z_ND = R_ND[2]

            # distance between plane of Z=Z_ND and the point of mu decay
            self.Dzmu = self.Z_ND - self.zmu

            # v angles
            # cosine of the azimuthal angle wrt to the Z direction (beam straight)
            self.ctheta_numu = Cfv.get_cosTheta(self.pnumu)
            self.ctheta_nue = Cfv.get_cosTheta(self.pnue)
            self.ctheta_e = Cfv.get_cosTheta(self.pe)

            self.ttheta_numu = np.sqrt(1.0 - self.ctheta_numu**2) / self.ctheta_numu
            self.ttheta_nue = np.sqrt(1.0 - self.ctheta_nue**2) / self.ctheta_nue
            self.ttheta_e = np.sqrt(1.0 - self.ctheta_e**2) / self.ctheta_e

            # polar angle with Y=0 plane being the ground/plane of the racetrack.
            self.phi_numu = np.arctan2(self.pnumu[:, 2], self.pnumu[:, 1])
            self.phi_nue = np.arctan2(self.pnue[:, 2], self.pnue[:, 1])
            self.phi_e = np.arctan2(self.pe[:, 2], self.pe[:, 1])

            # X and Y coordinates of the intersection between pnu and the plane Z=Z_ND
            self.Y_intersec_numu = (self.ttheta_numu * self.Dzmu) * np.sin(
                self.phi_numu
            )
            self.X_intersec_numu = (self.ttheta_numu * self.Dzmu) * np.cos(
                self.phi_numu
            )
            self.Y_intersec_nue = (self.ttheta_nue * self.Dzmu) * np.sin(self.phi_nue)
            self.X_intersec_nue = (self.ttheta_nue * self.Dzmu) * np.cos(self.phi_nue)

            # radius of the ring defined by the intersection of the neutrino momentum and the plane Z=Z_ND
            self.dR_numu = (
                np.sqrt(1.0 - self.ctheta_numu**2) / self.ctheta_numu * self.Dzmu
            )
            self.dR_nue = (
                np.sqrt(1.0 - self.ctheta_nue**2) / self.ctheta_nue * self.Dzmu
            )

    def flux_in_detector(
        self, DIM_ND=[3e2, 3e2, 3e2], NBINS=100, acceptance=False, circular=False
    ):

        if circular:

            # distances - we assume that the detector is at (0, 0, Racc)
            self.d_numu = np.sqrt(
                (self.int_points[:, 0, 0, 1] - 0) ** 2
                + (self.int_points[:, 0, 2, 1] - self.Racc) ** 2
            )
            self.d_nue = np.sqrt(
                (self.int_points[:, 0, 0, 2] - 0) ** 2
                + (self.int_points[:, 0, 2, 2] - self.Racc) ** 2
            )

            # masks
            self.mask_numu = (self.d_numu < self.Rdet) & (self.d_numu > self.Rhole)

            self.mask_nue = (self.d_nue < self.Rdet) & (self.d_nue > self.Rhole)

            # additional mask: if delta < pi, emitted particle cannot reach detector
            not_accepted = []
            for i, d in enumerate(self.delta):
                if d < np.pi:
                    not_accepted.append(i)

            for j, i in enumerate(not_accepted):
                self.mask_numu[i] = 0
                self.mask_nue[i] = 0

            # distance done within detector; Racc is radius of accelerator, Rdet is radius of detector, Rhole is radius of hole.
            """planes for detector: z_top = Racc + Rdet
                                    z_bottom = Racc - Rdet
                                    x_front = Rdet
                                    x_back = Rdet
            """
            self.heights = np.empty(
                [self.sample_size, 3]
            )  # first is e, then numu, then nue
            self.cases = np.empty([self.sample_size, 3])
            for i in range(self.sample_size):
                for j in range(3):

                    if (
                        (self.int_points[i, 1, 1, j] < self.det_height)
                        & (self.int_points[i, 1, 1, j] > 0)
                        & (self.int_points[i, 1, 0, j] < self.Rdet)
                        & (self.int_points[i, 1, 0, j] > -1 * self.Rdet)
                    ):
                        self.cases[i, j] = 1

                    elif (
                        (self.int_points[i, 2, 0, j] < self.Rdet)
                        & (self.int_points[i, 2, 0, j] > -1 * self.Rdet)
                        & (self.int_points[i, 2, 2, j] < self.Rdet + self.Racc)
                        & (self.int_points[i, 2, 2, j] > self.Racc - self.Rdet)
                    ):
                        self.cases[i, j] = 2

                    elif (
                        (self.int_points[i, 3, 1, j] < self.det_height)
                        & (self.int_points[i, 3, 1, j] > 0)
                        & (self.int_points[i, 3, 2, j] < self.Rdet + self.Racc)
                        & (self.int_points[i, 3, 2, j] > self.Racc - self.Rdet)
                    ):
                        self.cases[i, j] = 3

                    elif (
                        (self.int_points[i, 4, 1, j] < self.det_height)
                        & (self.int_points[i, 4, 1, j] > 0)
                        & (self.int_points[i, 4, 2, j] < self.Rdet + self.Racc)
                        & (self.int_points[i, 4, 2, j] > self.Racc - self.Rdet)
                    ):
                        self.cases[i, j] = 4

                    else:
                        self.cases[i, j] = 0

                    # case = self.cases[i,j]
                    self.heights[i, j] = D3distance(
                        self.int_points[i, int(self.cases[i, j]), :, j],
                        self.int_points[i, 0, :, j],
                    )

            for p in self.heights[:, 1][self.mask_numu]:
                assert p < np.sqrt(8 * self.Rdet**2 + self.det_height**2)
            for p in self.heights[:, 2][self.mask_nue]:
                assert p < np.sqrt(8 * self.Rdet**2 + self.det_height**2)

        else:
            self.mask_numu = np.array(
                (np.abs(self.Y_intersec_numu - self.R_ND[1]) < DIM_ND[1])
                & (np.abs(self.X_intersec_numu - self.R_ND[0]) < DIM_ND[0])
            )
            self.mask_nue = np.array(
                (np.abs(self.Y_intersec_nue - self.R_ND[1]) < DIM_ND[1])
                & (np.abs(self.X_intersec_nue - self.R_ND[0]) < DIM_ND[0])
            )

        self.nue_eff_ND = np.sum(self.w[self.mask_nue]) / np.sum(self.w)
        self.numu_eff_ND = np.sum(self.w[self.mask_numu]) / np.sum(self.w)

        if acceptance:
            return self.nue_eff_ND, self.numu_eff_ND

        print(
            "Detector acceptance: {:.2e} for nue and {:.2e} for numu".format(
                self.nue_eff_ND, self.numu_eff_ND
            )
        )

        if self.nue_eff_ND > 0 and self.numu_eff_ND > 0:
            self.wnue_ND = self.w[self.mask_nue]
            self.wnumu_ND = self.w[self.mask_numu]
            self.Enue_ND, self.flux_nue_ND_p = get_flux(
                self.Enue[self.mask_nue], self.wnue_ND, NBINS
            )
            self.Enumu_ND, self.flux_numu_ND_p = get_flux(
                self.Enumu[self.mask_numu], self.wnumu_ND, NBINS
            )

            # area of detector
            if circular:
                self.area = np.pi * (self.Rdet**2 - self.Rhole**2)
            else:
                self.area = DIM_ND[1] * DIM_ND[2]

            self.flux_nue_ND = (
                self.flux_nue_ND_p
                / np.sum(self.flux_nue_ND_p)
                * self.N_mu
                / self.area
                / (self.Enue_ND[1] - self.Enue_ND[0])
                * self.nue_eff_ND
            )
            self.flux_numu_ND = (
                self.flux_numu_ND_p
                / np.sum(self.flux_numu_ND_p)
                * self.N_mu
                / self.area
                / (self.Enumu_ND[1] - self.Enumu_ND[0])
                * self.numu_eff_ND
            )

            return self.flux_nue_ND_p, self.flux_numu_ND_p
        else:
            print("No flux in the detector")
            return 0, 0


def D3distance(position, intersection):
    return np.sqrt(
        (position[0] - intersection[0]) ** 2
        + (position[1] - intersection[1]) ** 2
        + (position[2] - intersection[2]) ** 2
    )


class MC_events(object):
    def __init__(
        self,
        model="LOmudecay_unpol",
        Mparent=const.m_mu,
        Mdaughter=const.m_e,
        pmin=0.010,
        pmax=10.0,
        theta_min=0.0,  # rad
        theta_max=1.0,  # rad
        include_beamdiv=True,
        helicity=-1,
        mnu1=0,
        mnu2=0,
        beam_p0=3.8,  # GeV
        beam_dpop=0.10,  # sigma / p
        beam_theta0=0.0,  # rad #what is this again
        beam_dtheta=0.005,  # rad
        NINT=10,
        NINT_warmup=10,
        NEVAL=1e5,
        NEVAL_warmup=1e4,
        save_mem=True,
    ):

        # set target properties
        self.pmin = pmin
        self.pmax = pmax
        self.theta_min = theta_min
        self.theta_max = theta_max
        self.Mparent = Mparent
        self.Mdaughter = Mdaughter
        self.mnu1 = mnu1
        self.mnu2 = mnu2
        self.helicity = helicity
        self.model = model
        self.include_beamdiv = include_beamdiv

        self.beam_p0 = beam_p0
        self.beam_dpop = beam_dpop
        self.beam_theta0 = beam_theta0
        self.beam_dtheta = beam_dtheta

        self.NINT = NINT
        self.NINT_warmup = NINT_warmup
        self.NEVAL = NEVAL
        self.NEVAL_warmup = NEVAL_warmup

        self.save_mem = save_mem

    def get_MC_events(self):
        if (
            self.model == "LOmudecay_unpol"
        ):  # need to look into this decay; try different ones?
            # BATCH SAMPLE INTEGRAN OF INTEREST
            DIM = 6  # dim of phase space
            batch_f = integrands.LOmudecay_unpol(dim=DIM, MC_case=self)
            integ = vg.Integrator(DIM * [[0.0, 1.0]])

            _ = run_vegas(
                batch_f,
                integ,
                NINT_warmup=self.NINT_warmup,
                NEVAL_warmup=self.NEVAL_warmup,
                NINT=self.NINT,
                NEVAL=self.NEVAL,
            )

        #########################
        # Get the int variables and weights
        samples, weights = get_samples(integ, batch_f)

        four_momenta = integrands.get_momenta_from_vegas_samples(
            vsamples=samples, MC_case=self
        )

        # SAVE ALL EVENTS AS A PANDAS DATAFRAME
        df_gen = create_df_from_vegas(four_momenta=four_momenta)

        # add weights to it
        df_gen["w_flux"] = weights["diff_rate"]
        df_gen["w_decay_rate"] = weights["diff_decay_rate"]

        if self.save_mem:
            del weights
            del four_momenta
            del samples
            del integ
            del batch_f
            gc.collect()

        return df_gen


def create_df_from_vegas(four_momenta, sparse=0):
    particles = list(four_momenta.keys())

    if sparse >= 2:  # keep visible, and parent momenta -- Enu to be added later
        particles = [
            x
            for x in particles
            if "target" not in x and "recoils" not in x and "daughter" not in x
        ]
    if sparse == 4:  # keep only visible momenta
        particles = [x for x in particles if "w_decay" not in x]

    columns_index = pd.MultiIndex.from_product([particles, ["0", "1", "2", "3"]])

    df_gen = pd.DataFrame(
        np.hstack([four_momenta[p] for p in particles]), columns=columns_index
    )

    # differential weights
    for column in df_gen:
        if "w_" in str(column):
            df_gen[column, ""] = df_gen[column]
    return df_gen
