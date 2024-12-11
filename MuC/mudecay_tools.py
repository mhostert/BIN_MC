# from memory_profiler import profile

import numpy as np
import pandas as pd
import vegas as vg
from collections import OrderedDict
from scipy.special import spence

# DarkNews modules
from DarkNews import Cfourvec as Cfv
from DarkNews.MC import run_vegas
from DarkNews.MC import get_samples

from MuC import const


def gauss_pdf(x, x0, sigma):
    if sigma == 0:
        return np.ones_like(x)
    else:
        return np.exp(-((x - x0) ** 2) / (2 * sigma**2)) / sigma / np.sqrt(2 * np.pi)


def Fnue0(x):
    return 6 * x**2 * (1 - x)


def Fnumu0(x):
    return x**2 * (3 - 2 * x)


def Jnue0(x):
    return 6 * x**2 * (1 - x)


def Jnumu0(x):
    return x**2 * (1 - 2 * x)


def L_Spence(x):
    """
    L(x) = - int_0^x log(1 - t)/t dt
    """
    return spence(1 - x)


def k_radiative(x):
    return 2 * L_Spence(x) + 2 * np.pi**2 / 3 + np.log(1 - x) ** 2


def Fnumu1(x):
    return (
        Fnumu0(x) * k_radiative(x)
        + 1 / 6 * (41 - 36 * x + 42 * x**2 - 16 * x**3) * np.log(1 - x)
        + 1 / 12 * x * (82 - 153 * x + 86 * x**2)
    )


def Jnumu1(x):
    term1 = Jnumu0(x) * k_radiative(x)
    term2 = (1 / 6) * (11 - 36 * x + 14 * x**2 - 16 * x**3 - 4 / x) * np.log(1 - x)
    term3 = (1 / 12) * (-8 + 18 * x - 103 * x**2 + 78 * x**3)
    return term1 + term2 + term3


<<<<<<< HEAD
def Fnue1(x):
    term1 = Fnue0(x) * k_radiative(x)
    term2 = (1 - x) * (
        (5 + 8 * x + 8 * x**2) * np.log(1 - x) + (1 / 2) * x * (10 - 19 * x)
    )
    return term1 + term2
=======
def get_particles(design, N_evals=1e5):
    """Monte Carlo generation of muon decays"""
    mdb = MuonDecay(N_mu=design["Nmu_per_bunch"])
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929


def Jnue1(x):
    term1 = Jnue0(x) * k_radiative(x)
    term2 = (1 - x) * (
        (-3 + 12 * x + 8 * x**2 + 4 / x) * np.log(1 - x)
        + (1 / 2) * (8 - 2 * x - 15 * x**2)
    )
    return term1 + term2


def mudecay_matrix_element_sqr(
    x_nu, costheta, muon_polarization, muon_charge, nuflavor, NLO=True
):
    if "nue" in nuflavor:
        return (
            Fnue0(x_nu)
            - muon_charge * muon_polarization * Jnue0(x_nu) * costheta
            - NLO
            * const.alphaQED
            / 2
            / np.pi
            * (Fnue1(x_nu) - Jnue1(x_nu) * muon_charge * muon_polarization * costheta)
        )

    elif "numu" in nuflavor:
        return (
            Fnumu0(x_nu)
            - muon_charge * muon_polarization * Jnumu0(x_nu) * costheta
            - NLO
            * const.alphaQED
            / 2
            / np.pi
            * (Fnumu1(x_nu) - Jnumu1(x_nu) * muon_charge * muon_polarization * costheta)
        )
    else:
        raise ValueError(f"nuflavor {nuflavor} not recognized.")


class NLOmudecay_pol(vg.BatchIntegrand):
    def __init__(self, dim, MC_case):
        self.dim = dim
        self.MC_case = MC_case

        # Find the normalization factor
        self.norm = {}
        self.norm["diff_rate"] = 1
        self.norm["diff_decay_rate"] = 1
        # normalize integrand with an initial throw
        _throw = self.__call__(
            np.random.rand(self.dim, 2000), np.ones((self.dim, 2000))
        )
        for key, val in _throw.items():
            self.norm[key] = np.mean(val)
            # cannot normalize zero integrand
            if self.norm[key] == 0.0:
                print(f"Warning: mean of integrand {key} is zero. Vegas may break.")
                self.norm[key] = 1

    def __call__(self, x, jac):

        MC_case = self.MC_case
        ######################################
        # mu -> e nu nu
        i_var = 0  # 0 is for the momentum pmu
        pmu = MC_case.pmin + (MC_case.pmax - MC_case.pmin) * x[:, i_var]  # GeV
        if np.min(pmu) <= 0:
            raise ValueError("Found a negative or 0 momentum for muons, stopping.")

        i_var += 1  # 1 is for costheta
        costheta = (2.0) * x[:, i_var] - 1.0

        i_var += 1  # 2 is for x_nu = 2 Enu/ m_mu
        x_nu = x[:, i_var]  # 0 to 1

        dgamma = mudecay_matrix_element_sqr(
            x_nu=x_nu,
            costheta=costheta,
            muon_polarization=MC_case.muon_polarization,
            muon_charge=MC_case.muon_charge,
            nuflavor=MC_case.nuflavor,
            NLO=MC_case.NLO,
        )

        # hypercube jacobian (vegas hypercube --> physical limits) transformation
        dgamma *= const.Gf**2 * const.m_mu**5 / 192 / np.pi**3

        dgamma *= 2  # d costheta
        dgamma *= 1  # d x_nu

        ##########################
        # beam energy spread
        drate = dgamma * gauss_pdf(
            pmu, MC_case.beam_p0, MC_case.beam_p0 * MC_case.beam_dpop
        )

        ##############################################
        # return all differential quantities of interest
        self.int_dic = OrderedDict()
        self.int_dic["diff_rate"] = drate
        self.int_dic["diff_decay_rate"] = dgamma

        ##############################################
        # normalization
        self.int_dic["diff_rate"] /= self.norm["diff_rate"]
        self.int_dic["diff_decay_rate"] /= self.norm["diff_decay_rate"]

        return self.int_dic


def get_momenta_from_vegas_samples_x_costheta(vsamples, MC_case):
    """
    Construct the four momenta of all particles in the upscattering+decay process from the
    vegas weights.

    Args:
            vsamples (np.ndarray): integration samples obtained from vegas
                            as hypercube coordinates. Always in the interval [0,1].

            MC_case (DarkNews.MC.MC_events): the main Monte-Carlo class of DarkNews

    Returns:
            dict: each key corresponds to a set of four momenta for a given particle involved,
                    so the values are 2D np.ndarrays with each row a different event and each column a different
                    four momentum component. Contains also the weights.
    """

    four_momenta = {}

    ########################
    # decay
    # energy of projectile
    pP = (MC_case.pmax - MC_case.pmin) * np.array(vsamples[0]) + MC_case.pmin
    EP = np.sqrt(MC_case.Mparent**2 + pP**2)
    masses_decay = {
        "m1": MC_case.Mparent,
        "m2": MC_case.Mdaughter,
        "m3": MC_case.mnu1,
        "m4": MC_case.mnu2,
    }

    # parent particle boost parameters
    boost_parent = {
        "EP_LAB": EP,
        "costP_LAB": np.ones((np.size(EP),)),
        "phiP_LAB": np.zeros((np.size(EP),)),
    }

    decay_samples = {"costheta": 2 * vsamples[1] - 1, "x_nu": vsamples[2]}

    # Pmu Pnu Pnubar Pe
    (
        PmuLAB_decay,
        PnuLAB_decay,
    ) = three_body_decay_x_costheta(decay_samples, boost=boost_parent, **masses_decay)

    four_momenta["P_decay_mu"] = PmuLAB_decay
    four_momenta["P_decay_nu"] = PnuLAB_decay

    return four_momenta


# Three body decay
def three_body_decay_x_costheta(
    samples, boost=False, m1=1, m2=0, m3=0, m4=0, rng=np.random.random
):

    if not samples:
        raise ValueError("Error! No samples were passed to three_body_decay.")
    else:
        # get sample size of the first item
        sample_size = np.shape(list(samples.values())[0])[0]

    # Mandelstam t = m23^2

    x_nu = samples["x_nu"]
    costheta = samples["costheta"]

    M = m1
    Enu = M / 2 * x_nu

    # p1
    PmuCM = Cfv.build_fourvec(
        m1 * np.ones(sample_size),
        np.zeros(sample_size),
        np.ones(sample_size),
        np.zeros(sample_size),
    )

    # p2
    phinu = Cfv.random_generator(sample_size, 0, 2 * np.pi)
    PnuCM = Cfv.build_fourvec(Enu, Enu, costheta, phinu)

    # four-momenta in the LAB frame
    if boost:
        Emu_LAB = boost["EP_LAB"]
        costmu_LAB = boost["costP_LAB"]
        phimu_LAB = boost["phiP_LAB"]

<<<<<<< HEAD
        # Outgoing neutrino
        PmuLAB_decay = Cfv.Tinv(
            PmuCM, -np.sqrt(Emu_LAB**2 - m1**2) / Emu_LAB, costmu_LAB, phimu_LAB
        )
        # Outgoing neutrino
        PnuLAB_decay = Cfv.Tinv(
            PnuCM, -np.sqrt(Emu_LAB**2 - m1**2) / Emu_LAB, costmu_LAB, phimu_LAB
=======
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
        Rpm=0.0,  # balance of + and - helicities
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
            _description_, by default 0.0
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

        if abs(Rpm) > 1:
            raise ValueError("Rpm must be between -1 and 1")

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
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929
        )

        return PmuLAB_decay, PnuLAB_decay

<<<<<<< HEAD
    # four-momenta in the parent particle rest frame
    else:
        return PmuCM, PnuCM
=======
        df_plus = h_plus.get_MC_events()  # gets all events for +1 helicity
        df_minus = h_minus.get_MC_events()  # same but for all -1 helicity

        for w in ["w_flux", "w_decay_rate"]:
            df_plus[w] *= (1 + Rpm) / 2
            df_minus[w] *= (1 - Rpm) / 2

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
        # else:
        #     # for linear
        #     # X,Y,Z coordinates of mu decay
        #     self.xmu = self.rvec_mu[:, 0]
        #     self.ymu = self.rvec_mu[:, 1]
        #     self.zmu = self.rvec_mu[:, 2]

        #     # distance to the detector
        #     self.X_ND = R_ND[0]
        #     self.Y_ND = R_ND[1]
        #     self.Z_ND = R_ND[2]

        #     # distance between plane of Z=Z_ND and the point of mu decay
        #     self.Dzmu = self.Z_ND - self.zmu

        #     # v angles
        #     # cosine of the azimuthal angle wrt to the Z direction (beam straight)
        #     self.ctheta_numu = Cfv.get_cosTheta(self.pnumu)
        #     self.ctheta_nue = Cfv.get_cosTheta(self.pnue)
        #     self.ctheta_e = Cfv.get_cosTheta(self.pe)

        #     self.ttheta_numu = np.sqrt(1.0 - self.ctheta_numu**2) / self.ctheta_numu
        #     self.ttheta_nue = np.sqrt(1.0 - self.ctheta_nue**2) / self.ctheta_nue
        #     self.ttheta_e = np.sqrt(1.0 - self.ctheta_e**2) / self.ctheta_e

        #     # polar angle with Y=0 plane being the ground/plane of the racetrack.
        #     self.phi_numu = np.arctan2(self.pnumu[:, 2], self.pnumu[:, 1])
        #     self.phi_nue = np.arctan2(self.pnue[:, 2], self.pnue[:, 1])
        #     self.phi_e = np.arctan2(self.pe[:, 2], self.pe[:, 1])

        #     # X and Y coordinates of the intersection between pnu and the plane Z=Z_ND
        #     self.Y_intersec_numu = (self.ttheta_numu * self.Dzmu) * np.sin(
        #         self.phi_numu
        #     )
        #     self.X_intersec_numu = (self.ttheta_numu * self.Dzmu) * np.cos(
        #         self.phi_numu
        #     )
        #     self.Y_intersec_nue = (self.ttheta_nue * self.Dzmu) * np.sin(self.phi_nue)
        #     self.X_intersec_nue = (self.ttheta_nue * self.Dzmu) * np.cos(self.phi_nue)

        #     # radius of the ring defined by the intersection of the neutrino momentum and the plane Z=Z_ND
        #     self.dR_numu = (
        #         np.sqrt(1.0 - self.ctheta_numu**2) / self.ctheta_numu * self.Dzmu
        #     )
        #     self.dR_nue = (
        #         np.sqrt(1.0 - self.ctheta_nue**2) / self.ctheta_nue * self.Dzmu
        #     )

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
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929


# class LOmudecay_pol_invmass_method(vg.BatchIntegrand):
#     def __init__(self, dim, MC_case):
#         self.dim = dim
#         self.MC_case = MC_case

#         # Find the normalization factor
#         self.norm = {}
#         self.norm["diff_rate"] = 1
#         self.norm["diff_decay_rate"] = 1
#         # normalize integrand with an initial throw
#         _throw = self.__call__(
#             np.random.rand(self.dim, 20000), np.ones((self.dim, 20000))
#         )
#         for key, val in _throw.items():
#             self.norm[key] = np.mean(val)
#             # cannot normalize zero integrand
#             if self.norm[key] == 0.0:
#                 print(f"Warning: mean of integrand {key} is zero. Vegas may break.")
#                 self.norm[key] = 1

#     def __call__(self, x, jac):

#         MC_case = self.MC_case
#         ######################################
#         # mu -> e nu nu
#         i_var = 0  # 0 is for the momentum pmu
#         pmu = MC_case.pmin + (MC_case.pmax - MC_case.pmin) * x[:, i_var]  # GeV
#         if np.min(pmu) <= 0:
#             raise ValueError("Found a negative or 0 momentum for muons, stopping.")

#         # Emu = np.sqrt(MC_case.Mparent**2 + pmu**2)
#         i_var += 1  # 1 is for the angle thetamu
#         thetamu = (
#             MC_case.theta_min + (MC_case.theta_max - MC_case.theta_min) * x[:, i_var]
#         )
#         i_var += 1  # 2 is for phase space theta_max and theta_min

#         m1 = MC_case.Mparent
#         m2 = MC_case.Mdaughter
#         m3 = MC_case.mnu1
#         m4 = MC_case.mnu2
#         masses = np.array([m1, m2, m3, m4])

#         # limits
#         tmax = phase_space.three_body_tmax(*masses)
#         tmin = phase_space.three_body_tmin(*masses)
#         t = (tmax - tmin) * x[:, i_var] + tmin
#         i_var += 1

#         umax = phase_space.three_body_umax(*masses, t)
#         umin = phase_space.three_body_umin(*masses, t)
#         u = (umax - umin) * x[:, i_var] + umin
#         i_var += 1

#         # v = np.sum(masses**2) - u - t

#         c3 = (2.0) * x[:, i_var] - 1.0
#         i_var += 1
#         phi34 = (2.0 * np.pi) * x[:, i_var]
#         i_var += 1

#         cphi34 = np.cos(phi34)
#         E3CM_decay = (m1**2 + m3**2 - u) / 2.0 / m1
#         E4CM_decay = (m1**2 + m4**2 - t) / 2.0 / m1
#         k3CM = const.kallen_sqrt(m1**2, u, m3**2) / 2.0 / m1
#         k4CM = const.kallen_sqrt(m1**2, t, m4**2) / 2.0 / m1
#         c34 = (t + u - m2**2 - m1**2 + 2 * E3CM_decay * E4CM_decay) / (2 * k3CM * k4CM)
#         #
#         c4 = c3 * c34 - np.sqrt(1.0 - c3 * c3) * np.sqrt(1.0 - c34 * c34) * cphi34

#         x2 = m2 / m1
#         x23 = t / m1**2
#         dgamma = (
#             # const.Gf**2
#             # * const.m_mu
#             # / 256
#             # / np.pi**4
#             # * (1.0 - u / m1**2)
#             # * u
#             # / m1**2
#             const.Gf**2
#             * const.m_mu**5
#             / 256
#             / np.pi**5
#             * (1.0 - x23)
#             * (1 + c4 * MC_case.muon_polarization * MC_case.muon_charge)
#             * (x23 - x2**2)
#         )

#         # hypercube jacobian (vegas hypercube --> physical limits) transformation
#         dgamma *= tmax - tmin
#         dgamma *= umax - umin
#         dgamma *= 2  # c3
#         dgamma *= 2 * np.pi  # phi34

#         ##########################
#         # beam energy spread
#         drate = dgamma * gauss_pdf(
#             pmu, MC_case.beam_p0, MC_case.beam_p0 * MC_case.beam_dpop
#         )

#         # # beam angular spread
#         drate *= gauss_pdf(thetamu, MC_case.beam_theta0, MC_case.beam_dtheta)

#         ##############################################
#         # return all differential quantities of interest
#         self.int_dic = OrderedDict()
#         self.int_dic["diff_rate"] = drate
#         self.int_dic["diff_decay_rate"] = dgamma

#         ##############################################
#         # normalization
#         self.int_dic["diff_rate"] /= self.norm["diff_rate"]
#         self.int_dic["diff_decay_rate"] /= self.norm["diff_decay_rate"]

#         return self.int_dic


# def get_momenta_from_vegas_samples(vsamples, MC_case):
#     """
#     Construct the four momenta of all particles in the upscattering+decay process from the
#     vegas weights.

#     Args:
#             vsamples (np.ndarray): integration samples obtained from vegas
#                             as hypercube coordinates. Always in the interval [0,1].

#             MC_case (DarkNews.MC.MC_events): the main Monte-Carlo class of DarkNews

#     Returns:
#             dict: each key corresponds to a set of four momenta for a given particle involved,
#                     so the values are 2D np.ndarrays with each row a different event and each column a different
#                     four momentum component. Contains also the weights.
#     """

#     four_momenta = {}

#     ########################
#     # decay
#     # energy of projectile
#     pP = (MC_case.pmax - MC_case.pmin) * np.array(vsamples[0]) + MC_case.pmin
#     EP = np.sqrt(MC_case.Mparent**2 + pP**2)
#     thetaP = (MC_case.theta_max - MC_case.theta_min) * np.array(
#         vsamples[1]
#     ) + MC_case.theta_min
#     masses_decay = {
#         "m1": MC_case.Mparent,
#         "m2": MC_case.Mdaughter,
#         "m3": MC_case.mnu1,
#         "m4": MC_case.mnu2,
#     }

#     # parent particle boost parameters
#     boost_parent = {
#         "EP_LAB": EP,
#         "costP_LAB": np.cos(thetaP),
#         "phiP_LAB": Cfv.random_generator(np.size(EP), 0.0, 2 * np.pi),
#     }

#     decay_samples = {
#         "unit_t": vsamples[2],
#         "unit_u": vsamples[3],
#         "unit_c3": vsamples[4],
#         "unit_phi34": vsamples[5],
#     }

#     # Pmu Pnu Pnubar Pe
#     (
#         P1LAB_decay,
#         P2LAB_decay,
#         P3LAB_decay,
#         P4LAB_decay,
#     ) = phase_space.three_body_decay(decay_samples, boost=boost_parent, **masses_decay)

#     four_momenta["P_decay_mu"] = P1LAB_decay
#     four_momenta["P_decay_e"] = P2LAB_decay
#     four_momenta["P_decay_nu_e"] = P3LAB_decay
#     four_momenta["P_decay_nu_mu"] = P4LAB_decay

#     return four_momenta


class MC_events(object):
    def __init__(
        self,
        model="NLOmudecay_pol",
        Mparent=const.m_mu,
        Mdaughter=const.m_e,
        pmin=0.010,
        pmax=10.0,
        include_beamdiv=True,
        nuflavor="mu",
        NLO=True,
        muon_polarization=-1,
        muon_charge=+1,
        mnu1=0,
        mnu2=0,
        beam_p0=3.8,  # GeV
        beam_dpop=0.10,  # sigma / p
        NINT=10,
        NINT_warmup=10,
        NEVAL=1e5,
        NEVAL_warmup=1e4,
        save_mem=True,
    ):

        # set target properties
        self.pmin = pmin
        self.pmax = pmax
        self.Mparent = Mparent
        self.Mdaughter = Mdaughter
        self.mnu1 = mnu1
        self.mnu2 = mnu2
        self.NLO = NLO
        self.muon_polarization = muon_polarization
        self.nuflavor = nuflavor
        self.muon_charge = muon_charge
        self.model = model
        self.include_beamdiv = include_beamdiv

        self.beam_p0 = beam_p0
        self.beam_dpop = beam_dpop

        self.NINT = NINT
        self.NINT_warmup = NINT_warmup
        self.NEVAL = NEVAL
        self.NEVAL_warmup = NEVAL_warmup

        self.save_mem = save_mem

    def get_MC_events(self):
        if self.model == "NLOmudecay_pol":

            DIM = 3  # dim of phase space
            batch_f = NLOmudecay_pol(dim=DIM, MC_case=self)
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

        four_momenta = get_momenta_from_vegas_samples_x_costheta(
            vsamples=samples, MC_case=self
        )

        # SAVE ALL EVENTS AS A PANDAS DATAFRAME
        df_gen = create_df_from_vegas(four_momenta=four_momenta)

        # NOTE: Also saving CoM variables in case we want to reweigh events later
        df_gen["costheta_CM"] = 2 * samples[1] - 1
        df_gen["x_CM"] = samples[2]

        # add weights to it
        df_gen["w_flux"] = weights["diff_rate"] * batch_f.norm["diff_rate"]
        df_gen["w_decay_rate"] = (
            weights["diff_decay_rate"] * batch_f.norm["diff_decay_rate"]
        )

        if self.save_mem:
            del weights
            del four_momenta
            del samples
            del integ
            del batch_f

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
