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


def Fnue1(x):
    term1 = Fnue0(x) * k_radiative(x)
    term2 = (1 - x) * (
        (5 + 8 * x + 8 * x**2) * np.log(1 - x) + (1 / 2) * x * (10 - 19 * x)
    )
    return term1 + term2


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

        # Outgoing neutrino
        PmuLAB_decay = Cfv.Tinv(
            PmuCM, -np.sqrt(Emu_LAB**2 - m1**2) / Emu_LAB, costmu_LAB, phimu_LAB
        )
        # Outgoing neutrino
        PnuLAB_decay = Cfv.Tinv(
            PnuCM, -np.sqrt(Emu_LAB**2 - m1**2) / Emu_LAB, costmu_LAB, phimu_LAB
        )

        return PmuLAB_decay, PnuLAB_decay

    # four-momenta in the parent particle rest frame
    else:
        return PmuCM, PnuCM


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
