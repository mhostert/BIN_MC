import numpy as np
import vegas as vg
from collections import OrderedDict

from DarkNews import const
from DarkNews import phase_space
from DarkNews import Cfourvec as Cfv


def gauss_pdf(x, x0, sigma):
    if sigma == 0:
        return np.ones_like(x)
    else:
        return np.exp(-((x - x0) ** 2) / (2 * sigma**2)) / sigma / np.sqrt(2 * np.pi)


class LOmudecay_unpol(vg.BatchIntegrand):
    def __init__(self, dim, MC_case):
        self.dim = dim
        self.MC_case = MC_case

        # Find the normalization factor
        self.norm = {}
        self.norm["diff_rate"] = 1
        self.norm["diff_decay_rate"] = 1
        # normalize integrand with an initial throw
        _throw = self.__call__(
            np.random.rand(self.dim, 20000), np.ones((self.dim, 20000))
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

        # Emu = np.sqrt(MC_case.Mparent**2 + pmu**2)
        i_var += 1  # 1 is for the angle thetamu
        thetamu = (
            MC_case.theta_min + (MC_case.theta_max - MC_case.theta_min) * x[:, i_var]
        )
        i_var += 1  # 2 is for phase space theta_max and theta_min

        m1 = MC_case.Mparent
        m2 = MC_case.Mdaughter
        m3 = MC_case.mnu1
        m4 = MC_case.mnu2
        masses = np.array([m1, m2, m3, m4])

        # limits
        tmax = phase_space.three_body_tmax(*masses)
        tmin = phase_space.three_body_tmin(*masses)
        t = (tmax - tmin) * x[:, i_var] + tmin
        i_var += 1

        umax = phase_space.three_body_umax(*masses, t)
        umin = phase_space.three_body_umin(*masses, t)
        u = (umax - umin) * x[:, i_var] + umin
        i_var += 1

        # v = np.sum(masses**2) - u - t

        c3 = (2.0) * x[:, i_var] - 1.0
        i_var += 1
        phi34 = (2.0 * np.pi) * x[:, i_var]
        i_var += 1

        cphi34 = np.cos(phi34)
        E3CM_decay = (m1**2 + m3**2 - u) / 2.0 / m1
        E4CM_decay = (m1**2 + m4**2 - t) / 2.0 / m1
        k3CM = const.kallen_sqrt(m1**2, u, m3**2) / 2.0 / m1
        k4CM = const.kallen_sqrt(m1**2, t, m4**2) / 2.0 / m1
        c34 = (t + u - m2**2 - m1**2 + 2 * E3CM_decay * E4CM_decay) / (2 * k3CM * k4CM)
        #
        c4 = c3 * c34 - np.sqrt(1.0 - c3 * c3) * np.sqrt(1.0 - c34 * c34) * cphi34

        x2 = m2 / m1
        x23 = t / m1**2
        dgamma = (
            # const.Gf**2
            # * const.m_mu
            # / 256
            # / np.pi**4
            # * (1.0 - u / m1**2)
            # * u
            # / m1**2
            const.Gf**2
            * const.m_mu**5
            / 256
            / np.pi**5
            * (1.0 - x23)
            * (1 + c4 * MC_case.helicity)
            * (x23 - x2**2)
        )

        # hypercube jacobian (vegas hypercube --> physical limits) transformation
        dgamma *= tmax - tmin
        dgamma *= umax - umin
        dgamma *= 2  # c3
        dgamma *= 2 * np.pi  # phi34

        ##########################
        # beam energy spread
        drate = dgamma * gauss_pdf(
            pmu, MC_case.beam_p0, MC_case.beam_p0 * MC_case.beam_dpop
        )

        # # beam angular spread
        drate *= gauss_pdf(thetamu, MC_case.beam_theta0, MC_case.beam_dtheta)

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


def get_momenta_from_vegas_samples(vsamples, MC_case):
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
    thetaP = (MC_case.theta_max - MC_case.theta_min) * np.array(
        vsamples[1]
    ) + MC_case.theta_min
    masses_decay = {
        "m1": MC_case.Mparent,
        "m2": MC_case.Mdaughter,
        "m3": MC_case.mnu1,
        "m4": MC_case.mnu2,
    }

    # parent particle boost parameters
    boost_parent = {
        "EP_LAB": EP,
        "costP_LAB": np.cos(thetaP),
        "phiP_LAB": Cfv.random_generator(np.size(EP), 0.0, 2 * np.pi),
    }

    decay_samples = {
        "unit_t": vsamples[2],
        "unit_u": vsamples[3],
        "unit_c3": vsamples[4],
        "unit_phi34": vsamples[5],
    }

    # Pmu Pnu Pnubar Pe
    (
        P1LAB_decay,
        P2LAB_decay,
        P3LAB_decay,
        P4LAB_decay,
    ) = phase_space.three_body_decay(decay_samples, boost=boost_parent, **masses_decay)

    four_momenta["P_decay_mu"] = P1LAB_decay
    four_momenta["P_decay_e"] = P2LAB_decay
    four_momenta["P_decay_nu_e"] = P3LAB_decay
    four_momenta["P_decay_nu_mu"] = P4LAB_decay

    return four_momenta
