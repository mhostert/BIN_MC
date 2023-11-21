import numpy as np
import random
import pandas as pd
from scipy.stats import truncexpon
from particle import literals as lp


from DarkNews import const
from DarkNews import Cfourvec as Cfv

from nuflux import MC


def get_flux(x, w, nbins):
    TMIN = np.min(x)
    TMAX = np.max(x)
    hist1 = np.histogram(x, weights=w, bins=nbins, density=False, range=(TMIN, TMAX))

    ans0 = hist1[1][:nbins]
    ans1 = hist1[0]  # /(ans0[1]-ans0[0])
    comb1 = ans1 / np.sum(ans1)
    return ans0, comb1


class MuonDecay(object):
    def __init__(self, N_mu=1e18):
        self.N_mu = N_mu  # number of muons to be decayed

    #######################################
    # Use vegas to simulate mu decays
    def simulate_decays(
        self,
        Rpm=0.5,  # fraction that are plus helicity
        model="LOmudecay_unpol",
        pmin=0.010,  # GeV
        pmax=10.0,  # GeV
        tmin=0.0,  # rad
        tmax=1.0,  # rad
        beam_p0=3.8,  # GeV
        beam_dpop=0.10,  #
        beam_theta0=0.0,  # rad
        beam_dtheta=0.005,  # rad
        NINT=10,
        NINT_warmup=10,
        NEVAL=1e5,
        NEVAL_warmup=1e4,
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
        tmin : float, optional
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
        h_plus = MC.MC_events(
            model=model,
            Mparent=Mparent,
            Mdaughter=Mdaughter,
            helicity=+1,  # For all h_plus muon events
            tmin=tmin,  # rad
            tmax=tmax,  # rad
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

        h_minus = MC.MC_events(
            model=model,
            Mparent=Mparent,
            Mdaughter=Mdaughter,
            helicity=-1,
            tmin=tmin,  # rad
            tmax=tmax,  # rad
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

        self.df_gen = pd.concat([df_plus, df_minus], axis=0).reset_index(
            drop=True
        )  # adds all helicities together

        return self.df_gen

    #######################################
    # auxiliary func to decay particle along its trajectory
    def decay_position(self):
        Emu = self.pmu[:, 0]
        pmu = np.sqrt(Emu**2 - const.m_mu**2)
        gammabetak = pmu / const.m_mu
        lmudecay = gammabetak * const.c_LIGHT * (lp.mu_minus.lifetime * 1e-9)  # cm
        cmu1 = Cfv.get_cosTheta(self.pmu)

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

    def propagate_to_detector(
        self,
        ZBEAMEND=250e2,  # cm
        ZBEAMEXIT=0,  # cm
        BEAM_HALFSIZE=12,  # cm
        R_ND=[0, 0, 50e2],  # cm
        smearing=True,
        include_beamdiv=False,
        truncate_exp=True,
    ):
        self.ZBEAMEND = ZBEAMEND  # cm
        self.ZBEAMEXIT = ZBEAMEXIT  # cm
        self.BEAM_HALFSIZE = BEAM_HALFSIZE  # cm
        self.R_ND = R_ND  # cm

        self.include_beamdiv = include_beamdiv
        self.truncate_exp = truncate_exp

        self.pmu = self.df_gen["P_decay_mu"].to_numpy()  # all muon decaying momenta
        self.pe = self.df_gen["P_decay_e"].to_numpy()  # all emitted electrons momenta
        self.pnue = self.df_gen[
            "P_decay_nu_e"
        ].to_numpy()  # all emitted electron neutrinos momenta
        self.pnumu = self.df_gen[
            "P_decay_nu_mu"
        ].to_numpy()  # all emitted muonic neutrinos momenta
        self.w = self.df_gen["w_flux"].to_numpy()  # flux of emitted W; momenta of W?

        self.sample_size = np.size(self.pmu[:, 0])

        # Energies
        self.Emu = self.pmu[:, 0]
        self.Enumu = self.pnumu[:, 0]
        self.Enue = self.pnue[:, 0]
        self.Ee = self.pe[:, 0]

        # generate position of the mu decay along straight
        self.decay_position()

        # geometry

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
        self.Y_intersec_numu = (self.ttheta_numu * self.Dzmu) * np.sin(self.phi_numu)
        self.X_intersec_numu = (self.ttheta_numu * self.Dzmu) * np.cos(self.phi_numu)
        self.Y_intersec_nue = (self.ttheta_nue * self.Dzmu) * np.sin(self.phi_nue)
        self.X_intersec_nue = (self.ttheta_nue * self.Dzmu) * np.cos(self.phi_nue)

        # radius of the ring defined by the intersection of the neutrino momentum and the plane Z=Z_ND
        self.dR_numu = (
            np.sqrt(1.0 - self.ctheta_numu**2) / self.ctheta_numu * self.Dzmu
        )
        self.dR_nue = np.sqrt(1.0 - self.ctheta_nue**2) / self.ctheta_nue * self.Dzmu

    def flux_in_detector(self, DIM_ND=[3e2, 3e2, 3e2], NBINS=100, acceptance=False):
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

        if acceptance == True:
            return self.nue_eff_ND, self.numu_eff_ND

        print(
            "Detector acceptance:",
            self.nue_eff_ND,
            " for nue, and",
            self.numu_eff_ND,
            "for numu.",
        )

        if self.nue_eff_ND > 0 and self.numu_eff_ND > 0:
            self.wnue_ND = self.w[self.mask_nue]
            self.wnumu_ND = self.w[self.mask_numu]
            self.Enue_ND, self.flux_nue_ND = get_flux(
                self.Enue[self.mask_nue], self.wnue_ND, NBINS
            )
            self.Enumu_ND, self.flux_numu_ND = get_flux(
                self.Enumu[self.mask_numu], self.wnumu_ND, NBINS
            )

            # area of detector
            self.area = DIM_ND[1] * DIM_ND[2]
            self.flux_nue_ND = (
                self.flux_nue_ND
                / np.sum(self.flux_nue_ND)
                * self.N_mu
                / self.area
                / (self.Enue_ND[1] - self.Enue_ND[0])
                * self.nue_eff_ND
            )
            self.flux_numu_ND = (
                self.flux_numu_ND
                / np.sum(self.flux_numu_ND)
                * self.N_mu
                / self.area
                / (self.Enumu_ND[1] - self.Enumu_ND[0])
                * self.numu_eff_ND
            )

            return self.flux_nue_ND, self.flux_numu_ND
        else:
            print("No flux in the detector")
            return 0, 0
