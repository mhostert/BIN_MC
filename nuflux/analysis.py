import numpy as np
import random
from scipy.stats import truncexpon
from DarkNews import Cfourvec as Cfv
from DarkNews import const

from nuflux import analysis
from nuflux import fluxMC as MC


class neutrino_flux(object):
    def __init__(
        self,
        ZBEAMEND=250e2,  # cm
        ZBEAMEXIT=0,  # cm
        BEAM_HALFSIZE=12,  # cm
        include_beamdiv=True,
        truncate_exp=True,
        R_ND=[0, 0, 50e2],  # cm
    ):
        self.ZBEAMEND = ZBEAMEND  # cm
        self.ZBEAMEXIT = ZBEAMEXIT  # cm
        self.BEAM_HALFSIZE = BEAM_HALFSIZE  # cm
        self.include_beamdiv = include_beamdiv
        self.truncate_exp = truncate_exp

    #######################################
    # Use vegas to simulate mu decays
    def simulate_decays(
        self,
        Rpm=0.5,
        model="LOmudecay_unpol",
        pmin=0.200,
        pmax=10.0,
        beam_dpop=3.8,  # GeV
        beam_sigmaE=0.10,  # %
        beam_theta0=0.0,  # rad
        beam_dtheta=0.005,  # rad
        NINT=10,
        NINT_warmup=10,
        NEVAL=1e5,
        NEVAL_warmup=1e4,
    ):
        Mparent = const.m_mu
        Mdaughter = const.m_e
        # muon helicities h #
        h_plus = MC.MC_events(
            model=model,
            Mparent=Mparent,
            Mdaughter=Mdaughter,
            helicity=+1,
            ptmax=0.001,
            pmin=pmin,
            pmax=pmax,
            beam_dpop=beam_dpop,  # GeV
            beam_sigmaE=beam_sigmaE,  # %
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
            ptmax=0.001,
            pmin=pmin,
            pmax=pmax,
            beam_dpop=beam_dpop,  # GeV
            beam_sigmaE=beam_sigmaE,
            beam_theta0=beam_theta0,
            beam_dtheta=beam_dtheta,
            NINT=NINT,
            NINT_warmup=NINT_warmup,
            NEVAL=NEVAL,
            NEVAL_warmup=NEVAL_warmup,
        )

        bag_plus = h_plus.get_MC_events()
        bag_minus = h_minus.get_MC_events()

        self.dic = MC.Combine_MC_output(
            bag_plus,
            bag_minus,
            Ifactor_case0=Rpm,
            Ifactor_case1=(1 - Rpm),
            case0_flag=0,
            case1_flag=1,
        )
        ##########
        # get observables
        true = analysis.observables(
            self.dic, smearing=True, include_pt=True, truncate_exp=True
        )
        return true

    def decay_position(self):
        Emu = self.pmu[:, 0]
        pmu = np.sqrt(Emu**2 - const.m_mu**2)
        gammabetak = pmu / const.m_mu
        lmudecay = gammabetak * const.c_LIGHT * const.Tmuon * 100  # cm
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

    #######################################
    # propagate the flux to the ND
    def propagate_to_ND(
        self,
        smearing=True,
    ):
        self.pmu = self.df_gen["P_decay_mu"].values
        self.pe = self.df_gen["P_decay_e"].values
        self.pnue = self.df_gen["P_decay_nu_e"].values
        self.pnumu = self.df_gen["P_decay_nu_mu"].values

        self.sample_size = np.size(self.pmu[:, 0])

        ############################################################################
        # generate position of the mu decay along straight
        self.decay_position()

        #############################################################################
        # geometry

        # X,Y,Z coordinates of mu decay
        self.xmu = self.rvec_mu[:, 0]
        self.ymu = self.rvec_mu[:, 1]
        self.zmu = self.rvec_mu[:, 2]

        # distance to the detector
        self.X_ND = self.R_ND[0]
        self.Y_ND = self.R_ND[1]
        self.Z_ND = self.R_ND[2]

        # distance between plane of Z=Z_ND and the point of mu decay
        self.Dzmu = self.Z_ND - self.zmu

        ##########
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

        ########
        # radius of the ring defined by the intersection of the neutrino momentum and the plane Z=Z_ND
        self.dR_numu = (
            np.sqrt(1.0 - self.ctheta_numu**2) / self.ctheta_numu * self.Dzmu
        )
        self.dR_nue = np.sqrt(1.0 - self.ctheta_nue**2) / self.ctheta_nue * self.Dzmu


# ##############################################
# # SMEAR ENERGY AND X,Y POSITION MEASUREMENTS
# def det_smear_E_and_pos(E,r):
# 	sigma_E =  np.sqrt((1.935e-2*np.sqrt(E))*(1.935e-2*np.sqrt(E)) + (5.954e-3*E)*(5.954e-3*E))
# 	sigma_xy = np.sqrt(0.2148*0.2148 + (0.3663/np.sqrt(E))*(0.3663/np.sqrt(E)))

# 	E_r = Cfv.random_normal(E,sigma_E)
# 	x_r = Cfv.random_normal(r[:,0], sigma_xy)
# 	y_r = Cfv.random_normal(r[:,1], sigma_xy)
# 	r_r = np.append(x_r,[y_r,r[:,2]]).reshape(3,np.size(sigma_E)).T
# 	return E_r,r_r
