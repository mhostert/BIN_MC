import numpy as np
import matplotlib.pyplot as plt
import DarkNews as dn
from scipy.stats import binom
from memory_profiler import profile
from nuflux import fluxMC
import sys
import os
cfp = os.path.dirname(os.path.abspath(__file__))
td = os.path.join(cfp,'detector_geometries')
sys.path.append(td)
import useful_data
import gc

def show_events(RA = 1e6, DH = 1e3,DD = [3e2, 0], p_min = 0, p_max = 3e3, p0 = 1e3, dpop = 0.1, Nmu = 1e18, rho = 1.5, mn = 939e6 *1.6e-19/(3e8)**2 * 10**3, bins = 50, plots = True):
    
    mdb = fluxMC.MuonDecay()

    df = mdb.simulate_decays(
                            pmin = p_min, #minimum momentum
                            pmax = p_max,  #maximum momentum
                            beam_p0 = p0, # average momentum (GeV)
                            beam_dpop = dpop,# little beam spread
                            Rpm=0.5, #Fraction of total muons that are plus (or minus?)
                            NINT=10, #for MC integration
                            NINT_warmup=10,
                            NEVAL=1e5,
                            NEVAL_warmup=1e4,
                            )

    _= mdb.propagate_to_detector(
                            Racc = RA,
                            circular = True,
                            det_height = DH, #10 m in cm
                            Ddetector = DD) #arbitrary dimensions for circ detect SHEET at x = 0, y=0, z = Racc ; [0] is radius of detector and [1] is the hollow hole in the middle 
    
    _ = mdb.flux_in_detector(circular = True, NBINS=bins)

    print("If decay uniform, then acceptance = {:.2e}".format(1 / 2 /np.pi * np.arccos(RA / (RA + DD[0]))))

    # Cross sections
    sigmanue, sigmanuebar, sigmanumu, sigmanumubar = useful_data.cs_interp()

    # Approximation for event counts
    Ntargets= rho /mn
    deltaE_nue = mdb.Enue_ND[1] -mdb.Enue_ND[0]
    deltaE_numu = mdb.Enumu_ND[1] - mdb.Enumu_ND[0]
    epsilon_nue = mdb.flux_nue_ND_p/np.sum(mdb.flux_nue_ND_p) / deltaE_nue
    epsilon_numu = mdb.flux_numu_ND_p/np.sum(mdb.flux_numu_ND_p) / deltaE_numu
    dNnue_dE = epsilon_nue * Ntargets*Nmu * mdb.nue_eff_ND * DH * sigmanue(mdb.Enue_ND)
    dNnumu_dE= epsilon_numu * Ntargets *DH*Nmu * mdb.numu_eff_ND * sigmanumubar(mdb.Enumu_ND)
    app_count_nue = round(sum(dNnue_dE) * (mdb.Enue_ND[1] - mdb.Enue_ND[0]))
    app_count_numu = round(sum(dNnumu_dE) * (mdb.Enumu_ND[1] - mdb.Enumu_ND[0]))

    # Actual number event from ray tracing
    na_numu = np.sum(mdb.mask_numu)
    na_nue = np.sum(mdb.mask_nue)
    factors_numu = np.full((1, na_numu),round(Nmu / mdb.sample_size))
    factors_nue = np.full((1, na_nue), round(Nmu / mdb.sample_size))
    earg_numu = - 1 * Ntargets * mdb.heights[:, 1][mdb.mask_numu] * sigmanumubar(mdb.Enumu)[mdb.mask_numu]
    earg_nue = - 1 * Ntargets * mdb.heights[:, 1][mdb.mask_nue] * sigmanue(mdb.Enue)[mdb.mask_nue]
    p_ints_numu = 1 - np.exp(earg_numu)
    p_ints_nue = 1 - np.exp(earg_nue)
    count_nue = binom.rvs(factors_nue, p_ints_nue)
    print("nue interaction count (ray tracing): {:.3g}, vs. approximate: {:.3g}".format(sum(count_nue), app_count_nue))
    count_numu = binom.rvs(factors_numu, p_ints_numu)
    print("numu interaction count (ray tracing): {:.3g}, vs. approximate {:.3g}".format(sum(count_numu), app_count_numu))
    print("\n\n\n")

    # Plots!
    if plots:
        bin_size_mu = np.max(mdb.Enumu[mdb.mask_numu]) / bins
        bin_size_e = np.max(mdb.Enue[mdb.mask_nue]) / bins
        e_bins = np.empty((2, bins))
        for i in range(bins):
            mask_nue = [(mdb.Enue[mdb.mask_nue] > i *bin_size_e) & (mdb.Enue[mdb.mask_nue] <= (i+1) *bin_size_e)]
            e_bins[0, i] = sum(count_nue[mask_nue])
            mask_numu = [(mdb.Enumu[mdb.mask_numu] > i *bin_size_mu) & (mdb.Enumu[mdb.mask_numu] <= (i+1) *bin_size_mu)]
            e_bins[1, i] = sum(count_numu[mask_numu])
        
        plt.figure(1)
        fig = plt.figure(figsize=(15,9))
        plt.title("Interaction counts distribution")
        plt.step([i*bin_size_e for i in range(50)], e_bins[0,:], label =r'$\nu_{e}$', color = 'blue')
        plt.step([i*bin_size_mu for i in range(50)], e_bins[1,:], label =r'$\nu_{\mu}$', color = 'darkorange')
        plt.ylabel("Counts")
        plt.xlabel(r"$E_{\nu}/$GeV")
        plt.legend()

        plt.figure(2)
        fig = plt.figure(figsize=(15,9))
        plt.title(r"Interaction counts distribution ($\nu_{e}$)")
        plt.step([i*bin_size_e for i in range(50)], e_bins[0,:], label =r'$\nu_{e}$', color = 'blue')
        plt.step(mdb.Enue_ND,dNnue_dE*(mdb.Enue_ND[1] - mdb.Enue_ND[0]), label=r'$\nu_{e}$ approximate', linestyle = '--', color = 'blue')
        plt.ylabel("Counts")
        plt.xlabel(r"$E_{\nu}/$GeV")
        plt.legend()

        plt.figure(3)
        fig = plt.figure(figsize=(15,9))
        plt.title(r"Interaction counts distribution ($\nu_{\mu}$)")
        plt.step([i*bin_size_mu for i in range(50)], e_bins[1,:], label =r'$\nu_{\mu}$', color = 'darkorange')
        plt.step(mdb.Enumu_ND,dNnumu_dE*(mdb.Enumu_ND[1] - mdb.Enumu_ND[0]), label=r'$\nu_{\mu}$ approximate', linestyle = '--', color = 'darkorange')
        plt.ylabel("Counts")
        plt.xlabel(r"$E_{\nu}/$GeV")
        plt.legend()

        plt.figure(4)
        fig = plt.figure(figsize=(15,9))
        plt.title(r"Neutrino in detector energy distribution ($\nu$)")
        plt.ylabel("Counts")
        plt.xlabel(r"$E_{\nu}/$GeV")
        plt.hist(mdb.Enue[mdb.mask_nue], weights = np.full((len(mdb.Enue[mdb.mask_nue])), round(Nmu / mdb.sample_size)), histtype ='step', label =r'$\nu_{e}$', color = 'blue')
        plt.hist(mdb.Enumu[mdb.mask_numu], weights = np.full((len(mdb.Enumu[mdb.mask_numu])), round(Nmu / mdb.sample_size)), histtype = 'step', label =r'$\nu_{\mu}$', color = 'darkorange')
        plt.legend()

    return mdb


#@profile
def get_particles(parameters="mutristan_small", N_evals = 1e5):

    param_set = useful_data.parameters[parameters]
    mdb = fluxMC.MuonDecay(N_mu =param_set['Nmu'] )

    mdb.simulate_decays(
                            pmin = param_set["pmin"], #minimum momentum
                            pmax = param_set["pmax"],  #maximum momentum
                            beam_p0 = param_set["beam_p0"], # average momentum (GeV)
                            beam_dpop = param_set["beam_dpop"],# little beam spread
                            beam_dtheta = param_set["beam_dtheta"],
                            Rpm=0.5, #Fraction of total muons that are plus (or minus?)
                            NINT=10, #for MC integration
                            NINT_warmup=10,
                            NEVAL=N_evals,
<<<<<<< HEAD
                            NEVAL_warmup=N_evals/10,
                            luc = True
=======
                            NEVAL_warmup=N_evals/10
>>>>>>> 502e55835199c64fe1919866c5dab403db34d53d
                            )

    mdb.propagate_to_detector(
                            Racc = param_set["Racc"],
                            circular = param_set["circular"],
                            get_int=False) 
    
    
    R = mdb.Racc
    w = np.copy(mdb.w)
    sample_size = mdb.sample_size
    Enumu = np.copy(mdb.Enumu)
    Enue = np.copy(mdb.Enue)
    N_mu=  mdb.N_mu
    pnumu_ar = np.copy(mdb.pnumu)
    pnue_ar = np.copy(mdb.pnue)
    pos_at = np.copy(mdb.pos_at)
    del mdb
    gc.collect()
    return R, w, sample_size, Enumu, Enue, N_mu, pnumu_ar, pnue_ar, pos_at