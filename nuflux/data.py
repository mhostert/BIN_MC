import sys, os
cfp = os.path.dirname(os.path.abspath(__file__))
td = os.path.join(cfp,'detector_geometries')
sys.path.append(td)

import gc
import numpy as np
import matplotlib.pyplot as plt
from memory_profiler import profile

import DarkNews as dn
import useful_data
from nuflux import fluxMC


def get_particles(parameters="mutristan_small", N_evals = 1e5):
    '''Monte Carlo generation of muon decays'''
    param_set = useful_data.parameters[parameters]
    mdb = fluxMC.MuonDecay(N_mu =param_set['Nmu'] )
    
    mdb.simulate_decays(pmin = param_set["pmin"], #minimum momentum
                        pmax = param_set["pmax"],  #maximum momentum
                        beam_p0 = param_set["beam_p0"], # average momentum (GeV)
                        beam_dpop = param_set["beam_dpop"],# little beam spread
                        beam_dtheta = param_set["beam_dtheta"],
                        Rpm=0.5, #Fraction of total muons that are plus (or minus?)
                        NINT=10, #for MC integration
                        NINT_warmup=10,
                        NEVAL=N_evals,
                        NEVAL_warmup=N_evals/10,
                        luc = True)

    mdb.propagate_to_detector(C = param_set["C"],
                              circular = param_set["circular"],
                              get_int=False) 
    
    #summary quantities for det simulation
    C = mdb.C
    w = mdb.w
    sample_size = mdb.sample_size
    N_mu=  mdb.N_mu
    pnumu = mdb.pnumu
    pnue = mdb.pnue
    pos = mdb.rvec_mu
    Emu = mdb.pmu[:,0]
    
    #free memory
    del mdb
    gc.collect()
    
    return C, w, sample_size, N_mu, pnumu, pnue, pos, param_set['name'], Emu