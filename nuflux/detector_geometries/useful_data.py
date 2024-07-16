import numpy as np

from scipy import interpolate

from helpers import material, subs, comp, unif


def cs_interp():
    '''Returns necessary cross neutrino cross sections from xsecs'''
    log10E,sigmaeo,sigmamuo,_,sigmaebaro,sigmamubaro,_ = np.genfromtxt('xsecs/XCC.dat',unpack=True)
    exs = 10**(log10E)
    exs=np.append(exs, 1e4)
    
    #adding values linearily up to 10 TeV
    sigmae = np.append(sigmaeo, 7.046434082801659e-1)
    sigmamu = np.append(sigmamuo, 7.046434082801659e-1)
    sigmaebar = np.append(sigmaebaro, 3.758631195493489e-1)
    sigmamubar = np.append(sigmamubaro, 3.758631195493489e-1)
    
    #generating interpolators
    sigmanue = interpolate.interp1d(exs,sigmae*exs*1e-38,bounds_error=False,fill_value=0.0)
    sigmanuebar = interpolate.interp1d(exs,sigmaebar*exs*1e-38,bounds_error=False,fill_value=0.0)
    sigmanumu = interpolate.interp1d(exs,sigmamu*exs*1e-38,bounds_error=False,fill_value=0.0)
    sigmanumubar = interpolate.interp1d(exs,sigmamubar*exs*1e-38,bounds_error=False,fill_value=0.0)
    
    return sigmanue, sigmanuebar, sigmanumu, sigmanumubar

#density in g/cm**3; atomic mass in g/mol
Si = subs(2.329, 28.0855, 28)
WSi2 = subs(9.3, 240.01, 240)
Fe = subs(7.874, 55.845,56)
Al = subs(2.7,26.981539, 27)
W = subs(19.3, 183.84, 184)
Cu = subs(8.96, 63.546, 64)
PS = subs(1.05, 104.1, 104)

#from CLICdet paper
hcal_CLICdet = comp([[Fe, 20/26.5],[Al, 0.7/26.5],[Cu, 0.1 / 26.5],[PS, 3/26.5]]) 
ecal_CLICdet = comp([[W, 1.9 / 5.05],[Cu, 2.3/5.05],[Si, 0.5/5.05]])

#from online
EARTH_DENSITY = unif(5.51).N

#Parameter sets; from muTristan (Hamada et al., 2022), and MuCoL (Intl Muon Coll Collab, 2023)
mutristan_small = {"name": "muTRISTAN (s)","beam_p0": 1e3, "pmax": 3e3, "pmin":0,"C": 3e5, "beam_dpop": 0.1, "beam_dtheta": 1.12e-4,"circular": False, "Nmu": 3.6e-9 / (1.6e-19) * (1 - 1/np.e) *40* 365.25 * 24* 3600 * 50*0.7}
mutristan_large = {"name": "muTRISTAN (L)","beam_p0": 3e3, "pmax": 9e3, "pmin":0,"C": 9e5, "beam_dpop": 0.1, "beam_dtheta": 7e-5,"circular": False, "Nmu": 3.6e-9 / (1.6e-19) * (1 - 1/np.e) *40* 365.25 * 24* 3600 * 50*0.7}
mucol_s1 = {"name": "MuCoL (s1)","beam_p0": 2.5e3, "pmax": 7.5e3, "pmin":0,"C": 4.5e5, "beam_dpop": 0.1, "beam_dtheta": 3.52e-4,"circular": False, "Nmu": 4.9e9  * 4.5e3* 1.2e7}
mucol_s2 = {"name": "MuCoL (s2)","beam_p0": 5e3, "pmax": 15e3, "pmin":0,"C": 10e5, "beam_dpop": 0.1, "beam_dtheta": 5.88e-4,"circular": False, "Nmu": 1.8e9 * 10e3 * 1.2e7}
parameters = {"mutristan_small": mutristan_small, "mutristan_large": mutristan_large, 'mucol_s1' : mucol_s1, 'mucol_s2': mucol_s2}
