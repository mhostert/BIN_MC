#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
from nuflux import detgeo, data
from memory_profiler import profile


# In[2]:


mdb = data.get_particles("mutristan_small")


# In[3]:


h = 645
R = 3e5/2/np.pi
Lc = 2*  np.sqrt((-1*h**2 + np.sqrt(h**4 + 4*h**2 * R**2))/2)


# In[4]:


L = Lc*4
cc = detgeo.get_quantities(mdb)
_ = cc.straight_segment_at_detector(L)


# In[5]:


geom = "approximate_muon_detector_2"
particle = 'both'
sim,sim2 = detgeo.SimulateDetector(cc, geom, particle).run()


# In[ ]:




