#!/usr/bin/env python
# coding: utf-8

# In[9]:


import numpy as np
from scipy.integrate import odeint
import pylab as pl


# In[10]:


lcap = 6e-9
num  = 50000
R0   = (1 + abs(np.random.randn(num))) * 1e-9
cinf = 55.33
cs   = 5.53e-2
Vm   = 3.29e-5
D    = 3.01e-18
k    = 7.97e-10
N0   = 8.04e21
beta = 4 * np.pi * N0 * R0**3 / (3 * Vm)
Da1  = D/(k * R0)


# In[11]:


delta_C = cinf - cs * np.exp(lcap/R0)
minarg  = np.argmin(R0**2 / delta_C)
t0      = (R0**2 / (Vm * D * delta_C))[minarg]


# In[12]:


def dR_dt(R, t, Da):
    # Identifies particles which are less than 
    ind  = ((R * R0) / 1e-9 < 1)
    nind = np.invert(ind)

    N       = nind.size
    csolute = cinf - np.sum(beta[nind] * R[nind]**3) / N
    
    dRdt       = np.zeros(num)
    dRdt[nind] =   ((R0[minarg] / R0[nind])**2 / delta_C[minarg])                  * (csolute - cs * np.exp(lcap / (R0[nind] * R[nind])))/(Da[nind] + R[nind])
    
    return dRdt


# In[13]:


t = np.append(np.linspace(0, 20, 401)[:-1], np.linspace(20, 1050, 301))


# In[ ]:


Rinit = np.ones(num)
sol1  = odeint(dR_dt, Rinit, t, args = (Da1,), rtol = 1e-8, full_output = 1)


# In[7]:


print(sol1[1]['message'])


# In[10]:


# for time_index, t0 in enumerate(t):
#     print(time_index)
#     pl.hist(R0 * sol1[0][time_index] / 1e-9, 50)
#     pl.xlim([0.8, 7])
#     pl.ylim([0, 4000])    
#     pl.xlabel(r'Size(in nm)')
#     pl.ylabel(r'$N$')
#     pl.title('Time = %3.2f'%(t[time_index]) + r'$\tau$')
#     pl.savefig('images/%04d'%(time_index) + '.png', dpi = 50)
#     pl.clf()


# In[8]:


pl.hist(R0 * sol1[0][0] / 1e-9, 50)
pl.hist(R0 * sol1[0][-1] / 1e-9, 50)
pl.show()


# In[12]:


# surf_energy = np.zeros(t.size)
# for time_index, t0 in enumerate(t):
#     surf_energy[time_index] = ((R0 * sol1[0][time_index] / 1e-9)**2).mean()


# # In[13]:


# pl.plot(t, surf_energy)


# # In[14]:


# sol1[0][-1].min()


# # In[ ]:




