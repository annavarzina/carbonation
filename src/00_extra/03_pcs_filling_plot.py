# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 17:52:33 2020

@author: avarzina
"""

import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(precision = 3)

cs = np.array([0.7, 0.5, 0.3])
ps = np.array([0.005, 0.01, 0.05])

poros = np.array([[0.041442492, 0.072463019, 0.31893285],
                  [0.057594873, 0.102263581, 0.410802807],
                  [0.093609196, 0.16709336, 0.566013554]])
poros = np.around(poros, decimals=3)
fig, ax = plt.subplots()
#ax.grid( color="w", linestyle='-', linewidth=3)
#ax.tick_params(which="minor", bottom=False, left=False)
ax.set_xticks(np.arange(len(ps)))
ax.set_yticks(np.arange(len(cs)))
ax.set_xticklabels(ps)
ax.set_yticklabels(cs)

im = ax.imshow(poros,cmap='PuBu_r', alpha = 0.7)
for i in range(len(cs)):
    for j in range(len(ps)):
        text = ax.text(j, i, poros[i, j],
                       ha="center", va="center", color="black", fontsize=12)
        
fig.tight_layout()
plt.xlabel(r"Pore size ($\mu m$)", fontsize=14)
plt.ylabel(r"Crystal size ($\mu m$)", fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.tick_params(labelsize=12)
#cbar.set_label(r'Threshold porosity $\theta$', fontsize=16)
plt.show()

#%%
ie = np.array([0.1, 0.5, 1.0])
poros_ie = np.array([0.00486247, 0.024193441, 0.048108211])

plt.plot(ie, poros_ie, 'x-')
plt.xlabel(r"Interfacial tension energy $\gamma$  ($J/m^2$)", fontsize=14)
plt.ylabel(r"Voxel residual porosity", fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()

#%% POre radius vs SR
line = np.array(['-', '--', '-.', ':', '-'])
ie = np.array([0.1, 0.5, 1.0]) # J/m2
radius = np.arange(1e-3, 1e-1, 1e-4) * 1e-6
def satratio(int_energy, rad):
    mvol = 3.69e-5 #m3/mol
    angle =  np.pi 
    gas_c = 8.314 #J/mol/K
    temp_kelv = 298.3 #K
    omega = np.exp(-1*mvol*np.cos(angle)* 2 * int_energy / (gas_c * temp_kelv * rad ))
    return omega

plt.figure(figsize = (10,5))
for i in range(0, len(ie)):
    plt.plot(radius*1e+6, satratio(ie[i], radius), label = r'$\gamma$=' + str(ie[i]) + r' ($J/m^2$)', ls = line[i])
plt.legend(fontsize=12)
plt.xscale("log")
plt.yscale("log")
plt.tick_params(axis='both', which='major', labelsize=12)
plt.ylabel(r"$\Omega$", fontsize=14)
plt.xlabel(r"Pore size ($\mu m$)", fontsize=14)

#%% Upfilling porosity
cs = np.array([0.7, 0.5, 0.3])
ps = np.array([0.005, 0.01, 0.05])
d = ps*2
theta = np.zeros((3,3))
for i in range(len(cs)):
    for j in range(len(ps)):
        Vcp =( cs[i] + d[j])**3
        Vm = cs[i]**3
        theta[i,j] = (Vcp - Vm)/Vcp
        
print(theta)
poros = np.around(theta, decimals=2)
fig, ax = plt.subplots()
#ax.grid( color="w", linestyle='-', linewidth=3)
#ax.tick_params(which="minor", bottom=False, left=False)
ax.set_xticks(np.arange(len(ps)))
ax.set_yticks(np.arange(len(cs)))
ax.set_xticklabels(ps)
ax.set_yticklabels(cs)

im = ax.imshow(poros,cmap='PuBu_r', alpha = 0.7)
for i in range(len(cs)):
    for j in range(len(ps)):
        text = ax.text(j, i, poros[i, j],
                       ha="center", va="center", color="black", fontsize=12)
        
fig.tight_layout()
plt.xlabel(r"Pore size ($\mu m$)", fontsize=14)
plt.ylabel(r"Crystal size ($\mu m$)", fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.tick_params(labelsize=12)
#cbar.set_label(r'Threshold porosity $\theta$', fontsize=16)
plt.show()
#%%
cs2 = np.array([0.7])
ps2 = np.array([0.005, 0.01, 0.025, 0.05])
d2 = ps2*2
theta2 = np.zeros((len(cs2),len(ps2)))
for i in range(len(cs2)):
    for j in range(len(ps2)):
        Vcp =( cs2[i] + d2[j])**3
        Vm = cs2[i]**3
        theta2[i,j] = (Vcp - Vm)/Vcp
        
print(theta2)
