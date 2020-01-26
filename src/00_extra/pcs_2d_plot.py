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
                       ha="center", va="center", color="black")
        
fig.tight_layout()
plt.xlabel(r"Pore size ($\mu m$)")
plt.ylabel(r"Crystal size ($\mu m$)")
ax.figure.colorbar(im, ax=ax)
plt.show()

#%%
ie = np.array([0.1, 0.5, 1.0])
poros_ie = np.array([0.00486247, 0.024193441, 0.048108211])

plt.plot(ie, poros_ie, 'x-')
plt.xlabel(r"Interfacial tension energy  ($J/m^2$)")
plt.ylabel(r"Porosity in a voxel ")
plt.show()