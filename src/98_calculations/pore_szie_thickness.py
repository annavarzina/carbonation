# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 13:56:36 2020

@author: avarzina
"""
import numpy as np
import matplotlib.pylab as plt


linetype = np.array(['-', '--', '-.', ':'])

ps_0005  = np.array([0, 0.166666667, 0.655277778, 1.514722222, 2.735277778, 4.386666667])
ps_001  = np.array([0, 0.154444444, 0.603055556, 1.394166667, 2.518055556, 4.111666667])
ps_0025  = np.array([0, 0.123333333, 0.475833333, 1.097777778, 1.982777778, 3.417777778])
ps_005  = np.array([0, 0.0875, 0.327222222, 0.754444444, 1.361666667, 2.575833333])

#1	0.166666667	0.154444444	0.123333333	0.0875
#2	0.655277778	0.603055556	0.475833333	0.327222222
#3	1.514722222	1.394166667	1.097777778	0.754444444
#4	2.735277778	2.518055556	1.982777778	1.361666667
#5	4.386666667	4.111666667	3.417777778	2.575833333

x = np.arange(1, len(ps_0005)+1) - 1
#%%

plt.figure(figsize=(8,4), dpi = 500)
plt.plot(ps_0005, x, label = "pore size 0.005", ls = linetype[0])
plt.plot(ps_001, x, label = "pore size 0.01", ls = linetype[1])
plt.plot(ps_0025, x, label = "pore size 0.025", ls = linetype[2])
plt.plot(ps_005, x, label = "pore size 0.05", ls = linetype[3])
plt.xlabel('Time (h)', fontsize=14)
plt.ylabel(r'Thickness ($\mu m$)', fontsize=14)
plt.legend(fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()