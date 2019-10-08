# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 14:00:32 2019

@author: avarzina
"""


import numpy as np 
import math
import matplotlib.pylab as plt

cs = 900 #ppm
x = 2*1e-3 #m
D = 1e-10 #m2/s
t = 600*3600 #s

c = cs*(1 - math.erf(x/(2*np.sqrt(D*t))))
print(c)
