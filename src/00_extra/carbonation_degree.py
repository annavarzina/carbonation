# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
#%% DATA

alpha = np.array([0.,0.3,0.55,1.05,1.25,1.55])
time = np.array([0.,2.,4.5,10.,12.,20.]) #days
time_new = np.arange(0.,100.,1.)
#%% PREDICTION 1
#a = 2.1231*(1-exp(-t*0.07188))
def degree_predict(t, c0, c1):
    a = c0*(1 - np.exp(-t*c1))
    return a


c, cov = curve_fit(degree_predict, time, alpha)
print(c)

d_opt = degree_predict(time, c[0],c[1])
d_opt2=degree_predict(time_new, c[0],c[1])
print('R2: ', r2_score(d_opt, alpha))


plt.loglog(time, d_opt, 'r.', label = "Result")
plt.loglog(time_new, d_opt2, 'b-', label = "Fitting")
plt.xlabel('Time (days)')
plt.ylabel('Degree of carbonation (%)')
plt.legend()
plt.show()

plt.plot(time, d_opt, 'r.', label = "Result")
plt.plot(time_new, d_opt2, 'b-', label = "Fitting")
plt.xlabel('Time (days)')
plt.ylabel('Degree of carbonation (%)')
#plt.yscale("log")
plt.legend()
plt.show()
#%%

time_new = np.arange(0.,10.,0.01)
a=degree_predict(time_new, c[0],c[1])
plt.plot(time_new, a, 'b-', label = "Fitting")
plt.xlabel('Time (days)')
plt.ylabel('Degree of carbonation (%)')
#plt.yscale("log")
plt.legend()
plt.show()

time_new = np.arange(0.,1.e-2,1.e-5) 
a=degree_predict(time_new, c[0],c[1])

mass0 = 22.977 
mass = (a-1)*mass0

plt.plot(time_new*24.*3600, (mass - mass0), 'b-')
plt.xlabel('Time (seconds)')
plt.ylabel('Mass increase (mg)')
plt.legend()
plt.show()

print((mass[-1]- mass0)*10e-3)
#%% MASS Predict
def d_predict(t, c0, c1):
    d = c0*(1 - np.exp(-t*c1))
    return d

mass0 = 22.977 *1e-3 #g
dmass = alpha*mass0/100 #delta mass per 1 um2
xmol = dmass/26.812 #constant = rho_CC/C_CC - rho_CH/C_CH

sa  = 29.19e+6
d = xmol*0.0331*1e+15/sa
c, cov = curve_fit(d_predict, time, d)
print(c)

dnew = d_predict(time, c[0],c[1])

f = 0.476*1e-12 #surface area?

plt.plot(time, d, 'r.', label = "Result")
plt.plot(time, dnew, 'b-', label = "Fitting")
plt.xlabel('Time (days)')
plt.ylabel('Depth (um)')
plt.legend()
plt.show()


#%%

time_new = np.arange(0.,1.,0.001) 
dm2 = d_predict(time_new, c[0],c[1])



plt.plot(time_new*24, dm2, 'b-', label = "Fitting")
plt.xlabel('Time (hour)')
plt.ylabel('Depth (um)')
plt.legend()
plt.show()

#%%

time_new = np.arange(0.,100.,1.)
dm2 = d_predict(time_new, c[0],c[1])

plt.plot(time, d, 'r.', label = "Result")
plt.plot(time_new, dm2, 'b-', label = "Fitting")
plt.xlabel('Time (day)')
plt.ylabel('Depth (um)')
plt.legend()
plt.show()

#%%

plt.plot(time, d, 'b.')
#plt.plot(time, d, 'b-')
plt.xlabel('Time (days)')
plt.ylabel('Depth (um)')
plt.legend()
plt.show()


plt.plot(time*24, d, 'b-')
plt.xlabel('Time (hours)')
plt.ylabel('Depth (um)')
plt.xlim(0,16.5)
plt.ylim(0,1)
plt.legend()
plt.show()

