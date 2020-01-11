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
#plt.loglog(time_new, d_opt2, 'b-', label = "Fitting")
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
#%% depth Predict
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
plt.figure(figsize = (8,4))
plt.plot(time, d, 'b.')
#plt.loglog(time_new, d_opt2, 'b-', label = "Fitting")
plt.xlabel('Time (days)')
plt.ylabel('Degree of carbonation (%)')
plt.show()
#%% depth Predict
'''
def t_predict(d, c0, c1):
    t = -c0*np.log(1.+d/c1)
    return t

mass0 = 22.977 *1e-3 #g
dmass = alpha*mass0/100 #delta mass per 1 um2
xmol = dmass/26.812 #constant = rho_CC/C_CC - rho_CH/C_CH

sa  = 29.19e+6
d = xmol*0.0331*1e+15/sa
c, cov = curve_fit(t_predict, d, time)
print(c)
depths = np.arange(0.,25.,1.) 

tnew = t_predict(depths, c[0],c[1])

f = 0.476*1e-12 #surface area?

plt.plot(d, time, 'r.', label = "Result")
plt.plot(depths, tnew, 'b-', label = "Fitting")
plt.ylabel('Time (days)')
plt.xlabel('Depth (um)')
plt.legend()
plt.show()
'''

#%%

time_new = np.arange(0.,20.,0.01) 
dm2 = d_predict(time_new, c[0],c[1])

plt.figure(figsize=(8,4))
#plt.plot(time, d, 'r.', label = "Result")
plt.plot(time_new, dm2, 'b-', label = "Fitting")
plt.xlabel('Time (day)')
plt.ylabel('Depth (um)')
plt.legend()
plt.show()

depths = np.arange(0.,17.,1.)
k = 0
times = []
for i in range(1, len(time_new)):
    if dm2[i-1] <= depths[k] and dm2[i] > depths[k]:
        times.append(time_new[i])
        k +=1
        if k == len(depths)-1:
            break
#%%
t11 = np.array([0.8026666666666116,
 2.4080000000003117,
 7.224000000006197,
 14.19000000001471,
 23.99399999994149,
 36.20599999982626,
 50.8833333330211,
 68.08333333303288,
 87.66266666701618,
 109.6213333344118,
 134.76200000191304,
 162.62600000283803,
 193.81533333720674,
 227.38400000498777,
 263.4466666715817])


x = np.arange(1, len(t11)+1)

s = t11[14]/times[15]
simtimes11 = [i/s for i in t11]

plt.figure(figsize=(8,4))
plt.plot(time, d, 'r.', label = "Experiment")
#plt.plot(times, np.arange(0.,16.,1.), 'b-', label = "Fitting")
plt.plot(simtimes11, x, 'g-', label = "Simulation")
plt.xlabel('Time (day)')
plt.ylabel('Depth (um)')
plt.legend()
plt.show()

#%%
t12 = np.array([23.827999999943057, 59.48366666627328, 109.03900000105914, 172.1486666698208,
     248.5536666723572, 338.4266666612844, 441.50866664712777, 604.0743333228716, 
     839.3326667041333])
    
    
x = np.arange(1, len(t12)+1)

s12 = t12[8]/t11[8]*s
simtimes12 = [i/s12 for i in t12]

plt.figure(figsize=(8,4))
plt.plot(time, d, 'r.', label = "Experiment")
plt.plot(simtimes11, np.arange(1, len(t11)+1), 'g-', label = "Simulation 11")
plt.plot(simtimes12, np.arange(1, len(t12)+1), 'b-', label = "Simulation 12")
plt.xlabel('Time (day)')
plt.ylabel('Depth (um)')
plt.legend()
plt.show()
#%%
t13 = np.array([
 90.73633333378488,
 334.45533332849647,
 733.3153333492014])
    
    
x = np.arange(1, len(t13)+1)

s13 = t13[2]/t12[2]*s12
simtimes13 = [i/s13 for i in t13]

plt.figure(figsize=(8,4))
plt.plot(time, d, 'r.', label = "Experiment")
plt.plot(simtimes11, np.arange(1, len(t11)+1), 'g-', label = "Simulation 11")
plt.plot(simtimes12, np.arange(1, len(t12)+1), 'b-', label = "Simulation 12")
plt.plot(simtimes13, np.arange(1, len(t13)+1), 'r-', label = "Simulation 13")
plt.xlabel('Time (day)')
plt.ylabel('Depth (um)')
plt.legend()
plt.show()
#%% 


