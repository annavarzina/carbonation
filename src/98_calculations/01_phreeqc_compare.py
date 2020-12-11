# -*- coding: utf-8 -*-

import IPhreeqcPy
import numpy as np
import matplotlib.pylab as plt
import time

def phrqc_string_kin(steps, k):
    phrqc_input = []    
    
    phrqc_input.append('solution\t1')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    phrqc_input.append('\tCa\t0\n')
    phrqc_input.append('kinetics\t1')
    phrqc_input.append('Portlandite')
    phrqc_input.append('\t-m\t1.e+0')
    phrqc_input.append('\t-m0\t1.e+0')
    phrqc_input.append('\t-tol\t1e-8')
    phrqc_input.append('\t-steps '+steps) # seconds
    phrqc_input.append('\t-runge_kutta 6')

    phrqc_input.append('INCREMENTAL_REACTIONS true')
    phrqc_input.append('RATES')
    phrqc_input.append('\tPortlandite')    
    phrqc_input.append('-start')
    phrqc_input.append('10\tsi_p = SI("Portlandite")') #define saturation index (log(Q/K))for rate equation
    phrqc_input.append('20\tQ_K = 10^si_p') #calculate Q/K
    phrqc_input.append('30\tk = '+ str(k))#1/s, kinetic rate constant 
    phrqc_input.append('40\tr = k * (1-Q_K)')#r is kinetic rate in mol/l.s, equals to k*(1-Q/K)
    phrqc_input.append('50\tmoles = r * TIME')
    phrqc_input.append('60\tSAVE moles')
    phrqc_input.append('-end')
    
    phrqc_input.append('SELECTED_OUTPUT 1')
    phrqc_input.append('\t-reset false')
    phrqc_input.append('\t-time false')
    phrqc_input.append('\t-high_precision true')
    phrqc_input.append('\t-solution true')
    phrqc_input.append('\t-pH false')
    phrqc_input.append('\t-pe false')
    phrqc_input.append('\t-charge_balance false')
    phrqc_input.append('\t-alkalinity false')
    phrqc_input.append('\t-ionic_strength false')
    phrqc_input.append('\t-percent_error false')
    
    phrqc_input.append('USER_PUNCH')
    phrqc_input.append('\t-headings\tCa\tPortlandite')
    phrqc_input.append('\t-start')
    phrqc_input.append('\t10\tpunch\ttot("Ca")')
    phrqc_input.append('\t20\tpunch\ttot("Portlandite")')
    phrqc_input.append('\t30\tpunch')
    phrqc_input.append('\t-end')
    phrqc_input.append('END') 

        
    return '\n'.join(phrqc_input)

def phrqc_string_mix(steps, fraction):
    phrqc_input = []    
    
    phrqc_input.append('solution\t1')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    phrqc_input.append('\tCa\t0\n')
    phrqc_input.append('equilibrium_phases\t1')
    phrqc_input.append('\tPortlandite\t0\t0')
    phrqc_input.append('end\n')
    
    phrqc_input.append('solution\t2')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    phrqc_input.append('\tCa\t0\n')
    phrqc_input.append('equilibrium_phases\t2')
    phrqc_input.append('\tPortlandite\t0\t1')
    phrqc_input.append('end\n')
    
    for i in np.arange(0, steps):
        phrqc_input.append('mix\t2')
        phrqc_input.append('1\t' + str(fraction))
        phrqc_input.append('use solution 1')
        phrqc_input.append('use equilibrium_phase 2')
        phrqc_input.append('save equilibrium_phase 2')
        phrqc_input.append('save solution 2')
        phrqc_input.append('end\n')
        
        phrqc_input.append('SELECTED_OUTPUT 1')
        phrqc_input.append('\t-reset false')
        phrqc_input.append('\t-time false')
        phrqc_input.append('\t-high_precision true')
        phrqc_input.append('\t-solution true')
        phrqc_input.append('\t-pH false')
        phrqc_input.append('\t-pe false')
        phrqc_input.append('\t-charge_balance false')
        phrqc_input.append('\t-alkalinity false')
        phrqc_input.append('\t-ionic_strength false')
        phrqc_input.append('\t-percent_error false')
        
        phrqc_input.append('USER_PUNCH')
        phrqc_input.append('\t-headings\tCa\tPortlandite')
        phrqc_input.append('\t-start')
        phrqc_input.append('\t10\tpunch\ttot("Ca")')
        phrqc_input.append('\t20\tpunch\ttot("Portlandite")')
        phrqc_input.append('\t30\tpunch')
        phrqc_input.append('\t-end')
        phrqc_input.append('END') 
                
        phrqc_input.append('mix\t3')
        phrqc_input.append('2\t1')
        phrqc_input.append('1\t' + str(1.-fraction))
        phrqc_input.append('use solution 2')
        phrqc_input.append('use equilibrium_phase 1')
        phrqc_input.append('save solution 1')
        phrqc_input.append('end\n')
        
    return '\n'.join(phrqc_input)

#%%

n = 1000
steps=''
for i in range(n):
    steps = steps + ' ' +str(1)
# fractions = array([1.   , 0.7  , 0.5  , 0.3  , 0.2  , 0.1  , 
#       0.07 , 0.05 , 0.03 , 0.02 , 0.01 , 0.007,
#       0.005, 0.003, 0.002, 0.001])
#k = array([1.614e-04, 1.612e-04, 1.608e-04, 1.590e-04, 1.564e-04, 1.481e-04,
#       1.412e-04, 1.327e-04, 1.159e-04, 9.989e-05, 7.037e-05, 5.606e-05,
#       4.404e-05, 2.934e-05, 2.070e-05, 1.098e-05])
k = 5e-5
f = 0.004
scale = 80
f = k * scale
print(f)
it=time.time() 
ps_k = phrqc_string_kin(steps, k)
IPhreeqc = IPhreeqcPy.IPhreeqc()
IPhreeqc.LoadDatabase('C:\Anaconda2\lib\site-packages\databases\cemdata18.dat')
IPhreeqc.RunString(ps_k)
so_k = IPhreeqc.GetSelectedOutputArray()
simulation_time = time.time()-it
print(simulation_time)

it=time.time() 
ps_m = phrqc_string_mix(n, f)
IPhreeqc = IPhreeqcPy.IPhreeqc()
IPhreeqc.LoadDatabase('C:\Anaconda2\lib\site-packages\databases\cemdata18.dat')
IPhreeqc.RunString(ps_m)
so_m = IPhreeqc.GetSelectedOutputArray()
simulation_time = time.time()-it
print(simulation_time)
#%%
ca_m = []
for i in range(len(so_m)):
    if so_m[i][0]==3:
        ca_m.append(so_m[i][1])

ca_k = []
for i in range(len(so_k)):
    if so_k[i][0]==1:
        ca_k.append(so_k[i][1])
plt.figure()
plt.plot(ca_m, label = "mix")
plt.plot(ca_k, label = "kin")
plt.xlabel('time (s)')
plt.ylabel('Ca (mol/l)')
plt.legend()
#plt.plot(camix)