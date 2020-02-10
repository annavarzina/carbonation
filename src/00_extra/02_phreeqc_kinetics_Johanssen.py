# -*- coding: utf-8 -*-

import IPhreeqcPy
import numpy as np
import matplotlib.pylab as plt
import time

def phrqc_string(steps, k):
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
    phrqc_input.append('\t-runge_kutta 1')

    phrqc_input.append('INCREMENTAL_REACTIONS true')
    phrqc_input.append('RATES')
    phrqc_input.append('\tPortlandite')    
    phrqc_input.append('-start')
    phrqc_input.append('10\tsi_p = SI("Portlandite")') #define saturation index (log(Q/K))for rate equation
    phrqc_input.append('20\tQ_K = 10^si_p') #calculate Q/K
    phrqc_input.append('30\tk = '+ str(k))#1/s, kinetic rate constant 
    phrqc_input.append('40\tr = k * (1-Q_K)')#r is kinetic rate in mol/l.s, equals to k*(1-Q/K)
    phrqc_input.append('41\tarea = 1')# cm2  
    phrqc_input.append('50\tmoles = r * TIME * area')
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
#%%
n = 2000
steps=''
for i in range(0,n,100):
    steps = steps + ' ' +str(n)

k = 1.*1e+5# mmol/l/s/cm2
k = k/1e-3# mol/l/s/cm2

it=time.time() 
ps = phrqc_string(steps, k)
IPhreeqc = IPhreeqcPy.IPhreeqc()
#IPhreeqc.LoadDatabase('C:\Anaconda2\lib\site-packages\databases\phreeqc.dat')
IPhreeqc.LoadDatabase('C:\Anaconda2\lib\site-packages\databases\cemdata07.dat')
IPhreeqc.RunString(ps)
#print(IPhreeqc.GetSelectedOutputArray())
so = IPhreeqc.GetSelectedOutputArray()
simulation_time = time.time()-it
print(simulation_time)
#%%
ca = []
for i in range(len(so)):
    if so[i][0]==1:
        ca.append(so[i][1])

plt.plot(ca)
#plt.plot(camix)