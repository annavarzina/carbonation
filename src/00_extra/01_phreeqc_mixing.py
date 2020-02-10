# -*- coding: utf-8 -*-

import IPhreeqcPy
import numpy as np
import matplotlib.pylab as plt
import time

def phrqc_string(steps, fraction):
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
steps = 2000
fraction = 0.004
it=time.time() 
ps = phrqc_string(steps, fraction)
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
    if so[i][0]==3:
        ca.append(so[i][1])

plt.plot(ca)