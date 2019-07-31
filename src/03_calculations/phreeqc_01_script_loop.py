# -*- coding: utf-8 -*-

import IPhreeqcPy
import numpy as np
import time

def phrqc_string(pco2, ca):
    n = ca.size
    phrqc_input = []    
    
    for i in np.arange(0, n):
        phrqc_input.append('solution\t' + str(i+1))
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        phrqc_input.append('\tC\t1\tCO2(g)\t-' + str(pco2) )
        phrqc_input.append('\tCa\t' + str(ca[i]) + '\n')
        
    
    phrqc_input.append('SELECTED_OUTPUT')
    phrqc_input.append('\t-reset false')
    phrqc_input.append('\t-time false')
    phrqc_input.append('\t-high_precision true')
    phrqc_input.append('\t-solution false')
    phrqc_input.append('\t-pH false')
    phrqc_input.append('\t-pe false')
    phrqc_input.append('\t-charge_balance false')
    phrqc_input.append('\t-alkalinity false')
    phrqc_input.append('\t-ionic_strength false')
    phrqc_input.append('\t-percent_error false')
    
    phrqc_input.append('USER_PUNCH')
    phrqc_input.append('\t-headings\tC\tCa')
    phrqc_input.append('\t-start')
    phrqc_input.append('\t10\tpunch\ttot("C")')
    phrqc_input.append('\t20\tpunch\ttot("Ca")')
    phrqc_input.append('\t30\tpunch')
    phrqc_input.append('\t-end')
    phrqc_input.append('END')
    return '\n'.join(phrqc_input)
    
co2 = 3.41# np.arange(1,4,0.05)#np.array([3.41, 3.0, 2.0, 1.0])
ca = np.arange(0.01,0.05,0.01)#co2.size)#np.array([0.01, 0.01, 0.01, 0.01])
it=time.time() 
ps = phrqc_string(co2, ca)
IPhreeqc = IPhreeqcPy.IPhreeqc()
#IPhreeqc.LoadDatabase('C:\Anaconda2\lib\site-packages\databases\phreeqc.dat')
IPhreeqc.LoadDatabase('C:\Anaconda2\lib\site-packages\databases\cemdata07.dat')
IPhreeqc.RunString(ps)
#print(IPhreeqc.GetSelectedOutputArray())
so = IPhreeqc.GetSelectedOutputArray()
simulation_time = time.time()-it
print(simulation_time)

print([row[0] for row in so]) #C
print([row[1] for row in so]) #Ca
