# -*- coding: utf-8 -*-
import matplotlib.pylab as plt
from  mixing import PhreeqcMixing
from kinetics import PhreeqcKinetics

#%% PARAMETERS
database = 'C:\Anaconda2\lib\site-packages\databases\cemdata18.dat'
n = 2000 # time
krate = 5e-5
s = 80#scale factor
fraction = krate * s 
print('Kinetic rate = ' + str(krate))
print('Mixing fraction = ' + str(fraction))

#%% RUN
pm = PhreeqcMixing(n, fraction, database)
pm.run_phreeqc()

pk = PhreeqcKinetics(n, krate, database)
pk.run_phreeqc()

print('Kinetic rate simulation time = ' + str(pk.simulation_time))
print('Mixing fraction simulation time = ' + str(pm.simulation_time))

#%% PLOT
ca_m = []
for i in range(len(pm.selected_output)):
    if pm.selected_output[i][0]==3:
        ca_m.append(pm.selected_output[i][1])

ca_k = []
for i in range(len(pk.selected_output)):
    if pk.selected_output[i][0]==1:
        ca_k.append(pk.selected_output[i][1])
plt.figure()
plt.plot(ca_m, label = "mix")
plt.plot(ca_k, label = "kin")
plt.xlabel('time (s)')
plt.ylabel('Ca (mol/l)')
plt.legend()