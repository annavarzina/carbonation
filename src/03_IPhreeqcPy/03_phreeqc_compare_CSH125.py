# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pylab as plt
from  mixing import PhreeqcMixing
from kinetics import PhreeqcKinetics

class PhreeqcMixingCSH(PhreeqcMixing):
    def __init__(self, n, fraction, csh, database):
        self.phase = csh['name']
        self.csh = csh
        self.steps = n
        self.fraction = fraction
        self.database = database
        self.phrqc_input = []
        self.selected_output = []
        self.phrqc_string = ''
        self.simulation_time = 0
    
    def generate_phrqc_string(self):
        self.phases()
        self.solution_1()
        self.solution_2()        
        for i in np.arange(0, self.steps):
            self.mix_2()
            self.selected_output_1()
            self.user_punch()
            self.mix_3()            
        self.phrqc_string = '\n'.join(self.phrqc_input)
        
    def phases(self):
        phrqc_input = []   
        # CSH stochiometry
        s = self.csh['stochiometry']
        h = s['H+']
        h2o =  s['H2O'] +  s['H+'] - s['Ca']
        sign1 = '+'
        if h < 0: 
            sign1 = '-'
            h *= -1
        sign2 = '+'
        if h2o < 0: 
            sign2 = '-'
            h2o *= -1
        # input 
        phrqc_input.append('PHASES')  
        phrqc_input.append(self.phase)  
        phrqc_input.append('\t(CaO)' + str(s['Ca']) +'(SiO2)'+ str(s['Si']) + \
                           '(H2O)' + str(s['H2O']) + ' ' + sign1 + ' ' + str(h) + 'H+ = ' + \
                           str(s['Ca']) + 'Ca+2 + ' + str(s['Si']) + 'SiO2 '  + sign2 +\
                           ' ' + str(h2o) + ' H2O') 
        #phrqc_input.append('\t-Vm\t' + str(csh['vm']) ) 
        phrqc_input.append('\t-log_K\t' + str(self.csh['log_k']) + '\n')
        self.phrqc_input +=  phrqc_input
        
    
    def user_punch(self):
        phrqc_input = []         
        phrqc_input.append('USER_PUNCH')
        phrqc_input.append('\t-headings\tCa\t' + self.phase)
        phrqc_input.append('\t-start')
        phrqc_input.append('\t10\tpunch\ttot("Ca")')
        phrqc_input.append('\t15\tpunch\ttot("Si")')
        phrqc_input.append('\t20\tpunch\ttot("' + self.phase + '")')
        phrqc_input.append('\t30\tpunch')
        phrqc_input.append('\t-end')
        phrqc_input.append('END')         
        self.phrqc_input +=  phrqc_input
        
class PhreeqcKineticsCSH(PhreeqcKinetics):
    def __init__(self, n, kinrate, csh, database, h = 1):
        self.phase = csh['name']
        self.csh = csh
        self.generate_steps(n, h)
        self.krate = kinrate
        self.database = database
        self.phrqc_input = []
        self.selected_output = []
        self.phrqc_string = ''
        self.simulation_time = 0
        
    def generate_steps(self, n, h):
        steps=''
        for i in range(n/h):
            steps = steps + ' ' + str(h)
        self.steps = steps
        
    def generate_phrqc_string(self):
        self.phases()
        self.solution_1()
        self.kinetics_1()
        self.rates_1()
        self.selected_output_1()
        self.user_punch()        
        self.phrqc_string = '\n'.join(self.phrqc_input)
        
    def phases(self):
        phrqc_input = []   
        # CSH stochiometry
        s = self.csh['stochiometry']
        h = s['H+']
        h2o =  s['H2O'] +  s['H+'] - s['Ca']
        sign1 = '+'
        if h < 0: 
            sign1 = '-'
            h *= -1
        sign2 = '+'
        if h2o < 0: 
            sign2 = '-'
            h2o *= -1
        # input 
        phrqc_input.append('PHASES')  
        phrqc_input.append(self.phase)  
        phrqc_input.append('\t(CaO)' + str(s['Ca']) +'(SiO2)'+ str(s['Si']) + \
                           '(H2O)' + str(s['H2O']) + ' ' + sign1 + ' ' + str(h) + 'H+ = ' + \
                           str(s['Ca']) + 'Ca+2 + ' + str(s['Si']) + 'SiO2 '  + sign2 +\
                           ' ' + str(h2o) + ' H2O') 
        #phrqc_input.append('\t-Vm\t' + str(csh['vm']) ) 
        phrqc_input.append('\t-log_K\t' + str(self.csh['log_k']) + '\n')
        self.phrqc_input +=  phrqc_input
        
    def kinetics_1(self):
        phrqc_input = []                        
        phrqc_input.append('KINETICS\t1')
        phrqc_input.append('\t' + self.phase)
        phrqc_input.append('\t-m\t1.e+0')
        phrqc_input.append('\t-m0\t1.e+0')
        phrqc_input.append('\t-tol\t1e-8')
        phrqc_input.append('\t-steps ' + self.steps) # seconds
        phrqc_input.append('\t-runge_kutta 2')
        phrqc_input.append('INCREMENTAL_REACTIONS true')        
        self.phrqc_input +=  phrqc_input
        
        
    def user_punch(self):
        phrqc_input = []         
        phrqc_input.append('USER_PUNCH')
        phrqc_input.append('\t-headings\tCa\t' + self.phase)
        phrqc_input.append('\t-start')
        phrqc_input.append('\t10\tpunch\ttot("Ca")')
        phrqc_input.append('\t15\tpunch\ttot("Si")')
        phrqc_input.append('\t20\tpunch\ttot("' + self.phase + '")')
        phrqc_input.append('\t30\tpunch')
        phrqc_input.append('\t-end')
        phrqc_input.append('END')         
        self.phrqc_input +=  phrqc_input
    
#%% PARAMETERS
database = 'C:\Anaconda2\lib\site-packages\databases\cemdata18.dat'
h = 1000 # step for kinetic rate so that it converges

csh = {'name':'CSH', 'stochiometry':{'Ca':1.25, 'Si':1.0, 'H2O':3.42, 'H+':2.5}, 'log_k':19.873,}

n = 20000000 # time should be ~1.7 hours for mixing
krate = 10**(-9.98) #-9.98
s = 1200#scale factor
fraction = krate * s 
print('Kinetic rate = ' + str(krate))
print('Mixing fraction = ' + str(fraction))

#%% RUN
pm = PhreeqcMixingCSH(n, fraction, csh, database)
pm.run_phreeqc()

pk = PhreeqcKineticsCSH(n, krate, csh, database, h)
pk.run_phreeqc()

print('Kinetic rate simulation time = ' + str(pk.simulation_time))
print('Mixing fraction simulation time = ' + str(pm.simulation_time))

#%% PLOT

t = range(1, n+1)
t = [i/3600. for i in t]

ca_m = []
si_m = []
for i in range(len(pm.selected_output)):
    if pm.selected_output[i][0]==3:
        ca_m.append(pm.selected_output[i][1])
        si_m.append(pm.selected_output[i][2])

ca_k = []
si_k = []
for i in range(len(pk.selected_output)):
    if pk.selected_output[i][0]==1:
        ca_k.append(pk.selected_output[i][1])
        si_k.append(pk.selected_output[i][2])
     
plt.figure()
plt.plot(t, ca_m, label = "mix")
plt.plot(t[0::h], ca_k[1:], label = "kin")
plt.xlabel('time (h)')
plt.ylabel('Ca (mol/l)')
plt.legend()

plt.figure()
plt.plot(t, si_m, label = "mix")
plt.plot(t[0::h], si_k[1:], label = "kin")
plt.xlabel('time (h)')
plt.ylabel('Si (mol/l)')
plt.legend()
plt.show()
