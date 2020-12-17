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
h = 1 # step for kinetic rate so that it converges
krate = 10**(-6) #-10.99
#s = 1600#scale factor
ca_si = np.array([1.67, 1.65, 1.60, 1.55, 1.50, 1.45, 1.40, 1.35, 1.30, 1.25,
                  1.20, 1.15, 1.10, 1.05, 1.00, 0.95, 0.90, 0.85, 0.83])
h2o = np.array([4.34, 4.3, 4.19, 4.08, 3.97, 3.86, 3.75, 3.64, 3.53, 3.42, 3.31, 
                3.19, 3.08, 2.97, 2.86, 2.75, 2.64, 2.53, 2.49])
h_plus = np.array([3.34, 3.3, 3.2, 3.1, 3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 
                   2.3, 2.2, 2.1, 2.0, 1.9, 1.8, 1.7, 1.66]) 
log_k = np.array([29.133, 28.724, 27.572, 26.443, 25.328, 24.222, 23.124, 
                  22.034, 20.95, 19.873, 18.801, 17.736, 16.678, 15.627, 
                  14.583, 13.55, 12.529, 11.531, 11.15])
s = np.array([865, 880, 925,  970, 1015, 1060, 1105, 1150, 1195, 1240, 1285, 1330,
              1375, 1420, 1465, 1510, 1555, 1600, 1615])
n = np.array(range(10000, 5000, -250))

#%% LOOP
#for i in range(0, len(ca_si)):
for j in [0,9, 18]: # range(16, 18):
    print('Ca/Si %.2f' %ca_si[j])    
    csh = {'name':'CSH', 'stochiometry':{'Ca':ca_si[j], 'Si':1.0, 
                                         'H2O':h2o[j], 'H+':h_plus[j]}, 
           'log_k':log_k[j]}
    fraction = krate * s[j] 
    print('Kinetic rate = ' + str(krate))
    print('Mixing fraction = ' + str(fraction))
    #%% RUN
    
    pm = PhreeqcMixingCSH(n[j], fraction, csh, database)
    pm.run_phreeqc()
    
    pk = PhreeqcKineticsCSH(n[j], krate, csh, database, h)
    pk.run_phreeqc()
    
    print('Kinetic rate simulation time = ' + str(pk.simulation_time))
    print('Mixing fraction simulation time = ' + str(pm.simulation_time))

#%% PLOT    
    t = range(1, n[j]+1)
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
    
    '''     
    plt.figure()
    plt.plot(t, ca_m, label = "mix")
    plt.plot(t[0::h], ca_k[1:], label = "kin")
    plt.xlabel('time (h)')
    plt.ylabel('Ca (mol/l)')
    plt.title('Ca concentration change for Ca/Si %.2f' %ca_si[j])    
    plt.legend()
    
    plt.figure()
    plt.plot(t, si_m, label = "mix")
    plt.plot(t[0::h], si_k[1:], label = "kin")
    plt.xlabel('time (h)')
    plt.ylabel('Si (mol/l)')
    plt.title('Si concentration change for Ca/Si %.2f' %ca_si[j])    
    plt.legend()
    plt.show()
    '''
#%% Predict model
d = np.array([-8.40, -8.47, -8.53, -8.53, -8.56, -8.67, -8.89, -9.20, -9.58,
              -9.98, -10.36, -10.67, -10.88, -10.97, -10.95, -10.86, -10.79, 
              -10.87, -10.99])
#d = np.array([-9.01,-9.12,-9.28,-9.30,-9.07,-9.06,-9.28,-9.45,-9.79,-10.49,-10.70,
#              -10.45,-10.4,-10.9,-10.95,-10.95,-10.9,-10.9,-10.99])
sigma = (10**d) * s
plt.plot(ca_si, sigma)
plt.yscale("log")
plt.show()
''' #or
d = np.array([-8.40, -8.43, -8.47, -8.51, -8.56, -8.67, -8.89, -9.20, -9.58,
              -9.98, -10.36, -10.67, -10.88, -10.92, -10.95, -10.96, -10.97, 
              -10.98, -10.99])
plt.plot(ca_si, d)
plt.show()
'''
#%%linear regression model
from sklearn.linear_model import LinearRegression

X = ca_si.reshape(-1,1) # put your dates in here
y = sigma # put your kwh in here

model = LinearRegression()
model.fit(X, y)

X_predict = np.array([1.52654]).reshape(-1,1)  # put the dates of which you want to predict kwh here
y_predict = model.predict(X_predict)

#%%
'''
for j in [0,9]: # range(16, 18):
    h = int(10**(-6-d[j]))*10
    krate = 10**d[j]
    print('Ca/Si %.2f' %ca_si[j])    
    csh = {'name':'CSH', 'stochiometry':{'Ca':ca_si[j], 'Si':1.0, 
                                         'H2O':h2o[j], 'H+':h_plus[j]}, 
           'log_k':log_k[j]}
    print('Kinetic rate = ' + str(krate))

    pk = PhreeqcKineticsCSH(int(n[j]*h/10), krate, csh, database, h)
    pk.run_phreeqc()
    
    print('Kinetic rate simulation time = ' + str(pk.simulation_time))

 
    t = range(1, int(n[j]*h/10)+1)
    t = [i/3600. for i in t]
    ca_k = []
    si_k = []
    for i in range(len(pk.selected_output)):
        if pk.selected_output[i][0]==1:
            ca_k.append(pk.selected_output[i][1])
            si_k.append(pk.selected_output[i][2])
    
    plt.figure()
    plt.plot(t[0::h], ca_k[1:], label = "kin")
    plt.xlabel('time (h)')
    plt.ylabel('Ca (mol/l)')
    plt.title('Ca concentration change for Ca/Si %.2f' %ca_si[j])    
    plt.legend()
        
    plt.figure()
    plt.plot(t[0::h], si_k[1:], label = "kin")
    plt.xlabel('time (h)')
    plt.ylabel('Si (mol/l)')
    plt.title('Si concentration change for Ca/Si %.2f' %ca_si[j])    
    plt.legend()
    plt.show()
'''