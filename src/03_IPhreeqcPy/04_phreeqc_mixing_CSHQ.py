import numpy as np
import matplotlib.pylab as plt
from  mixing import PhreeqcMixing
from kinetics import PhreeqcKinetics

class PhreeqcMixingCSHQ(PhreeqcMixing):
    def __init__(self, n, fraction, database):
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
        #print(self.phrqc_string)
        
    def phases(self):
        phrqc_input = []   
        phrqc_input.append('PHASES')  
        phrqc_input.append('CSHQ_TobH')  
        phrqc_input.append('\t(CaO)0.6667(SiO2)1(H2O)1.5 + 1.3334H+ = 0.6667Ca+2 + 2.1667H2O + SiO2') 
        phrqc_input.append('\t-log_K\t8.286642') 
        phrqc_input.append('\t-analytical_expression\t-12.519254 0 2163.381583 5.476331 0 0 0') 
        phrqc_input.append('CSHQ_TobD') 
        phrqc_input.append('\t((CaO)1.25(SiO2)1(H2O)2.75)0.6667 + 1.66675H+ = 0.833375Ca+2 + 2.6668H2O + 0.6667SiO2') 
        phrqc_input.append('\t-log_K\t13.655314') 
        phrqc_input.append('\t-analytical_expression\t-10.916344 0 3959.367696 4.563888 0 0 0') 
        phrqc_input.append('CSHQ_JenH') 
        phrqc_input.append('\t(CaO)1.3333(SiO2)1(H2O)2.1667 + 2.6666H+ = 1.3333Ca+2 + 3.5H2O + SiO2') 
        phrqc_input.append('\t-log_K\t22.179305') 
        phrqc_input.append('\t-analytical_expression\t-17.10944 0 6470.553982 7.107847 0 0 0') 
        phrqc_input.append('CSHQ_JenD') 
        phrqc_input.append('\t(CaO)1.5(SiO2)0.6667(H2O)2.5 + 3H+ = 1.5Ca+2 + 4H2O + 0.6667SiO2') 
        phrqc_input.append('\t-log_K\t28.730362')  
        phrqc_input.append('\t-analytical_expression\t-15.591756 0 8609.739692 6.24251 0 0 0') 
        phrqc_input.append('\n') 
        self.phrqc_input +=  phrqc_input
        
    
    def solution_1(self):
        phrqc_input = []         
        phrqc_input.append('SOLUTION\t1')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        phrqc_input.append('\tCa\t0')         
        phrqc_input.append('\tSi\t0\n')    
        phrqc_input.append('SOLID_SOLUTIONS\t1')
        phrqc_input.append('Tob_jen_ss')
        phrqc_input.append('\t-comp\tCSHQ_TobH\t0')
        phrqc_input.append('\t-comp\tCSHQ_TobD\t0')
        phrqc_input.append('\t-comp\tCSHQ_JenH\t0')
        phrqc_input.append('\t-comp\tCSHQ_JenD\t0')
        phrqc_input.append('END\n')        
        self.phrqc_input +=  phrqc_input
        
    def solution_2(self):
        phrqc_input = []         
        phrqc_input.append('SOLUTION\t2')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        phrqc_input.append('\tCa\t0\n')     
        phrqc_input.append('\tSi\t0\n')     
        phrqc_input.append('SOLID_SOLUTIONS\t2')
        phrqc_input.append('Tob_jen_ss')
        phrqc_input.append('\t-comp\tCSHQ_TobH\t0.1041')
        phrqc_input.append('\t-comp\tCSHQ_TobD\t2.5050')
        phrqc_input.append('\t-comp\tCSHQ_JenH\t2.1555')
        phrqc_input.append('\t-comp\tCSHQ_JenD\t3.2623')
        phrqc_input.append('END\n')        
        self.phrqc_input +=  phrqc_input
        
    def mix_2(self):
        phrqc_input = [] 
                
        phrqc_input.append('MIX\t2')
        phrqc_input.append('1\t' + str(self.fraction))
        phrqc_input.append('USE solution 1')
        phrqc_input.append('USE SOLID_SOLUTIONS 2')
        phrqc_input.append('SAVE SOLID_SOLUTIONS 2')
        phrqc_input.append('SAVE solution 2')
        phrqc_input.append('END\n')
        
        self.phrqc_input +=  phrqc_input
       
    def mix_3(self):
        phrqc_input = []         
        phrqc_input.append('MIX\t3')
        phrqc_input.append('2\t1')
        phrqc_input.append('1\t' + str(1.-self.fraction))
        phrqc_input.append('USE solution 2')
        phrqc_input.append('USE SOLID_SOLUTIONS 1')
        phrqc_input.append('SAVE solution 1')
        phrqc_input.append('END\n')        
        self.phrqc_input +=  phrqc_input
        
    def selected_output_1(self):
        phrqc_input = []          
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
        self.phrqc_input +=  phrqc_input
    
    def user_punch(self):
        phrqc_input = []         
        phrqc_input.append('USER_PUNCH')
        phrqc_input.append('\t-headings\tCa\tSi\tCSHQ_TobH\tCSHQ_TobD\tCSHQ_JenH\tCSHQ_JenD')
        phrqc_input.append('\t-start')
        phrqc_input.append('\t10\tpunch\ttot("Ca")')
        phrqc_input.append('\t15\tpunch\ttot("Si")')
        phrqc_input.append('\t20\tpunch\ttot("CSHQ_TobH")')
        phrqc_input.append('\t30\tpunch\ttot("CSHQ_TobD")')
        phrqc_input.append('\t40\tpunch\ttot("CSHQ_JenH")')
        phrqc_input.append('\t50\tpunch\ttot("CSHQ_JenD")')
        phrqc_input.append('\t60\tpunch')
        phrqc_input.append('\t-end')
        phrqc_input.append('END')         
        self.phrqc_input +=  phrqc_input
    
        
class PhreeqcKineticsCSHQ(PhreeqcKinetics):
    def __init__(self, n, kinrate, database, h = 1):
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
        phrqc_input.append('PHASES')  
        phrqc_input.append('CSHQ_TobH')  
        phrqc_input.append('\t(CaO)0.6667(SiO2)1(H2O)1.5 + 1.3334H+ = 0.6667Ca+2 + 2.1667H2O + SiO2') 
        phrqc_input.append('\t-log_K\t8.286642') 
        phrqc_input.append('\t-analytical_expression\t-12.519254 0 2163.381583 5.476331 0 0 0') 
        phrqc_input.append('CSHQ_TobD') 
        phrqc_input.append('\t((CaO)1.25(SiO2)1(H2O)2.75)0.6667 + 1.66675H+ = 0.833375Ca+2 + 2.6668H2O + 0.6667SiO2') 
        phrqc_input.append('\t-log_K\t13.655314') 
        phrqc_input.append('\t-analytical_expression\t-10.916344 0 3959.367696 4.563888 0 0 0') 
        phrqc_input.append('CSHQ_JenH') 
        phrqc_input.append('\t(CaO)1.3333(SiO2)1(H2O)2.1667 + 2.6666H+ = 1.3333Ca+2 + 3.5H2O + SiO2') 
        phrqc_input.append('\t-log_K\t22.179305') 
        phrqc_input.append('\t-analytical_expression\t-17.10944 0 6470.553982 7.107847 0 0 0') 
        phrqc_input.append('CSHQ_JenD') 
        phrqc_input.append('\t(CaO)1.5(SiO2)0.6667(H2O)2.5 + 3H+ = 1.5Ca+2 + 4H2O + 0.6667SiO2') 
        phrqc_input.append('\t-log_K\t28.730362')  
        phrqc_input.append('\t-analytical_expression\t-15.591756 0 8609.739692 6.24251 0 0 0') 
        phrqc_input.append('\n') 
        self.phrqc_input +=  phrqc_input
        
    def solution_1(self):
        phrqc_input = []                
        phrqc_input.append('SOLUTION\t1')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        phrqc_input.append('\tCa\t0')  
        phrqc_input.append('\tSi\t0\n')       
        self.phrqc_input +=  phrqc_input
        
    def kinetics_1(self):
        phrqc_input = []                        
        phrqc_input.append('KINETICS\t1')
        phrqc_input.append('\tCSHQ_TobH')
        phrqc_input.append('\t-m\t0.1041')
        phrqc_input.append('\t-m0\t0.1041')
        phrqc_input.append('\tCSHQ_TobD')
        phrqc_input.append('\t-m\t2.5050')
        phrqc_input.append('\t-m0\t2.5050')
        phrqc_input.append('\tCSHQ_JenH')
        phrqc_input.append('\t-m\t2.1555')
        phrqc_input.append('\t-m0\t2.1555')
        phrqc_input.append('\tCSHQ_JenD')
        phrqc_input.append('\t-m\t3.2623')
        phrqc_input.append('\t-m0\t3.2623')
        phrqc_input.append('\t-tol\t1e-8')
        phrqc_input.append('\t-steps ' + self.steps) # seconds
        phrqc_input.append('\t-runge_kutta 2')
        phrqc_input.append('INCREMENTAL_REACTIONS true')        
        self.phrqc_input +=  phrqc_input
        
    def rates_1(self):
        phrqc_input = []                
        phrqc_input.append('RATES')
        phrqc_input.append('\tCSHQ_TobH')    
        phrqc_input.append('-start')
        phrqc_input.append('10\tsi_p = SI("CSHQ_TobH")') #define saturation index (log(Q/K))for rate equation
        phrqc_input.append('20\tQ_K = 10^si_p') #calculate Q/K
        phrqc_input.append('30\tk = '+ str(self.krate/125))#mol/l/s/m2, kinetic rate constant 
        phrqc_input.append('40\tr = k * (1-Q_K)')#r is kinetic rate in mol/l.s, equals to k*(1-Q/K)
        phrqc_input.append('50\tmoles = r * TIME')
        phrqc_input.append('60\tSAVE moles')
        phrqc_input.append('-end')   
        phrqc_input.append('\tCSHQ_TobD')    
        phrqc_input.append('-start')
        phrqc_input.append('10\tsi_p = SI("CSHQ_TobD")') #define saturation index (log(Q/K))for rate equation
        phrqc_input.append('20\tQ_K = 10^si_p') #calculate Q/K
        phrqc_input.append('30\tk = '+ str(self.krate/25))#mol/l/s/m2, kinetic rate constant 
        phrqc_input.append('40\tr = k * (1-Q_K)')#r is kinetic rate in mol/l.s, equals to k*(1-Q/K)
        phrqc_input.append('50\tmoles = r * TIME')
        phrqc_input.append('60\tSAVE moles')
        phrqc_input.append('-end')      
        phrqc_input.append('\tCSHQ_JenH')    
        phrqc_input.append('-start')
        phrqc_input.append('10\tsi_p = SI("CSHQ_JenH")') #define saturation index (log(Q/K))for rate equation
        phrqc_input.append('20\tQ_K = 10^si_p') #calculate Q/K
        phrqc_input.append('30\tk = '+ str(self.krate/5))#mol/l/s/m2, kinetic rate constant 
        phrqc_input.append('40\tr = k * (1-Q_K)')#r is kinetic rate in mol/l.s, equals to k*(1-Q/K)
        phrqc_input.append('50\tmoles = r * TIME')
        phrqc_input.append('60\tSAVE moles')
        phrqc_input.append('-end') 
        phrqc_input.append('\tCSHQ_JenD')    
        phrqc_input.append('-start')
        phrqc_input.append('10\tsi_p = SI("CSHQ_JenD")') #define saturation index (log(Q/K))for rate equation
        phrqc_input.append('20\tQ_K = 10^si_p') #calculate Q/K
        phrqc_input.append('30\tk = '+ str(self.krate))#mol/l/s/m2, kinetic rate constant 
        phrqc_input.append('40\tr = k * (1-Q_K)')#r is kinetic rate in mol/l.s, equals to k*(1-Q/K)
        phrqc_input.append('50\tmoles = r * TIME')
        phrqc_input.append('60\tSAVE moles')
        phrqc_input.append('-end')    
        self.phrqc_input +=  phrqc_input
        
     
    def selected_output_1(self):
        phrqc_input = []          
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
        self.phrqc_input +=  phrqc_input
    
    def user_punch(self):
        phrqc_input = []         
        phrqc_input.append('USER_PUNCH')
        phrqc_input.append('\t-headings\tCa\tSi\tCSHQ_TobH\tCSHQ_TobD\tCSHQ_JenH\tCSHQ_JenD')
        phrqc_input.append('\t-start')
        phrqc_input.append('\t10\tpunch\ttot("Ca")')
        phrqc_input.append('\t15\tpunch\ttot("Si")')
        phrqc_input.append('\t20\tpunch\ttot("CSHQ_TobH")')
        phrqc_input.append('\t30\tpunch\ttot("CSHQ_TobD")')
        phrqc_input.append('\t40\tpunch\ttot("CSHQ_JenH")')
        phrqc_input.append('\t50\tpunch\ttot("CSHQ_JenD")')
        phrqc_input.append('\t60\tpunch')
        phrqc_input.append('\t-end')
        phrqc_input.append('END')         
        self.phrqc_input +=  phrqc_input
        
        #%% PARAMETERS
database = 'C:\Anaconda2\lib\site-packages\databases\cemdata18.dat'

h = 10

n = 100000 # time should be ~10 minutes
krate = 10**(-5) #1e-7
s = 80#scale factor
fraction = 0.0001#krate * s 
print('Kinetic rate = ' + str(krate))
print('Mixing fraction = ' + str(fraction))

#%% RUN
pm = PhreeqcMixingCSHQ(n, fraction, database)
pm.run_phreeqc()
print('Mixing fraction simulation time = ' + str(pm.simulation_time))

pk = PhreeqcKineticsCSHQ(n, krate, database, h)
pk.run_phreeqc()
print('Kinetic rate simulation time = ' + str(pk.simulation_time))

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
#%%    
plt.figure()
plt.plot(t, ca_m, label = "mix")
plt.plot(t, ca_k[1:], label = "kin")#[1::h]
plt.xlabel('time (h)')
plt.ylabel('Ca (mol/l)')
plt.legend()

plt.figure()
plt.plot(t, si_m, label = "mix")
plt.plot(t, si_k[1:], label = "kin")#[1::h]
plt.xlabel('time (h)')
plt.ylabel('Si (mol/l)')
plt.legend()
plt.show()

#%% csh phases
'''
tobh = []
tobd = []
jenh = []
jend = []
for i in range(len(pm.selected_output)):
    if pm.selected_output[i][0]==1:
        tobh.append(pm.selected_output[i][3])
        tobd.append(pm.selected_output[i][4])
        jenh.append(pm.selected_output[i][5])
        jend.append(pm.selected_output[i][6])
        
plt.figure()
plt.plot(tobh, label = "tobh")
plt.plot(tobd, label = "tobd")
plt.plot(jenh, label = "jenh")
plt.plot(jend, label = "jend")
plt.xlabel('time (h)')
plt.ylabel('CSH (mol/dm3)')
plt.legend()
'''