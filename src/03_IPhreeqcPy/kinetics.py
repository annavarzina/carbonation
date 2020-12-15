# -*- coding: utf-8 -*-
#import numpy as np
import time
import IPhreeqcPy

class PhreeqcKinetics():
    def __init__(self, n, kinrate, database, phase_name = 'Portlandite'):
        self.phase = phase_name
        self.generate_steps(n)
        self.krate = kinrate
        self.database = database
        self.phrqc_input = []
        self.selected_output = []
        self.phrqc_string = ''
        self.simulation_time = 0
        
    def generate_steps(self,n):
        steps=''
        for i in range(n):
            steps = steps + ' 1'
        self.steps = steps
        
    def run_phreeqc(self):        
        it=time.time() 
        self.generate_phrqc_string()
        IPhreeqc = IPhreeqcPy.IPhreeqc()
        IPhreeqc.LoadDatabase(self.database)
        IPhreeqc.RunString(self.phrqc_string)
        self.selected_output = IPhreeqc.GetSelectedOutputArray()
        self.simulation_time = time.time()-it
        
    def generate_phrqc_string(self):
        self.solution_1()
        self.kinetics_1()
        self.rates_1()
        self.selected_output_1()
        self.user_punch()        
        self.phrqc_string = '\n'.join(self.phrqc_input)
    
    
    def solution_1(self):
        phrqc_input = []                
        phrqc_input.append('SOLUTION\t1')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        phrqc_input.append('\tCa\t0\n')        
        self.phrqc_input +=  phrqc_input
        
    def kinetics_1(self):
        phrqc_input = []                        
        phrqc_input.append('KINETICS\t1')
        phrqc_input.append('\t' + self.phase)
        phrqc_input.append('\t-m\t1.e+0')
        phrqc_input.append('\t-m0\t1.e+0')
        phrqc_input.append('\t-tol\t1e-8')
        phrqc_input.append('\t-steps ' + self.steps) # seconds
        phrqc_input.append('\t-runge_kutta 6')
        phrqc_input.append('INCREMENTAL_REACTIONS true')        
        self.phrqc_input +=  phrqc_input
        
    def rates_1(self):
        phrqc_input = []                
        phrqc_input.append('RATES')
        phrqc_input.append('\t' + self.phase)    
        phrqc_input.append('-start')
        phrqc_input.append('10\tsi_p = SI("' + self.phase + '")') #define saturation index (log(Q/K))for rate equation
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
        phrqc_input.append('\t-headings\tCa\t' + self.phase)
        phrqc_input.append('\t-start')
        phrqc_input.append('\t10\tpunch\ttot("Ca")')
        phrqc_input.append('\t20\tpunch\ttot("' + self.phase + '")')
        phrqc_input.append('\t30\tpunch')
        phrqc_input.append('\t-end')
        phrqc_input.append('END')                 
        self.phrqc_input +=  phrqc_input
        
    def save_phrqc_input(self, root_dir, name):
        with open(root_dir +'\\phreeqc_input\\' + name + '.phrq', 'w') as f:
            for item in self.phrqc_input:
                f.write("%s\n" % item)