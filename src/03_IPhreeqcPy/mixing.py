# -*- coding: utf-8 -*-
import numpy as np
import time
import IPhreeqcPy

class PhreeqcMixing():
    def __init__(self, n, fraction, database, phase_name = 'Portlandite'):
        self.phase = phase_name
        self.steps = n
        self.fraction = fraction
        self.database = database
        self.phrqc_input = []
        self.selected_output = []
        self.phrqc_string = ''
        self.simulation_time = 0
        
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
        self.solution_2()        
        for i in np.arange(0, self.steps):
            self.mix_2()
            self.selected_output_1()
            self.user_punch()
            self.mix_3()            
        self.phrqc_string = '\n'.join(self.phrqc_input)
    
    def solution_1(self):
        phrqc_input = []         
        phrqc_input.append('SOLUTION\t1')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        phrqc_input.append('\tCa\t0\n')        
        phrqc_input.append('EQUILIBRIUM_PHASES\t1')
        phrqc_input.append('\t' + self.phase + '\t0\t0')
        phrqc_input.append('END\n')        
        self.phrqc_input +=  phrqc_input
        
    def solution_2(self):
        phrqc_input = []         
        phrqc_input.append('SOLUTION\t2')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        phrqc_input.append('\tCa\t0\n')        
        phrqc_input.append('EQUILIBRIUM_PHASES\t2')
        phrqc_input.append('\t' + self.phase + '\t0\t1')
        phrqc_input.append('END\n')        
        self.phrqc_input +=  phrqc_input
        
    def mix_2(self):
        phrqc_input = [] 
                
        phrqc_input.append('MIX\t2')
        phrqc_input.append('1\t' + str(self.fraction))
        phrqc_input.append('USE solution 1')
        phrqc_input.append('USE equilibrium_phases 2')
        phrqc_input.append('SAVE equilibrium_phases 2')
        phrqc_input.append('SAVE solution 2')
        phrqc_input.append('END\n')
        
        self.phrqc_input +=  phrqc_input
       
    def mix_3(self):
        phrqc_input = []         
        phrqc_input.append('MIX\t3')
        phrqc_input.append('2\t1')
        phrqc_input.append('1\t' + str(1.-self.fraction))
        phrqc_input.append('USE solution 2')
        phrqc_input.append('USE equilibrium_phases 1')
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