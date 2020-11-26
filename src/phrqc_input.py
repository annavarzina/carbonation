# -*- coding: utf-8 -*-
class PhreeqcInput():
    
    def __init__(self, input_concentration):
        self.c = input_concentration
        self.phrqc_input = []
        self.make_phrqc_input()
        
    def make_phrqc_input(self):
        '''
        Example:
        p = {'c_bc':{'type':'conc', 'value': 0.01}, 
             'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
             'c_liq':{'type':'eq', 'value': 'calcite'},
             'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
             'ca_liq':{'type':'eq', 'value': 'calcite'} }
        '''
        # default is for CH carbonation
        #self.phrqc_input = [] 
        self.phrqc_boundary_voxel(self.c['c_bc'])
        self.phrqc_liquid_voxel(self.c['c_liq'], self.c['ca_liq'])
        self.phrqc_multilevel_voxel(self.c['c_mlvl'], self.c['ca_mlvl'])        
        self.phrqc_solid_voxel()
        return(self.phrqc_input)
    
    def save_phrqc_input(self, root_dir, name):
        with open(root_dir +'\\phreeqc_input\\' + name + '.phrq', 'w') as f:
            for item in self.phrqc_input:
                f.write("%s\n" % item)
    
    def phrqc_boundary_voxel(self,c):
        phrqc_input = [] 
        phrqc_input.append('#boundary_solution')    
        phrqc_input.append('SOLUTION\t100001')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(c['type'] == 'conc'):
            phrqc_input.append('\tC\t' + str(c['value']) + '\n')
        elif(c['type'] == 'pco2'):
            phrqc_input.append('\tC\t1\tCO2(g)\t-' + str(c['value']) + '\n')
        phrqc_input.append('EQUILIBRIUM_PHASES\t100001\n')
        self.phrqc_input +=  phrqc_input
    
    def phrqc_liquid_voxel(self, c, ca):
        phrqc_input = [] 
        phrqc_input.append('#solution_liquid')    
        phrqc_input.append('SOLUTION\t100002')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(c['type'] == 'conc'):
            phrqc_input.append('\tC\t' + str(c['value']))
        elif(c['type'] == 'eq'):
            phrqc_input.append('\tC\t1\t' + str(c['value']))
        else:
            pass        
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100002')
        phrqc_input.append('portlandite\t0\t0')
        phrqc_input.append('calcite\t0\t0\n')
        self.phrqc_input +=  phrqc_input

    def phrqc_multilevel_voxel(self, c, ca):
        phrqc_input = [] 
        phrqc_input.append('#solution_multilevel')    
        phrqc_input.append('SOLUTION\t100003')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(c['type'] == 'conc'):
            phrqc_input.append('\tC\t' + str(c['value']))
        elif(c['type'] == 'eq'):
            phrqc_input.append('\tC\t1\t' + str(c['value']))
        else:
            pass              
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100003')
        phrqc_input.append('portlandite\t0\t1')
        phrqc_input.append('calcite\t0\t0\n')
        self.phrqc_input +=  phrqc_input

    def phrqc_solid_voxel(self):
        phrqc_input = [] 
        phrqc_input.append('#solution_solid')    
        phrqc_input.append('SOLUTION\t100005')
        phrqc_input.append('\t-water\t1\n')
        self.phrqc_input +=  phrqc_input
            