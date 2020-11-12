# -*- coding: utf-8 -*-

def set_phrqc_input(p, ptype ='CSH'):
	# TODO add initial Si concentrations
    # TODO CSH parameters
    '''
    Example:
    p = {'c_bc':{'type':'conc', 'value': 0.01}, 
         'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
         'c_liq':{'type':'eq', 'value': 'calcite'},
         'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
         'ca_liq':{'type':'eq', 'value': 'calcite'} }
    '''
    phrqc_input = [] 
    if(ptype=='CSH'):
        phrqc_input.append('PHASES')  
        phrqc_input.append('CSHQ_TobH')  
        phrqc_input.append('\t(CaO)0.66666667(SiO2)1(H2O)1.5 = 0.66666667Ca++ + 1 SiO(OH)3- + 0.33333334OH- -0.16666667 H2O') 
        phrqc_input.append('\tlog_K -6.190832') 
        phrqc_input.append('CSHQ_TobD') 
        phrqc_input.append('\t(CaO)0.8333333333(SiO2)0.6666666667(H2O)1.8333333333 = 0.8333333333 Ca++ + 0.6666666667 SiO(OH)3- + 0.99999999990 OH- + 0.3333333333 H2O') 
        phrqc_input.append('\tlog_K -6.8995533') 
        phrqc_input.append('CSHQ_JenH') 
        phrqc_input.append('\t(CaO)1.3333333333(SiO2)1(H2O)2.1666666667 = 1.3333333333 Ca++ + 1 SiO(OH)3- + 1.6666666667 OH- -0.1666666667 H2O') 
        phrqc_input.append('\tlog_K -10.96765') 
        phrqc_input.append('CSHQ_JenD') 
        phrqc_input.append('\t(CaO)1.5(SiO2)0.6666666667(H2O)2.5 = 1.5 Ca++ + 0.6666666667 SiO(OH)3- + 2.3333333333 OH- + 0.3333333333 H2O') 
        phrqc_input.append('\tlog_K -10.47635') 
        phrqc_input.append('knobs') 
        phrqc_input.append('\t-iterations 8000') 
    phrqc_input += set_phrqc_bc(p['ca_bc'], ptype)
    phrqc_input += set_phrqc_liquid(p['ca_liq'], ptype)
    phrqc_input += set_phrqc_mlvl(p['ca_mlvl'], ptype)    

    phrqc_input += set_phrqc_solid()
    return phrqc_input
    
def set_phrqc_bc(ca, ptype):
    phrqc_input = [] 
    phrqc_input.append('#boundary_solution')    
    phrqc_input.append('SOLUTION\t100001')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    if(ca['type'] == 'conc'):
        phrqc_input.append('\tCa\t' + str(ca['value']) + '\n')
    else:
        pass
    phrqc_input.append('EQUILIBRIUM_PHASES\t100001\n')
    return phrqc_input
    
def set_phrqc_liquid(ca, ptype):
    phrqc_input = [] 
    phrqc_input.append('#solution_liquid')    
    phrqc_input.append('SOLUTION\t100002')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    if(ca['type'] == 'conc'):
        phrqc_input.append('\tCa\t' + str(ca['value']))
    elif(ca['type'] == 'eq'):
        phrqc_input.append('\tCa\t1\t' + str(ca['value']))
    else:
        pass        
    phrqc_input.append('EQUILIBRIUM_PHASES\t100002')
    phrqc_input.append('portlandite\t0\t0')
    return phrqc_input

def set_phrqc_mlvl(ca, ptype):
    phrqc_input = [] 
    phrqc_input.append('#solution_multilevel')    
    phrqc_input.append('SOLUTION\t100003')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    if(ca['type'] == 'conc'):
        phrqc_input.append('\tCa\t' + str(ca['value']))
    elif(ca['type'] == 'eq'):
        phrqc_input.append('\tCa\t1\t' + str(ca['value']))
    else:
        pass        
    phrqc_input.append('EQUILIBRIUM_PHASES\t100003')
    phrqc_input.append('portlandite\t0\t1')
    if(ptype=='CSH'):
        phrqc_input.append('#solution_csh_multilevel') 
        phrqc_input.append('SOLUTION\t100004') 
        phrqc_input.append('\t-water\t0.448230266981165') 
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\tpH\t12\tcharge') 
        phrqc_input.append('\tCa\t1.955e-002') 
        phrqc_input.append('\tSi\t3.018e-005') 
        phrqc_input.append('SOLID_SOLUTIONS\t100004') 
        phrqc_input.append('Tob_jen_ss') 
        phrqc_input.append('\t-comp\tCSHQ_TobH\t0.1041') 
        phrqc_input.append('\t-comp\tCSHQ_TobD\t2.5050') 
        phrqc_input.append('\t-comp\tCSHQ_JenH\t2.1555') 
        phrqc_input.append('\t-comp\tCSHQ_JenD\t3.2623') 
        phrqc_input.append('EQUILIBRIUM_PHASES\t100004')
        phrqc_input.append('portlandite\t0\t0')
    return phrqc_input

def set_phrqc_solid():
    phrqc_input = [] 
    phrqc_input.append('#solution_solid')    
    phrqc_input.append('SOLUTION\t100005')
    phrqc_input.append('\t-water\t1\n')
    return phrqc_input