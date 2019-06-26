# -*- coding: utf-8 -*-
'''
The example of portlandite dissolution and calcite precipitation 
when CO2 flows in from the left.
'''

#%% Python modules
from __future__ import division  #using floating everywhere
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../..'))
import matplotlib.pylab as plt
import numpy as np
import time
from copy import deepcopy
np.set_printoptions(precision=5, threshold=np.inf)

#%% Custom modules
sys.path.append('C:\\Users\\avarzina\\Documents\\Python_Projects\\yantra') # change the path
import yantra
import modules.add_ons.cell_type as ct # change the path to cell_type file
import func as fn
#%%
Ts =10.

params = fn.get_params_list()

#%% 
res_si_05 = fn.load_obj('C:/Users/avarzina/Documents/Python_Projects/yantra/in_progress/25_refactoring_and_si_based_porosity/output_19-05-26/' + 'test' +'_SI-0.5' +'_results') 
res_si_05D = fn.load_obj('output/' + 'test' +'_results') 

#%%
    
plt.figure(figsize=(8,4))
plt.plot(res_si_05['time'], res_si_05['portlandite'], ls='--', label = 'Archie')    
plt.plot(res_si_05D['time'], res_si_05D['portlandite'], ls='--', label = 'Fixed D')    
plt.title('Portlandite')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
plt.plot(res_si_05['time'], res_si_05['calcite'], ls='--', label = 'Archie')    
plt.plot(res_si_05D['time'], res_si_05D['calcite'], ls='--', label = 'Fixed D')  
plt.title('Calcite')
plt.ylabel('CC mass (*1e-12) [mol]')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
plt.plot(res_si_05['time'], res_si_05['avg_poros'], ls='--', label = 'Archie')    
plt.plot(res_si_05D['time'], res_si_05D['avg_poros'], ls='--', label = 'Fixed D')  
plt.title('Porosity')
plt.ylabel('Porosity [-]')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 


#%% 
res_si_50000 = fn.load_obj('output/' + 'test_50000_pores' +'_results') 
res_si_1000 = fn.load_obj('output/' + 'test_1000_pores' +'_results') 
res_si_5000 = fn.load_obj('output/' + 'test_5000_pores' +'_results') 
res_si_10000 = fn.load_obj('output/' + 'test_10000_pores' +'_results') 
#res_si_10000_acc = fn.load_obj('output/' + 'test_10000_pores_10C' +'_results') 
#%%
    
plt.figure(figsize=(8,4))
plt.plot(res_si_5000['time'], res_si_5000['portlandite'], ls='--', label = '5000')    
plt.plot(res_si_10000['time'], res_si_10000['portlandite'], ls='--', label = '10000')    
plt.plot(res_si_50000['time'], res_si_50000['portlandite'], ls='--', label = '50000')    
plt.plot(res_si_1000['time'], res_si_1000['portlandite'], ls='--', label = '1000')    
#plt.plot(res_si_10000_acc['time'], res_si_10000_acc['portlandite'], ls='--', label = '10000 acc')    
plt.title('Portlandite')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
plt.plot(res_si_5000['time'], res_si_5000['calcite'], ls='--', label = '5000')    
plt.plot(res_si_10000['time'], res_si_10000['calcite'], ls='--', label = '10000') 
plt.plot(res_si_50000['time'], res_si_50000['calcite'], ls='--', label = '50000')    
plt.plot(res_si_1000['time'], res_si_1000['calcite'], ls='--', label = '1000') 
#plt.plot(res_si_10000_acc['time'], res_si_10000_acc['calcite'], ls='--', label = '10000 acc')      
plt.title('Calcite')
plt.ylabel('CC mass (*1e-12) [mol]')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
plt.plot(res_si_5000['time'], res_si_5000['avg_poros'], ls='--', label = '5000')    
plt.plot(res_si_10000['time'], res_si_10000['avg_poros'], ls='--', label = '10000')  
plt.plot(res_si_50000['time'], res_si_50000['avg_poros'], ls='--', label = '50000')    
plt.plot(res_si_1000['time'], res_si_1000['avg_poros'], ls='--', label = '1000')  
#plt.plot(res_si_10000_acc['time'], res_si_10000_acc['avg_poros'], ls='--', label = '10000 acc')    
plt.title('Porosity')
plt.ylabel('Porosity [-]')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 


plt.figure(figsize=(8,4))
plt.plot(res_si_5000['time'], res_si_5000['Ca'], ls='--', label = '5000')    
plt.plot(res_si_10000['time'], res_si_10000['Ca'], ls='--', label = '10000')
plt.plot(res_si_50000['time'], res_si_50000['Ca'], ls='--', label = '50000')    
plt.plot(res_si_1000['time'], res_si_1000['Ca'], ls='--', label = '1000')  
#plt.plot(res_si_10000_acc['time'], res_si_10000_acc['Ca'], ls='--', label = '10000 acc')    
plt.title('Ca')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 


plt.figure(figsize=(8,4))
plt.plot(res_si_5000['time'], res_si_5000['C'], ls='--', label = '5000')    
plt.plot(res_si_10000['time'], res_si_10000['C'], ls='--', label = '10000')  
plt.plot(res_si_50000['time'], res_si_50000['C'], ls='--', label = '50000')    
plt.plot(res_si_1000['time'], res_si_1000['C'], ls='--', label = '10000') 
#plt.plot(res_si_10000_acc['time'], res_si_10000_acc['C'], ls='--', label = '10000 acc')    
plt.title('C')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 
