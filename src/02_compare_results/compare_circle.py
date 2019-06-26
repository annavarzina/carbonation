

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
res_nat = fn.load_obj('output_circle/' + 'nat' +'_results') 
res_acc_10 = fn.load_obj('output_circle/' + 'acc10' +'_results') 

#%%
    
plt.figure(figsize=(8,4))
plt.plot(res_nat['time'], res_nat['portlandite'], ls='--', label = 'nat')    
plt.plot(res_acc_10['time'], res_acc_10['portlandite'], ls='--', label = 'acc 10%')     
plt.title('Portlandite')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
plt.plot(res_nat['time'], res_nat['calcite'], ls='--', label = 'nat')    
plt.plot(res_acc_10['time'], res_acc_10['calcite'], ls='--', label = 'acc 10%')     
plt.title('Calcite')
plt.ylabel('CC mass (*1e-12) [mol]')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
plt.plot(res_nat['time'], res_nat['avg_poros'], ls='--', label = 'nat')    
plt.plot(res_acc_10['time'], res_acc_10['avg_poros'], ls='--', label = 'acc 10%')    
plt.title('Porosity')
plt.ylabel('Porosity [-]')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 


plt.figure(figsize=(8,4))
plt.plot(res_nat['time'], res_nat['Ca'], ls='--', label = 'nat')    
plt.plot(res_acc_10['time'], res_acc_10['Ca'], ls='--', label = 'acc 10%')  
plt.title('Ca')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 


plt.figure(figsize=(8,4))
plt.plot(res_nat['time'], res_nat['C'], ls='--', label = 'nat')  
plt.plot(res_acc_10['time'], res_acc_10['C'], ls='--', label = 'acc 10%')   
plt.title('C')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 
