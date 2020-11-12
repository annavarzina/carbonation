# -*- coding: utf-8 -*-
'''
Compare the results for different PCO2
'''
#%% MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(root_dir)
sys.path.append(src_dir)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import misc_func as fn
import plot_func as cf
#%% SETTINGS
Ts =100.
fname = 'pco2'
fpath = root_dir+'\\results\\output\\07_pco2\\compare\\'
fn.make_output_dir(fpath)
#names = np.array([ '02_pco2_0', '02_pco2_1', '02_pco2_2', '02_pco2_3', '01_reference'])
#names = np.array(['02_pco2_1_p005', '03_pco2_2_p005','04_pco2_3_p005', '05_pco2_34_p005'])
#names = np.array(['02_pco2_1_p005', '03_pco2_2_p005', '04_pco2_3_p005', '05_pco2_34_p005'])
#names = np.array(['01_p005_c34m_ll1', '03_p005_c252m_ll1', '06_p005_c152m_ll1', '05_p005_c1m_ll1'])
#label = np.array(['0.04%', '0.3%', '3%', '10%'])
#names = np.array(['01_p005_c34', '03_p005_c252', '06_p005_c152', '05_p005_c1'])
#names = np.array(['01_p005_c34m', '02_p005_c3m', '03_p005_c252m', '04_p005_c2m', '06_p005_c152m', '05_p005_c1m'])
#names = np.array(['01_p005_c34m', '03_p005_c252m', '06_p005_c152m', '05_p005_c1m'])
#names = np.array(['08_p0_D13', '08_p052_D13', '08_p1_D13', '08_p152_D13', '08_p2_D13', '08_p34_D13'])
#names = np.array(['08_p052_D13', '08_p1_D13', '08_p152_D13', '08_p252_D13', '08_p34_D13'])
#names = np.array(['08_p1_D13', '08_p152_D13', '08_p252_D13', '08_p34_D13'])

#names = np.array('01_pco2_0.03', '01_pco2_0.1', '01_pco2_0.3', '01_pco2_1.0', '01_pco2_3.0','01_pco2_10.0', '01_pco2_30.0','01_pco2_100.0')
#label = np.array(['0.03%', '0.1%','0.3%', '1%','3%', '10%', '30%', '100%'])

names = np.array(['01_pco2_0.03', '01_pco2_0.3', '01_pco2_3.0','01_pco2_10.0', '01_pco2_30.0'])
label = np.array(['0.04%', '0.3%', '3%', '10%', '30%'])
linetype = np.array(['-', '--', '-.', ':', '-', '--','-', '--', '-.', ':', '-', '--'])

scale = 50

#%%
porosity = {}
solid = {}
de = {}
    
for i in range(0, len(names)): 
    porosity[names[i]] = np.load(root_dir+'\\results\\output\\07_pco2\\' + names[i] + '\\' + 'poros.npy')[1,1:-1]
    de[names[i]] = np.load(root_dir+'\\results\\output\\07_pco2\\' + names[i] + '\\' + 'De.npy')[1,1:-1]
    idx = np.concatenate((np.where(de[names[i]] == 1e-15)[0], np.where(de[names[i]] > 1e-10)[0]), axis = None)
    solid[names[i]] = 1 - porosity[names[i]]
    solid[names[i]][idx] *= -1

fig = plt.figure(figsize=(8,4))
axes = [fig.add_subplot(5, 1, i+1) for i in range(5)]
#fig, axes = plt.subplots(5, 1, constrained_layout=True)
for i in range(0, len(names)): 
    axes[i].imshow(np.expand_dims(solid[names[i]], axis=0),cmap=plt.get_cmap('RdBu'), vmin=-1.0, vmax=1.0)
    
for ax in plt.gcf().axes:
    try:
        ax.label_outer()
    except:
        pass
    
for ax in axes:
    ax.set_yticks([])
    #ax.set_aspect('equal')
    
plt.sca(axes[4])
#plt.xticks([0, 5, 10, 15, 20, 25],[1, 6, 11, 16, 21, 26])
plt.xticks([0, 5, 10, 15, 20, 25],[2, 7, 12, 17, 22, 27])
plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()
plt.show()   

#%%
np.savetxt('solid.txt', 
           np.transpose([solid[names[0]], solid[names[1]], solid[names[2]], solid[names[3]], solid[names[4]]]), 
           )
#%%
s = solid[names[4]]
s[np.where(s<0)] = 0
plt.figure(figsize = (8,4))
plt.imshow(np.expand_dims(solid[names[4]], axis=0),cmap=plt.get_cmap('Blues'), vmin=0, vmax=1.0) #Reds
plt.colorbar()
plt.show()   
#%%
    
plt.figure(figsize = (8,4))
for i in range(0, len(names)): 
    pH = np.load(root_dir+'\\results\\output\\07_pco2\\' + names[i] + '\\' + 'pH.npy')
    plt.plot(pH[1,1:-2], label=label[i])
plt.xlabel(r'Distance ($\mu$m)', fontsize = 14)
plt.ylabel('pH', fontsize = 14)
plt.legend(loc= "right", fontsize = 12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()
