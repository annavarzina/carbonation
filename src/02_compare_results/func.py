# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt

#%% FUNCTIONS
def plot_results(res, names, xfield, yfield, labels, ltype,
                 title, xlab, ylab, fpath, fname, ftitle,
                 fsize = (8,4)):
    plt.figure(figsize=fsize)
    for i in range(0, len(names)):
        plt.plot(res[names[i]][xfield], res[names[i]][yfield], 
                 ls=ltype[i], label = labels[i])      
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.legend()
    plt.savefig(fpath + fname + ftitle)
    plt.show()
   
def get_CH_dissolution_time(res, T):
    c = np.array(res['portlandite']) <=0
    if np.any(c):
        t = np.array(res['time'])[c]
        return str(np.min(t)) + ' s'
    else:
        return '>' + str(T) + ' s'  
    
    
def get_porosity_val(res, t, dt):
    
    c = np.logical_and(np.array(res['time']) >= t, 
                   np.array(res['time']) < t+dt)
    if np.any(c):
        p = np.array(res['avg_poros'])[c]
        return p[0]
    else:
        return 'Error' #raise error
  
def get_val_at_time(res, field, t, dt):
    
    c = np.logical_and(np.array(res['time']) >= t, 
                   np.array(res['time']) < t+dt)
    if np.any(c):
        p = np.array(res[field])[c]
        return p[0]
    else:
        return 'Error' #raise error
    
def get_rate(res, dt):
    ch = np.array(res)
    dCH = np.zeros(ch.shape,np.float)
    dCH[0:-1] = np.diff(ch)/dt
    return dCH
