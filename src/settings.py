import os
import json
import pickle
import numpy as np
import matplotlib.pylab as plt
import cell_type as ct 

from yantra._pyevtk.hl  import pointsToVTK
from yantra._pyevtk.hl  import gridToVTK
from yantra._pyevtk.hl  import imageToVTK

#%% SETTINGS
class Settings():
    @staticmethod
    def set_mvols(scale=50):
        return([])
        
    @staticmethod
    def get_max_pqty(mvol):
        max_pqty = map(lambda x: 1./x, mvol)
        return(max_pqty)
        
    @staticmethod
    def set_init_pqty(mvol, poros, scale): 
        '''
        mvol - molar volumes of all minerals in the model
        poros - initial porosity of all minerals in the model 
        '''
        pass
        
    @staticmethod
    def get_pqty(init_pq, domain): 
        return([])
    
    @staticmethod
    def set_labels(domain):
        slabels = np.zeros(domain.nodetype.shape)
        return(slabels)
    
    @staticmethod
    def get_porosity(domain, pqty, mvol):
        poros = np.zeros(domain.nodetype.shape)
        return(poros)
        
    @staticmethod
    def get_app_tort(domain,degree, porosity):
        app_tort = np.zeros(domain.nodetype.shape)
        app_tort = porosity ** degree
        return(app_tort)
    
    @staticmethod
    def set_domain_params(D, mvol, pqty, poros, app_tort, slabels, input_file):
        return({})
        
    @staticmethod
    def reset_portlandite_nodes(domain):
        domain.nodetype[domain.nodetype == ct.Type.MULTILEVEL_CH] = ct.Type.MULTILEVEL
        return(domain)
        
    
    @staticmethod
    def set_solver_params(tfact = None, smart_thres = 1e-8, 
                          cphi_fact = 1./3., cphi = 0):
        sp={}
        sp['collision_model']= 'trt' #'diff_vel' #
        sp['magic_para']=1.0/4.0
        sp['cphi_fact']= cphi_fact
        if cphi>0:        
            sp['cphi']= cphi
        
        sp['phrqc_flags'] = {}
        sp['phrqc_flags']['smart_run']=True
        sp['phrqc_smart_run_tol']=smart_thres
        if(tfact):
            sp['tfactbased']=1
            sp['tfact']= tfact
        return sp
    
    @staticmethod
    def set_bc_params(bc_slabels):
        bcp = {'solution_labels':bc_slabels, 
               'top':['flux', 0.0],
               'bottom':['flux', 0.0],
               'left':['flux', 0.0],
               'right':['flux', 0.0],}
        return bcp

    @staticmethod
    def save_settings(settings, bc_params, solver_params, path, name):
        def write_txt(settings, bc_params, solver_params, path, name):
            with open(path + name +'_settings.txt', 'w') as file:
                file.write(json.dumps(settings))
                file.write('\n\n')
                file.write(json.dumps(bc_params))
                file.write('\n\n')
                file.write(json.dumps(solver_params))
            file.close()
        try:
            write_txt(settings, bc_params, solver_params, path, name)
            #os.mkdir(dirName)   
        except IOError:
            os.mkdir(path)
            write_txt(settings, bc_params, solver_params, path, name)
            
#%% RESULTS
class Results((object)):
    
    def __init__(self, nodes=[]):
        self.results = {}
        self.init_results(nodes)
        
    def init_results(self, nodes=[]):
        '''
        pavg_list - list of average values to save 
        nodes - list of parameters per node to save 
        '''
        results = {}     
        params = [] 
        params += ['portlandite', 'Ca', 'O', 'H'] # always save these parameters
        results={name: [] for name in params}       
        results['params'] = params
        pavg_list = ['tot_vol', 'mean_D_eff', 'mean_poros', 
                          'dissolution_CH', 'nodes_CH', 'active_nodes', 
                          'mean_pH', 'dCH', 'dt']  
        results.update({name: [] for name in pavg_list})
        results['pavg_list'] = pavg_list
        if nodes: 
            pointparamslist = []
            pointparamslist += params
            pointparamslist += ['vol', 'poros', 'pH', 'De', 'vol_CH'] 
            l = [[m+ ' '+n for n in [str(n) for n in nodes]] for m in [str(m) for m in pointparamslist]]
            par_points = [item for sublist in l for item in sublist] #sum(l, [])
            results.update({name: [] for name in par_points}) 
            results['pointparamslist'] = pointparamslist
        results['nodes'] = nodes    
        results['time'] = []   
        self.results = results
    
    def append_results(self, rtmodel, step = 1e+2):
        if (rtmodel.iters%step == 0):
            self.results['time'].append(rtmodel.time)
            for num, phase in enumerate(rtmodel.solid.diffusive_phase_list, start=1):
                self.results[phase].append(np.sum(rtmodel.solid._diffusive_phaseqty[num-1]))        
            for num, comp in enumerate(rtmodel.fluid.components, start=1):
                self.results[comp].append(np.mean(getattr(rtmodel.fluid, comp)._c*getattr(rtmodel.fluid, comp).poros))  
            
            #general
            self.results['tot_vol'].append(self.__class__.total_mineral_volume(rtmodel))
            self.results['mean_D_eff'].append(self.__class__.mean_effective_D(rtmodel))
            self.results['mean_poros'].append(self.__class__.mean_porosity(rtmodel))
            self.results['dissolution_CH'].append(self.__class__.CH_dissolution(rtmodel))
            self.results['nodes_CH'].append(self.__class__.CH_nodes(rtmodel))
            self.results['active_nodes'].append(self.__class__.active_nodes(rtmodel))
            self.results['dt'].append(self.__class__.get_dt(rtmodel))
            self.results['mean_pH'].append(self.__class__.mean_pH(rtmodel))
            self.results['dCH'].append(self.__class__.delta_CH(rtmodel))
            if (rtmodel.iters ==0): 
                self.results['dCH'].append(0) 
            else:
                self.results['dCH'].append(self.__class__.CH_dissolution(rtmodel))
            
            # nodes
            if self.results['nodes']: # points
                for p in self.results['nodes']:
                    self.results['portlandite'+' ' + str(p)].append(rtmodel.solid.portlandite.c[p])
                    self.results['Ca'+' ' + str(p)].append(rtmodel.fluid.Ca._c[p]+\
                                 rtmodel.fluid.Ca._ss[p]/rtmodel.fluid.Ca.poros[p])
                    self.results['H'+' ' + str(p)].append(rtmodel.fluid.H._c[p]+\
                                 rtmodel.fluid.H._ss[p]/rtmodel.fluid.H.poros[p])
                    self.results['O'+' ' + str(p)].append(rtmodel.fluid.O._c[p]+\
                                 rtmodel.fluid.O._ss[p]/rtmodel.fluid.O.poros[p])
                    self.results['poros'+' ' + str(p)].append(rtmodel.solid.poros[p])
                    self.results['vol'+' ' + str(p)].append(rtmodel.solid.vol[p])
                    self.results['De'+' ' + str(p)].append(rtmodel.fluid.H.De[p])
                    self.results['vol_CH'+' ' + str(p)].append(rtmodel.solid.portlandite.vol[p])
                    self.results['pH'+' ' + str(p)].append(rtmodel.phrqc.selected_output()['pH'][p])
                
    def filter_results(self, path, name, length = 1e+4):
        l = len(self.results['time'])
        r = l/length
        filtered_results ={}
        if (r <= 1):
            filtered_results = self.results
        else:
            filt = int(r)
            filtered_results = {k: v[0::filt] for k, v in self.results.items()}
        #save_obj(filtered_results, path + str(name) +'_results')
        return filtered_results
    
    @staticmethod
    def total_mineral_volume(rt):
        '''
        Total mineral volume
        '''
        #vol = np.sum(rt.solid.vol)
        not_boundary = (rt.phrqc.boundcells !=1)
        vol = np.sum(rt.solid.vol[not_boundary])
        return(vol)
        
    @staticmethod
    def mean_effective_D(rt):
        '''
        Avarage effective diffusivity calculated by Archie's formula
        '''
        De = rt.fluid.H.De
        not_boundary = (rt.phrqc.boundcells !=1)  
        return np.mean(De[not_boundary])
    
    @staticmethod
    def mean_porosity(rt):
        '''
        Average poroisty in the domain
        '''
        #poros = np.mean(rt.solid.poros)
        not_boundary = (rt.phrqc.boundcells !=1)
        poros = np.mean(rt.solid.poros[not_boundary])
        return(poros)
            
    @staticmethod
    def delta_CH(rt):
        '''
        Differentiation of Portladite
        '''
        s = np.sum(rt.phrqc.dphases['portlandite'])#*getattr(rt.solid, 'portlandite').mvol)
        return(s)
    
    @staticmethod
    def CH_dissolution(rt):
        '''
        N cells that are dissolving
        '''
        s = np.sum(rt.phrqc.dphases['portlandite']<0)
        return(s)
         
    @staticmethod
    def CH_nodes(rt):
        '''
        N cells containing portlandite
        '''
        s = np.sum(rt.solid.portlandite.c>0)
        return(s)
          
    @staticmethod      
    def mean_pH(rt):
        '''
        N cells containing calcite
        '''
        nx = rt.fluid.Ca.nx -1
        ny = rt.fluid.Ca.ny -1
        s = np.mean(rt.phrqc.selected_output()['pH'][1:ny,1:nx])
        return(s)
    
    @staticmethod
    def active_nodes(rt):
        '''
        N active cells in phreeqc
        '''
        a = rt.phrqc.nactive
        return(a)
        
    @staticmethod
    def get_dt(rt):
        return(rt.dt)
                      
            
    #%% PLOT FUNCTIONS
      
    def get_titles(self):
        title ={}
        title = {'Ca': 'Calcium',
                 'O': 'Oxygen',
                 'H': 'Hydrogen',
                 'portlandite':'Portlandite',
                 'De':  'Effective diffusivity',
                 'poros': 'Porosity',
                 'vol_CH': 'Portlandite volume',
                 'vol':'Mineral volume',
                 'dt':'Time step',
                 'pH': 'pH',
                 'active_nodes':'Number of active nodes (PHREEQC)',
                 'mean_poros':'Average porosity',
                 'mean_D_eff':'Average effective diffusivity',
                 'mean_pH': 'Average pH ',
                 'dCH': 'Portlandite differentiation',
                 'dissolution_CH':'Number of dissolving nodes',
                 'nodes_CH':'Number of Portlandite nodes',
                 'tot_vol':'Total mineral volume'
                 }
        return(title)
         
    def get_ylabs(self):
        ylab = {}
        ylab = {'Ca':'Ca (*1e-12) [mol]',
                'O':'O (*1e-12) [mol]',
                'H':'H (*1e-12) [mol]',
                'portlandite': 'CH (*1e-12) [mol]',
                'De':  'D_eff [m^2/s]',
                'poros': 'Porosity [-]',
                'vol_CH': 'Portlandite volume [um^3]',
                'vol':'Mineral volume [um^3]',
                'dt':'dt [s]',
                'pH': 'pH [-]',                
                'active_nodes':'Nodes [-]',
                'mean_poros':'Average porosity [-]',
                'mean_D_eff':'Average D_eff [m^2/s]',
                'mean_pH': 'Average  pH [-]',              
                'dCH':'dCH [mol]',
                'dissolution_CH':'Nodes [-]',
                'nodes_CH':'Nodes [-]',
                'tot_vol': 'Total mineral volume [um^3]',
                }        
        return(ylab)
        
    def plot_species(self, names=[], fsize=(8,4)):
        ylab = self.get_ylabs()
        title = self.get_titles()
        if not names:
            names = self.results['params']
        for k in names:        
            plt.figure(figsize=fsize)
            plt.plot(self.results['time'], self.results[k])
            plt.legend()
            plt.title(title[k])
            plt.ylabel(ylab[k])
            plt.xlabel('Time [s]')
            plt.show()
            
    def plot_avg(self, fsize=(8,4)):
        ylab = self.get_ylabs()
        title = self.get_titles()    
        names = self.results['pavg_list']
        for k in names:        
            plt.figure(figsize=fsize)
            plt.plot(self.results['time'], self.results[k])
            plt.legend()
            plt.title(title[k])
            plt.ylabel(ylab[k])
            plt.xlabel('Time [s]')
            plt.show()
            
    def plot_nodes(self, names=[], fsize=(8,4)):
        ylab = self.get_ylabs()
        title = self.get_titles()
        if not self.results['nodes']:
            pass    
        else:
            if not names:
                names = self.results['pointparamslist']
            for k in names:
                plt.figure(figsize=fsize)        
                for p in self.results['nodes']:  
                    plt.plot(self.results['time'], self.results[k+' '+str(p)], label = str(p))
                plt.legend()
                plt.title(title[k] + ' in nodes ' + str(self.results['nodes']))
                plt.ylabel(ylab[k])
                plt.xlabel('Time [s]')
                plt.show()
    
    def plot_fields(self, rt, names={}, fsize = (8,4)):
        #nl = {name: limit}
        
        def plot(field, title, ylab, limit=[], size=fsize):
            plt.figure() 
            if not limit:
                plt.imshow(field)
            else:
                plt.imshow(field, vmin = limit[0], vmax = limit[1])
            
            clb = plt.colorbar()
            clb.ax.get_yaxis().labelpad = 15
            clb.ax.set_ylabel(ylab, rotation=270)
            
            plt.title(title)
            plt.show()
        
        fields = {'portlandite': rt.solid.portlandite.c[1:-1,1:-1],
                  'Ca': rt.fluid.Ca.c[1:-1,0:-1],
                  'H':  rt.fluid.H.c[1:-1,0:-1],
                  'O':  rt.fluid.O.c[1:-1,0:-1],
                  'poros': rt.solid.poros[1:-1,0:-1]}
            
        if not names:
            names = fields.keys()
        
        ylab = self.get_ylabs()
        title = self.get_titles()
        for k in names: 
            plot(fields[k], title[k], ylab[k])
        
    #%% SAVE
    
    def  make_output_dir(path):
        try:
            os.mkdir(path)
        except WindowsError:
            pass
                    
                
    def save_obj(obj, name ):
        with open(name + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
            f.close()
    
    def load_obj(name ):
        with open(name + '.pkl', 'rb') as f:
            return pickle.load(f)
        
    def save_pickle(rt, t, path, name):
        Results.save_obj(rt.phases, path + str(name) +'_nodetype_' + str(t))  
        Results.save_obj(rt.solid._poros, path + str(name) +'_porosity_' + str(t)) 
        
    def save_vti(rt, t, path, name):
        
        nx = rt.fluid.Ca.nx -1
        ny = rt.fluid.Ca.ny -1     
        filename = path + str(name) +'_all_' + str(t) 
        
        cCa = rt.fluid.Ca.c[1:ny,1:nx,np.newaxis]
        cH = rt.fluid.H.c[1:ny,1:nx,np.newaxis]
        cO = rt.fluid.O.c[1:ny,1:nx,np.newaxis]
        cCH = rt.solid.portlandite.c[1:ny,1:nx,np.newaxis]
        porosity = rt.solid._poros[1:ny,1:nx,np.newaxis]
        
        imageToVTK(filename, cellData = {"Ca":cCa,
                                         "O":cO,
                                         "H":cH,
                                         "Portlandite":cCH,
                                         "porosity":porosity,
                                         })
        
    def save_vts(domain, t, path, name):
        filename = path + str(name) +'_nodetype_' + str(t) #"./fname"
        x,y = domain.meshgrid()
        x = x[:,:,np.newaxis]
        y = y[:,:,np.newaxis]
        z = np.zeros(x.shape)
        p = domain.nodetype[:,:,np.newaxis]
        pointsToVTK(filename, x, y, z, data = {"Nodetype" : p})
        gridToVTK(filename, x, y, z, pointData = {"Nodetype" : p})
        
    
    #%% PRINT FUNCTIONS
    @staticmethod
    def print_time(st, rt):
        
        print('==========================')
        hours = int(st/3600)
        minutes = int((st/3600 - hours)*60)
        seconds = int(((st/3600 - hours)*60 - minutes)*60)
        print ('time :%s' %(rt.time))
        print ('iterations :%s' %(rt.iters))
        print ('real time taken for simulation:%s hours %s min %s sec' %(hours, minutes, seconds))
        print('==========================') 
       
    @staticmethod
    def print_nodes(rt, nodes, names=[]):    
        title = Results.get_titles()
        
        fields = {'portlandite': rt.solid.portlandite.c,
                  'Ca': rt.fluid.Ca.c,
                  'H':  rt.fluid.H.c,
                  'O':  rt.fluid.O.c,
                  'poros': rt.solid._poros,
                  'vol_CH': rt.solid.portlandite.vol}
        if not names:
            names = fields.keys()        
        for k in names: 
            print("%s: %s" % (title[k], fields[k][[i[0] for i in nodes], [i[1] for i in nodes]]))
    
    @staticmethod
    def print_profiles(rt):    
        print('Ca +ss %s' %str(np.array(rt.fluid.Ca.c[1,:]) + \
               np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
        print('H +ss %s' %str(np.array(rt.fluid.H.c[1,:]) + \
               np.array(rt.fluid.H._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
        print('O +ss %s' %str(np.array(rt.fluid.O.c[1,:]) + \
               np.array(rt.fluid.O._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
        print('CH %s' %str(np.array(rt.solid.portlandite.c[1,:])))
        print('dCH %s' %str(np.array(rt.phrqc.dphases['portlandite'][1,:])))
        print('Vol %s' %str(np.array(rt.solid.vol[1,:])))
        print('D %s' %str(np.array(rt.fluid.C.De[1,:])))
        print('SI %s' %str(np.array(rt.solid.target_SI[1,:])))
        print('pH %s' %str(np.array(rt.phrqc.selected_output()['pH'][1,:])))
        print('poros %s' %str(np.array(rt.solid.poros[1,:])))
        print('phrqc poros %s' %str(np.array(np.array(rt.phrqc.poros[1,:]))))