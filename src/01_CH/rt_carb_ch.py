# -*- coding: utf-8 -*-

import sys, os
src_dir = os.path.dirname(__file__)
sys.path.append(src_dir)
import numpy as np
import matplotlib.pylab as plt
from copy import deepcopy
import yantra
import cell_type as ct # change the path to cell_type file
import defaults as df
from rt_carb import CarbonationRT
from phrqc_input import PhreeqcInput
from settings import Settings, Results
from yantra._pyevtk.hl  import imageToVTK

class PhreeqcInputCH(PhreeqcInput):
    pass
    

class CarbonationCH(CarbonationRT):  
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        super(CarbonationCH, self).__init__(eqn,domain,domain_params,bc_params,solver_params, settings)
        self.ptype ='CH'
     
    def update_bc(self, settings):
        for key in self.phrqc.boundary_solution_labels:  
            self.fluid.Ca.bc[key] = ['flux', 0.0]
            self.fluid.Ca._bc[key+'bc'] = 'flux'            
            if (settings['bc']['type'] == 'pco2'):                
                self.fluid.C.bc[key] = ['flux', 0.0]
                self.fluid.C._bc[key+'bc'] = 'flux'
                self.fluid.H.bc[key] = ['flux', 0.0]
                self.fluid.H._bc[key+'bc'] = 'flux'
                self.fluid.O.bc[key] = ['flux', 0.0]
                self.fluid.O._bc[key+'bc'] = 'flux'
    
    def update_init_phrqc(self):
        self.phrqc.boundcells = self.update_boundary_cells()
        if 'portlandite' in self.solid.diffusive_phase_list:
            self.phrqc.init_port = self.solid.portlandite.c>0
        if 'calcite' in self.solid.diffusive_phase_list:
            self.phrqc.is_calc = self.solid.calcite.c>0
        self.phrqc._target_SI = np.zeros(self.solid.shape)  
        self.phrqc.nodetype = self.nodetype 
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty) 
            
            
    def update_source(self):                
        c=self.fluid.get_attr('c')
        phaseqty=self.solid.phaseqty_for_phrqc()
        phase_list = self.solid.diffusive_phase_list
        for num, phase in enumerate(phase_list, start=1):
            phaseqty[phase] = deepcopy(self.solid._diffusive_phaseqty[num-1])  
        if(self.settings['dissolution']=='subgrid'):
            phaseqty['portlandite'][np.where(self.solid.border)]=0
            #phaseqty['portlandite'] = np.zeros(phaseqty['portlandite'].shape)
        self.phrqc.modify_solid_phases(phaseqty)        
        ss = self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)  
        if (self.settings['dissolution']=='multilevel'): 
            pqty=self.solid.update(self.phrqc.dphases)  
        else:
            if self.iters >1 :pqty=self.solid.update(self.phrqc.dphases)  
            self.update_solid_params() 
            ss=self.update_border_solution(c,ss)
                  
        self.update_solid_params()  
        self.update_target_SI() 
        self.update_phrqc()  
        ss = self.update_no_flux(ss)
        self.fluid.set_attr('ss',ss)  
        
    def update_nodetype(self):
        '''
        find neighbous of the multilevel cells
        '''
        prev_nodetype = deepcopy(self.nodetype) 
        
        self.solid.is_port = (self.solid.portlandite.c>0)
        self.solid.is_calc = (self.solid.calcite.c>0)
        self.solid.is_solid = (self.solid.nodetype == ct.Type.SOLID)
        self.solid.is_critical = (self.solid.nodetype == ct.Type.SOLID)
        is_port = self.solid.is_port
        is_calc = self.solid.is_calc
        #is_clinker = self.solid.nodetype == ct.Type.CLINKER
        is_solid = self.solid.is_solid
        is_critical = np.zeros(np.shape(is_port))
        is_interface = np.zeros(np.shape(is_port))
        #calc_c = self.phrqc.solid_phase_conc['calcite']
        is_liquid = (~is_port)
        if(self.iters>1):
            is_critical = (self.solid.pcs.pore_size <= 
                           self.solid.pcs.threshold_pore_size) & is_calc & (~is_port)
            is_liquid =  (~is_critical)&(~is_port)&(~is_solid)&(~is_calc)#&((prev_nodetype==-1)|(prev_nodetype==-2))
            is_interface = (~is_critical)&(~is_port)&(~is_solid)&(~is_liquid)
        self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
            ct.Type.MULTILEVEL * is_critical + \
            ct.Type.MULTILEVEL * is_port +\
            ct.Type.INTERFACE * is_interface + \
            ct.Type.SOLID * is_solid
        self.solid.is_liquid = is_liquid
        self.solid.is_interface = is_interface
        self.solid.is_critical = is_critical 
        yantra._solvers.update2d.reassign_mlvl(self.solid.nodetype) 
        self.solid.nodetype[prev_nodetype == ct.Type.SOLID] = ct.Type.SOLID
        yantra._solvers.update2d.reassign_mlvl_solid(self.solid.nodetype) 
        #self.solid.prev_calc_c = deepcopy(self.phrqc.solid_phase_conc['calcite']) 
        self.nodetype = self.solid.nodetype
        self.update_border()
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)  
        self.solid.prev_nodetype = prev_nodetype
        
        
    def update_diffusivity(self):
        Dref = self.Dref        
        cc = self.settings['diffusivity']['CC']
        ch = self.settings['diffusivity']['CH']
        
        is_border = self.solid.border
        is_port = (self.solid.portlandite.c >0) & (~is_border)
        is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
        is_liquid = np.logical_and(~is_port, ~is_calc)
        is_liquid = np.logical_and(is_liquid, ~is_border)
        
        D_archie = Dref * self.solid.poros * self.solid.app_tort 
        D_CC = D_archie
        D_CH = D_archie
        if(cc[0] == 'const'):
            D_CC = cc[1]*np.ones(D_CC.shape)
        elif(cc[0] == 'inverse'):
            mineral = self.solid.calcite.c * self.solid.calcite.mvol
            #D_CC =np.nan_to_num(1./((1-mineral)/Dref/self.solid.poros/self.solid.app_tort + mineral/cc[1]), Dref)
            D_CC =np.nan_to_num(1./((1-mineral)/Dref + mineral/cc[1]), Dref)
        
        if(ch[0] == 'const'):
            D_CH = ch[1]*np.ones(D_CH.shape)
        elif(ch[0] == 'inverse'):            
            mineral = self.solid.vol
            D_CH =np.nan_to_num(1./((1-mineral)/Dref/self.solid.poros/self.solid.app_tort + mineral/ch[1]), Dref)
        
        De = D_CH*is_port + D_CC*is_calc + D_CC*is_border + Dref*is_liquid 
        De = De/self.solid.poros/self.solid.app_tort 
            
        self.fluid.set_attr('D0',De,component_dict=False)
        self.fluid.set_attr('Deref',np.max(De),component_dict=False)
        self.fluid.call("_set_relaxation_params")  
        
    def update_no_flux(self, ss):        
        ss['Ca'] = ss['Ca']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1) 
        return(ss)
        
    def update_border(self):
        is_port = (self.solid.portlandite.c>0)
        is_mineral = is_port | (self.solid.nodetype==ct.Type.SOLID)
        val = -16
        temp = deepcopy(self.solid.nodetype) 
        temp = temp*(~is_mineral) + val*is_mineral        
        rolled = np.roll(temp, -1, axis = 1)/4+np.roll(temp, 1, axis = 1)/4 +\
              np.roll(temp, -1, axis = 0)/4+np.roll(temp, 1, axis = 0)/4
        border = (rolled!=val)&is_port
        self.solid.border = border 
            
    def update_border_solution(self,c,ss):
        #poros = self.solid.poros#
        phrqc_poros = self.solid.poros#self.phrqc.selected_output()['poros'] #
        fraction = self.settings['subgrid']['fraction']
        result = {}
        by = np.where(self.solid.border)[0]
        bx = np.where(self.solid.border)[1]
        bf = np.where(self.solid.border.flatten())[0]
        for i in np.arange(0, np.sum(self.solid.border)):          
            if fraction >0:    
                cell_m = bf[i]+1                
                result  = self.update_equilibrium(result, cell_m,  
                                          self.solid.portlandite.c[by[i], bx[i]],
                                          phrqc_poros[by[i], bx[i]],
                                          fraction)
                ssnew={}
                for name in self.phrqc.components:
                    ssnew[name] = (result[str(cell_m)][name]-c[name][by[i], bx[i]])/self.dt
                    ssnew[name] *= phrqc_poros[by[i], bx[i]]
                    ss[name][by[i], bx[i]] = ssnew[name]
        
        for i in np.arange(0, np.sum(self.solid.border)):
            self.solid.portlandite.c[by[i], bx[i]] = result[str(bf[i]+1)]['portlandite_m']
            #self.solid.calcite.c[by[i], bx[i]] = result[str(bf[i]+1)]['calcite_m']
        return ss
    
    def update_equilibrium(self, result, n_ch, m_ch, porosity, fraction=1):
        ncell = 123456
        fract = fraction/porosity
        if fract>1: fract =1
        #print(fract)
        modify_str = []
        modify_str.append("EQUILIBRIUM_PHASES %i" % ncell)
        modify_str.append("Portlandite 0 %.20e dissolve only" %(m_ch))   
        modify_str.append("END") 
        modify_str.append('MIX %i' % ncell)         
        modify_str.append('%i %.20e' %( n_ch, fract)) #modify_str.append('%i 1' %n_int)  
        modify_str.append('SAVE solution %i' % ncell)  
        modify_str.append("END") 
        modify_str.append('USE solution %i' % ncell)  
        modify_str.append('USE equilibrium_phase %i' % ncell)
        modify_str.append('SAVE solution %i' % ncell)   
        modify_str.append("END") 
        modify_str.append("END") 
        modify_str ='\n'.join(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str) 
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        #print(modify_str)
        #print(output)
        
        port = output[2][10] 
        modify_str = [] 
        modify_str.append("EQUILIBRIUM_PHASES %i" %n_ch)
        modify_str.append("Portlandite 0 %.20e" %(0)) # dissolve only
        modify_str.append("END") 
        modify_str.append('USE equilibrium_phase %i' %n_ch)      
        modify_str.append('MIX %i' % ncell)      
        modify_str.append('%i 1' % ncell)   
        modify_str.append('%i %.20e' %(n_ch, 1.-fract)) #modify_str.append('%i 0' %n_int)   
        modify_str.append('SAVE solution %i' %(n_ch))  
        modify_str.append('SAVE equilibrium_phase %i' %(n_ch))  
        modify_str.append("END") 
        modify_str ='\n'.join(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str)  
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        #print(modify_str)
        #print(output)
        comp = {}        
        comp['portlandite_m'] = port
        comp['C'] = output[1][6]
        comp['Ca'] = output[1][7]
        comp['H'] = (output[1][8] - self.phrqc.H_norm)
        comp['O'] = (output[1][9] - self.phrqc.O_norm)
        result[str(n_ch)] = comp        
        return(result)     
        
    def update_volume(self):
        vol = np.zeros(self.solid.shape)
        phase_list = self.solid.diffusive_phase_list
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            vol += val.c * self.solid.mvol[num-1] 
        self.solid.vol = vol
        self.solid.portlandite.vol = self.volume_CH()
        self.solid.calcite.vol = self.volume_CC() 
        self.solid.dissolving_mineral_vol = self.solid.portlandite.vol
        
    def volume_CH(self):
        CH_vol = self.solid.portlandite.c * self.solid.portlandite.mvol
        return CH_vol
    
    def volume_CC(self):
        CC_vol = self.solid.calcite.c * self.solid.calcite.mvol
        return CC_vol
       
    def free_volume(self):
        v = self.solid.voxel_vol - self.solid.vol_ch
        return v

    
    def csh_conc(self):    
        csh = self.solid.CSHQ_TobH.c+ self.solid.CSHQ_TobD.c +\
                self.solid.CSHQ_JenH.c + self.solid.CSHQ_JenD.c
        return csh 
    
#%% Additional functions
class SettingsCH(Settings):
    
    @staticmethod
    def set_mvols(scale=50):
        '''
        idx 0 - portlandite,
        idx 1 - calcite
        '''
        mvolCH = df.MOLAR_VOLUME['CH']*scale
        mvolCC = df.MOLAR_VOLUME['CC']*scale
        return  [mvolCH, mvolCC]
    
    @staticmethod
    def set_init_pqty(mvol, poros, scale=50): 
        max_pqty = Settings.get_max_pqty(mvol)
        init_conc = [0.0]*len(mvol)
        for i in range(len(mvol)):
            init_conc[i] = (1 - poros[i]) * max_pqty[i]
        return init_conc
        
    @staticmethod
    def get_pqty(init_pqty, domain):
        pqty_CH = init_pqty[0] * (domain.nodetype == ct.Type.MULTILEVEL) + \
            init_pqty[0] * (domain.nodetype == ct.Type.MULTILEVEL_CH)
        pqty_CC = init_pqty[1] * np.ones(domain.nodetype.shape)
        return([pqty_CH, pqty_CC] )
        
    @staticmethod
    def set_labels(domain):
        slabels = np.zeros(domain.nodetype.shape)
        slabels = 100002 * (domain.nodetype == ct.Type.LIQUID) + \
                  100002 * (domain.nodetype == ct.Type.INTERFACE) + \
                  100003 * (domain.nodetype == ct.Type.MULTILEVEL)  + \
                  100005 * (domain.nodetype == ct.Type.SOLID)
        return(slabels)
    
    @staticmethod
    def get_porosity(domain, pqty, mvol):
        poros = np.zeros(domain.nodetype.shape)
        poros = (1- pqty[0]*mvol[0])*(domain.nodetype == ct.Type.MULTILEVEL) + \
                1*(domain.nodetype != ct.Type.MULTILEVEL)
        return(poros)
    
    @staticmethod
    def set_domain_params(D, mvol, pqty, poros, app_tort, slabels, 
                          input_file, database = 'cemdata18.dat' ):
        dp={}
        
        dp['D0'] = D                     
        dp['voxel_vol']=1
        dp['mvol']= mvol        
        dp['poros'] = poros                 
        dp['app_tort'] = app_tort         
        dp['solution_labels'] =slabels
        dp['database'] = database
        
        dp['phrqc_input_file']=input_file#'CH_CC_Ceq.phrq'
        dp['eq_names'] = ['portlandite', 'calcite']
        dp['solid_phases']={'portlandite': {'type':'diffusive','mvol':mvol[0],'c':pqty[0]},
                            'calcite':     {'type':'diffusive','mvol':mvol[1],'c':pqty[1]}}  
        return dp


class ResultsCH(Results):
    
    def init_results(self, nodes=[]):
        '''
        pavg_list - list of average values to save 
        nodes - list of parameters per node to save 
        '''
        results = {}     
        params = [] 
        params += ['portlandite', 'calcite', 'Ca', 'C', 'O', 'H'] # always save these parameters
        results={name: [] for name in params}       
        results['params'] = params
        pavg_list = ['tot_vol', 'mean_D_eff', 'mean_poros', 
                          'dissolution_CH', 'nodes_CH', 'active_nodes', 
                          'mean_pH', 'dCH', 'dt',
                          'precipitation_CC', 'nodes_CC', 'dCC']  
        results.update({name: [] for name in pavg_list})
        results['pavg_list'] = pavg_list
        if nodes: 
            pointparamslist = []
            pointparamslist += params
            pointparamslist += ['vol', 'poros', 'pH', 'De', 'vol_CH', 'vol_CC'] 
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
            self.results['nodes_CH'].append(self.__class__.CH_nodes(rtmodel))
            self.results['nodes_CC'].append(self.__class__.CC_nodes(rtmodel))
            self.results['active_nodes'].append(self.__class__.active_nodes(rtmodel))
            self.results['dt'].append(self.__class__.get_dt(rtmodel))
            self.results['mean_pH'].append(self.__class__.mean_pH(rtmodel))
            self.results['dissolution_CH'].append(self.__class__.CH_dissolution(rtmodel))
            self.results['precipitation_CC'].append(self.__class__.CC_precipitation(rtmodel))
            if (rtmodel.iters ==0): 
                self.results['dCH'].append(0) 
                self.results['dCC'].append(0) 
            else:
                self.results['dCH'].append(self.__class__.delta_CH(rtmodel))
                self.results['dCC'].append(self.__class__.delta_CC(rtmodel))
            
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
                    self.results['calcite'+' ' + str(p)].append(rtmodel.solid.calcite.c[p])
                    self.results['C'+' ' + str(p)].append(rtmodel.fluid.C._c[p]+\
                                 rtmodel.fluid.C._ss[p]/rtmodel.fluid.C.poros[p])
                    self.results['vol_CC'+' ' + str(p)].append(rtmodel.solid.calcite.vol[p])
    @staticmethod
    def delta_CC(rt):
        '''
        Differentiation of Portladite
        '''
        s = np.sum(rt.phrqc.dphases['calcite'])#*getattr(rt.solid, 'calcite').mvol)
        return(s)
    
    @staticmethod
    def CC_precipitation(rt):
        '''
        N cells that are dissolving
        '''
        s = np.sum(rt.phrqc.dphases['calcite']>0)
        return(s)
         
    @staticmethod
    def CC_nodes(rt):
        '''
        N cells containing portlandite
        '''
        s = np.sum(rt.solid.calcite.c>0)
        return(s)
    #%% PRINT and PLOT
    
    @staticmethod
    def print_profiles(rt):    
        Results.print_profiles(rt)
        print('C +ss %s' %str(np.array(rt.fluid.C.c[1,:]) + \
               np.array(rt.fluid.C._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
        print('CC %s' %str(np.array(rt.solid.calcite.c[1,:])))
        print('dCC %s' %str(np.array(rt.phrqc.dphases['calcite'][1,:])))
        
    def get_titles(self):
        title = {'Ca': 'Calcium',
                 'C': 'Carbon',
                 'O': 'Oxygen',
                 'H': 'Hydrogen',
                 'portlandite':'Portlandite',
                 'calcite': 'Calcite',
                 'De':  'Effective diffusivity',
                 'poros': 'Porosity',
                 'vol_CH': 'Portlandite volume',
                 'vol_CC': 'Calcite volume',
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
                 'tot_vol':'Total mineral volume',
                 'dCC': 'Calcite differentiation',
                 'precipitation_CC':'Number of precipitating nodes',
                 'nodes_CC':'Number of Calcite nodes',
                 }
        return(title)
         
    def get_ylabs(self):
        ylab = {'Ca':'Ca (*1e-12) [mol]',
                'C':'C (*1e-12) [mol]',
                'O':'O (*1e-12) [mol]',
                'H':'H (*1e-12) [mol]',
                'portlandite': 'CH (*1e-12) [mol]',
                'calcite': 'CC (*1e-12) [mol]',
                'De':  'D_eff [m^2/s]',
                'poros': 'Porosity [-]',
                'vol_CH': 'Portlandite volume [um^3]',
                'vol_CC': 'Calcite volume [um^3]',
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
                'dCC':'dCC [mol]',
                'precipitation_CC':'Nodes [-]',
                'nodes_CC':'Nodes [-]',
                }
        return(ylab)
    
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
                  'calcite': rt.solid.calcite.c[1:-1,1:-1],
                  'Ca': rt.fluid.Ca.c[1:-1,0:-1],
                  'C': rt.fluid.C.c[1:-1,0:-1],
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
    def save_vti(rt, phases, t, path, name):
        
        nx = rt.fluid.Ca.nx -1
        ny = rt.fluid.Ca.ny -1     
        filename = path + str(name) +'_all_' + str(t) 
        
        cCa = rt.fluid.Ca.c[1:ny,1:nx,np.newaxis]
        cH = rt.fluid.H.c[1:ny,1:nx,np.newaxis]
        cO = rt.fluid.O.c[1:ny,1:nx,np.newaxis]
        cC = rt.fluid.C.c[1:ny,1:nx,np.newaxis]
        cCC = rt.solid.calcite.c[1:ny,1:nx,np.newaxis]
        cCH = rt.solid.portlandite.c[1:ny,1:nx,np.newaxis]
        porosity = rt.solid._poros[1:ny,1:nx,np.newaxis]
        
        imageToVTK(filename, cellData = {"Ca":cCa,
                                         "C":cC,
                                         "O":cO,
                                         "H":cH,
                                         "Calcite":cCC,
                                         "Portlandite":cCH,
                                         "porosity":porosity,
                                         })