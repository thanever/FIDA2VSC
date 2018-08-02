# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 13:34:06 2016

@author: dell
"""
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import sys
import time

class NetworkData:

    def __init__(self, path_casedata):
        self.casedata_0 = loadmat(path_casedata)
        self._proc_casedata()   
        self.bus_num = self.casedata['bus'].shape[0]
        self.branch_num = self.casedata['branch'].shape[0]
        self._meter_info()    
        self._bra_bus_ind()
        self._make_y()
        self._runpf()
        self._get_x_true()
        self._get_z()
        self._R_W()
        
    def _proc_casedata(self):
        self.casedata = dict()
        
        self.casedata['bl_ac'] = self.casedata_0['bl_ac']
        self.casedata['bl_dc'] = self.casedata_0['bl_dc']
        self.casedata['bus_mc'] = self.casedata_0['bus_mc']
        self.casedata['branch_mc'] = self.casedata_0['branch_mc']
        self.casedata['busdc_mc'] = self.casedata_0['busdc_mc']
        self.casedata['Ybus'] = self.casedata_0['Ybus']
        self.casedata['Yf'] = self.casedata_0['Yf']
        self.casedata['Yt'] = self.casedata_0['Yt']
        self.casedata['baseMVA'] = self.casedata_0['baseMVA']
        self.casedata['bus'] = self.casedata_0['bus']
        self.casedata['gen'] = self.casedata_0['gen']
        self.casedata['branch'] = self.casedata_0['branch']
        self.casedata['baseMVAac'] = self.casedata_0['baseMVAac']
        self.casedata['baseMVAdc'] = self.casedata_0['baseMVAdc']
        self.casedata['pol'] = self.casedata_0['pol']
        self.casedata['busdc'] = self.casedata_0['busdc']
        self.casedata['convdc'] = self.casedata_0['convdc']
        self.casedata['branchdc'] = self.casedata_0['branchdc']

    def _make_y(self):
        """compute the bus admittance matrix and branch admittance, and also the some parameters of the dc"""
        self.bases = float(self.casedata['baseMVA'][0,0])
        self.bases_ac = float(self.casedata['baseMVAac'][0,0])
        self.bases_dc = float(self.casedata['baseMVAdc'][0,0])
        self.baseu_ac = self.casedata['convdc'][:,11]
        self.baseu_dc = self.casedata['busdc'][:,5]
        # bus admittance 
        self.Ybus = self.casedata['Ybus']
        # branch admittance | Ybra = [bus_f, bus_t, y_f, y_ft, y_t] |
        Yf = self.casedata['Yf']
        Yt = self.casedata['Yt']
        self.Ybra = np.c_[self.casedata['branch'][:,[0,1]],Yf.sum(axis=1), -Yf[np.arange(self.branch_num), map(int, self.bra_bus_all_ind[:,1])],Yt.sum(axis=1)]
        
        self.g_sc = (1/(self.casedata['convdc'][:,6] + self.casedata['convdc'][:,7]*1j + self.casedata['convdc'][:,9] + self.casedata['convdc'][:,10]*1j)).real
        self.b_sc = (1/(self.casedata['convdc'][:,6] + self.casedata['convdc'][:,7]*1j + self.casedata['convdc'][:,9] + self.casedata['convdc'][:,10]*1j)).imag
        self.r_dc = self.casedata['branchdc'][0,2]
        self.lossa = self.casedata['convdc'][:,16]/self.bases_ac
        self.lossb = self.casedata['convdc'][:,17]/self.baseu_ac
        self.losscr = self.casedata['convdc'][:,18]/((self.baseu_ac**2)/self.bases_ac)
        self.lossci = self.casedata['convdc'][:,19]/((self.baseu_ac**2)/self.bases_ac)
    
    def _find_k(self,k):
        """k-> ybra, b, i, j"""
        ybra = self.Ybra[k]
        #b = np.where(np.logical_and(self.bra_bus[:,0]==self.casedata['branch'][k,0], self.bra_bus[:,1]==self.casedata['branch'][k,1]))[0][0]
        i = np.where(self.casedata['bus'][:,0]==self.casedata['branch'][k,0])[0][0]
        j = np.where(self.casedata['bus'][:,0]==self.casedata['branch'][k,1])[0][0]
        return ybra, i, j
        
    def _bra_bus_ind(self):
        """ """
        #VSC-HVDC bus i and bus j, index number
        self.busdc_i = np.where(self.casedata['bus'][:,0]==self.casedata['busdc'][0,1])[0][0]  # bus i, corresponding with the first dcbus in buddc
        self.busdc_j = np.where(self.casedata['bus'][:,0]==self.casedata['busdc'][1,1])[0][0]  # bus j, corresponding with the second dcbus in buddc
        self.busdc_ij = np.array([self.busdc_i, self.busdc_j])
        # branch_bus - remove the shunt branch
        self.bra_bus_all = self.casedata['branch'][:,[0, 1]]
        bra_ind_same = list()
        for i in range(self.branch_num):
               for j in range(0, i):
                      if (self.bra_bus_all[i,:]==self.bra_bus_all[j,:]).all():
                             bra_ind_same.append(i)
        self.bra_bus = np.delete(self.bra_bus_all,bra_ind_same,axis=0)
        # index change
        self.bra_bus_ind = self.bra_bus.copy()
        self.bra_bus_all_ind = self.bra_bus_all.copy()
        for i in range(self.bra_bus_ind.shape[0]):
               for j in range(self.bra_bus_ind.shape[1]):
                      self.bra_bus_ind[i,j] = np.where(self.casedata['bus'][:,0]==self.bra_bus_ind[i,j])[0][0]
        for i in range(self.bra_bus_all_ind.shape[0]):
               for j in range(self.bra_bus_all_ind.shape[1]):
                      self.bra_bus_all_ind[i,j] = np.where(self.casedata['bus'][:,0]==self.bra_bus_all_ind[i,j])[0][0]

    def _runpf(self):
        self.casedata['bus_inj'] = np.c_[self.casedata['bus'][:,0],-self.casedata['bus'][:,2],-self.casedata['bus'][:,3]]
        for i in range(len(self.casedata['gen'])):
               i_gen = np.where(self.casedata['bus_inj'][:,0]==self.casedata['gen'][i,0])[0][0]  #index of gen in pf['bus']
               self.casedata['bus_inj'][i_gen, 1] = self.casedata['bus_inj'][i_gen, 1] + self.casedata['gen'][i,1]
               self.casedata['bus_inj'][i_gen, 2] = self.casedata['bus_inj'][i_gen, 2] + self.casedata['gen'][i,2]    

    def _meter_info(self):
        self.bus_mv = self.casedata['bus_mc'][:,1].nonzero()[0]
        self.bus_mp = self.casedata['bus_mc'][:,3].nonzero()[0]
        self.bus_mq = self.casedata['bus_mc'][:,5].nonzero()[0]
        self.branch_mfp = self.casedata['branch_mc'][:,2].nonzero()[0]
        self.branch_mfq = self.casedata['branch_mc'][:,4].nonzero()[0]
        self.branch_mtp = self.casedata['branch_mc'][:,6].nonzero()[0]
        self.branch_mtq = self.casedata['branch_mc'][:,8].nonzero()[0]
     
        self.busdc_mps    = self.casedata['busdc_mc'][:,2].nonzero()[0] 
        self.busdc_mqs    = self.casedata['busdc_mc'][:,4].nonzero()[0] 
        self.busdc_mpc    = self.casedata['busdc_mc'][:,6].nonzero()[0]  
        self.busdc_mqc    = self.casedata['busdc_mc'][:,8].nonzero()[0] 
        self.busdc_mpcinj = self.casedata['busdc_mc'][:,10].nonzero()[0]  
        self.busdc_mudc   = self.casedata['busdc_mc'][:,12].nonzero()[0] 
        self.busdc_midc   = self.casedata['busdc_mc'][:,14].nonzero()[0] 
        self.busdc_mm     = self.casedata['busdc_mc'][:,16].nonzero()[0] 
        self.busdc_mtheta = self.casedata['busdc_mc'][:,18].nonzero()[0]

        self.meter_num = self.bus_mv.shape[0]+self.bus_mp.shape[0]+self.bus_mq.shape[0]+self.branch_mfp.shape[0]\
        +self.branch_mfq.shape[0]+self.branch_mtp.shape[0]+self.branch_mtq.shape[0]+self.busdc_mps.shape[0]+self.busdc_mqs.shape[0]\
        +self.busdc_mpc.shape[0]+self.busdc_mqc.shape[0]\
        +self.busdc_mpcinj.shape[0]+self.busdc_mudc.shape[0]+self.busdc_midc.shape[0]+self.busdc_mm.shape[0]\
        +self.busdc_mtheta.shape[0]

    def _h(self):
           
        x_u = self.x_true[0:self.bus_num]  
        x_theta = np.r_[[0], self.x_true[self.bus_num:2*self.bus_num-1]]  #include reference bus
        x_uc = self.x_true[2*self.bus_num-1:2*self.bus_num+1]  
        x_thetac = self.x_true[2*self.bus_num+1:2*self.bus_num+3] 
        x_udci = self.x_true[2*self.bus_num+3:2*self.bus_num+4] 
        x_idci = self.x_true[2*self.bus_num+4:2*self.bus_num+5] 
        h = np.array([])
        # measurenment function of buses
        for i in self.bus_mv:
            h = np.append(h, x_u[i])
        for i in self.bus_mp:
            if i in self.busdc_ij:
                   i_dc = np.where(self.busdc_ij==i)[0][0] # the index in busdc of the bus i, = 0 or 1
                   ps = -x_u[i]**2 * self.g_sc[i_dc] + x_u[i] * x_uc[i_dc] * (self.g_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]) + \
                   self.b_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]))
            else:
                   ps = 0   
            h = np.append(h, (x_u[i] * x_u * (self.Ybus.real[i,:]*np.cos(x_theta[i]-x_theta) + \
            self.Ybus.imag[i,:]*np.sin(x_theta[i]-x_theta))).sum() - ps)
        for i in self.bus_mq:
            if i in self.busdc_ij:
                   i_dc = np.where(self.busdc_ij==i)[0][0] # the index in busdc of the bus i, = 0 or 1
                   qs = x_u[i]**2 * self.b_sc[i_dc] + x_u[i] * x_uc[i_dc] * (self.g_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]) - \
                   self.b_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]))
            else:
                   qs = 0 
            h = np.append(h, (x_u[i] * x_u * (self.Ybus.real[i,:]*np.sin(x_theta[i]-x_theta) - \
            self.Ybus.imag[i,:]*np.cos(x_theta[i]-x_theta))).sum() - qs)        
        # measurement function of branches
        for k in self.branch_mfp:
            ybra,i,j = self._find_k(k)  
            h = np.append(h, (ybra[2].real+ybra[3].real)*x_u[i]**2 - ybra[3].real*x_u[i]*x_u[j]*np.cos(x_theta[i]-x_theta[j])\
            - ybra[3].imag*x_u[i]*x_u[j]*np.sin(x_theta[i]-x_theta[j]))
        for k in self.branch_mfq:
            ybra,i,j = self._find_k(k)         
            h = np.append(h, -(ybra[2].imag+ybra[3].imag)*x_u[i]**2 + ybra[3].imag*x_u[i]*x_u[j]*np.cos(x_theta[i]-x_theta[j])\
            - ybra[3].real*x_u[i]*x_u[j]*np.sin(x_theta[i]-x_theta[j]))
        for k in self.branch_mtp:
            ybra,i,j = self._find_k(k)            
            h = np.append(h, (ybra[4].real+ybra[3].real)*x_u[j]**2 - ybra[3].real*x_u[j]*x_u[i]*np.cos(x_theta[j]-x_theta[i])\
            - ybra[3].imag*x_u[j]*x_u[i]*np.sin(x_theta[j]-x_theta[i]))
        for k in self.branch_mtq:
            ybra,i,j = self._find_k(k)           
            h = np.append(h, -(ybra[4].imag+ybra[3].imag)*x_u[j]**2 + ybra[3].imag*x_u[j]*x_u[i]*np.cos(x_theta[j]-x_theta[i])\
            - ybra[3].real*x_u[j]*x_u[i]*np.sin(x_theta[j]-x_theta[i]))
            
        # measurement function of dc system
        for i in self.busdc_mps:
            i_ac = self.busdc_ij[i]  #the index in bus of the busdc i
            h = np.append(h,-x_u[i_ac]**2 * self.g_sc[i] + x_u[i_ac] * x_uc[i] * (self.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) + \
            self.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i])))
        for i in self.busdc_mqs:
            i_ac = self.busdc_ij[i]  #the index in bus of the busdc i
            h = np.append(h, x_u[i_ac]**2 * self.b_sc[i] + x_u[i_ac] * x_uc[i] * (self.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) - \
            self.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i])))
        for i in self.busdc_mpc:
            i_ac = self.busdc_ij[i]
            h = np.append(h, x_uc[i]**2 * self.g_sc[i] - x_u[i_ac] * x_uc[i] * (self.g_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac]) + \
            self.b_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac])))
        for i in self.busdc_mqc:
            i_ac = self.busdc_ij[i]
            h = np.append(h, -x_uc[i]**2 * self.b_sc[i] - x_u[i_ac] * x_uc[i] * (self.g_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac]) - \
            self.b_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac])))
            
        for i in self.busdc_mpcinj:
            i_ac = self.busdc_ij[i]
            pc = x_uc[i]**2 * self.g_sc[i] - x_u[i_ac] * x_uc[i] * (self.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) - \
            self.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))
            ic = (((self.g_sc[i]**2 + self.b_sc[i]**2) * (x_u[i_ac]**2 + x_uc[i]**2 - 2*x_u[i_ac]*x_uc[i]*np.cos(x_thetac[i]-x_theta[i_ac])))**0.5)/(3**0.5)
            if i==0:
                   ploss = self.lossa[i] + self.lossb[i] * ic + self.losscr[i] * ic**2
                   pdc   = x_udci * x_idci
            elif i==1:
                   ploss = self.lossa[i] + self.lossb[i] * ic + self.lossci[i] * ic**2
                   pdc   = (x_udci - x_idci * self.r_dc) * (-x_idci)
            h = np.append(h, pc + ploss + pdc)
     
        for i in self.busdc_mudc:
            if i==0:
                   h = np.append(h, x_udci)
            elif i==1:
                   h = np.append(h, x_udci - x_idci * self.r_dc)
        for i in self.busdc_midc:
            if i==0:
                   h = np.append(h, x_idci)
            elif i==1:
                   h = np.append(h, -x_idci)
        for i in self.busdc_mm:
            i_ac = self.busdc_ij[i]
            if i==0:
                   h = np.append(h, 2**0.5 * x_uc[i]/x_udci)
            elif i==1:
                   h = np.append(h, 2**0.5 * x_uc[i]/(x_udci - x_idci * self.r_dc))
        for i in self.busdc_mtheta:
            i_ac = self.busdc_ij[i]
            h = np.append(h, x_thetac[i] - x_theta[i_ac])
        return h
                         
    def _get_z(self):
            
        """get the values of measurements z"""
        self.z_true_1 = np.r_[self.casedata['bus'][:,7][self.bus_mv],self.casedata['bus_inj'][:,1][self.bus_mp]/self.bases,\
        self.casedata['bus_inj'][:,2][self.bus_mq]/self.bases,self.casedata['branch'][:,13][self.branch_mfp]/self.bases,\
        self.casedata['branch'][:,14][self.branch_mfq]/self.bases,self.casedata['branch'][:,15][self.branch_mtp]/self.bases,\
        self.casedata['branch'][:,16][self.branch_mtq]/self.bases,\
        self.casedata['convdc'][:,3][self.busdc_mps]/self.bases_ac,self.casedata['convdc'][:,4][self.busdc_mqs]/self.bases_ac,\
        self.casedata['convdc'][:,26][self.busdc_mpc]/self.bases_ac,self.casedata['convdc'][:,27][self.busdc_mqc]/self.bases_ac,\
        0*self.casedata['convdc'][:,0][self.busdc_mpcinj],\
        self.casedata['busdc'][:,4][self.busdc_mudc],\
        -(self.casedata['busdc'][:,3][self.busdc_midc]/self.bases_dc)/self.casedata['busdc'][:,4][self.busdc_midc],\
        2**0.5 * self.casedata['convdc'][:,24][self.busdc_mm]/self.casedata['busdc'][:,4][self.busdc_mm],\
        (self.casedata['convdc'][:,25][self.busdc_mtheta] - self.casedata['bus'][:,8][self.busdc_ij[self.busdc_mtheta]])*np.pi/180.0]

        self.z_true_2 = self._h()
        self.z_true = self.z_true_2
        self.sigma = np.r_[self.casedata['bus_mc'][:,2][self.bus_mv],self.casedata['bus_mc'][:,4][self.bus_mp],\
        self.casedata['bus_mc'][:,6][self.bus_mq],self.casedata['branch_mc'][:,3][self.branch_mfp],\
        self.casedata['branch_mc'][:,5][self.branch_mfq],self.casedata['branch_mc'][:,7][self.branch_mtp],\
        self.casedata['branch_mc'][:,9][self.branch_mtq],\
        self.casedata['busdc_mc'][:,3][self.busdc_mps], self.casedata['busdc_mc'][:,5][self.busdc_mqs], self.casedata['busdc_mc'][:,7][self.busdc_mpc],\
        self.casedata['busdc_mc'][:,9][self.busdc_mqc], self.casedata['busdc_mc'][:,11][self.busdc_mpcinj],\
        self.casedata['busdc_mc'][:,13][self.busdc_mudc], self.casedata['busdc_mc'][:,15][self.busdc_midc], self.casedata['busdc_mc'][:,17][self.busdc_mm],\
        self.casedata['busdc_mc'][:,19][self.busdc_mtheta]]
        
        self.z = np.random.normal(self.z_true, self.sigma)

    def _R_W(self):
        """W = R^-1, R = diag{sigma^2...}""" 
        self.R = np.diag(self.sigma**2)
        self.W = np.diag(1/(self.sigma**2))

    def _get_x_true(self):
        """ get the true vlaues of x_u and x_theta"""
        self.x_u_true = self.casedata['bus'][:,7]
        self.x_theta_true = (self.casedata['bus'][:,8])*np.pi/180.0
        self.x_uc_true = self.casedata['convdc'][:,24]
        self.x_thetac_true =  (self.casedata['convdc'][:,25])*np.pi/180.0
        self.x_udci_true = self.casedata['busdc'][0,4]
        self.x_idci_true = -(self.casedata['busdc'][0,3]/self.bases_dc)/self.casedata['busdc'][0,4]
        self.x_true = np.r_[self.x_u_true, self.x_theta_true[1:], self.x_uc_true, self.x_thetac_true, self.x_udci_true, self.x_idci_true]  
   
    def get_z(self):
        self.z = np.random.normal(self.z_true, self.sigma)

  
class WlsSe(NetworkData):
    
    def __init__(self, nd):
        """ """
        self.nd = nd

    def _x_init(self):
        """initilize values of state variables x"""
        
        self.x_u = np.ones(self.nd.bus_num)
        self.x_theta = np.zeros(self.nd.bus_num)
        #take the angle of the first bus as reference angle
        self.x_theta[0] = 0
        self.x_uc = np.ones(2)
        self.x_thetac = np.zeros(2) + 0.1
        self.x_udci = np.ones(1)
        self.x_idci = np.ones(1)
        self.x = np.r_[self.x_u, self.x_theta[1:], self.x_uc, self.x_thetac, self.x_udci, self.x_idci]
        '''
        self.x_u = self.nd.x_u_true
        self.x_theta = self.nd.x_theta_true
        #take the angle of the first bus as reference angle
        self.x_theta[0] = 0
        self.x_uc = self.nd.x_uc_true
        self.x_thetac = self.nd.x_thetac_true
        self.x_udci = self.nd.x_udci_true
        self.x_idci = self.nd.x_idci_true
        self.x = np.r_[self.x_u, self.x_theta[1:], self.x_uc, self.x_thetac, self.x_udci, self.x_idci]        
        '''

    def _find_k(self,k):
        """k -> ybra,i,j  & k,i,j are both index of array, not man-made"""
        ybra = self.nd.Ybra[k]
        i =np.where(self.nd.casedata['bus'][:,0]==self.nd.casedata['branch'][k,0])[0][0]
        j =np.where(self.nd.casedata['bus'][:,0]==self.nd.casedata['branch'][k,1])[0][0]
        return ybra, i, j
        
    def _h(self):
        """compute the values of measurement function at the point of [x_u, x_theta]"""
        x_u, x_theta, x_uc, x_thetac, x_udci, x_idci = self.x_u, self.x_theta, self.x_uc, self.x_thetac, self.x_udci, self.x_idci
        h = np.array([])
        # measurenment function of buses
        for i in self.nd.bus_mv:
            h = np.append(h, x_u[i])
        for i in self.nd.bus_mp:
            if i in self.nd.busdc_ij:
                   i_dc = np.where(self.nd.busdc_ij==i)[0][0] # the index in busdc of the bus i, = 0 or 1
                   ps = -x_u[i]**2 * self.nd.g_sc[i_dc] + x_u[i] * x_uc[i_dc] * (self.nd.g_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]) + \
                   self.nd.b_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]))
            else:
                   ps = 0   
            h = np.append(h, (x_u[i] * x_u * (self.nd.Ybus.real[i,:]*np.cos(x_theta[i]-x_theta) + \
            self.nd.Ybus.imag[i,:]*np.sin(x_theta[i]-x_theta))).sum() - ps)
        for i in self.nd.bus_mq:
            if i in self.nd.busdc_ij:
                   i_dc = np.where(self.nd.busdc_ij==i)[0][0] # the index in busdc of the bus i, = 0 or 1
                   qs = x_u[i]**2 * self.nd.b_sc[i_dc] + x_u[i] * x_uc[i_dc] * (self.nd.g_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]) - \
                   self.nd.b_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]))
            else:
                   qs = 0 
            h = np.append(h, (x_u[i] * x_u * (self.nd.Ybus.real[i,:]*np.sin(x_theta[i]-x_theta) - \
            self.nd.Ybus.imag[i,:]*np.cos(x_theta[i]-x_theta))).sum() - qs)        
        # measurement function of branches
        for k in self.nd.branch_mfp:
            ybra,i,j = self._find_k(k)  
            h = np.append(h, (ybra[2].real+ybra[3].real)*x_u[i]**2 - ybra[3].real*x_u[i]*x_u[j]*np.cos(x_theta[i]-x_theta[j])\
            - ybra[3].imag*x_u[i]*x_u[j]*np.sin(x_theta[i]-x_theta[j]))
        for k in self.nd.branch_mfq:
            ybra,i,j = self._find_k(k)            
            h = np.append(h, -(ybra[2].imag+ybra[3].imag)*x_u[i]**2 + ybra[3].imag*x_u[i]*x_u[j]*np.cos(x_theta[i]-x_theta[j])\
            - ybra[3].real*x_u[i]*x_u[j]*np.sin(x_theta[i]-x_theta[j]))
        for k in self.nd.branch_mtp:
            ybra,i,j = self._find_k(k)            
            h = np.append(h, (ybra[4].real+ybra[3].real)*x_u[j]**2 - ybra[3].real*x_u[j]*x_u[i]*np.cos(x_theta[j]-x_theta[i])\
            - ybra[3].imag*x_u[j]*x_u[i]*np.sin(x_theta[j]-x_theta[i]))
        for k in self.nd.branch_mtq:
            ybra,i,j = self._find_k(k)            
            h = np.append(h, -(ybra[4].imag+ybra[3].imag)*x_u[j]**2 + ybra[3].imag*x_u[j]*x_u[i]*np.cos(x_theta[j]-x_theta[i])\
            - ybra[3].real*x_u[j]*x_u[i]*np.sin(x_theta[j]-x_theta[i]))
            
        # measurement function of dc system
        for i in self.nd.busdc_mps:
            i_ac = self.nd.busdc_ij[i]  #the index in bus of the busdc i
            h = np.append(h,-x_u[i_ac]**2 * self.nd.g_sc[i] + x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) + \
            self.nd.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i])))
        for i in self.nd.busdc_mqs:
            i_ac = self.nd.busdc_ij[i]  #the index in bus of the busdc i
            h = np.append(h, x_u[i_ac]**2 * self.nd.b_sc[i] + x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) - \
            self.nd.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i])))
        for i in self.nd.busdc_mpc:
            i_ac = self.nd.busdc_ij[i]
            h = np.append(h, x_uc[i]**2 * self.nd.g_sc[i] - x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac]) + \
            self.nd.b_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac])))
        for i in self.nd.busdc_mqc:
            i_ac = self.nd.busdc_ij[i]
            h = np.append(h, -x_uc[i]**2 * self.nd.b_sc[i] - x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac]) - \
            self.nd.b_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac])))
        for i in self.nd.busdc_mpcinj:
            i_ac = self.nd.busdc_ij[i]
            pc = x_uc[i]**2 * self.nd.g_sc[i] - x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac]) + \
            self.nd.b_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac]))
            ic = (((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (x_u[i_ac]**2 + x_uc[i]**2 - 2*x_u[i_ac]*x_uc[i]*np.cos(x_thetac[i]-x_theta[i_ac])))**0.5)/(3**0.5)
            if i==0:
                   ploss = self.nd.lossa[i] + self.nd.lossb[i] * ic + self.nd.losscr[i] * ic**2
                   pdc   = x_udci * x_idci
            elif i==1:
                   ploss = self.nd.lossa[i] + self.nd.lossb[i] * ic + self.nd.lossci[i] * ic**2
                   pdc   = (x_udci - x_idci * self.nd.r_dc) * (-x_idci)
            h = np.append(h, pc + ploss + pdc)
        for i in self.nd.busdc_mudc:
            if i==0:
                   h = np.append(h, x_udci)
            elif i==1:
                   h = np.append(h, x_udci - x_idci * self.nd.r_dc)
        for i in self.nd.busdc_midc:
            if i==0:
                   h = np.append(h, x_idci)
            elif i==1:
                   h = np.append(h, -x_idci)
        for i in self.nd.busdc_mm:
            i_ac = self.nd.busdc_ij[i]
            if i==0:
                   h = np.append(h, 2**0.5 * x_uc[i]/x_udci)
            elif i==1:
                   h = np.append(h, 2**0.5 * x_uc[i]/(x_udci - x_idci * self.nd.r_dc))
        for i in self.nd.busdc_mtheta:
            i_ac = self.nd.busdc_ij[i]
            h = np.append(h, x_thetac[i] - x_theta[i_ac])
        self.h = h

    def _H(self):
        """compute the Jacobia matrix"""  
        x_u, x_theta, x_uc, x_thetac, x_udci, x_idci = self.x_u, self.x_theta, self.x_uc, self.x_thetac, self.x_udci, self.x_idci
        H_u     = np.zeros([self.nd.meter_num, self.nd.bus_num])
        H_theta = np.zeros([self.nd.meter_num, self.nd.bus_num])
        H_uc    = np.zeros([self.nd.meter_num, 2])
        H_thetac= np.zeros([self.nd.meter_num, 2])
        H_udci  = np.zeros([self.nd.meter_num, 1])
        H_idci  = np.zeros([self.nd.meter_num, 1])
        c = 0
        # dV/dv dV/dtheta
        for i in self.nd.bus_mv:
            H_u[c,i] = 1
            c = c+1
        # dPi/dv dPi/dtheta
        for i in self.nd.bus_mp:
            if i in self.nd.busdc_ij:
                   i_dc = np.where(self.nd.busdc_ij==i)[0][0] # the index in busdc of the bus i, = 0 or 1 
                   h_u_i = -2 * x_u[i] * self.nd.g_sc[i_dc] +  x_uc[i_dc] * (self.nd.g_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]) + \
                   self.nd.b_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]))
                   
                   h_theta_i = x_u[i] * x_uc[i_dc] * (-self.nd.g_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]) + \
                   self.nd.b_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]))
                   
                   H_uc[c, i_dc] = x_u[i] * (self.nd.g_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]) + \
                   self.nd.b_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]))
                   
                   H_thetac[c, i_dc] = x_u[i] * x_uc[i_dc] * (self.nd.g_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]) - \
                   self.nd.b_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]))      
            else:
                   h_u_i = 0
                   h_theta_i = 0  
            H_u[c,:] = x_u[i] * (self.nd.Ybus.real[i,:]*np.cos(x_theta[i]-x_theta) + self.nd.Ybus.imag[i,:]*np.sin(x_theta[i]-x_theta))
            
            H_u[c,i] = (self.nd.Ybus.real[i,i]*x_u[i]**2 + (x_u[i] * x_u * (self.nd.Ybus.real[i,:]*np.cos(x_theta[i]-x_theta) + \
            self.nd.Ybus.imag[i,:]*np.sin(x_theta[i]-x_theta))).sum())/x_u[i] - h_u_i

            H_theta[c,:] = x_u[i] * x_u * (self.nd.Ybus.real[i,:]*np.sin(x_theta[i]-x_theta) - self.nd.Ybus.imag[i,:]*np.cos(x_theta[i]-x_theta))
            
            H_theta[c,i] = -self.nd.Ybus.imag[i,i]*x_u[i]**2 - (x_u[i] * x_u * (self.nd.Ybus.real[i,:]*np.sin(x_theta[i]-x_theta) - \
            self.nd.Ybus.imag[i,:]*np.cos(x_theta[i]-x_theta))).sum() - h_theta_i
            c = c+1
        # dQi/dv dQi/dtheta
        for i in self.nd.bus_mp:
               
            if i in self.nd.busdc_ij:
                   i_dc = np.where(self.nd.busdc_ij==i)[0][0] # the index in busdc of the bus i, = 0 or 1 
                   h_u_i = 2 * x_u[i] * self.nd.b_sc[i_dc] + x_uc[i_dc] * (self.nd.g_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]) - \
                   self.nd.b_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]))
                   
                   h_theta_i = x_u[i] * x_uc[i_dc] * (self.nd.g_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]) + \
                   self.nd.b_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]))
                   
                   H_uc[c, i_dc] = x_u[i] * (self.nd.g_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]) - \
                   self.nd.b_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]))
                   
                   H_thetac[c, i_dc] = x_u[i] * x_uc[i_dc] * (-self.nd.g_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]) - \
                   self.nd.b_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]))
            else:
                   h_u_i = 0
                   h_theta_i = 0                 
                
            H_u[c,:] = x_u[i] * (self.nd.Ybus.real[i,:]*np.sin(x_theta[i]-x_theta) - self.nd.Ybus.imag[i,:]*np.cos(x_theta[i]-x_theta))
            
            H_u[c,i] = (-self.nd.Ybus.imag[i,i]*x_u[i]**2 + (x_u[i] * x_u * (self.nd.Ybus.real[i,:]*np.sin(x_theta[i]-x_theta) - \
            self.nd.Ybus.imag[i,:]*np.cos(x_theta[i]-x_theta))).sum())/x_u[i] - h_u_i

            H_theta[c,:] = -x_u[i] * x_u * (self.nd.Ybus.real[i,:]*np.cos(x_theta[i]-x_theta) + self.nd.Ybus.imag[i,:]*np.sin(x_theta[i]-x_theta))
            
            H_theta[c,i] = -self.nd.Ybus.real[i,i]*x_u[i]**2 + (x_u[i] * x_u * (self.nd.Ybus.real[i,:]*np.cos(x_theta[i]-x_theta) + \
            self.nd.Ybus.imag[i,:]*np.sin(x_theta[i]-x_theta))).sum() - h_theta_i
            c = c+1
        # dPij/dv dPij/dtheta
        for k in self.nd.branch_mfp:
            ybra,i,j = self._find_k(k)
            H_u[c,i] = 2*x_u[i]*(ybra[2].real+ybra[3].real) - x_u[j]*ybra[3].real*np.cos(x_theta[i]-x_theta[j]) - \
            x_u[j]*ybra[3].imag*np.sin(x_theta[i]-x_theta[j])
            H_theta[c,i] = x_u[i]*x_u[j]*(ybra[3].real*np.sin(x_theta[i]-x_theta[j])-ybra[3].imag*np.cos(x_theta[i]-x_theta[j]))
            H_u[c,j] = -x_u[i]*(ybra[3].real*np.cos(x_theta[i]-x_theta[j])+ybra[3].imag*np.sin(x_theta[i]-x_theta[j]))
            H_theta[c,j] = -x_u[i]*x_u[j]*(ybra[3].real*np.sin(x_theta[i]-x_theta[j])-ybra[3].imag*np.cos(x_theta[i]-x_theta[j]))
            c = c+1
        # dQij/dv dQij/dtheta
        for k in self.nd.branch_mfq:
            ybra,i,j = self._find_k(k)
            H_u[c,i] = -2*x_u[i]*(ybra[2].imag+ybra[3].imag) - x_u[j]*ybra[3].real*np.sin(x_theta[i]-x_theta[j]) + \
            x_u[j]*ybra[3].imag*np.cos(x_theta[i]-x_theta[j])
            H_theta[c,i] = -x_u[i]*x_u[j]*(ybra[3].real*np.cos(x_theta[i]-x_theta[j])+ybra[3].imag*np.sin(x_theta[i]-x_theta[j]))
            H_u[c,j] = -x_u[i]*(ybra[3].real*np.sin(x_theta[i]-x_theta[j])-ybra[3].imag*np.cos(x_theta[i]-x_theta[j]))
            H_theta[c,j] = x_u[i]*x_u[j]*(ybra[3].real*np.cos(x_theta[i]-x_theta[j])+ybra[3].imag*np.sin(x_theta[i]-x_theta[j]))
            c = c+1
        # dPji/dv dPji/dtheta
        for k in self.nd.branch_mtp:
            ybra,i,j = self._find_k(k) 
            H_u[c,i] = x_u[j]*(-ybra[3].real*np.cos(x_theta[i]-x_theta[j])+ybra[3].imag*np.sin(x_theta[i]-x_theta[j]))
            H_theta[c,i] = x_u[i]*x_u[j]*(ybra[3].real*np.sin(x_theta[i]-x_theta[j])+ybra[3].imag*np.cos(x_theta[i]-x_theta[j]))
            H_u[c,j] = 2*x_u[j]*(ybra[4].real+ybra[3].real) - x_u[i]*ybra[3].real*np.cos(x_theta[i]-x_theta[j]) + \
            x_u[i]*ybra[3].imag*np.sin(x_theta[i]-x_theta[j])
            H_theta[c,j] = -x_u[i]*x_u[j]*(ybra[3].real*np.sin(x_theta[i]-x_theta[j])+ybra[3].imag*np.cos(x_theta[i]-x_theta[j]))
            c = c+1
        # dQji/dv dQji/dtheta
        for k in self.nd.branch_mtq:
            ybra,i,j = self._find_k(k) 
            H_u[c,i] = x_u[j]*(ybra[3].real*np.sin(x_theta[i]-x_theta[j])+ybra[3].imag*np.cos(x_theta[i]-x_theta[j]))
            H_theta[c,i] = x_u[i]*x_u[j]*(ybra[3].real*np.cos(x_theta[i]-x_theta[j])-ybra[3].imag*np.sin(x_theta[i]-x_theta[j]))
            H_u[c,j] = -2*x_u[j]*(ybra[4].imag+ybra[3].imag) + x_u[i]*ybra[3].real*np.sin(x_theta[i]-x_theta[j]) + \
            x_u[i]*ybra[3].imag*np.cos(x_theta[i]-x_theta[j])
            H_theta[c,j] = -x_u[i]*x_u[j]*(ybra[3].real*np.cos(x_theta[i]-x_theta[j])-ybra[3].imag*np.sin(x_theta[i]-x_theta[j]))
            c = c+1
            
            
        for i in self.nd.busdc_mps:
            i_ac = self.nd.busdc_ij[i]  #the index in bus of the busdc i
            H_u[c, i_ac] = -2 * x_u[i_ac] * self.nd.g_sc[i] +  x_uc[i] * (self.nd.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) + \
            self.nd.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))
            H_theta[c, i_ac] = x_u[i_ac] * x_uc[i] * (-self.nd.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) + \
            self.nd.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]))
            H_uc[c, i] = x_u[i_ac] * (self.nd.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) + \
            self.nd.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))
            H_thetac[c, i] = x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) - \
            self.nd.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]))
            c = c + 1

        for i in self.nd.busdc_mqs:
            i_ac = self.nd.busdc_ij[i]  #the index in bus of the busdc i
            H_u[c, i_ac] = 2 * x_u[i_ac] * self.nd.b_sc[i] +  x_uc[i] * (self.nd.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) - \
            self.nd.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]))
            H_theta[c, i_ac] = x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.cos(x_theta[i_ac] - x_thetac[i]) + \
            self.nd.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))
            H_uc[c, i] = x_u[i_ac] * (self.nd.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) - \
            self.nd.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]))
            H_thetac[c, i] = x_u[i_ac] * x_uc[i] * (-self.nd.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) - \
            self.nd.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))
            c = c + 1            
        
        for i in self.nd.busdc_mpc:
            i_ac = self.nd.busdc_ij[i]  #the index in bus of the busdc i
            H_u[c, i_ac] = -x_uc[i] * (self.nd.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) + \
            self.nd.b_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac]))
            H_theta[c, i_ac] = -x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac]) - \
            self.nd.b_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac]))
            H_uc[c, i] = 2 * x_uc[i] * self.nd.g_sc[i] -  x_u[i_ac] * (self.nd.g_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac]) + \
            self.nd.b_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac]))
            H_thetac[c, i] = -x_u[i_ac] * x_uc[i] * (-self.nd.g_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac]) + \
            self.nd.b_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac]))
            c = c + 1 
            
        for i in self.nd.busdc_mqc:
            i_ac = self.nd.busdc_ij[i]  #the index in bus of the busdc i
            H_u[c, i_ac] = x_uc[i] * (self.nd.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) + \
            self.nd.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]))
            
            H_theta[c, i_ac] = x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) - \
            self.nd.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))
            
            H_uc[c, i] = -2 * x_uc[i] * self.nd.b_sc[i] -  x_u[i_ac] * (self.nd.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) + \
            self.nd.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]))
            
            H_thetac[c, i] = x_u[i_ac] * x_uc[i] * (-self.nd.g_sc[i] * np.cos(x_theta[i_ac] - x_thetac[i]) + \
            self.nd.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))
            c = c + 1 
            
        for i in self.nd.busdc_mpcinj:
            i_ac = self.nd.busdc_ij[i]
            #the following four items are from 'for i in self.nd.busdc_mpc:'
            h_u_pc = -x_uc[i] * (self.nd.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) - \
            self.nd.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))

            h_theta_pc = -x_u[i_ac] * x_uc[i] * (-self.nd.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) - \
            self.nd.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]))
            
            h_uc_pc = 2 * x_uc[i] * self.nd.g_sc[i] -  x_u[i_ac] * (self.nd.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) - \
            self.nd.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))
            
            h_thetac_pc = -x_u[i_ac] * x_uc[i] * (self.nd.g_sc[i] * np.sin(x_theta[i_ac] - x_thetac[i]) + \
            self.nd.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]))
            
            icc = ((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (x_u[i_ac]**2 + x_uc[i]**2 - 2*x_u[i_ac]*x_uc[i]*np.cos(x_thetac[i]-x_theta[i_ac])))
            if i==0:
                   lossc = self.nd.losscr[i]
            elif i==1:
                   lossc = self.nd.lossci[i]
            h_u_ploss = 3**(-0.5) * self.nd.lossb[i] * 0.5 * icc**(-0.5) *\
            ((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (2*x_u[i_ac] - 2*x_u[i_ac]*np.cos(x_thetac[i]-x_theta[i_ac]))) + \
            3**(-1) * lossc * ((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (2*x_u[i_ac] - 2*x_u[i_ac]*np.cos(x_thetac[i]-x_theta[i_ac])))
            
            h_theta_ploss = 3**(-0.5) *self.nd.lossb[i] * 0.5 * icc**(-0.5) *\
            ((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (-2*x_u[i_ac]*x_uc[i]*np.sin(x_thetac[i]-x_theta[i_ac]))) + \
            3**(-1) * lossc * ((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (- 2*x_u[i_ac]*x_uc[i]*np.sin(x_thetac[i]-x_theta[i_ac])))
            
            h_uc_ploss = 3**(-0.5) *self.nd.lossb[i] * 0.5 * icc**(-0.5) *\
            ((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (2*x_uc[i] - 2*x_uc[i]*np.cos(x_thetac[i]-x_theta[i_ac]))) + \
            3**(-1) *lossc * ((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (2*x_uc[i] - 2*x_uc[i]*np.cos(x_thetac[i]-x_theta[i_ac])))
            
            h_thetac_ploss = 3**(-0.5) *self.nd.lossb[i] * 0.5 * icc**(-0.5) *\
            ((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (2*x_u[i_ac]*x_uc[i]*np.sin(x_thetac[i]-x_theta[i_ac]))) + \
            3**(-1) *lossc * ((self.nd.g_sc[i]**2 + self.nd.b_sc[i]**2) * (2*x_u[i_ac]*x_uc[i]*np.sin(x_thetac[i]-x_theta[i_ac])))
            
            H_u[c, i_ac] = h_u_pc + h_u_ploss
            H_theta[c, i_ac] = h_theta_pc + h_theta_ploss
            H_uc[c, i] = h_uc_pc + h_uc_ploss
            H_thetac[c, i] = h_thetac_pc + h_thetac_ploss
            
            if i==0:
                   H_udci[c, 0] = x_idci
                   H_idci[c, 0] = x_udci
            elif i==1:
                   H_udci[c, 0] = -x_idci
                   H_idci[c, 0] = -x_udci + 2 * self.nd.r_dc * x_idci             
            c = c+1  
        
        for i in self.nd.busdc_mudc:
            if i==0:
                   H_udci[c, 0] = 1
            elif i==1:
                   H_udci[c, 0] = 1
                   H_idci[c, 0] = -self.nd.r_dc
            c = c+1
        
        for i in self.nd.busdc_midc:
            if i==0:
                   H_idci[c, 0] = 1
            elif i==1:
                   H_idci[c, 0] = -1
            c = c+1
        
        for i in self.nd.busdc_mm:
            if i==0:
                   H_uc[c, i] = (2**0.5)/x_udci
                   H_udci[c,0] = -(2**0.5) * x_uc[i] * (x_udci**(-2))
            elif i==1:
                   H_uc[c, i] = (2**0.5)/(x_udci - x_idci * self.nd.r_dc)
                   H_udci[c,0] = -(2**0.5) * x_uc[i] * ((x_udci - x_idci * self.nd.r_dc)**(-2))
                   H_idci[c,0] = (2**0.5) * x_uc[i] * ((x_udci - x_idci * self.nd.r_dc)**(-2)) * self.nd.r_dc
            c = c+1
        
        for i in  self.nd.busdc_mtheta:
            i_ac = self.nd.busdc_ij[i]
            H_theta[c, i_ac] = -1
            H_thetac[c, i] = 1
            c =c+1
            
        H_theta = np.delete(H_theta,0,axis=1)
        H = np.c_[H_u,H_theta, H_uc, H_thetac, H_udci, H_idci]
        self.H = H
 

        
    def _G_K(self):
        """G = H' * W * H"""       
        self.G = np.dot(np.dot(self.H.T, self.nd.W), self.H)
        self.K = np.dot(self.H, np.linalg.solve(self.G, np.dot(self.H.T, self.nd.W)))
        
    def _x_delta(self):
        """ x_delta"""   
        self.x_delta = np.linalg.solve(self.G, np.dot(np.dot(self.H.T, self.nd.W), self.nd.z-self.h))
        
    def _x_updt(self):
        """ x  = x + delta"""
        self.x_u = self.x_u + self.x_delta[0:self.nd.bus_num]  
        self.x_theta = self.x_theta + np.r_[[0], self.x_delta[self.nd.bus_num:2*self.nd.bus_num-1]]  #include reference bus
        self.x_uc = self.x_uc + self.x_delta[2*self.nd.bus_num-1:2*self.nd.bus_num+1]  
        self.x_thetac = self.x_thetac + self.x_delta[2*self.nd.bus_num+1:2*self.nd.bus_num+3] 
        self.x_udci = self.x_udci + self.x_delta[2*self.nd.bus_num+3:2*self.nd.bus_num+4] 
        self.x_idci = self.x_idci + self.x_delta[2*self.nd.bus_num+4:2*self.nd.bus_num+5] 

        self.x = np.r_[self.x_u, self.x_theta[1:], self.x_uc, self.x_thetac, self.x_udci, self.x_idci]  #not include reference bus
    def _bddi(self, flag):
        """ 
        flag == true : plot the figure, else not.
        """
        self.r = self.nd.z - self.h #residual
        self.S = np.eye(self.nd.meter_num) - self.K  #residual sensitivity matrix
        self.Omega = np.dot(self.S, self.nd.R)  #residual covariance matrix
        self.r_N = np.abs(self.r)/np.sqrt(np.diag(self.Omega)) #normalized residual
        # 2-norm of residuals
        self.r_norm2 = np.linalg.norm(self.r)
        # largest normalized residual
        self.r_N_max = self.r_N.max()

#        if flag == True:
#               print "-"*18, "r_N", "-"*18
#               fig = plt.figure()
#               fig.set_figheight(10)
#               fig.set_figwidth(10)
#               ax1 = fig.add_subplot(211)
#               ax1.vlines(range(len(self.r_N)),[0],self.r_N)
#               ax1.plot(range(len(self.r_N)), self.r_N, 'o', markersize=4)
#               ax1.grid(axis='y')
#               plt.show()
#               print  "-"*40

    def wls(self):
        """
        main program of wls, excluding bddi
        """
        self._x_init()
        c = 0
        self.intera = list()
        while True:
               self._H()
               self._h()
               self._G_K()
               self._x_delta()
               #print c,'\t',abs(self.x_delta).max()
               self.intera.append(abs(self.x_delta).max())
               if abs(self.x_delta).max() < 1e-5:
                      break
               self._x_updt()
               c = c+1
        self._bddi(False)

        

        