# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 21:42:11 2016

@author: dell
"""

import pyomo.environ as pe
from pyomo.opt import SolverFactory
import numpy as np
import itertools
import copy
#from se import NetworkData, WlsSe
#from ssa import StaticSA

class Attack:
       def __init__(self, s_c, conv_i, r1, r2):
              """ """
              self.s_c = s_c
              self.conv_i = conv_i  # the index of converter attacked, [0,1...]
              self.r1 = r1
              self.r2 = r2
              self._w()
       def _w(self):
              """ """
              w = np.r_[self.s_c.nd.casedata['bus_mc'][:,1],self.s_c.nd.casedata['bus_mc'][:,3],self.s_c.nd.casedata['bus_mc'][:,5],
              self.s_c.nd.casedata['branch_mc'][:,2],self.s_c.nd.casedata['branch_mc'][:,4],self.s_c.nd.casedata['branch_mc'][:,6],
              self.s_c.nd.casedata['branch_mc'][:,8],self.s_c.nd.casedata['busdc_mc'][:,2],self.s_c.nd.casedata['busdc_mc'][:,4],
              self.s_c.nd.casedata['busdc_mc'][:,6],self.s_c.nd.casedata['busdc_mc'][:,8],self.s_c.nd.casedata['busdc_mc'][:,10],
              self.s_c.nd.casedata['busdc_mc'][:,12],self.s_c.nd.casedata['busdc_mc'][:,14],self.s_c.nd.casedata['busdc_mc'][:,16],
              self.s_c.nd.casedata['busdc_mc'][:,18]]
              w = w[np.where(w!=0)] #remove 0
              w = 2 - w #replace 1 by 1, 2 by 0; and 0 means feigned measurements.
              return w
       
       def _IF(self):
              """ """
              if self.k==1:     
                     ifi = np.array(list(set(range(self.s_c.x.size)).difference(set([self.s_c.nd.busdc_ij[self.conv_i], self.s_c.nd.busdc_ij[self.conv_i]\
                                         + self.s_c.nd.bus_num - 1, 2*self.s_c.nd.bus_num - 1 + self.conv_i, 2*self.s_c.nd.bus_num - 1 + self.conv_i + 2]))))
              else:
                     u = (1-self.v[self.k-1, self.n[self.k - 1]]) + self.w  # v=1 & w = 0 <-> u=0
                     i_z = np.where(u==0)[0]
                     H_sub = (abs(self.s_c.H[i_z])).sum(axis = 0)
                     i_x = np.where(H_sub==0)[0]
                     ifi = i_x
                     for i in range(1, self.k):
                            ifi = np.array(list(set(ifi).intersection(set(self.IF[i]))))
              self.IF[self.k] = ifi
              
       def _IV(self):
              """ """
              if self.k==1:
                     ivi = np.array([self.s_c.nd.busdc_ij[self.conv_i], self.s_c.nd.busdc_ij[self.conv_i] + self.s_c.nd.bus_num - 1, \
                                     2*self.s_c.nd.bus_num - 1 + self.conv_i, 2*self.s_c.nd.bus_num - 1 + self.conv_i + 2])
                           # index corresponding with ac voltage amplitude, ac phase angle, dc side voltage amplitude and dc side phase angle.
              else:
                     u = (1-self.v[self.k-1, self.n[self.k - 1]]) + self.w  # v=1 & w = 0 <-> u=0
                     i_z = np.where(u==0)[0]
                     H_sub = (abs(self.s_c.H[i_z])).sum(axis = 0)
                     i_x = np.where(H_sub!=0)[0]
                     ivi = i_x
                     for i in range(1, self.k):
                            ivi = np.array(list(set(ivi).difference(set(self.IV[i]))))
              self.IV[self.k] = ivi
       
       def _B(self):
              """ all members in the bool space"""
              dim = len(self.IV[self.k])
              k = 0
              b = [int(i) for i in list(bin(k)[2:])]
              b = [0] * (dim - len(b)) + b
              nb = np.array([b])
              while True:
                     k = k + 1
                     b = [int(i) for i in list(bin(k)[2:])]
                     b = [0] * (dim - len(b)) + b
                     if len(b) <= dim:
                            nb = np.r_[nb, np.array([b])]
                     else:
                            break
              nb = nb[0:-1]  #remove [1,1,...]
              nb = nb[::-1]
              self.B[self.k] = nb  
              
       def _v(self):
              """ """
              if self.k==1:
                     ind = self.IV[self.k][np.nonzero(1 - self.B[self.k][self.n[self.k]-1])[0]] #column index of H and corresponding di=0
                     H_sub = self.s_c.H[:,ind]
                     v = abs(H_sub).sum(axis=1)
                     v[np.nonzero(v)]=1
              else:
                     ind = self.IV[self.k][np.nonzero(1 - self.B[self.k][self.n[self.k]-1])[0]] #column index of H and corresponding di=0
                     H_sub = self.s_c.H[:,ind]
                     v = abs(H_sub).sum(axis=1)
                     v[np.nonzero(v)]=1
                     i_v_0 = np.array(list(set(np.where(self.v[self.k-1, self.n[self.k-1]]==1)[0]).intersection(set(np.where(self.w==0)[0]))))
                     v[i_v_0] = 0
                     
              self.v[self.k, self.n[self.k]] = v
              
       def _Z_L(self):
              """ """
              v_sum = 0
              for i in range(1,self.k+1):
                     v_sum = v_sum + self.v[i,self.n[i]]
              w_v = v_sum * self.w
              Z = np.nonzero(w_v)[0]
              L = Z.size
              self.Z[self.t] = Z
              self.L[self.t] = L
       
       def _equa(self, i_equa):
              """ """
              ind_x = range(self.s_c.x.size)  # list of index of x
              x_initial = {key: self.s_c.x[key] for key in ind_x}  #initial values of x
              
              model = pe.ConcreteModel()
              model.x = pe.Var(ind_x, initialize=x_initial)
              
              x_u = [model.x[i] for i in range(0, self.s_c.nd.bus_num)]
              x_theta = [model.x[i] for i in range(self.s_c.nd.bus_num, 2*self.s_c.nd.bus_num-1)]
              x_theta.insert(0,0)  
              x_uc = [model.x[i] for i in range(2*self.s_c.nd.bus_num-1, 2*self.s_c.nd.bus_num+1)]
              x_thetac = [model.x[i] for i in range(2*self.s_c.nd.bus_num+1, 2*self.s_c.nd.bus_num+3)]
              x_udci = [model.x[i] for i in range(2*self.s_c.nd.bus_num+3, 2*self.s_c.nd.bus_num+4)][0]
              x_idci = [model.x[i] for i in range(2*self.s_c.nd.bus_num+4, 2*self.s_c.nd.bus_num+5)][0]               

              model.obj = pe.Objective( expr = sum((model.x[i] - self.s_c.x[i])**2 for i in ind_x) )  # without sqrt
              model.con = pe.ConstraintList()
              
              #constraints about variables themselves
              Pi = 2*np.pi
              for i in range(0, self.s_c.nd.bus_num):
                     model.con.add(0.9 <= model.x[i] <= 1.1)  #x_u
              for i in range(self.s_c.nd.bus_num, 2*self.s_c.nd.bus_num-1):
                     model.con.add(-Pi <= model.x[i] <= Pi)  #x_theta                     
              for i in range(2*self.s_c.nd.bus_num-1, 2*self.s_c.nd.bus_num+1):
                     model.con.add(0.9 <= model.x[i] <= 1.1)  #x_uc
              for i in range(2*self.s_c.nd.bus_num+1, 2*self.s_c.nd.bus_num+3):
                     model.con.add(-Pi <= model.x[i] <= Pi)  #x_thetac              
              for i in range(2*self.s_c.nd.bus_num+3, 2*self.s_c.nd.bus_num+4):
                     model.con.add(0.8 <= model.x[i] <= 1.2)  #x_udci              
              for i in range(2*self.s_c.nd.bus_num+4, 2*self.s_c.nd.bus_num+5):
                     model.con.add( (model.x[i])**2 <= 1.2**2)  #x_idci
    
              # constraints g_z'(x)
              conv_i_ac = self.s_c.nd.busdc_ij[self.conv_i]  ##the index in bus of the busdc i
              pc_conv_i = x_uc[self.conv_i]**2 * self.s_c.nd.g_sc[self.conv_i] - x_u[conv_i_ac] * x_uc[self.conv_i] * (self.s_c.nd.g_sc[self.conv_i] * pe.cos(x_thetac[self.conv_i]-x_theta[conv_i_ac]) + \
                          self.s_c.nd.b_sc[self.conv_i] * pe.sin(x_thetac[self.conv_i]-x_theta[conv_i_ac]))    
              qc_conv_i = -x_uc[self.conv_i]**2 * self.s_c.nd.b_sc[self.conv_i] - x_u[conv_i_ac] * x_uc[self.conv_i] * (self.s_c.nd.g_sc[self.conv_i] * pe.sin(x_thetac[self.conv_i]-x_theta[conv_i_ac]) - \
                          self.s_c.nd.b_sc[self.conv_i] * pe.cos(x_thetac[self.conv_i]-x_theta[conv_i_ac]))
              
              model.con.add( pc_conv_i**2 + qc_conv_i**2 <= (self.r1 * x_u[conv_i_ac] * self.s_c.nd.casedata['convdc'][:,14][self.conv_i])**2 )
              model.con.add( (pc_conv_i + x_u[conv_i_ac]**2 * self.s_c.nd.g_sc[self.conv_i] )**2 + (qc_conv_i - x_u[conv_i_ac]**2 * self.s_c.nd.b_sc[self.conv_i] )**2 <= \
                            (self.r2 * x_u[conv_i_ac] * self.s_c.nd.casedata['convdc'][:,12][self.conv_i] * (self.s_c.nd.g_sc[self.conv_i]**2 + self.s_c.nd.b_sc[self.conv_i]**2)**0.5)**2 )
              
              # constraints x - x_c = 0
              for i in self.IF[self.k]:
                     model.con.add( model.x[i] - self.s_c.x[i] == 0)
              
              # constraints d(x - x_c) = 0
              for k in range(1, self.k + 1):
                     for (i,j) in itertools.izip(self.IV[k], self.B[k][self.n[k] - 1]):
                            if j == 1:
                                    model.con.add( model.x[i] - self.s_c.x[i] == 0)
             
              # constraints z_i = h_i(x)
              if i_equa == 2 or i_equa == 3:
                     v_sum = self.v[1, self.n[1]]
                     if i_equa == 2:
                            for i in range(2, self.k+1):
                                   v_sum = v_sum + self.v[i, self.n[i]]
                     elif i_equa == 3:
                            for i in range(2, self.k):
                                   v_sum = v_sum + self.v[i, self.n[i]]
                     v_sum[np.nonzero(v_sum)[0]] = 1
                     v_sum_w = (1 - v_sum) + self.w
                     self.ind_z_virtual = np.where(v_sum_w==0)[0]   # array made up of i, z_i = h_i(x)
                     
                     for i_z in self.ind_z_virtual:
                            type_iz, i_m = self._iz_info(i_z)
                            
                            if type_iz == 'mv':
                                   i = self.s_c.nd.bus_mv[i_m]
                                   cons = x_u[i]
       
                            elif type_iz == 'mp':
                                   i = self.s_c.nd.bus_mp[i_m]
                                   if i in self.s_c.nd.busdc_ij:
                                          i_dc = np.where(self.s_c.nd.busdc_ij==i)[0][0]
                                          ps = -x_u[i]**2 * self.s_c.nd.g_sc[i_dc] + x_u[i] * x_uc[i_dc] * (self.s_c.nd.g_sc[i_dc] * pe.cos(x_theta[i]-x_thetac[i_dc]) + \
                                          self.s_c.nd.b_sc[i_dc] * pe.sin(x_theta[i]-x_thetac[i_dc]))
                                   else:
                                          ps = 0
                                   cons = sum( (x_u[i] * x_u[j] * (self.s_c.nd.Ybus.real[i,j]*pe.cos(x_theta[i]-x_theta[j]) + self.s_c.nd.Ybus.imag[i,j]*pe.sin(x_theta[i]-x_theta[j]))) \
                                              for j in range(self.s_c.nd.bus_num) ) - ps
                                              
                            elif type_iz == 'mq':
                                   i = self.s_c.nd.bus_mq[i_m]
                                   if i in self.s_c.nd.busdc_ij:
                                          i_dc = np.where(self.s_c.nd.busdc_ij==i)[0][0]
                                          qs = x_u[i]**2 * self.s_c.nd.b_sc[i_dc] + x_u[i] * x_uc[i_dc] * (self.s_c.nd.g_sc[i_dc] * pe.sin(x_theta[i]-x_thetac[i_dc]) - \
                                          self.s_c.nd.b_sc[i_dc] * pe.cos(x_theta[i]-x_thetac[i_dc]))
                                   else:
                                          qs = 0
                                   cons = sum( (x_u[i] * x_u[j] * (self.s_c.nd.Ybus.real[i,j]*pe.sin(x_theta[i]-x_theta[j]) - self.s_c.nd.Ybus.imag[i,j]*pe.cos(x_theta[i]-x_theta[j]))) \
                                              for j in range(self.s_c.nd.bus_num) ) - qs
                                              
                            elif type_iz == 'mfp':
                                   k = self.s_c.nd.bus_mfp[i_m]
                                   ybra,i,j = self._find_k(k)
                                   cons = (ybra[2].real+ybra[3].real)*x_u[i]**2 - ybra[3].real*x_u[i]*x_u[j]*pe.cos(x_theta[i]-x_theta[j])\
                                           - ybra[3].imag*x_u[i]*x_u[j]*pe.sin(x_theta[i]-x_theta[j])
                                           
                            elif type_iz == 'mfq':
                                   k = self.s_c.nd.bus_mfq[i_m]
                                   ybra,i,j = self._find_k(k)
                                   cons = -(ybra[2].imag+ybra[3].imag)*x_u[i]**2 + ybra[3].imag*x_u[i]*x_u[j]*pe.cos(x_theta[i]-x_theta[j])\
                                          - ybra[3].real*x_u[i]*x_u[j]*pe.sin(x_theta[i]-x_theta[j])
                                              
                            elif type_iz == 'mtp':
                                   k = self.s_c.nd.bus_mtp[i_m]
                                   ybra,i,j = self._find_k(k)
                                   cons = (ybra[4].real+ybra[3].real)*x_u[j]**2 - ybra[3].real*x_u[j]*x_u[i]*pe.cos(x_theta[j]-x_theta[i])\
                                          - ybra[3].imag*x_u[j]*x_u[i]*pe.sin(x_theta[j]-x_theta[i])                                  
                                              
                            elif type_iz == 'mtq':
                                   k = self.s_c.nd.bus_mtq[i_m]
                                   ybra,i,j = self._find_k(k)
                                   cons = -(ybra[4].imag+ybra[3].imag)*x_u[j]**2 + ybra[3].imag*x_u[j]*x_u[i]*pe.cos(x_theta[j]-x_theta[i])\
                                        - ybra[3].real*x_u[j]*x_u[i]*pe.sin(x_theta[j]-x_theta[i])        
                                        
                            elif type_iz == 'mps':
                                   i = self.s_c.nd.busdc_mps[i_m]
                                   i_ac = self.s_c.nd.busdc_ij[i]
                                   cons = -x_u[i_ac]**2 * self.s_c.nd.g_sc[i] + x_u[i_ac] * x_uc[i] * (self.s_c.nd.g_sc[i] * pe.cos(x_theta[i_ac]-x_thetac[i]) + \
                                          self.s_c.nd.b_sc[i] * pe.sin(x_theta[i_ac]-x_thetac[i]))
       
                            elif type_iz == 'mqs':
                                   i = self.s_c.nd.busdc_mqs[i_m]
                                   i_ac = self.s_c.nd.busdc_ij[i]
                                   cons = x_u[i_ac]**2 * self.s_c.nd.b_sc[i] + x_u[i_ac] * x_uc[i] * (self.s_c.nd.g_sc[i] * pe.sin(x_theta[i_ac]-x_thetac[i]) - \
                                          self.s_c.nd.b_sc[i] * pe.cos(x_theta[i_ac]-x_thetac[i]))
       
                            elif type_iz == 'mpc':
                                   i = self.s_c.nd.busdc_mpc[i_m]
                                   i_ac = self.s_c.nd.busdc_ij[i]
                                   cons = x_uc[i]**2 * self.s_c.nd.g_sc[i] - x_u[i_ac] * x_uc[i] * (self.s_c.nd.g_sc[i] * pe.cos(x_thetac[i]-x_theta[i_ac]) + \
                                          self.s_c.nd.b_sc[i] * pe.sin(x_thetac[i]-x_theta[i_ac]))
       
                            elif type_iz == 'mqc':
                                   i = self.s_c.nd.busdc_mqc[i_m]
                                   i_ac = self.s_c.nd.busdc_ij[i]
                                   cons = -x_uc[i]**2 * self.s_c.nd.b_sc[i] - x_u[i_ac] * x_uc[i] * (self.s_c.nd.g_sc[i] * pe.sin(x_thetac[i]-x_theta[i_ac]) - \
                                          self.s_c.nd.b_sc[i] * pe.cos(x_thetac[i]-x_theta[i_ac]))
       
                            elif type_iz == 'mpcinj':
                                   i = self.s_c.nd.busdc_mpcinj[i_m]
                                   i_ac = self.s_c.nd.busdc_ij[i]
                                   pc = x_uc[i]**2 * self.s_c.nd.g_sc[i] - x_u[i_ac] * x_uc[i] * (self.s_c.nd.g_sc[i] * pe.cos(x_theta[i_ac]-x_thetac[i]) - \
                                   self.s_c.nd.b_sc[i] * pe.sin(x_theta[i_ac]-x_thetac[i]))
                                   ic = (((self.s_c.nd.g_sc[i]**2 + self.s_c.nd.b_sc[i]**2) * (x_u[i_ac]**2 + x_uc[i]**2 - 2*x_u[i_ac]*x_uc[i]*pe.cos(x_thetac[i]-x_theta[i_ac])))**0.5)/(3**0.5)
                                   if i==0:
                                          ploss = self.s_c.nd.lossa[i] + self.s_c.nd.lossb[i] * ic + self.s_c.nd.losscr[i] * ic**2
                                          pdc   = x_udci * x_idci
                                   elif i==1:
                                          ploss = self.s_c.nd.lossa[i] + self.s_c.nd.lossb[i] * ic + self.s_c.nd.lossci[i] * ic**2
                                          pdc   = (x_udci - x_idci * self.s_c.nd.r_dc) * (-x_idci)
                                   cons = pc + ploss + pdc
       
                            elif type_iz == 'mudc':
                                   i = self.s_c.nd.busdc_mudc[i_m]
                                   if i == 0:
                                          cons = x_udci
                                   elif i == 1:
                                          cons = x_udci - x_idci * self.s_c.nd.r_dc
                   
                            elif type_iz == 'midc':
                                   i = self.s_c.nd.busdc_midc[i_m]
                                   if i == 0:
                                          cons = x_idci
                                   elif i == 1:
                                          cons = -x_idci      
                   
                            elif type_iz == 'mm':
                                   i = self.s_c.nd.busdc_mm[i_m]
                                   i_ac = self.s_c.nd.busdc_ij[i]
                                   if i==0:
                                          cons = 2**0.5 * x_uc[i]/x_udci
                                   elif i==1:
                                          cons = 2**0.5 * x_uc[i]/(x_udci - x_idci * self.s_c.nd.r_dc)
                                          
                            elif type_iz == 'mtheta':
                                  i = self.s_c.nd.busdc_mtheta[i_m]
                                  i_ac = self.s_c.nd.busdc_ij[i]
                                  cons = x_thetac[i] - x_theta[i_ac]
       
                            model.con.add( cons == self.s_c.nd.z_true[i_z])

              instance = model
              opt = SolverFactory('ipopt')
              try:
                     result  = opt.solve(instance)
                     a = result['Solver'][0]['Termination condition'].key
                     x_a = [model.x[ii]() for ii in ind_x]
                     x_a = np.array(x_a)
                     if a == 'infeasible':
                            self.x_a = np.array([])
                     else:
                            self.x_a = x_a
              except IOError:
                     self.x_a = np.array([])
              except ValueError:
                     self.x_a = np.array([])
                     
       def _iz_info(self, i_z):
              """According to i_z, the index of z, compute the measurement type of i_z, and the corresponding index of i_z in **_m** """
              m_iz = dict()
              a = list()
              m_iz['mv'] = range(0,               self.s_c.nd.bus_mv.size)
              a.extend(m_iz['mv'])
              m_iz['mp'] = range(len(a), len(a) + self.s_c.nd.bus_mp.size)
              a.extend(m_iz['mp'])
              m_iz['mq'] = range(len(a), len(a) + self.s_c.nd.bus_mq.size)
              a.extend(m_iz['mq'])              
              m_iz['mfp'] = range(len(a), len(a) + self.s_c.nd.branch_mfp.size)
              a.extend(m_iz['mfp'])              
              m_iz['mfq'] = range(len(a), len(a) + self.s_c.nd.branch_mfq.size)
              a.extend(m_iz['mfq'])              
              m_iz['mtp'] = range(len(a), len(a) + self.s_c.nd.branch_mtp.size)
              a.extend(m_iz['mtp'])              
              m_iz['mtq'] = range(len(a), len(a) + self.s_c.nd.branch_mtq.size)
              a.extend(m_iz['mtq'])              
              m_iz['mps'] = range(len(a), len(a) + self.s_c.nd.busdc_mps.size)
              a.extend(m_iz['mps'])
              m_iz['mqs'] = range(len(a), len(a) + self.s_c.nd.busdc_mqs.size)
              a.extend(m_iz['mqs'])
              m_iz['mpc'] = range(len(a), len(a) + self.s_c.nd.busdc_mpc.size)
              a.extend(m_iz['mpc'])
              m_iz['mqc'] = range(len(a), len(a) + self.s_c.nd.busdc_mqc.size)
              a.extend(m_iz['mqc'])
              m_iz['mpcinj'] = range(len(a), len(a) + self.s_c.nd.busdc_mpcinj.size)
              a.extend(m_iz['mpcinj'])
              m_iz['mudc'] = range(len(a), len(a) + self.s_c.nd.busdc_mudc.size)
              a.extend(m_iz['mudc'])              
              m_iz['midc'] = range(len(a), len(a) + self.s_c.nd.busdc_midc.size)
              a.extend(m_iz['midc'])
              m_iz['mm'] = range(len(a), len(a) + self.s_c.nd.busdc_mm.size)
              a.extend(m_iz['mm'])
              m_iz['mtheta'] = range(len(a), len(a) + self.s_c.nd.busdc_mtheta.size)
              a.extend(m_iz['mtheta'])

              for key in m_iz.keys():
                     if i_z in m_iz[key]:
                            type_iz = key
                            i_m = m_iz[key].index(i_z)
                            break
              return type_iz, i_m

       def _find_k(self,k):
              """k-> ybra, b, i, j"""  
              ybra = self.s_c.nd.Ybra[k]
              i =np.where(self.s_c.nd.casedata['bus'][:,0]==self.s_c.nd.casedata['branch'][k,0])[0][0]
              j =np.where(self.s_c.nd.casedata['bus'][:,0]==self.s_c.nd.casedata['branch'][k,1])[0][0]
              return ybra, i, j
               
       def _f_sub1(self):
              """ """
              self.M = 3000
              self.k, self.t = 1, 1
              self.w = self._w()
              self.n = dict()
              self.n[self.k] = 1
              self.l_min = self.M
              self.IF = dict()
              self.IV = dict()
              self.B = dict()
              self.v = dict()
              self.x_at = dict()
              self.L = dict()
              self.Z = dict()          
              self._IF()
              self._IV()
              self._B()
              self._v()
              
       
       def _f_sub2(self, flag3):
              """ """
              if flag3 == 0:
                     self.L[self.t] = self.M

              elif flag3 == 1:
                     self.x_at[self.t] = self.x_a
                     self._Z_L()
                     
              self.t = self.t + 1
              self.l_min = min(self.L.values())
              while True:
                     
                     if self.k > 1 and self.n[self.k] == 2**(self.IV[self.k].size) - 1:
                            self.k = self.k - 1
                            continue
                            self.n[self.k] = self.n[self.k] + 1
                            self._v()
                            flag2 = 0
                     elif self.k == 1 and self.n[self.k] == 2**(self.IV[self.k].size) - 1:
                            flag2 = 1 #finish
                            break
                     else:
                            self.n[self.k] = self.n[self.k] + 1
                            self._v()
                            flag2 = 0
                            
                     v_sum = self.v[1,self.n[1]]
                     for i in range(2,self.k+1):
                            v_sum = v_sum + self.v[self.k,self.n[self.k]]
                     w_v = v_sum * self.w
                     Z = np.nonzero(w_v)[0]
                     L = Z.size   
                     print self.k, self.n[self.k]
                     if L > self.l_min:
                            continue
                     else:
                            break
                            
              return flag2
              
       def _f_sub3(self, flag1):
              """ """
              label = 0
              if flag1 == 1:
                     self._equa(1)  # solve the equations 1
                     if self.x_a.size == 0:
                            flag3 = 0
                     else:
                            if (self.v[self.k, self.n[self.k]] * self.w).sum() == self.v[self.k, self.n[self.k]].sum():
                                   flag3 = 1
                            else:
                                   self._equa(2)
                                   if self.x_a.size != 0:
                                          flag3 = 1
                                   else:
                                          label = 1
                                          self.k = self.k + 1
                                          self.n[self.k] = 1
                                          self._IF()  # k > 1
                                          self._IV()  # k > 1
                                          self._B()  # k > 1
                 
              if label == 1 or flag1 == 0:
                     while True:
                            if label == 1:
                                   self._v()
                            self._equa(3)                     
                     
                            if self.x_a.size == 0:
                                   flag3 = 0
                                   break
                            else:
                                   if (self.v[self.k, self.n[self.k]] * self.w).sum() == self.v[self.k, self.n[self.k]].sum():
                                          flag3 = 1
                                          break
                                   else:
                                          self._equa(2)
                                          if self.x_a.size != 0:
                                                 flag3 = 1
                                                 break
                                          else:
                                                 label = 1
                                                 self.k = self.k + 1
                                                 self.n[self.k] = 1
                                                 self._IF()  # k > 1
                                                 self._IV()  # k > 1
                                                 self._B()  # k > 1
                                                 continue
              return flag3
              
       def _f_sub4(self):
              """ """
              if min(self.L.values()) == self.M:
                     self.z_a = np.array([])
                     self.success = 0
              else:
                     self.success = 1
                     L_value_n = np.array(self.L.values())
                     L_key_n = np.array(self.L.keys())
                     self.t_sub = L_key_n[np.where(L_value_n==L_value_n.min())[0]]
                     self.x_n2 = {key: np.linalg.norm(self.x_at[key]- self.s_c.x) for key in self.t_sub } 
                     x_n2_value = np.array(self.x_n2.values())
                     x_n2_key = np.array(self.x_n2.keys())
                     self.t_a = x_n2_key[np.where(x_n2_value==x_n2_value.min())[0]][0]  #[0] - use the first t_a if t_a is more than one, but this situetion is very uncommon
                     
                     self.x_a_optimal = self.x_at[self.t_a]
                     self.Z_a = self.Z[self.t_a]

       def _nd_a(self):
               """# generate nd_a which is forming by replacing nd.z by z_a"""
               self.nd_a = copy.deepcopy(self.s_c.nd)
               
               x_u = self.x_a_optimal[0:self.nd_a.bus_num]  
               x_theta = np.r_[[0], self.x_a_optimal[self.nd_a.bus_num:2*self.nd_a.bus_num-1]]  #include reference bus
               x_uc = self.x_a_optimal[2*self.nd_a.bus_num-1:2*self.nd_a.bus_num+1]  
               x_thetac = self.x_a_optimal[2*self.nd_a.bus_num+1:2*self.nd_a.bus_num+3] 
               x_udci = self.x_a_optimal[2*self.nd_a.bus_num+3:2*self.nd_a.bus_num+4] 
               x_idci = self.x_a_optimal[2*self.nd_a.bus_num+4:2*self.nd_a.bus_num+5] 
               h = np.array([])
               # measurenment function of buses
               for i in self.nd_a.bus_mv:
                   h = np.append(h, x_u[i])
               for i in self.nd_a.bus_mp:
                   if i in self.nd_a.busdc_ij:
                          i_dc = np.where(self.nd_a.busdc_ij==i)[0][0] # the index in busdc of the bus i, = 0 or 1
                          ps = -x_u[i]**2 * self.nd_a.g_sc[i_dc] + x_u[i] * x_uc[i_dc] * (self.nd_a.g_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]) + \
                          self.nd_a.b_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]))
                   else:
                          ps = 0   
                   h = np.append(h, (x_u[i] * x_u * (self.nd_a.Ybus.real[i,:]*np.cos(x_theta[i]-x_theta) + \
                   self.nd_a.Ybus.imag[i,:]*np.sin(x_theta[i]-x_theta))).sum() - ps)
               for i in self.nd_a.bus_mq:
                   if i in self.nd_a.busdc_ij:
                          i_dc = np.where(self.nd_a.busdc_ij==i)[0][0] # the index in busdc of the bus i, = 0 or 1
                          qs = x_u[i]**2 * self.nd_a.b_sc[i_dc] + x_u[i] * x_uc[i_dc] * (self.nd_a.g_sc[i_dc] * np.sin(x_theta[i]-x_thetac[i_dc]) - \
                          self.nd_a.b_sc[i_dc] * np.cos(x_theta[i]-x_thetac[i_dc]))
                   else:
                          qs = 0 
                   h = np.append(h, (x_u[i] * x_u * (self.nd_a.Ybus.real[i,:]*np.sin(x_theta[i]-x_theta) - \
                   self.nd_a.Ybus.imag[i,:]*np.cos(x_theta[i]-x_theta))).sum() - qs)        
               # measurement function of branches
               for k in self.nd_a.branch_mfp:
                   ybra,i,j = self._find_k(k)  
                   h = np.append(h, (ybra[2].real+ybra[3].real)*x_u[i]**2 - ybra[3].real*x_u[i]*x_u[j]*np.cos(x_theta[i]-x_theta[j])\
                   - ybra[3].imag*x_u[i]*x_u[j]*np.sin(x_theta[i]-x_theta[j]))
               for k in self.nd_a.branch_mfq:
                   ybra,i,j = self._find_k(k)         
                   h = np.append(h, -(ybra[2].imag+ybra[3].imag)*x_u[i]**2 + ybra[3].imag*x_u[i]*x_u[j]*np.cos(x_theta[i]-x_theta[j])\
                   - ybra[3].real*x_u[i]*x_u[j]*np.sin(x_theta[i]-x_theta[j]))
               for k in self.nd_a.branch_mtp:
                   ybra,i,j = self._find_k(k)            
                   h = np.append(h, (ybra[4].real+ybra[3].real)*x_u[j]**2 - ybra[3].real*x_u[j]*x_u[i]*np.cos(x_theta[j]-x_theta[i])\
                   - ybra[3].imag*x_u[j]*x_u[i]*np.sin(x_theta[j]-x_theta[i]))
               for k in self.nd_a.branch_mtq:
                   ybra,i,j = self._find_k(k)           
                   h = np.append(h, -(ybra[4].imag+ybra[3].imag)*x_u[j]**2 + ybra[3].imag*x_u[j]*x_u[i]*np.cos(x_theta[j]-x_theta[i])\
                   - ybra[3].real*x_u[j]*x_u[i]*np.sin(x_theta[j]-x_theta[i]))
                   
               # measurement function of dc system
               for i in self.nd_a.busdc_mps:
                   i_ac = self.nd_a.busdc_ij[i]  #the index in bus of the busdc i
                   h = np.append(h,-x_u[i_ac]**2 * self.nd_a.g_sc[i] + x_u[i_ac] * x_uc[i] * (self.nd_a.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) + \
                   self.nd_a.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i])))
               for i in self.nd_a.busdc_mqs:
                   i_ac = self.nd_a.busdc_ij[i]  #the index in bus of the busdc i
                   h = np.append(h, x_u[i_ac]**2 * self.nd_a.b_sc[i] + x_u[i_ac] * x_uc[i] * (self.nd_a.g_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]) - \
                   self.nd_a.b_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i])))
               for i in self.nd_a.busdc_mpc:
                   i_ac = self.nd_a.busdc_ij[i]
                   h = np.append(h, x_uc[i]**2 * self.nd_a.g_sc[i] - x_u[i_ac] * x_uc[i] * (self.nd_a.g_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac]) + \
                   self.nd_a.b_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac])))
               for i in self.nd_a.busdc_mqc:
                   i_ac = self.nd_a.busdc_ij[i]
                   h = np.append(h, -x_uc[i]**2 * self.nd_a.b_sc[i] - x_u[i_ac] * x_uc[i] * (self.nd_a.g_sc[i] * np.sin(x_thetac[i]-x_theta[i_ac]) - \
                   self.nd_a.b_sc[i] * np.cos(x_thetac[i]-x_theta[i_ac])))
                   
               for i in self.nd_a.busdc_mpcinj:
                   i_ac = self.nd_a.busdc_ij[i]
                   pc = x_uc[i]**2 * self.nd_a.g_sc[i] - x_u[i_ac] * x_uc[i] * (self.nd_a.g_sc[i] * np.cos(x_theta[i_ac]-x_thetac[i]) - \
                   self.nd_a.b_sc[i] * np.sin(x_theta[i_ac]-x_thetac[i]))
                   ic = (((self.nd_a.g_sc[i]**2 + self.nd_a.b_sc[i]**2) * (x_u[i_ac]**2 + x_uc[i]**2 - 2*x_u[i_ac]*x_uc[i]*np.cos(x_thetac[i]-x_theta[i_ac])))**0.5)/(3**0.5)
                   if i==0:
                          ploss = self.nd_a.lossa[i] + self.nd_a.lossb[i] * ic + self.nd_a.losscr[i] * ic**2
                          pdc   = x_udci * x_idci
                   elif i==1:
                          ploss = self.nd_a.lossa[i] + self.nd_a.lossb[i] * ic + self.nd_a.lossci[i] * ic**2
                          pdc   = (x_udci - x_idci * self.nd_a.r_dc) * (-x_idci)
                   h = np.append(h, pc + ploss + pdc)
            
               for i in self.nd_a.busdc_mudc:
                   if i==0:
                          h = np.append(h, x_udci)
                   elif i==1:
                          h = np.append(h, x_udci - x_idci * self.nd_a.r_dc)
               for i in self.nd_a.busdc_midc:
                   if i==0:
                          h = np.append(h, x_idci)
                   elif i==1:
                          h = np.append(h, -x_idci)
               for i in self.nd_a.busdc_mm:
                   i_ac = self.nd_a.busdc_ij[i]
                   if i==0:
                          h = np.append(h, 2**0.5 * x_uc[i]/x_udci)
                   elif i==1:
                          h = np.append(h, 2**0.5 * x_uc[i]/(x_udci - x_idci * self.nd_a.r_dc))
               for i in self.nd_a.busdc_mtheta:
                   i_ac = self.nd_a.busdc_ij[i]
                   h = np.append(h, x_thetac[i] - x_theta[i_ac])
               self.h = h
               #self.Z_a = [ 19,  33, 122, 124, 126, 128, 136, 138]
               self.nd_a.z[self.Z_a] = h[self.Z_a] 

       def main(self):
              """ """
              self._f_sub1()
              while True:
                     #print self.k, self.n[self.k]
                     if self.k == 1:
                            flag1 = 1
                     else:
                            flag1 = 0
                     flag3 = self._f_sub3(flag1)
                     flag2 = self._f_sub2(flag3)
                     if flag2 == 0:
                            continue
                     else:
                            break
              self._f_sub4()
              if self.success == 1:
                     self._nd_a()
              else:
                     self.nd_a = copy.deepcopy(self.s_c.nd)
