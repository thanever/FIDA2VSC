# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 12:06:54 2016

@author: dell
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 15:48:43 2016

Static secuirty analysis of VSC-HVDC

@author: dell
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
myfont = matplotlib.font_manager.FontProperties(fname=r'C:/Windows/Fonts/simsun.ttc')
plt.rc('font', size=9,family='STIXGeneral')
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

class StaticSA:
       
       def __init__(self, w):
              
              self.w = w
              self.g_sc = w.nd.g_sc
              self.b_sc = w.nd.b_sc
              self.i_c_max = w.nd.casedata['convdc'][:,14]
              self.u_c_max = w.nd.casedata['convdc'][:,12]
            
              self.x_u = w.x_u
              self.x_theta = w.x_theta
              self.x_uc = w.x_uc
              self.x_thetac = w.x_thetac

              self.x_u = w.x[0:w.nd.bus_num]  
              self.x_theta = np.r_[[0], w.x[w.nd.bus_num:2*w.nd.bus_num-1]]  #include reference bus
              self.x_uc = w.x[2*w.nd.bus_num-1:2*w.nd.bus_num+1]  
              self.x_thetac = w.x[2*w.nd.bus_num+1:2*w.nd.bus_num+3] 
              self.x_udci = w.x[2*w.nd.bus_num+3:2*w.nd.bus_num+4] 
              self.x_idci = w.x[2*w.nd.bus_num+4:2*w.nd.bus_num+5] 

              self.busdc_ij = w.nd.busdc_ij
              
              # for paper plot , add in 2017.12.27
              self.path_f1 = os.path.join(os.getcwd(),'result//f1.png')
              self.path_f2 = os.path.join(os.getcwd(),'result//f2.png')
              self.path_f3 = os.path.join(os.getcwd(),'result//f3.png')
              self.path_r = [os.path.join(os.getcwd(),'result//'+str(int(i+1))+'.txt') for i in range(3)]

       def _pq_parameter(self):
              
              self.pc, self.qc, self.c_ic, self.r_ic, self.c_uc, self.r_uc = [np.array([])]*6
              for i in [0,1]:
                     i_ac = self.busdc_ij[i]

                     self.pc = np.append(self.pc, self.x_uc[i]**2 * self.g_sc[i] - self.x_u[i_ac] * self.x_uc[i] \
                                              * (self.g_sc[i] * np.cos(self.x_thetac[i]-self.x_theta[i_ac]) + \
                                              self.b_sc[i] * np.sin(self.x_thetac[i]-self.x_theta[i_ac])))
                     self.qc = np.append(self.qc, -self.x_uc[i]**2 * self.b_sc[i] - self.x_u[i_ac] * self.x_uc[i] \
                                              * (self.g_sc[i] * np.sin(self.x_thetac[i]-self.x_theta[i_ac]) - \
                                              self.b_sc[i] * np.cos(self.x_thetac[i]-self.x_theta[i_ac])))
                     self.c_ic = np.append(self.c_ic, [0, 0])
                     self.r_ic = np.append(self.r_ic, self.x_u[i_ac] * self.i_c_max[i])
                     self.c_uc = np.append(self.c_uc, [-self.x_u[i_ac]**2 * self.g_sc[i], self.x_u[i_ac]**2 * self.b_sc[i]])
                     self.r_uc = np.append(self.r_uc, self.x_u[i_ac] * self.u_c_max[i] * (self.g_sc[i]**2 + self.b_sc[i]**2)**0.5)             
       
       def _intersec(self, xx1,yy1,r1,xx2,yy2,r2):

              d = ((abs(xx1-xx2))**2 + ((yy1-yy2))**2)**0.5
              A = ((r2)**2 - (r1)**2 + (d)**2) / (2 * d)
              h = ((r2)**2 - (A)**2)**0.5
              
              x2 = xx2 + A * (xx1-xx2)/d
              y2 = yy2 + A * (yy1-yy2)/d
              
              x3 = x2 - h * ( yy1 - yy2 ) / d
              y3 = y2 + h * ( xx1 - xx2 ) / d
              
              x4 = x2 + h * (yy1 - yy2) / d
              y4 = y2 - h * (xx1 - xx2) / d
                                    
              if x3>x4:
                     x3,x4 = x4,x3
                     y3,y4 = y4,y3
              xb = np.linspace(x3,x4,100)
              yb1 = -(r1**2 - (xb-xx1)**2)**0.5 + yy1
              yb2 = (r2**2 - (xb-xx2)**2)**0.5 + yy2
              
              return xb, yb1,yb2           
              
       def pq_analysis(self):
              
              self._pq_parameter()
              
              self.angles_circle = [i*np.pi/180 for i in range(0,360)] 

              self.x_ic_0 = self.r_ic[0]*np.cos(self.angles_circle) + self.c_ic[0]
              self.y_ic_0 = self.r_ic[0]*np.sin(self.angles_circle) + self.c_ic[1]
              self.x_uc_0 = self.r_uc[0]*np.cos(self.angles_circle) + self.c_uc[0]
              self.y_uc_0 = self.r_uc[0]*np.sin(self.angles_circle) + self.c_uc[1]
   
              self.xb, self.yb1,self.yb2 = self._intersec(self.c_ic[0], self.c_ic[1], self.r_ic[0], self.c_uc[0], self.c_uc[1], self.r_uc[0])


              
       def compare(self):
              
              self.re = [np.loadtxt(self.path_r[i]) for i in range(3)]
              
              f1, axarr1 = plt.subplots(3,1)          
              f1.set_figheight(6.5625)
              f1.set_figwidth(5.25)
              
              for i in range(3):
                     axarr1[i].plot(range(100),self.re[i][:,0],'o',markersize = 3,label=u'$r_{\mathrm{max,c}}^{N}$:FDI攻击前的$r_{\mathrm{max}}^{N}$')
                     if i!=2:
                            axarr1[i].plot(range(100),self.re[i][:,1],'rx',markersize = 3,label=u'$r_{\mathrm{max,a}}^{N}$:FDI攻击后$r_{\mathrm{max}}^{N}$')
                     else:
                            axarr1[i].plot(range(100),self.re[i][:,1],'rx',markersize = 3,label=u'$r_{\mathrm{max,a}}^{N}$:FDI攻击后$r_{\mathrm{max}}^{N}$')
                            
                     axarr1[i].vlines(range(100),[0],self.re[i][:,3],label='$\Delta r_{\mathrm{max}}^{N}:r_{\mathrm{max,a}}^{N}-r_{\mathrm{max,c}}^{N}$')
                     axarr1[i].set_ylim(-1.0,3.5)
                     axarr1[i].set_ylabel('$r_{\mathrm{max}}^{N}$')

              axarr1[0].text(50,3.75,'$r_1=1.0,r_2=1.0$',horizontalalignment='center', verticalalignment='center')
              axarr1[1].text(50,3.75,'$r_1=0.9,r_2=0.9$',horizontalalignment='center', verticalalignment='center')
              axarr1[2].text(50,3.75,'$r_1=0.85,r_2=0.85$',horizontalalignment='center', verticalalignment='center')
              axarr1[2].text(50,-1.7,u'试验次数',horizontalalignment='center', verticalalignment='center',fontproperties=myfont)
              axarr1[i].legend(prop=myfont, fontsize=10.5, loc = [-0.01,-0.55],numpoints=1,labelspacing=0.01,frameon =False,ncol=2)
              plt.setp([a.get_xticklabels() for a in axarr1[[0,1]]], visible=False)
              f1.savefig(self.path_f2,dpi = 300, transparent=False, bbox_inches='tight')
 