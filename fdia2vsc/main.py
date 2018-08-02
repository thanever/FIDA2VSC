# -*- coding: utf-8 -*-
"""
Created on Sat May 27 09:16:10 2017

@author: T.Han
"""
import os
import sys
from se import NetworkData, WlsSe
from time import sleep

package_1 = os.path.abspath(os.path.join(os.path.pardir, "pyui"))
os.sys.path.append(package_1)

import mainwindow
from PyQt4 import QtGui, QtCore
from cache_data import cache_data
import numpy as np
from pq import StaticSA
from ak import Attack


class FDIApp(QtGui.QMainWindow, mainwindow.Ui_MainWindow):
    def __init__(self, parent=None):

        QtGui.QMainWindow.__init__(self, parent)
        mainwindow.Ui_MainWindow.__init__(self)
        self.setupUi(self)


        self.matplotlibwidget_1.axes.axis('off')
        self.matplotlibwidget_2.axes.axis('off')
        self.matplotlibwidget_3.axes.axis('off')
        self.matplotlibwidget_4.axes.axis('off')
        self.matplotlibwidget_5.axes.axis('off')
        self.matplotlibwidget_5.axes.axis('off')
        self.matplotlibwidget_6.axes.axis('off')
        
        #parameters
        self.path_casedata =  os.path.abspath(os.path.join(os.path.dirname("__file__"), "extra_data/data/pf_result"))
        self.path_batch = os.path.abspath(os.path.join(os.path.dirname("__file__"), "extra_data/data/batch.txt"))

        
        self.buttonBrowse_1.clicked.connect(self.buttonBrowse_1_click)
        self.buttonBrowse_2.clicked.connect(self.buttonBrowse_2_click)
        self.buttonBrowse_3.clicked.connect(self.buttonBrowse_3_click)
        self.buttonBrowse_4.clicked.connect(self.buttonBrowse_4_click)
        self.pushButton_1.clicked.connect(self.pushButton_1_click)  # load data
        self.pushButton_2.clicked.connect(self.pushButton_2_click)  # plot sing-line diagram
        self.pushButton_3.clicked.connect(self.pushButton_3_click)  # generate measurement data
        self.pushButton_4.clicked.connect(self.pushButton_4_click)  # state estimation
        self.pushButton_5.clicked.connect(self.pushButton_5_click)  # plot r_N
        self.pushButton_6.clicked.connect(self.pushButton_6_click)  # plot P-Q
        self.pushButton_7.clicked.connect(self.pushButton_7_click)  # launch attack
        self.pushButton_8.clicked.connect(self.pushButton_8_click)  # state estimation after being attacked
        self.pushButton_9.clicked.connect(self.pushButton_9_click)  # plot r_N after being attacked
        self.pushButton_10.clicked.connect(self.pushButton_10_click)  # plot P-Q
        self.pushButton_11.clicked.connect(self.pushButton_11_click)  # batch analysis

        
            
    def buttonBrowse_1_click(self):
        path_ac = QtGui.QFileDialog.getOpenFileName(self, "Open File...", None,u"交流系统数据 (*.txt);;All Files (*)")
        print path_ac[0]
        self.editFileName_1.setText(path_ac) 

    def buttonBrowse_2_click(self):
        path_ac = QtGui.QFileDialog.getOpenFileName(self, "Open File...", None,u"VSC-HVDC系统数据 (*.txt);;All Files (*)")
        self.editFileName_2.setText(path_ac) 
        
    def buttonBrowse_3_click(self):
        path_ac = QtGui.QFileDialog.getOpenFileName(self, "Open File...", None,u"量测配置数据 (*.txt);;All Files (*)")
        self.editFileName_3.setText(path_ac) 

    def buttonBrowse_4_click(self):
        path_ac = QtGui.QFileDialog.getOpenFileName(self, "Open File...", None,u"节点地理位置数据 (*.txt);;All Files (*)")
        self.editFileName_4.setText(path_ac) 
    
    def pushButton_1_click(self):
        cache_data(self.editFileName_1.text(), self.editFileName_2.text(), self.editFileName_3.text(), self.editFileName_4.text())
        self.nd = NetworkData(self.path_casedata)
        self.display_data()
        
    def display_data(self):
        self.tableWidget_1.setRowCount(self.nd.casedata['busdc'].shape[0])
        for i in  range(self.nd.casedata['busdc'].shape[0]):
            self.tableWidget_1.verticalHeaderItem(i)
            for j in  range(self.nd.casedata['busdc'].shape[1]):
                if j in [0,1,2]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['busdc'][i,j])))
                    self.tableWidget_1.setItem(i, j, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['busdc'][i,j],4)))
                    self.tableWidget_1.setItem(i, j, newItem)
        self.tableWidget_2.setRowCount(self.nd.casedata['convdc'].shape[0])
        for i in  range(self.nd.casedata['convdc'].shape[0]):
            self.tableWidget_2.verticalHeaderItem(i)
            for j in  range(self.nd.casedata['convdc'].shape[1]):
                if j in [0,1,2]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['convdc'][i,j])))
                    self.tableWidget_2.setItem(i, j, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['convdc'][i,j],4)))
                    self.tableWidget_2.setItem(i, j, newItem)        
        self.tableWidget_3.setRowCount(self.nd.casedata['branchdc'].shape[0])
        for i in  range(self.nd.casedata['branchdc'].shape[0]):
            self.tableWidget_3.verticalHeaderItem(i)
            newItem = QtGui.QTableWidgetItem(str(int(i+1)))
            self.tableWidget_3.setItem(i, 0, newItem) 
            for j in  range(self.nd.casedata['branchdc'].shape[1]):                
                if j in [0,1]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['branchdc'][i,j])))
                    self.tableWidget_3.setItem(i, j+1, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['branchdc'][i,j],4)))
                    self.tableWidget_3.setItem(i, j+1, newItem)
        
        self.lineEdit_1.setText(str(self.nd.casedata['baseMVAac'][0,0])+'MW')
        self.lineEdit_2.setText(str(self.nd.casedata['baseMVAdc'][0,0])+'MW')
        if self.nd.casedata['pol'] ==1: self.comboBox_1.setCurrentIndex(0)
        else: self.comboBox_1.setCurrentIndex(1)
        self.tableWidget_4.setRowCount(self.nd.casedata['bus'].shape[0])
        for i in  range(self.nd.casedata['bus'].shape[0]):
            self.tableWidget_4.verticalHeaderItem(i)
            for j in  range(self.nd.casedata['bus'].shape[1]):
                if j in [0,1,6,10]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['bus'][i,j])))
                    self.tableWidget_4.setItem(i, j, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['bus'][i,j],4)))
                    self.tableWidget_4.setItem(i, j, newItem)        
        self.tableWidget_5.setRowCount(self.nd.casedata['gen'].shape[0])
        for i in  range(self.nd.casedata['gen'].shape[0]):
            self.tableWidget_5.verticalHeaderItem(i)
            for j in  range(self.nd.casedata['gen'].shape[1]):
                if j in [0]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['gen'][i,j])))
                    self.tableWidget_5.setItem(i, j, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['gen'][i,j],4)))
                    self.tableWidget_5.setItem(i, j, newItem)    
        self.tableWidget_6.setRowCount(self.nd.casedata['branch'].shape[0])
        for i in  range(self.nd.casedata['branch'].shape[0]):
            self.tableWidget_6.verticalHeaderItem(i)
            newItem = QtGui.QTableWidgetItem(str(int(i+1)))
            self.tableWidget_6.setItem(i, 0, newItem) 
            for j in  range(self.nd.casedata['branch'].shape[1]):                
                if j in [0,1]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['branch'][i,j])))
                    self.tableWidget_6.setItem(i, j+1, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['branch'][i,j],4)))
                    self.tableWidget_6.setItem(i, j+1, newItem)
        self.tableWidget_7.setRowCount(self.nd.casedata['bus_mc'].shape[0])
        for i in  range(self.nd.casedata['bus_mc'].shape[0]):
            self.tableWidget_7.verticalHeaderItem(i)
            for j in  range(self.nd.casedata['bus_mc'].shape[1]):
                if j in [0,1,3,5]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['bus_mc'][i,j])))
                    self.tableWidget_7.setItem(i, j, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['bus_mc'][i,j],6)))
                    self.tableWidget_7.setItem(i, j, newItem)  
        self.tableWidget_8.setRowCount(self.nd.casedata['branch_mc'].shape[0])
        for i in  range(self.nd.casedata['branch_mc'].shape[0]):
            self.tableWidget_8.verticalHeaderItem(i)
            for j in  range(self.nd.casedata['branch_mc'].shape[1]):
                if j in [0,1,2,4,6,8]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['branch_mc'][i,j])))
                    self.tableWidget_8.setItem(i, j, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['branch_mc'][i,j],6)))
                    self.tableWidget_8.setItem(i, j, newItem) 
        self.tableWidget_9.setRowCount(self.nd.casedata['busdc_mc'].shape[0])
        for i in  range(self.nd.casedata['busdc_mc'].shape[0]):
            self.tableWidget_9.verticalHeaderItem(i)
            for j in  range(self.nd.casedata['busdc_mc'].shape[1]):
                if j in [0,1,2,4,6,8,10,12,14,16,18]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['busdc_mc'][i,j])))
                    self.tableWidget_9.setItem(i, j, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['busdc_mc'][i,j],6)))
                    self.tableWidget_9.setItem(i, j, newItem)  
        self.tableWidget_10.setRowCount(self.nd.casedata['bl_ac'].shape[0])
        for i in  range(self.nd.casedata['bl_ac'].shape[0]):
            self.tableWidget_10.verticalHeaderItem(i)
            for j in  range(self.nd.casedata['bl_ac'].shape[1]):
                if j in [0]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['bl_ac'][i,j])))
                    self.tableWidget_10.setItem(i, j, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['bl_ac'][i,j],6)))
                    self.tableWidget_10.setItem(i, j, newItem) 
        self.tableWidget_11.setRowCount(self.nd.casedata['bl_dc'].shape[0])
        for i in  range(self.nd.casedata['bl_dc'].shape[0]):
            self.tableWidget_11.verticalHeaderItem(i)
            for j in  range(self.nd.casedata['bl_dc'].shape[1]):
                if j in [0]:
                    newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['bl_dc'][i,j])))
                    self.tableWidget_11.setItem(i, j, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.casedata['bl_dc'][i,j],6)))
                    self.tableWidget_11.setItem(i, j, newItem) 
                    
                    
    def pushButton_2_click(self):
        
        bl_ac = self.nd.casedata['bl_ac']
        bl_dc = self.nd.casedata['bl_dc']
        branch_ac = self.nd.casedata['branch'][:,0:2]
        branch_acdc = self.nd.casedata['busdc'][:,0:2]
        branch_dc = self.nd.casedata['branchdc'][:,0:2]
        
        self.matplotlibwidget_1.axes.hold(True)
        for i in range(branch_ac.shape[0]):
            i_1 = np.where(bl_ac[:,0]==branch_ac[i,0])[0][0]
            i_2 = np.where(bl_ac[:,0]==branch_ac[i,1])[0][0]
            self.matplotlibwidget_1.axes.plot(bl_ac[[i_1,i_2],1],bl_ac[[i_1,i_2],2], c='steelblue',linewidth=2)
            
            
        for i in range(branch_acdc.shape[0]):
            i_ac = np.where(bl_ac[:,0]==branch_acdc[i,1])[0][0]
            i_dc = np.where(bl_dc[:,0]==branch_acdc[i,0])[0][0]
            self.matplotlibwidget_1.axes.plot([bl_ac[i_ac,1], bl_dc[i_dc,1]], [bl_ac[i_ac,2], bl_dc[i_dc,2]], c='steelblue',linewidth=2)

            
        for i in range(branch_dc.shape[0]):
            i_1 = np.where(bl_dc[:,0]==branch_dc[i,0])[0][0]
            i_2 = np.where(bl_dc[:,0]==branch_dc[i,1])[0][0]
            self.matplotlibwidget_1.axes.plot(bl_dc[[i_1,i_2],1],bl_dc[[i_1,i_2],2], c='steelblue',linewidth=2)

            
        self.matplotlibwidget_1.axes.scatter(bl_ac[:,1], bl_ac[:,2],c='steelblue',s=160,edgecolor='black',alpha=1,linewidth=2)

        self.matplotlibwidget_1.axes.scatter(bl_dc[:,1], bl_dc[:,2],marker='s',c='steelblue',s=160,edgecolor='black')
 
        for i in range(bl_ac.shape[0]):
            self.matplotlibwidget_1.axes.text(bl_ac[i,1], bl_ac[i,2],str(int(bl_ac[i,0])),va='center',ha='center')
 
        for i in range(bl_dc.shape[0]):
            self.matplotlibwidget_1.axes.text(bl_dc[i,1], bl_dc[i,2],str(int(bl_dc[i,0])),va='center',ha='center')
 
        self.matplotlibwidget_1.axes.set_xticks([])
        self.matplotlibwidget_1.axes.set_yticks([])
        self.matplotlibwidget_1.axes.axis('off')
        self.matplotlibwidget_1.axes.set_position([0,0,1,1])
        self.matplotlibwidget_1.draw() 
        
    def pushButton_3_click(self):
        self.nd.get_z()
        c = 0
        self.tableWidget_12.setRowCount(self.nd.casedata['bus_mc'].shape[0])
        for i in  range(self.nd.casedata['bus_mc'].shape[0]):
            self.tableWidget_12.verticalHeaderItem(i)
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['bus_mc'][i,0])))
            self.tableWidget_12.setItem(i, 0, newItem)  
            
        for j in [1,3,5]:
            for i in range(self.nd.casedata['bus_mc'].shape[0]):
                if self.nd.casedata['bus_mc'][i,j]==0:
                    newItem = QtGui.QTableWidgetItem('-')
                    self.tableWidget_12.setItem(i, (j+1)/2, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.z[c],4)))           
                    self.tableWidget_12.setItem(i, (j+1)/2, newItem)  
                    c = c+1
        self.tableWidget_13.setRowCount(self.nd.casedata['branch_mc'].shape[0])    
        for i in  range(self.nd.casedata['branch_mc'].shape[0]):
            self.tableWidget_13.verticalHeaderItem(i)
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['branch_mc'][i,0])))
            self.tableWidget_13.setItem(i, 0, newItem)  
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['branch_mc'][i,1])))
            self.tableWidget_13.setItem(i, 1, newItem)
            
        for j in [2,4,6,8]:
            for i in range(self.nd.casedata['branch_mc'].shape[0]):
                if self.nd.casedata['branch_mc'][i,j]==0:
                    newItem = QtGui.QTableWidgetItem('-')
                    self.tableWidget_13.setItem(i, j/2+1, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.z[c],4)))
                    self.tableWidget_13.setItem(i, j/2+1, newItem)  
                    c = c+1            
        self.tableWidget_14.setRowCount(self.nd.casedata['busdc_mc'].shape[0])
        for i in  range(self.nd.casedata['busdc_mc'].shape[0]):
            self.tableWidget_14.verticalHeaderItem(i)
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['busdc_mc'][i,0])))
            self.tableWidget_14.setItem(i, 0, newItem)  
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['busdc_mc'][i,1])))
            self.tableWidget_14.setItem(i, 1, newItem)
            
        for j in [2,4,6,8,10,12,14,16,18]:
            for i in range(self.nd.casedata['busdc_mc'].shape[0]):
                if self.nd.casedata['busdc_mc'][i,j]==0:
                    newItem = QtGui.QTableWidgetItem('-')
                    self.tableWidget_14.setItem(i, j/2+1, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(self.nd.z[c],4)))
                    self.tableWidget_14.setItem(i, j/2+1, newItem)  
                    c = c+1  
    def pushButton_4_click(self):
        
        self.s_c = WlsSe(self.nd)
        self.s_c.wls()
        
    def pushButton_5_click(self):
        
        r_N = self.s_c.r_N
        self.matplotlibwidget_2.axes.clear()
        self.matplotlibwidget_2.axes.hold(True)
        self.matplotlibwidget_2.axes.vlines(range(len(r_N)),[0],r_N,linewidth=0.5)
        self.matplotlibwidget_2.axes.plot(range(len(r_N)), r_N, 'o', markersize=2)
        self.matplotlibwidget_2.axes.set_xlabel('Measurement')
        self.matplotlibwidget_2.axes.set_ylabel('Normalized Residual')
        self.matplotlibwidget_2.axes.grid(axis = 'y')
        self.matplotlibwidget_2.axes.axis('on')
        self.matplotlibwidget_2.axes.set_position([0.1,0.15,0.88,0.83])
        self.matplotlibwidget_2.draw()
        
    def pushButton_6_click(self):
        
        self.s = StaticSA(self.s_c)
        self.s.pq_analysis()
        
        self.matplotlibwidget_4.axes.clear()
        self.matplotlibwidget_4.axes.hold(True)
        self.matplotlibwidget_4.axes.plot(self.s.x_ic_0*100, self.s.y_ic_0*100, self.s.x_uc_0*100, self.s.y_uc_0*100,'--')
        self.matplotlibwidget_4.axes.plot(self.s.pc[0] * 100, self.s.qc[0] * 100, '.',c='r')
        self.matplotlibwidget_4.axes.plot(self.s.pc[0] * 100, self.s.qc[0] * 100, 'x',c='black')
        self.matplotlibwidget_4.axes.fill_between(self.s.xb*100, self.s.yb1*100,self.s.yb2*100, alpha=0.5,color='gray')

        self.matplotlibwidget_4.axes.set_xlim(-150,150)
        self.matplotlibwidget_4.axes.set_ylim(-150,150)
        self.matplotlibwidget_4.axes.set_xticks([-150,-100,-50,0,50,100,150])
        self.matplotlibwidget_4.axes.set_yticks([-150,-100,-50,0,50,100,150])
        self.matplotlibwidget_4.axes.set_xlabel('$P_{s1}$ $\mathrm{(MW)}$',fontsize=8)
        self.matplotlibwidget_4.axes.set_ylabel('$Q_{s1}$ $\mathrm{(MVar)}$',fontsize=8)

        self.matplotlibwidget_4.axes.grid(True)  
        self.matplotlibwidget_4.axes.set_position([0.2,0.15,0.75,0.8])
        self.matplotlibwidget_4.draw()
        
    def pushButton_7_click(self):
        
        self.mark_vsc_attacked()

        # launch attack
        self.ak = Attack(self.s_c, self.a_i-1 ,self.doubleSpinBox_1.value(), self.doubleSpinBox_2.value())
        self.ak.main()
        
        # display measurement data after being attacked
        self.display_z_a()

    def mark_vsc_attacked(self):
        self.a_i = int(self.comboBox_2.currentText())
        self.matplotlibwidget_1.axes.scatter(self.nd.casedata['bl_dc'][self.a_i-1,1], self.nd.casedata['bl_dc'][self.a_i-1,2],c='red', marker='o',s=500,edgecolor='red',alpha=0.4)
        self.matplotlibwidget_1.draw()
        
    def display_z_a(self):

        z = self.nd.z
        z_a = self.ak.nd_a.z
 
        c = 0
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        self.tableWidget_15.setRowCount(self.nd.casedata['bus_mc'].shape[0])
        for i in  range(self.nd.casedata['bus_mc'].shape[0]):
            self.tableWidget_15.verticalHeaderItem(i)
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['bus_mc'][i,0])))
            self.tableWidget_15.setItem(i, 0, newItem)  
            
        for j in [1,3,5]:
            for i in range(self.nd.casedata['bus_mc'].shape[0]):
                if self.nd.casedata['bus_mc'][i,j]==0:
                    newItem = QtGui.QTableWidgetItem('-')
                    self.tableWidget_15.setItem(i, (j+1)/2, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(z_a[c],4)))
                    if z[c]!=z_a[c]:
                        newItem.setBackground(brush)
                    self.tableWidget_15.setItem(i, (j+1)/2, newItem)  
                    c = c+1
        self.tableWidget_16.setRowCount(self.nd.casedata['branch_mc'].shape[0])        
        for i in  range(self.nd.casedata['branch_mc'].shape[0]):
            self.tableWidget_16.verticalHeaderItem(i)
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['branch_mc'][i,0])))
            self.tableWidget_16.setItem(i, 0, newItem)  
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['branch_mc'][i,1])))
            self.tableWidget_16.setItem(i, 1, newItem)
            
        for j in [2,4,6,8]:
            for i in range(self.nd.casedata['branch_mc'].shape[0]):
                if self.nd.casedata['branch_mc'][i,j]==0:
                    newItem = QtGui.QTableWidgetItem('-')
                    self.tableWidget_16.setItem(i, j/2+1, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(z_a[c],4)))
                    if z[c]!=z_a[c]:
                        newItem.setBackground(brush)                  
                    self.tableWidget_16.setItem(i, j/2+1, newItem)  
                    c = c+1            
        self.tableWidget_17.setRowCount(self.nd.casedata['busdc_mc'].shape[0])
        for i in  range(self.nd.casedata['busdc_mc'].shape[0]):
            self.tableWidget_17.verticalHeaderItem(i)
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['busdc_mc'][i,0])))
            self.tableWidget_17.setItem(i, 0, newItem)  
            newItem = QtGui.QTableWidgetItem(str(int(self.nd.casedata['busdc_mc'][i,1])))
            self.tableWidget_17.setItem(i, 1, newItem)
            
        for j in [2,4,6,8,10,12,14,16,18]:
            for i in range(self.nd.casedata['busdc_mc'].shape[0]):
                if self.nd.casedata['busdc_mc'][i,j]==0:
                    newItem = QtGui.QTableWidgetItem('-')
                    self.tableWidget_17.setItem(i, j/2+1, newItem)  
                else:
                    newItem = QtGui.QTableWidgetItem(str(round(z_a[c],4)))
                    if z[c]!=z_a[c]:
                        newItem.setBackground(brush)
                    self.tableWidget_17.setItem(i, j/2+1, newItem)  
                    c = c+1  
                    
    def pushButton_8_click(self):
        self.ak.s_c_a = WlsSe(self.ak.nd_a)
        self.ak.s_c_a.wls()
        
    def pushButton_9_click(self):
        
        r_N = self.ak.s_c_a.r_N
        self.matplotlibwidget_3.axes.clear()
        self.matplotlibwidget_3.axes.hold(True)
        self.matplotlibwidget_3.axes.vlines(range(len(r_N)),[0],r_N,linewidth=0.5)
        self.matplotlibwidget_3.axes.plot(range(len(r_N)), r_N, 'o', markersize=2)
        self.matplotlibwidget_3.axes.set_xlabel('Measurement')
        self.matplotlibwidget_3.axes.set_ylabel('Normalized Residual')
        self.matplotlibwidget_3.axes.grid(axis = 'y')
        self.matplotlibwidget_3.axes.axis('on')
        self.matplotlibwidget_3.axes.set_position([0.1,0.15,0.88,0.83])
        self.matplotlibwidget_3.draw()

    def pushButton_10_click(self):
        self.s_a = StaticSA(self.ak.s_c_a)
        self.s_a.pq_analysis()
        
        self.matplotlibwidget_5.axes.clear()
        self.matplotlibwidget_5.axes.hold(True)
        self.matplotlibwidget_5.axes.plot(self.s_a.x_ic_0*100, self.s_a.y_ic_0*100, self.s_a.x_uc_0*100, self.s_a.y_uc_0*100,'--')
        self.matplotlibwidget_5.axes.plot(self.s_a.pc[0] * 100, self.s_a.qc[0] * 100, '.',c='r')
        self.matplotlibwidget_5.axes.plot(self.s_a.pc[0] * 100, self.s_a.qc[0] * 100, 'x',c='black')
        self.matplotlibwidget_5.axes.fill_between(self.s_a.xb*100, self.s_a.yb1*100,self.s_a.yb2*100, alpha=0.5,color='gray')

        self.matplotlibwidget_5.axes.set_xlim(-150,150)
        self.matplotlibwidget_5.axes.set_ylim(-150,150)
        self.matplotlibwidget_5.axes.set_xticks([-150,-100,-50,0,50,100,150])
        self.matplotlibwidget_5.axes.set_yticks([-150,-100,-50,0,50,100,150])
        self.matplotlibwidget_5.axes.set_xlabel('$P_{s1}$ $\mathrm{(MW)}$',fontsize=8)
        self.matplotlibwidget_5.axes.set_ylabel('$Q_{s1}$ $\mathrm{(MVar)}$',fontsize=8)

        self.matplotlibwidget_5.axes.grid(True)  
        self.matplotlibwidget_5.axes.set_position([0.2,0.15,0.75,0.8])
        self.matplotlibwidget_5.draw()
        
    def pushButton_11_click(self):
        
        self.mark_vsc_attacked()
        self.batch_main()
        self.batch_plot()
        self.batch_table_data()
        
    def batch_main(self):
        
        f = open(self.path_batch,'w')
        self.n_batch = self.doubleSpinBox_3.value()
        self.progressBar_1.setMaximum(self.n_batch)
        i = 0
        
        while True:               
               self.nd.get_z()
               s_c = WlsSe(self.nd)
               s_c.wls()
               if s_c.r_N_max <= 3:
        
                      ak = Attack(s_c, self.a_i-1 ,self.doubleSpinBox_1.value(), self.doubleSpinBox_2.value())
                      ak.main()
        
                      ak.s_a = WlsSe(ak.nd_a)
                      ak.s_a.wls()
               
                      f.write(str(s_c.r_N_max)+'\t'+str(ak.s_a.r_N_max)+'\t'+str(min(ak.L.values()))+'\t'+str(ak.s_a.r_N_max-s_c.r_N_max)+'\t')
                      for xx in ak.s_a.x:
                             f.write(str(xx)+'\t')
                      f.write('\n')
                      i =i+1
               self.progressBar_1.setValue(i)  
               if i==self.n_batch:
                      break                 
        f.close()
        
    def batch_plot(self):
        self.re = np.loadtxt(self.path_batch)
        
        ind = self.re[:,3].argsort()
        self.re_sort = self.re[ind,:]
        
        self.matplotlibwidget_6.axes.hold(True)
        self.matplotlibwidget_6.axes.axis('on') 

        self.matplotlibwidget_6.axes.plot(range(int(self.n_batch)), self.re_sort[:,0],'o',markersize = 3,label='$r_{\mathrm{max,c}}^{N}$')
        self.matplotlibwidget_6.axes.plot(range(int(self.n_batch)), self.re_sort[:,1],'rx',markersize = 3,label='$r_{\mathrm{max,a}}^{N}$')
        self.matplotlibwidget_6.axes.vlines(range(int(self.n_batch)),[0], self.re_sort[:,3],label='$\Delta r_{\mathrm{max}}^{N}:r_{\mathrm{max,a}}^{N}-r_{\mathrm{max,c}}^{N}$')
        
        self.matplotlibwidget_6.axes.set_ylim(-1.0,4)
        self.matplotlibwidget_6.axes.set_xlim(-1,self.n_batch)   
        self.matplotlibwidget_6.axes.set_ylabel('$r_{\mathrm{max}}^{N}$') 
        self.matplotlibwidget_6.axes.set_xlabel('Number of Trials')
        self.matplotlibwidget_6.axes.legend(fontsize=8, loc='upper right',numpoints=1,labelspacing=0.001,frameon =True,ncol=3)
        self.matplotlibwidget_6.axes.set_position([0.1,0.15,0.88,0.83])
        self.matplotlibwidget_6.draw()
        
    def batch_table_data(self):
        n_row = self.re.shape[0]
        self.tableWidget_18.setRowCount(n_row)
        for i in range(n_row):
            self.tableWidget_18.verticalHeaderItem(i)

            newItem = QtGui.QTableWidgetItem(str(i+1))
            self.tableWidget_18.setItem(i, 0, newItem)
        c = 0
        for [i,j,k] in self.re[:,[0,1,3]]:
            newItem_1 = QtGui.QTableWidgetItem(str(round(i,4)))
            newItem_2 = QtGui.QTableWidgetItem(str(round(j,4)))
            newItem_3 = QtGui.QTableWidgetItem(str(round(k,4)))
            self.tableWidget_18.setItem(c, 1, newItem_1) 
            self.tableWidget_18.setItem(c, 2, newItem_2)
            self.tableWidget_18.setItem(c, 3, newItem_3) 
            c=c+1

        
def main():
    app = QtGui.QApplication(sys.argv)
    form = FDIApp()
    form.show()
#    app.exec_()
    sys.exit(app.exec_()) 
 
    
    
if __name__ == '__main__':
    main()
    
    
  