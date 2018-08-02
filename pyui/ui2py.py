# -*- coding: utf-8 -*-
"""
Created on Sat May 27 10:08:29 2017

@author: T.Han
"""

import os
from string import replace


path_ui = os.path.abspath(os.path.join(os.path.pardir, "ui")) 
path_qrc = os.path.abspath(os.path.join(os.path.pardir, "qrc")) 
path_pyui = os.path.abspath(os.path.dirname("__file__"))
path_pyqrc = os.path.abspath(os.path.dirname("__file__"))
 

name_ui = ['mainwindow']
name_qrc = ['images']

for s in name_qrc:
    qrc = os.path.join(path_qrc,s+'.qrc')
    pyqrc = os.path.join(path_pyqrc,s+'_rc.py')
    cmd_s = 'pyrcc4 -o '+ pyqrc + ' '+ qrc
    os.system(cmd_s)


for s in name_ui:
    ui = os.path.join(path_ui,s+'.ui')
    pyui = os.path.join(path_pyui,s+'.py')
    cmd_s = 'pyuic4 -o '+ pyui + ' '+ ui
    
    os.system(cmd_s)
    
#E:/0-se/fdia2acdc-1.0.0/pyui/image/1.png