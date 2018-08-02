# -*- coding: utf-8 -*-
"""
Created on Thu Jun 01 15:50:31 2017

@author: T.Han
"""



from distutils.core import setup
import py2exe
import sys

sys.setrecursionlimit(5000)
#this allows to run it with a simple double click.
sys.argv.append('py2exe')
 
py2exe_options = {
        "includes": ["sip"],
        "dll_excludes": ["MSVCP90.dll",],
        "compressed": 1,
        "optimize": 2,
        "ascii": 0,
        "bundle_files": 3,
        }
 
setup(
      name = 'PyQt Demo',
      version = '1.0',
      windows = ['main.py',], 
      zipfile = None,
      options = {'py2exe': py2exe_options}
      )