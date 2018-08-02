# -*- coding: utf-8 -*-
"""
Created on Sun May 28 19:55:39 2017

@author: T.Han
"""
from os import path  
import shutil
from mlab.releases import latest_release as matlab
import scipy.io as sio
import numpy as np

def cache_data(path_ac, path_dc, path_mc, path_bl):
    # geenerate cache data
    path_ac_m = path.abspath(path.join(path.dirname("__file__"), "extra_data/data/case_ac.m"))
    path_dc_m = path.abspath(path.join(path.dirname("__file__"), "extra_data/data/case_dc.m"))
    path_mc_m = path.abspath(path.join(path.dirname("__file__"), "extra_data/data/case_mc.m"))
    path_bl_ac_m = path.abspath(path.join(path.dirname("__file__"), "extra_data/data/case_bl_ac.txt"))
    path_bl_dc_m = path.abspath(path.join(path.dirname("__file__"), "extra_data/data/case_bl_dc.txt"))
    path_pfresult =  path.abspath(path.join(path.dirname("__file__"), "extra_data/data/pf_result.mat"))
    # ac grid data
    ac = open(path_ac,'rt')
    ac_m = open(path_ac_m,'wt')
    ac_m.write('function mpc = case_ac\n')
    ac_m.write('mpc.version = \'2\';\n')
    while True:
        l = ac.readline()
        if l=='START-BASEMVA\n':
            ac_m.write('mpc.baseMVA = '+ac.readline())

        if l=='START-BUS\n':
            ac_m.write('mpc.bus = [\n')
            while True:
                l = ac.readline()
                if l == 'END-BUS\n': break
                ac_m.write(l)
            ac_m.write('];\n')
    
        if l=='START-GENERATOR\n':
            ac_m.write('mpc.gen = [\n')
            while True:
                l = ac.readline()
                if l == 'END-GENERATOR\n': break
                ac_m.write(l)
            ac_m.write('];\n')
        
        if l=='START-BRANCH\n':
            ac_m.write('mpc.branch = [\n')
            while True:
                l = ac.readline()
                if l == 'END-BRANCH': break
                ac_m.write(l)
            ac_m.write('];\n')
            break  
    ac_m.close() 
    ac.close()         
    # VSC-HVDC grid data 
    dc = open(path_dc,'rt')
    dc_m = open(path_dc_m,'wt')
    dc_m.write('function [baseMVAac, baseMVAdc, pol, busdc, convdc, branchdc]= case_dc\n')
    while True:
        l = dc.readline()
        if l=='START-BASEMVAac\n':
            dc_m.write('baseMVAac = '+dc.readline())

        if l=='START-BASEMVAdc\n':
            dc_m.write('baseMVAdc = '+dc.readline())

        if l=='START-POL\n':
            dc_m.write('pol='+dc.readline())

        if l=='START-BUS\n':
            dc_m.write('busdc = [\n')
            while True:
                l = dc.readline()
                if l == 'END-BUS\n': break
                dc_m.write(l)
            dc_m.write('];\n')

        if l=='START-CONVERTER\n':
            dc_m.write('convdc = [\n')
            while True:
                l = dc.readline()
                if l == 'END-CONVERTER\n': break
                dc_m.write(l)
            dc_m.write('];\n') 

        if l=='START-BRANCH\n':
            dc_m.write('branchdc = [\n')
            while True:
                l = dc.readline()
                if l == 'END-BRANCH': break
                dc_m.write(l)
            dc_m.write('];\n')
            break
    dc_m.close()     
    dc.close()
    # measurement configurement data
    mc = open(path_mc,'rt')
    mc_m = open(path_mc_m,'wt')
    mc_m.write('function mpc = case_mc\n')    
    while True:
        l = mc.readline()

        if l=='START-BUS\n':
            mc_m.write('mpc.bus_mc = [\n')
            while True:
                l = mc.readline()
                if l == 'END-BUS\n': break
                mc_m.write(l)
            mc_m.write('];\n')    

        if l=='START-BRANCH\n':
            mc_m.write('mpc.branch_mc = [\n')
            while True:
                l = mc.readline()
                if l == 'END-BRANCH\n': break
                mc_m.write(l)
            mc_m.write('];\n')    
        
        if l=='START-BUSDC\n':
            mc_m.write('mpc.busdc_mc = [\n')
            while True:
                l = mc.readline()
                if l == 'END-BUSDC': break
                mc_m.write(l)
            mc_m.write('];\n')      
            break
    mc_m.close()     
    mc.close()
    
#    busac_loc = np.array([])
#    busdc_loc = 
    bl = open(path_bl,'rt')
    bl_ac = open(path_bl_ac_m,'wt')
    bl_dc = open(path_bl_dc_m,'wt')
    while True:
        l = bl.readline()
        if l == 'START-ACBUS\n':
            while True:
                l = bl.readline()
                if l == 'END-ACBUS\n': break
                bl_ac.write(l)
        if l == 'START-DCBUS\n':
            while True:
                l = bl.readline()
                if l == 'END-DCBUS\n': break
                bl_dc.write(l)
            break
 
    bl.close()
    bl_ac.close()
    bl_dc.close()
    
    bl_ac = np.loadtxt(path_bl_ac_m)
    bl_dc = np.loadtxt(path_bl_dc_m)
    
    
    YY = matlab.makeYbus_py()
    Ybus = YY.Ybus_i + YY.Ybus_i *1j
    Yf = YY.Yf_i + YY.Yf_j *1j
    Yt = YY.Yt_i + YY.Yt_j *1j
    resultsacdc = matlab.runacdcpf_py('case_ac','case_dc')
    result_ac = resultsacdc.results_ac
    result_dc = resultsacdc.results_dc
    mc = matlab.case_mc()   
    pf_result = dict()
    pf_result['bl_ac'] = bl_ac
    pf_result['bl_dc'] = bl_dc
    pf_result['bus_mc'] = mc.bus_mc
    pf_result['branch_mc']  = mc.branch_mc
    pf_result['busdc_mc']  = mc.busdc_mc
    pf_result['Ybus']  = Ybus
    pf_result['Yf']  = Yf
    pf_result['Yt']  = Yt
    pf_result['baseMVA']  = result_ac.baseMVA
    pf_result['bus']  = result_ac.bus
    pf_result['gen']  = result_ac.gen
    pf_result['branch']  = result_ac.branch
    pf_result['baseMVAac']  = result_dc.baseMVAac
    pf_result['baseMVAdc']  = result_dc.baseMVAdc
    pf_result['pol']  = result_dc.pol
    pf_result['busdc']  = result_dc.busdc
    pf_result['convdc']  = result_dc.convdc
    pf_result['branchdc']  = result_dc.branchdc         
    sio.savemat(path_pfresult,pf_result)



 