#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:11:09 2022

@author: rachita
"""

import scipy
from scipy.optimize import minimize
from numpy import gradient, ndarray, diff, where, arange, argmin
import matplotlib.pyplot as plt
from larch.wxlib import plotlabels as plab
from ase import io
from ase.io import read, write
import wxmplot.interactive as wi
import numpy as np
from ase.io import read, write
import larch
from larch.io import read_athena, write_group, create_athena
from larch_plugins.xafs import feffdat
from larch import Group, isNamedClass
from larch.xafs import feffit, TransformGroup, FeffitDataSet, feffit_report, feffit_transform, pre_edge
from larch.fitting import guess, group2params, param_group
from larch import Group, Parameter
from larch.xafs import FeffPathGroup
from larch.xafs.feffit import (TransformGroup, FeffitDataSet, feffit, feffit_report)
from lmfit import Parameters, Parameter as param, Minimizer
from larch.xafs import feffrunner, feffpath, feff6l, feff8l
from wxmplot.interactive import plot
from larch.xafs import autobk
from larch.wxlib import plotlabels as plab
from larch.xafs import xftf
import scipy
from scipy.optimize import minimize
import os, shutil, sys
import glob
import scipy
from scipy.optimize import minimize
import pandas as pd

def read_experimental_data(filename,verbose=False, plot_expt = False):
    """
    Function to read in the experimental *.prj file
    """
    project = read_athena(filename)
    expt_data = {}
    for name, data in project._athena_groups.items():
        autobk(data.energy, data.mu, group=data, rbkg=1.0, kweight=3)
        expt_data[name] = data

    if verbose:
        for attr in dir(data):
            print(attr, type(getattr(data, attr)))
    if plot_expt:
        plt.plot(data.k, data.chi*data.k**2, label='$\chi$')
        plt.xlabel(r'$k\, ({\rm\AA}){-1}$')
        plt.ylabel(r'$k^2\chi, ({\rm\AA})^{-2}$')
        plt.legend()
    return expt_data

expt_csv = pd.read_csv('/Users/rachita/Library/CloudStorage/Box-Box/Larch_part3/site_fractions/_sim_data3/data0.10.csv')
expt_k = expt_csv['k']
expt_chi = expt_csv['chi']
expt_csv.k = expt_k
expt_csv.chi = expt_chi 

#expt_data = read_experimental_data('/Users/rachita/Library/CloudStorage/Box-Box/Larch_part3/site_fractions/_sim_data/Pt_MgO_Athena.prj')
#print(expt_data)
#Pd_MgO= expt_data['PtMgO400']
#Pt_MgO= expt_data['Pt14_50_50_smallNP_iso_10Kmerge16']
#print('Done reading the Athena Prj file or experimental data')
wd = os.getcwd()
strct_dirs = []

def read_structures ():
    wd_main = os.getcwd()
    folder_name = ['strct1', 'strct2']
    files_path = ['/Users/rachita/Library/CloudStorage/Box-Box/Larch_part3/site_fractions/_sim_data3/strct1/feff.inp',
                  '/Users/rachita/Library/CloudStorage/Box-Box/Larch_part3/site_fractions/_sim_data3/strct2/feff.inp']           
    for name, f in zip(folder_name, files_path):      
        if not os.path.exists(name):
            os.mkdir(name)
            shutil.copy(f, name)
            os.chdir(name)
            curr_dir = os.getcwd()
            strct_dirs.append(curr_dir)
            feff6l(feffinp='./feff.inp')
            os.chdir(wd_main)
               
    print(strct_dirs)       
    return (strct_dirs)
running_feff = read_structures()
###############
def exafs_fitting (initial_amp):
    
    """f1 =initial_amp
    f2=1-f1
    x_fract = [f1,f2]"""
    wd_main2 = os.getcwd()
    folder_name = ['strct1', 'strct2']
    main_list_pathargs = []
    for dirs, f_n in zip (strct_dirs, folder_name):
        curr_wd = os.chdir(dirs)
        if f_n == 'strct1':
            os.chdir(dirs)
            strct1_dat = glob.glob('files.dat')
            f = open('files.dat', "r")
            list_pathargs_1 = [] #
            list_nlegs = [] #
            paths = [] #
            dict_sig2s = {1:'s1_sig2_1',2:'s1_sig2_2'}
            #data = f.readlines()

           
            for line in f.readlines():
                i_line = line.strip()
                if 'feff' in i_line:
                    i_fname, i_sig2, i_amp, i_deg, i_nlegs, i_r_eff = i_line.split()
            
                    var_fname = i_fname
                    
                    """if float(i_r_eff) < 2.1: 
                        var_sig2_s1 = 's1_sig2_1' 
                    elif float(i_r_eff) > 2.1 and float(i_r_eff)<3.2: 
                        var_sig2_s1 = 's1_sig2_2' 
                    else:
                        pass"""
                    if float(i_r_eff) < 2.1: 
                        var_sig2_s1 = 's1_sig2_1' 
                    elif float(i_r_eff) > 2.1 and float(i_r_eff)<3.1: 
                        var_sig2_s1 = 's1_sig2_2' 
                    elif float(i_r_eff) > 3.1 and float(i_r_eff)<3.3: #and int(i_nlegs) == 2: 
                        var_sig2_s1 = 's1_sig2_3' 
                    elif float(i_r_eff) > 3.3 and float(i_r_eff)<4.2: #and int(i_nlegs)>2: 
                        var_sig2_s1 = 's1_sig2_4' 
                    else:# float(i_r_eff) > 4.5 and float(i_r_eff)<4.1: #and int(i_nlegs)>2: 
                        var_sig2_s1 = 's1_sig2_5'
                    #var_sig2_s1 = dict_sig2s[int(i_nlegs)]
                    var_degen_s1 = i_deg#'deg' #i_deg 
                    var_deltar_s1 = 's1_del_r*%s' % i_r_eff  #'del_r*reff' #del_r*%s' % i_r_eff 
                    #print(var_fname, var_sig2, var_degen, var_deltar)
                    
                    i_pathargs_1 = [var_fname, var_sig2_s1, var_degen_s1, var_deltar_s1]
                    list_pathargs_1.append(i_pathargs_1)
                    main_list_pathargs.append(list_pathargs_1)
                    #print(i_fname, i_nlegs, dict_sig2s[int(i_nlegs)])
                    list_nlegs.append(i_nlegs)
            f.close()
            num_paths = len(list_nlegs) #Total paths we care about 
            print('Total paths read = %s' % num_paths)
            os.chdir(wd_main2)
        elif f_n == 'strct2':
            os.chdir(dirs)
            strct2_dat = glob.glob('files.dat')
            f = open('files.dat', "r")
            list_pathargs_2 = [] #
            list_nlegs = [] #
            paths = [] #
            dict_sig2s = {1:'s2_sig2_1'}
            #data = f.readlines()

            
            for line in f.readlines():
                i_line = line.strip()
                if 'feff' in i_line:
                    i_fname, i_sig2, i_amp, i_deg, i_nlegs, i_r_eff = i_line.split()
            
                    var_fname = i_fname
                    
                    """if float(i_r_eff) < 2.8: 
                        var_sig2_s2 = 's2_sig2_1' 
                        var_cn1 = 's2_n1'
                    else:# float(i_r_eff) > 4.5 and float(i_r_eff)<4.1: #and int(i_nlegs)>2: 
                        pass"""
                    if float(i_r_eff) < 2.1: 
                        var_sig2_s2 = 's2_sig2_1' 
                    elif float(i_r_eff) > 2.1 and float(i_r_eff)<3.3: 
                        var_sig2_s2 = 's2_sig2_2' 
                    elif float(i_r_eff) > 3.3 and float(i_r_eff)<3.5: #and int(i_nlegs) == 2: 
                        var_sig2_s2 = 's2_sig2_3' 
                    elif float(i_r_eff) > 3.5 and float(i_r_eff)<4.0: #and int(i_nlegs)>2: 
                        var_sig2_s2 = 's2_sig2_4' 
                    elif float(i_r_eff) > 4.0 and float(i_r_eff)<4.5: #and int(i_nlegs)>2: 
                        var_sig2_s2 = 's2_sig2_5' 
                    
                    
                    else:# float(i_r_eff) > 4.5 and float(i_r_eff)<4.1: #and int(i_nlegs)>2: 
                        var_sig2_s2 = 's2_sig2_6'
                    #var_sig2_s2 = dict_sig2s[int(i_nlegs)]
                    
                    var_degen_s2 = i_deg #'deg' #i_deg # ******************************
                    
                    var_deltar_s2 = 's2_del_r*%s' % i_r_eff  #'del_r*reff' #del_r*%s' % i_r_eff 
                    #print(var_fname, var_sig2, var_degen, var_deltar)
                    i_pathargs_2 = [var_fname, var_sig2_s2, var_degen_s2, var_deltar_s2]
                    list_pathargs_2.append(i_pathargs_2)
                    main_list_pathargs.append(list_pathargs_2)
                    #print(i_fname, i_nlegs, dict_sig2s[int(i_nlegs)])
                    list_nlegs.append(i_nlegs)
            f.close()
            num_paths = len(list_nlegs) #Total paths we care about 
            print('Total paths read = %s' % num_paths)
            os.chdir(wd_main2)
    pars = param_group(s1_del_e0 = param(0.7, vary=True),
                       s2_del_e0 = param(0.7, vary=True),
                  #s2_n1 = param(9, min=1.0, max=12.0,  vary = True),
                  #s1_amp = param(0.5, vary=True, min=0.0, max=1.0),
                  #s2_amp = param(0.5, vary=True, min=0.0, max=1.0),
                  #s2_n2 = param(5, vary = True),
                  #s2_n3 = param(4, vary = True),
                  #s2_n4 = param(9, vary = True),
                  s1_sig2_1 = param(0.002, min=0.0, max=0.1, vary=True),
                  s1_sig2_2 = param(0.006, min=0.0, max=0.1, vary=True),
                  s1_sig2_3 = param(0.007, min=0.0, max=0.1, vary=True),
                  s1_sig2_4 = param(0.004, min=0.0, max=0.1, vary=True),
                  s1_sig2_5 = param(0.009, min=0.0, max=0.1, vary=True),
                  s2_sig2_1 = param(0.006, min=0.0, max=0.1, vary=True),
                  s2_sig2_2 = param(0.008, min=0.0, max=0.1, vary=True),
                  s2_sig2_3 = param(0.003, min=0.0, max=0.1, vary=True),
                  s2_sig2_4 = param(0.007, min=0.0, max=0.1, vary=True),
                  s2_sig2_5 = param(0.005, min=0.0, max=0.1, vary=True),
                  s2_sig2_6 = param(0.005, min=0.0, max=0.1, vary=True),
                  #s2_sig2_7 = param(0.005, min=0.0, max=0.1, vary=True),
                  s1_s02 = param(0.7, vary=True, min=0.0, max=1.0),
                  s2_s02 = param(0.2, vary=True, min=0.0, max=1.0),
                  #s1_s02 = 0.82,
                  s1_del_r = param(0.01, vary=True),
                  s2_del_r = param(0.02, vary=True))
    data_sets = []
    ## 
    for fld_name, pathargs, dirs in zip (folder_name, main_list_pathargs, strct_dirs):
        print (fld_name, dirs)
        if fld_name == 'strct1':
            os.chdir(dirs)
            for i, i_pathargs in enumerate(list_pathargs_1):
                var_fname, var_sig2_s1, var_degen_s1, var_deltar_s1 = i_pathargs
                pathargs = dict(e0='s1_del_e0', sigma2=var_sig2_s1, deltar=var_deltar_s1, s02= 's1_s02')
                #print (pathargs)
                var_full_fname = wd_main2+'/'+fld_name+'/'+var_fname
                #print(var_full_fname)
                paths.append(feffdat.feffpath(var_full_fname,**pathargs))
                #f_paths = int(x_fract[0])*paths
            print('********************^^^^^^^^^^^^^^^^^^^^***********************')
            #print (os.getcwd())
            os.chdir(wd_main2)
        elif fld_name == 'strct2':
            amp_2 = initial_amp-1 
            os.chdir(dirs)
            for i, i_pathargs in enumerate(list_pathargs_2):
                var_fname, var_sig2_s2, var_degen_s2, var_deltar_s2 = i_pathargs
                #pathargs = dict(e0='del_e0', sigma2=var_sig2_s2, deltar=var_deltar_s2, degen= 's2_n1*s02')
                #pathargs = dict(e0='s2_del_e0', sigma2=var_sig2_s2, deltar=var_deltar_s2)
                pathargs = dict(e0='s2_del_e0', sigma2=var_sig2_s2, deltar=var_deltar_s2, s02= 's2_s02')
                #print (pathargs)
                var_full_fname = wd_main2+'/'+fld_name+'/'+var_fname
                paths.append(feffdat.feffpath(var_full_fname,**pathargs))
                #f_paths = int(x_fract[1])*paths
            print('********************^^^^^^^^^^^^^^^^^^^^***********************')
            print (os.getcwd())
            os.chdir(wd_main2)
    ########################
    trans = TransformGroup(kmin =2.4, kmax=13.0, kw=3, dk=1, window = 'hanning', rmin=1.0, rmax=5.0)
    dset = FeffitDataSet(data=expt_csv, pathlist=paths, transform=trans)
    print ('running EXAFS fit... please be patient')
    out = feffit(pars, dset)
    report = feffit_report(out)
    file_name = 'report.txt'
    with open (file_name, 'w') as file:
        file.write(report)
    ############################           
    print(dset)
    print ('QuantEXAFS run successful, check report! ^-^')
    
    """x_frac = [f1,f2]
    answer= scipy.optimize.minimize(dset.data.r, x_frac)
    f1, f2 = answer.x
    print('The optimized fractions of sites are -')
    print(f1, f2)"""
    
    from scipy.interpolate import CubicSpline
    import similaritymeasures
    dataset = dset
    expt_x = dataset.data.r
    expt_y = dataset.data.chir_mag
    sim_x = dataset.model.r
    sim_y = dataset.model.chir_mag
    def compare_curves(sim_x = dataset.data.r, sim_y=dataset.data.chir_mag, expt_x=dataset.model.r, expt_y=dataset.model.chir_mag, x_range=None):
        """Returns Fréchet distance, area between curves, and RMSE
           for simulation and experimental spectra
    
           Input:
               sim_x = list of x values for simulation curve
    
               sim_y = list of y values for simulation curve
    
               expt_x = list of x values for experimental curve
    
               expt_y = list of y values for experimental curve
    
               x_range (optional) = [min_x, min_y], range of x values to use
        """
        if x_range is None: # uses full range if none supplied
            x_min = min(sim_x)
            x_max = max(sim_x)
        else:
            x_min, x_max = x_range
    
        if list(sim_x) != list(expt_x):
            cs = CubicSpline(expt_x, expt_y)
            expt_x = sim_x
            expt_y = cs(expt_x)
    
        sim = [[x, y] for x, y in zip(sim_x, sim_y) if x <= x_max and x >= x_min]
        expt = [[x, y] for x, y in zip(expt_x, expt_y) if x <= x_max and x >= x_min]
        sim, expt = np.array(sim), np.array(expt)
    
        frechet = similaritymeasures.frechet_dist(sim, expt)
        area = similaritymeasures.area_between_two_curves(sim, expt)
    
        mse = sum((sim[:, 1] - expt[:, 1])**2) / len(sim[:, 1])
        rmse = np.sqrt(mse)
        print ('*************************************************')
        
        return frechet, area, rmse
    compare = compare_curves(sim_x, sim_y, expt_x, expt_y, x_range=[1.0, 5.0])
    comp = f'frechet = {compare[0]}\narea = {compare[1]}\nRMSE = {compare[2]}\n'
    file_name = 'compare_curve.txt'
    with open (file_name, 'w') as file:
        file.write(comp)
    
    frechet, area, rmse = compare_curves(sim_x, sim_y, expt_x, expt_y, x_range=[1.0, 5.0])
    
    #creating seaborn plots of EXAFS expt data and model in k-space and R-space
    import seaborn as sns
    colors = sns.color_palette('muted')
    def plot_chifit(dataset, kmin=0, kmax=None, kweight=None, rmax=None,
                show_mag=True, show_real=True, show_imag=False,
                title='name', new=True, delay_draw=False, offset=1, win=1,
                _larch=None):
        if kweight is None:
            kweight = dataset.transform.kweight
        
        if isinstance(kweight, (list, tuple, ndarray)): kweight=kweight[0]
        data_chik  = dataset.data.chi * dataset.data.k**kweight
        model_chik = dataset.model.chi * dataset.model.k**kweight
        
        
        
            
        # k-weighted chi(k) in first plot window
        plt.figure()   
        plt.xlim([0,14])
        plt.ylim([-11, 20])
        plt.plot(dataset.data.k, (data_chik)+offset, color='black', label='data')     
        plt.plot(dataset.model.k, model_chik+offset, color=colors[0], label='fit')
        
        
               
        plt.xlabel(plab.k)
        plt.ylabel(plab.chikw.format(kweight))
        plt.legend(bbox_to_anchor=(1.2, 1.05))
        plt.savefig('k-space.svg')
        # plotting the real part and the magnitude of the fit
     
        #######################################################################################################            
        file_name = 'plot_data.txt'
        with open (file_name, 'w') as file:
            np.savetxt('plot_data', [dataset.data.r, dataset.data.chir_mag], header='expt data in r-space (real)', 
                   footer='second array is chir_mag')
        file_name = 'plot_model.txt'
        with open (file_name, 'w') as file:
            np.savetxt('plot_model', [dataset.model.r, dataset.model.chir_mag], header='model in r-space (real)', 
                   footer='second array is chir_mag')

        if show_imag:
            plt.plot(dataset.data.r, dataset.data.chir_im+offset, linestyle='dotted', color='black', label='Im|data|')
            plt.plot(dataset.model.r, dataset.model.chir_im+offset, linestyle='dashed', color='red', label='Im|fit|')
        
        plt.figure()
       
        
        plt.plot(dataset.data.r,  dataset.data.chir_mag+offset,
                 color='black',label='|data|')
        plt.plot(dataset.model.r, dataset.model.chir_mag+offset,
                  color=colors[0], label='|fit|', linewidth=2.0)
        plt.xlabel(plab.r)
        plt.ylabel(plab.chir.format(4))
        plt.legend()
        plt.xlim([0,5])
        plt.ylim([-18, 18])
        plt.plot(dataset.data.r, dataset.data.chir_im+offset, color='black', label='Im|data|')
        plt.plot(dataset.model.r, dataset.model.chir_im+offset,  color=colors[2], label='Im|fit|')
        plt.legend()
        plt.savefig('R-space.svg')
        
        
        
        
    plot_chifit(dataset, rmax = 3.5, show_mag=True, show_real=True)
    
    os.chdir(wd)
    print ('Report, R-space.svg and k-space.svg files generated')
    return 

     
    
    
    
    
