#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:11:09 2022

@author: rachita
"""
import os, sys
import numpy as np
from exafs_frac import read_structures, read_experimental_data
from exafs_frac import *
import scipy
from scipy.optimize import minimize


wd2 = os.getcwd()
#def_array = np.arange(start=0.1, stop=1, step=0.1)
def_array = np.array([1])
for my_amp in def_array:
    path = 'cc-all_paths-100-O2' #path = 'folder{}'.format('%.1f' % my_amp)
    if not os.path.exists(path):
        os.mkdir(path)
        os.chdir(path)
        read_structures ()
        data = exafs_fitting (my_amp)
        #data = exafs_fitting (my_amp)
    os.chdir(wd2)
    
