    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:41:07 2019

@author: fboucher
This function uploads output data from model and assings it to variables
Variables are named after output data (parameter, time)


"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

class obj_:
    pass

    def __repr__(self):
        return '<Datasets>'
a = obj_()
def load_data():
    """
    this function loads files from output and assings them to object a
    filename: shear1990_01_01_00_20_k0009 corresponds to shearYY_MM_DD_HH_MM
    object 'a' contains every file from output when called
    """    
    

    
    for filename in os.listdir(os.getcwd()):
        if not filename in ['read_data.py', 'info.00']:
            #file names are kept unchaged expeted for '.00'
            setattr(a, filename[:-3], np.loadtxt(filename))
#    print(a.filename)


def plot_dataLog(filename):
    
    data = getattr(a, filename[:-3])
    
    plt.imshow(data[1:-1,3:-3], origin='upper', norm=mpl.colors.LogNorm())
    plt.colorbar()
    """
    This function takes as an entry filename of output and plots it. 
    It works fine with shear/divergence/A/h/sigI/sigII some issues with sig1norm/sig2norm
    """

def plot_data(filename):
    
    data = getattr(a, filename[:-3])
    
    plt.imshow(data[1:-1,3:-3], origin='upper')
    plt.colorbar()    