##python libraries and functions used to perform fits

##libraries
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import colors as clrs
import matplotlib as mpl
import pandas as pd
%matplotlib inline

import glob
import os.path
import random
import csv
from collections import defaultdict
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

##functions to fit
def exponential_decay_function(x, a, b, c):
    return a * np.exp(b * x) + c

def logarithmic_decay_function(x, a, b, c):
    return a * np.log(b * x ) + c

##other functions
def r2_score(y_true, y_pred):
    ss_total = np.sum((y_true - np.mean(y_true)) ** 2)
    ss_residual = np.sum((y_true - y_pred) ** 2)
    return 1 - (ss_residual / ss_total)
    
##scipy optimizer fitting function structure: can modify to test functions of interest
def make_exp_fit_precursor(x_data, y_data, sample_name):
    x_data[0] = 0.01  
    a_range = np.linspace(0.001, 1, 10)
    b_range = np.linspace(-1., 4., 20)
    c_range = np.linspace(0.001, 1., 20)
    best_r2 = -np.inf  
    best_params = None
    best_fit = None

    for a in a_range:
        for b in b_range:
            for c in c_range:
                try:
                    
                    params, covariance = curve_fit(exponential_decay_function, x_data, y_data, 
                                                   p0=[a, b, c], bounds=([0, -np.inf, 0], [1, np.inf, 1]))
                    y_fit_exp = exponential_decay_function(x_data, *params)
                    r2_exp = r2_score(y_data, y_fit_exp)

                    
                    if r2_exp > best_r2:
                        best_r2 = r2_exp
                        best_params = params
                        best_fit = y_fit_exp
                        best_loss_function = 'exp'

                except Exception as e:
                    print(f"err: {e}")  
                    pass  
    
    return best_fit, best_params, best_r2

##scipy minimizer fitting function structure: can modify to test functions of interest
##def make_exp_fit_arag1(x_data, y_data):
    ##standard loss function
    def loss_function(params, x, y):
        a, b, c = params
        y_pred = exponential_decay_function(x, a, b, c)
        mse = np.mean((y - y_pred) ** 2)
        return mse
    def exponential_decay_function_arag(x, a, b, c):
        return a * np.exp(b * x) + c

    a_range = np.linspace(0.001, 1, 10)
    b_range = np.linspace(0.001, 0.2, 10)
    c_range = np.linspace(0.001, 0.2, 10)
    
    best_r2_exp = -np.inf
    best_params_exp = None
    
    for a in a_range:
        #print(a)
        for b in b_range:
            
            for c in c_range:
                try:

                    bounds = [(None, None), (None, None), (None, None)]  # No bounds for a and c

                    res = minimize(loss_function, x0=[a,b,c], args=(x_data, y_data), 
                    method='Nelder-Mead')
                    params = res.x
                    y_fit = exponential_decay_function(x_data, *params)
                    r2_exp = r2_score_manual(y_data, y_fit)
                    if r2_exp > best_r2_exp:
                        best_r2_exp = r2_exp
                        best_params_exp = params
                except Exception as e:
                    print(e)
                    pass  
    return y_data, best_params_exp, best_r2_exp
    
def make_exp_fit_arag2(x_data, y_data):
    def loss_function(params, x, y):
        ##non standard loss function, assign higher weights to points on the left side of the curve
        a, b, c = params
        y_pred = a * np.exp(-b * x) + c

        weights = np.exp(-x)
        mse = np.mean(weights * (y - y_pred) ** 2)
        return mse

    a_range = np.linspace(0.001, 1, 10)
    b_range = np.linspace(0.001, 0.2, 10)
    c_range = np.linspace(0.001, 0.2, 10)
    
    best_r2_exp = -np.inf
    best_params_exp = None
    
    for a in a_range:
        #print(a)
        for b in b_range:
            
            for c in c_range:
                try:

                    bounds = [(None, None), (None, None), (None, None)]

                    res = minimize(loss_function, x0=[a,b,c], args=(x_data, y_data), 
                    method='Nelder-Mead')
                    params = res.x
                    y_fit = exponential_decay_function(x_data, *params)
                    r2_exp = r2_score_manual(y_data, y_fit)
                    if r2_exp > best_r2_exp:
                        best_r2_exp = r2_exp
                        best_params_exp = params
                except Exception as e:
                    print(e)
                    pass
    
    return y_data, best_params_exp, best_r2_exp
