##All functions used in this work

##relevant libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import glob
import csv
from scipy.optimize import curve_fit
from scipy.optimize import minimize


##step 1: reading in information from .csv files, as a function of distance in micron or time in minutes
########################################################################################################
def get_info_spatial(file):
    data = {}  
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  
    
        for row in reader:
            distance = row[0]  
            values = row[1:]

            data[float(distance)] = [float(value) for value in values]
            
    original_keys = sorted(data.keys())

    depth = 4.  # 4 micron from the skeleton surface
    num_bins = len(original_keys) 
    bin_dict = {key: np.linspace(0, depth, num_bins)[int(key) - 1] for key in original_keys}   
    new_data = {bin_dict[key]: [value for value in values if value > 0.] for key, 
                  values in data.items()}
    
    means = {}  
    std_devs = {}
    for distance, values in new_data.items():
        mean_value = np.mean(values)
        means[distance] = mean_value
        
        std_deviation = np.std(values)
        std_devs[distance] = std_deviation
    return new_data, means, std_devs

def get_info_temporal(file):
    data = {}  
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  
    
        for row in reader:
            distance = row[0]  
            values = row[1:]

            data[float(distance)] = [float(value) for value in values]
            
    original_keys = sorted(data.keys())

    depth = 4.  # 4 micron from the skeleton surface
    num_bins = len(original_keys) 
    bin_dict = {key: np.linspace(0, depth, num_bins)[int(key) - 1] for key in original_keys}
    converted_bin_dict = {key: (value / growth_rate) * time_conversion for key, value in bin_dict.items()} 
    new_data = {converted_bin_dict[key]: [value for value in values if value > 0.] for key, 
                  values in data.items()}
    
    means = {}  
    std_devs = {}
    for distance, values in new_data.items():
        mean_value = np.mean(values)
        means[distance] = mean_value
        
        std_deviation = np.std(values)
        std_devs[distance] = std_deviation
    return new_data, means, std_devs

##step 2 (optional): combine information to produce meaned distributions as function of distance or time
########################################################################################################
def sum_dictionaries(num_files,one = None,two=None,three=None,
                      four=None,five=None):

    total_pixels = 0
    dicts = [one,two,three,four,five]
    sum_dict = defaultdict(list)

    for d in dicts:
        for key, values in d.items():
            total_pixels += len(values)
            if key in sum_dict:
                sum_dict[key] = [sum(x) for x in zip(sum_dict[key], values)]
            else:
                sum_dict[key] = values[:]
    
    sum_dict = dict(sum_dict)
    
    for key in sum_dict:
        sum_dict[key] = [val / num_files for val in sum_dict[key]]
        
    means = {}  
    std_devs = {}
    for distance, values in sum_dict.items():
        m_value = np.sum(values)/len(values)
        means[distance] = m_value
        std_deviation = np.std(values)
        std_devs[distance] = std_deviation
    return  sum_dict, means, std_devs

##step 3: fitting functional forms to data
#########################################################################################################
def exponential_decay_function(x, a, b, c):
    return a * np.exp(-b * x) + c

def logarithmic_decay_function(x, a, b, c):
    return a * np.log(b * x ) + c

def r2_score(y_true, y_pred):
    ss_total = np.sum((y_true - np.mean(y_true)) ** 2)
    ss_residual = np.sum((y_true - y_pred) ** 2)
    return 1 - (ss_residual / ss_total)

def best_fit(exp_r2, log_r2):
    exp_dif = np.abs(1 - exp_r2)
    log_dif = np.abs(1 - log_r2)
    best_model = None
    if exp_dif < log_dif:
        best_model = 'exp'
    else:
        best_model = 'log'
    return (f'exp R2: {np.round(exp_r2,4)}, ln R2: {np.round(log_r2,4)} --> best model: {best_model}')

def make_exp_fit_pre(x_data, y_data, sample_name):
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
                    print(f"An error occurred: {e}")  
                    pass  
    
    print(f'Best R^2 for exponential fit: {best_r2}')
    print(f'Best parameters for exponential fit: {best_params}')
    return best_fit, best_params, best_r2

def make_log_fit_pre(x_data, y_data, sample_name):
    x_data[0]=.01
    a_range = np.linspace(-1., -1e-5, 10)
    b_range = np.linspace(0.001, 4., 20)
    c_range = np.linspace(-1., 0., 20)
    best_r2 = -np.inf 
    best_params = None
    best_fit = None

    for a in a_range:
        for b in b_range:
            for c in c_range:
                try:
                    
                    params, covariance = curve_fit(logarithmic_decay_function, x_data, y_data, 
                                                   p0=[a, b, c], bounds=([-1, 1e-5, -1], [-1e-5, np.inf, 0]))
                    y_fit_log = logarithmic_decay_function(x_data, *params)
                    r2_log = r2_score(y_data, y_fit_log)

                    
                    if r2_log > best_r2:
                        best_r2 = r2_log
                        best_params = params
                        best_fit = y_fit_log
                        best_loss_function = 'log'

                except Exception as e:
                    print(f"An error occurred: {e}")  
                    pass  
    print(f'Best R^2 for logarithmic fit: {best_r2}')
    print(f'Best parameters for logarithmic fit: {best_params}')
    return best_fit, best_params, best_r2

def make_exp_fit_arag1(x_data, y_data):
    def loss_function(params, x, y):
        a, b, c = params
        y_pred = exponential_decay_function(x, a, b, c)
        mse = np.mean((y - y_pred) ** 2)
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
                    r2_exp = r2_score(y_data, y_fit)
                    if r2_exp > best_r2_exp:
                        best_r2_exp = r2_exp
                        best_params_exp = params
                except Exception as e:
                    print(e)
                    pass  
    
    print(f'Best R^2 for exponential fit: {best_r2_exp}')
    print(f'Best parameters for exponential fit: {best_params_exp}')
    y_data = exponential_decay_function(x_data,*best_params_exp)
    return y_data, best_params_exp, best_r2_exp

def make_exp_fit_arag2(x_data, y_data):
    def loss_function(params, x, y):
        a, b, c = params
        y_pred = a * np.exp(-b * x) + c
        # Assign higher weights to points on the left side of the curve
        weights = np.exp(-x)
        mse = np.mean(weights * (y - y_pred) ** 2)
        return mse

    a_range = np.linspace(0.001, 1, 10)
    b_range = np.linspace(0.001, 0.2, 10)
    c_range = np.linspace(0.001, 0.2, 10)
    
    best_r2_exp = -np.inf
    best_params_exp = None
    
    for a in a_range:
        for b in b_range:
            
            for c in c_range:
                try:

                    bounds = [(None, None), 
                              (None, None), 
                              (None, None)]  

                    res = minimize(loss_function, x0=[a,b,c], args=(x_data, y_data), 
                    method='Nelder-Mead')
                    params = res.x
                    y_fit = exponential_decay_function(x_data, *params)
                    r2_exp = r2_score(y_data, y_fit)
                    if r2_exp > best_r2_exp:
                        best_r2_exp = r2_exp
                        best_params_exp = params
                except Exception as e:
                    print(e)
                    pass  
    
    print(f'Best R^2 for exponential fit: {best_r2_exp}')
    print(f'Best parameters for exponential fit: {best_params_exp}')
    y_data = exponential_decay_function(x_data,*best_params_exp)
    return y_data, best_params_exp, best_r2_exp

##step 4: calculating 1/e length
########################################################################################################
def length_1e(params):
    return 1/params[1]