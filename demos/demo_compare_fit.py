##Determine best fit using R2 values

##import packages
import numpy as np
import matplotlib.pyplot as plt
import glob
import csv

from scipy.optimize import curve_fit
from scipy.optimize import minimize
plt.rcParams.update({
    'font.size': 26,  
    'text.color': 'black',  
    'xtick.labelsize': 26, 
    'ytick.labelsize': 26,  
    'legend.fontsize': 26, 
    'axes.labelsize': 26,  
})

time_conversion = 24 * 60 # 24 hours * 60 minutes
growth_rate = 210 #210 micron per day

##some commands to suppress warnings for initial parameters 
##that lead to large/small values in exponential or logarithmic functions
import warnings
warnings.filterwarnings("ignore", 
                        category=RuntimeWarning, 
                        message="overflow encountered in scalar divide")
warnings.filterwarnings("ignore", 
                        category=RuntimeWarning, 
                        message="overflow encountered in exp")

##define functions
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

def exponential_decay_function(x, a, b, c):
    return a * np.exp(-b * x) + c

def logarithmic_decay_function(x, a, b, c):
    return a * np.log(b * x ) + c

def r2_score(y_true, y_pred):
    ss_total = np.sum((y_true - np.mean(y_true)) ** 2)
    ss_residual = np.sum((y_true - y_pred) ** 2)
    return 1 - (ss_residual / ss_total)

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

def best_fit(exp_r2, log_r2):
    exp_dif = np.abs(1 - exp_r2)
    log_dif = np.abs(1 - log_r2)
    best_model = None
    if exp_dif < log_dif:
        best_model = 'exp'

    else:
        best_model = 'log'

    return (f'exp R2: {np.round(exp_r2,4)}, ln R2: {np.round(log_r2,4)} --> best model: {best_model}')





###call functions
def main():
    S49_0_data_x, S49_0_med_x, S49_0_std_x = get_info_spatial('data/S49_0.csv')
    S49_1_data_x, S49_1_med_x, S49_1_std_x = get_info_spatial('data/S49_1.csv')
    S49_2_data_x, S49_2_med_x, S49_2_std_x = get_info_spatial('data/S49_2.csv')
    S49_3_data_x, S49_3_med_x, S49_3_std_x = get_info_spatial('data/S49_3.csv')
    
    S49_0_x = np.array(list(S49_0_med_x.keys()))
    S49_0_y = np.array(list(S49_0_med_x.values()))
    S49_1_x = np.array(list(S49_1_med_x.keys()))
    S49_1_y = np.array(list(S49_1_med_x.values()))
    S49_2_x = np.array(list(S49_2_med_x.keys()))
    S49_2_y = np.array(list(S49_2_med_x.values()))
    S49_3_x = np.array(list(S49_3_med_x.keys()))
    S49_3_y = np.array(list(S49_3_med_x.values()))

    S49_0_fit,S49_0_params,S49_0_r2 = make_exp_fit_pre(S49_0_x,S49_0_y,'S49_0')
    S49_1_fit,S49_1_params,S49_1_r2 = make_exp_fit_pre(S49_1_x,S49_1_y,'S49_1')
    S49_2_fit,S49_2_params,S49_2_r2 = make_exp_fit_pre(S49_2_x,S49_2_y,'S49_2')
    S49_3_fit,S49_3_params,S49_3_r2 = make_exp_fit_pre(S49_3_x,S49_3_y,'S49_3')

    S49_0_fitl,S49_0_paramsl,S49_0_r2l = make_log_fit_pre(S49_0_x,S49_0_y,'S49_0')
    S49_1_fitl,S49_1_paramsl,S49_1_r2l = make_log_fit_pre(S49_1_x,S49_1_y,'S49_1')
    S49_2_fitl,S49_2_paramsl,S49_2_r2l = make_log_fit_pre(S49_2_x,S49_2_y,'S49_2')
    S49_3_fitl,S49_3_paramsl,S49_3_r2l = make_log_fit_pre(S49_3_x,S49_3_y,'S49_3')

    print(best_fit(S49_0_r2,S49_0_r2l))
    print(best_fit(S49_1_r2,S49_1_r2l))
    print(best_fit(S49_2_r2,S49_2_r2l))
    print(best_fit(S49_3_r2,S49_3_r2l))

if __name__ == "__main__":
    main()