'''
Import necessary packages
'''
import numpy as np
import scipy
import scipy.integrate
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
from matplotlib import rc
import pickle
from tqdm import tqdm

# from ip3_ca_ode import *
from ip3_ca_ode_cfg import *
import cfg
import os

save = 'figures/paper_plots/'

#set figure font sizes for readability
font = {'size' : 30,
       'family': 'serif',
       'sans-serif': ['Helvetica']}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

        
   
def run_spritz(period=10, conc=0.6, dur=1, num=10, noise=False, osc_type=1,
              fix_c_er=False):
    '''
    Run a spritz experiment - use a different solver step size for spritz versus
    rest time to significantly speed up the experiment
    period: second between spritzes
    conc: concentration
    dur: duration of spritz
    num: number of sprtizes1
    osc_type:
        1: exponential ramping
        2: flat pulses
    !Note that specific configurations of input_smoothing and oscillation on duration
    need to be used to get exponential_oscillation2 to work (check help details)
    '''
    
    max_step1 = dur * 0.1
    
        # print(max_step)
    final_t_f = (period + dur) * (num + 1)
    
    cfg.input_max = conc
    cfg.input_start = 0
    cfg.input_duration = 0
    cfg.oscillation_on_duration = dur
    cfg.oscillation_off_duration = period
    cfg.num_oscillations = num
    cfg.input_smoothing = dur/2
    
    max_step2 = period * 0.02
    
    cur_t = 0
    if osc_type == 1:
        input_type = 'exponential_oscillation2'
    elif osc_type == 2:
        input_type = 'oscillation'
        
    for i in range(num):
        #run spritz part
        if i == 0:
            run_experiment(input_type, t_f=dur, max_step=max_step1)
        else:
            continue_experiment(input_type, t_cont=dur, max_step=max_step1)
            
        #run rest part until next spike
        continue_experiment(input_type, t_cont=period, max_step=max_step2)
    
    


# def run_spritz(period=10, conc=0.6, dur=1, num=10, noise=False, osc_type=1,
#               fix_c_er=False):
#     '''
#     DEPRECATED
#     Run a spritz experiment
#     period: second between spritzes
#     conc: concentration
#     dur: duration of spritz
#     num: number of sprtizes1
#     osc_type:
#         1: exponential ramping
#         2: flat pulses
#     !Note that specific configurations of input_smoothing and oscillation on duration
#     need to be used to get exponential_oscillation2 to work (check help details)
#     '''
    
#     max_step = dur * 0.1
#     # print(max_step)
#     t_f = (period + dur) * (num + 1)
    
#     cfg.input_max = conc
#     cfg.input_start = 0
#     cfg.input_duration = 0
#     cfg.oscillation_on_duration = dur
#     cfg.oscillation_off_duration = period
#     cfg.num_oscillations = num
#     cfg.input_smoothing = dur/2
    
#     if osc_type == 1:
#         run_experiment('exponential_oscillation2', t_f=t_f, max_step=max_step, noise=noise,
#                       fix_c_er=fix_c_er)
#     if osc_type == 2:
#         run_experiment('oscillation', t_f=t_f, max_step=max_step, noise=noise,
#                       fix_c_er=fix_c_er)
    
    
    
def get_condition_bounds():
    """Perform a single long bath 1000s with concentration of 1 to get min and max
    ranges for each variable for grid search

    Returns:
        var_ranges: 8 tuples, one min and max value for each var
    """
    input_max = cfg.input_max
    input_duration = cfg.input_duration
    cfg.input_max = 1
    cfg.input_duration = 1000
    run_experiment('pulse', max_step=1)
    
    variables = ['c', 'c_tot', 'h', 'p', 'Gstar', 'Gd1', 'Gd2', 'lamb']
    var_ranges = []
    for v in variables:
        val = getattr(cfg, v)
        var_ranges.append((np.min(val), np.max(val)))
    
    return var_ranges
        
        
    
def test_single_spritz(init_conds, conc=1.3, t_record=7.5):
    """Record response to a single spritz of 1 second duration

    Args:
        init_conds (list): 8 values making up initial condition
        conc (float, optional): concentration of spritz. Defaults to 1.3.
        t_record (float, optional): when to record the variable conditions. Defaults to 7.5.

    Returns:
        _type_: _description_
    """
    temp_conds = cfg.all_init.copy()
    cfg.all_init = init_conds
    run_spritz(period=5, conc=conc, dur=1, num=1)
    cfg.all_init = temp_conds
    
    t_idx = np.argmax(cfg.t > 7.5)
    next_var_vals = [cfg.c[t_idx], cfg.c_tot[t_idx], cfg.h[t_idx], cfg.p[t_idx], 
                     cfg.Gstar[t_idx], cfg.Gd1[t_idx], cfg.Gd2[t_idx], cfg.lamb[t_idx]]
    
    return np.max(cfg.c), next_var_vals



def get_random_init(var_ranges):
    var_ranges = np.array(var_ranges)
    base = var_ranges[:, 0]
    rnges = var_ranges[:, 1] - var_ranges[:, 0]
    init = np.random.random(8) * rnges + base
    return init
