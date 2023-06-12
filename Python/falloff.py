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
import proplot as pplt

# from ip3_ca_ode import *
from ip3_ca_ode_cfg import *
import cfg
import os

save = 'figures/paper_plots/falloff/'

#set figure font sizes for readability
# font = {'size' : 30,
#        'family': 'serif',
#        'sans-serif': ['Helvetica']}
# matplotlib.rc('font', **font)
# matplotlib.rc('text', usetex=True)
color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

'''
Functions used in generating calcium falloff plots
'''

def report_max_ratios(times=[0, 100, 160, 300], var='c'):
    '''
    Report the max values and ratio of maxes for a variable in the time frames given
    Times should always be a list of 4 integers giving two time intervals to measure maxes from
    Returns:
        ratio of second max to first max (e.g. if second max is 50% of first, return 0.5)
        list of max values reported
        list of index location of max values (for plotting)
    '''
    indices = [np.argmax(cfg.t > t) for t in times]
    if indices[-1] == 0:
        indices[-1] = -1
    y = getattr(cfg, var)
    
    maxes = []
    pos = []
    for i in range(2):
        low = indices[i*2]
        high = indices[i*2 + 1]
        
        ymax = np.max(y[low:high])
        maxpos = np.argmax(y[low:high])
        
        pos.append(maxpos + low)
        maxes.append(ymax)
        
    ratio = maxes[1] / maxes[0]
    return ratio, maxes, pos

def slope_at_index(idx, var='c'):
    y = getattr(cfg, var)
    t2 = cfg.t[idx + 1]
    t1 = cfg.t[idx - 1]
    y2 = y[idx + 1]
    y1 = y[idx - 1]
    return (y2 - y1) / (t2 - t1)


def conditional_double_bath(glut=0.4, first=240, rest=310, mod=0, noise=False,
                           force_run=False):
    '''
    get experimental double bath c ratio from globals()['df'] if it is loaded
    otherwise run and save
    '''
    if 'df' in globals()['df']:
        df = globals()['df']
    else:
        df = pd.read_pickle('data/experiments/double_bath/dataframe')
        globals()['df'] = df
    
    df = df[df['glut'] == glut]
    df = df[df['first'] == first]
    df = df[df['rest'] == rest]
    df = df[df['mod'] == mod]
    df = df[df['noise'] == noise]
    
    kd1 = 0.02
    kd2 = 0.6
    
    #if no such experiment found, run it and save the results
    if len(df) == 0 or force_run:
        set_init()
        ratio_times = [10, 10+first, 10+first+rest, 70+first+rest]
        cfg.custom_input_times = [0, 10, 10+first, 10+first+rest, 70+first+rest]
        cfg.custom_input_vals = [0, glut, 0, glut, 0]
        cfg.kd1 = kd1 * mod
        cfg.kd2 = kd2 * mod
        sol = run_experiment('custom', t_f=80+first+rest, max_step=0.1)
        ratio, _, _ = report_max_ratios(times=ratio_times)
        
        df_row = pd.DataFrame([{
            'glut': glut,
            'first': first,
            'rest': rest,
            'mod': mod,
            'noise': noise,
            'ratio': ratio
        }], columns=['glut','first','rest','mod','noise','ratio'])
        df = pd.read_pickle('data/experiments/double_bath/dataframe')
        df = df.append(df_row, ignore_index=True)
        df.to_pickle('data/experiments/double_bath/dataframe')
        globals()['df'] = df
        return ratio
    else:
        ratio = df.iloc[0]['ratio']
        return ratio
    
    
def capture_c(t):
    '''
    get c at values of t
    '''
    idx = np.argmax(cfg.t >= t)
    c = cfg.c[idx]
    return c



'''Oscillation functions'''

def get_peaks():
    peaks = scipy.signal.find_peaks(cfg.c)[0]
    return_low_lim = 0.2 #how far does calcium have to return before we consider
                         #ourselves to be in standard oscillation range
    first_peak = peaks[0]
    # return_low = np.argmax(cfg.c[first_peak:] < return_low_lim)
    # if return_low == 0:
    #     osc_start = -1
    # else:
    #     osc_start = first_peak + return_low
    # peaks = np.append([first_peak], scipy.signal.find_peaks(cfg.c[osc_start:])[0] + osc_start)
    c_peaks = cfg.c[peaks]
    t_peaks = cfg.t[peaks]
    
    return t_peaks, c_peaks


import matplotlib.image as mpimg

falloff1 = np.array([    
    [0, 65, 66, 67, 68, 72, 75, 80, 87, 95, 100, 105, 107, 110, 120, 130, 140, 150, 160, 165, 170, 190, 230, 250],
    [0.02, 0.02, 0.3, 0.6, 0.7, 0.94, 0.95, 0.7, 0.55, 0.4, 0.35, 0.37, 0.35, 0.3, 0.33, 0.3, 0.25, 0.2, 0.12, 0.1, 0.09, 0.08, 0.06, 0.055]
])
falloff2 = np.array([
    [0, 65, 66, 67, 68, 69, 70, 71, 72, 76, 78, 85, 88, 92, 97, 104, 107, 130, 135, 140, 170, 200, 250],
    [0.02, 0.02, 0.3, 0.6, 0.7, 0.72, 0.7, 0.65, 0.6, 0.4, 0.35, 0.3, 0.31, 0.33, 0.22, 0.25, 0.2, 0.19, 0.15, 0.1, 0.12, 0.1, 0.08]
])
falloff3 = np.array([
    [0, 65, 66, 67, 68, 69, 72, 73, 77, 80, 83, 88, 94, 105, 140, 200, 250],
    [0.02, 0.02, 0.25, 0.35, 0.52, 0.45, 0.3, 0.21, 0.22, 0.15, 0.16, 0.12, 0.1, 0.09, 0.072, 0.055, 0.05]
])