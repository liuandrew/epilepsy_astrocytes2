'''
Import necessary packages
'''
import numpy as np
import pandas as pd
import scipy
import scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import pickle
import proplot as pplt
import random
import math

import cfg

#set figure font sizes for readability
# font = {'size' : 30,
#        'family': 'serif',
#        'sans-serif': ['Helvetica']}
# matplotlib.rc('font', **font)
# matplotlib.rc('text', usetex=True)
color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

colors = {
    'stable_ss': 'green',   #stable steady state
    'unstable_ss': 'black', #unstable steady state
    'stable_po': 'blue',    #stable periodic orbit
    'unstable_po': 'red'    #unstable periodic orbit
}

linestyles = {
    'stable_ss': '-',   #stable steady state
    'unstable_ss': '-', #unstable steady state
    'stable_po': '--',    #stable periodic orbit
    'unstable_po': '--'    #unstable periodic orbit    
}

pplt.rc.update({'font.size': 10})
c_label = r'[Ca$^{2+}$]$_{cyt}$'
c_er_label = r'[Ca$^{2+}$]$_{ER}$'


#-------
#IP3 double exponential values from Marsa's paper
#-------
compute_r_decay = lambda A, d_decay: -1 / d_decay * np.log(0.005 / A)
#This dict contains some preset IP3 curve parameters so we can easily set
#which curve we want later
ip3_curves = {
    #Single peak
    'SP': {
        'A': 0.2,
        'd_rise': 21,
        'r_rise': 0.002,
        'd_decay': 97,
        'r_decay': compute_r_decay(0.2, 97)
    },
    #Multi peak
    'MP': {
        'A': 0.2,
        'd_rise': 41,
        'r_rise': 0.15,
        'd_decay': 179,
        'r_decay': compute_r_decay(0.2, 97)
    },
    #Plateau
    'PL': {
        'A': 0.375,
        'd_rise': 36,
        'r_rise': 0.002,
        'd_decay': 120,
        'r_decay': compute_r_decay(0.375, 120)
    }
}



'''
#############################

ODE Functions

#############################
'''




'''
--------------------------
ODE equations
--------------------------
These functions are called by the numerical solver. Given a vector of values
x, they return the functions for xdot
'''

def all_ode_equations(t, x, glutamate_input_type='pulse', old_gpcr=False, Gstar_input_type=False,
                        fix_c_er=False, fix_c_er_func=None):
    '''
    This is the ODE for the full Glutamate/GPCR/IP3/Calcium system
    We will in effect call ip3_ca_ode_equations (with manual_ip3=False)
    as well as gpcr_ode_equations
    fix_c_er: if passed in a number value, c_er will be fixed to that value
    
    Please make sure manual_ip3 == False
    
    IP3/Ca:
    c = x[0]
    c_tot = x[1]
    h = x[2]
    p = x[3]
    
    GPCR:
    Gstar = x[4]
    Gd1 = x[5]
    Gd2 = x[6]
    lambda = x[7] #this is the downstream PKA type signal for Gd2
    '''
    xdot = np.zeros(len(x))
    
    #manual Gstar
    if(Gstar_input_type is not False): 
        Gstar = get_input_value(Gstar_input_type, t)
        x[4] = Gstar

    #fix c_er to certain value
    if(fix_c_er is not False):
        c_tot = (fix_c_er / cfg.gamma) + x[0]
        x[1] = c_tot
    elif(fix_c_er_func is not None):
        c_tot = (fix_c_er_func(t) / cfg.gamma) + x[0]
        x[1] = c_tot

    #Note here we are passing Gstar to ip3_ode_equations
    xdot[:4] = ip3_ca_ode_equations(t, x[:4], ip3_input_type=None,
                                    Gstar=x[4])
    if(old_gpcr):
        xdot[4:] = old_gpcr_ode_equations(t, x[4:], glutamate_input_type)
    else:
        xdot[4:] = gpcr_ode_equations(t, x[4:], glutamate_input_type)
    
    return xdot




#This old set of GPCR ode equations were those that Alla's team originally had
#but now we are using the one with a downstream signal molecule for Gd2
def old_gpcr_ode_equations(t, x, glutamate_input_type='pulse'):
    '''
    ODE equations following the GPCR equations given earlier
    x[0] = Gstar
    x[1] = Gd1
    x[2] = Gd2
    This returns x_dot, which is an array of the right hand sides of 
    each of the ODE equations in the same order as above
    '''
    Gstar = x[0]
    Gd1 = x[1]
    Gd2 = x[2]
    G = 1 - Gstar - Gd1 - Gd2
    
    #query our glutamate function for what level glutamate is input
    glut = get_input_value(glutamate_input_type, t)
    
    x_dot = np.zeros(4)
    
    x_dot[0] = cfg.kp*glut*G - cfg.km*Gstar - cfg.kd1*Gstar
    x_dot[1] = cfg.kd1*Gstar - cfg.kr1*Gd1
    x_dot[2] = cfg.kd2*Gstar*G - cfg.kr2*Gd2
    
    return x_dot





def gpcr_ode_equations(t, x, glutamate_input_type='pulse'):
    '''
    ODE equations following the GPCR equations given earlier
    x[0] = Gstar
    x[1] = Gd1
    x[2] = Gd2
    x[3] = lambda
    This returns x_dot, which is an array of the right hand sides of 
    each of the ODE equations in the same order as above
    
    input_type: given by ODE args, can be 'pulse' or 'step'
    gstar_to_d2: whether to have the G* -> Gd2 arrow on or off (True/False)
    sat_kinetics: whether lambda should act by saturating kinetics
    '''    
    Gstar = x[0]
    Gd1 = x[1]
    Gd2 = x[2]
    lamb = x[3]
    G = 1 - Gstar - Gd1 - Gd2
    
    #query our glutamate function for what level glutamate is input
    glut = get_input_value(glutamate_input_type, t)
    
    x_dot = np.zeros(4)
        
    x_dot[0] = cfg.kp*glut*G - cfg.km*Gstar - cfg.kd1*Gstar - cfg.kd2*Gstar*lamb
    x_dot[1] = cfg.kd1*Gstar - cfg.kr1*Gd1
    x_dot[2] = cfg.kd2*(Gstar+G)*lamb - cfg.kr2*Gd2
    x_dot[3] = cfg.klp*Gstar - cfg.klm*lamb
        
    return x_dot



def ip3_ca_ode_equations(t, x, ip3_input_type=None, Gstar=None):
    '''
    ODE equations detailing the calcium transient response to a
    given IP3 level
    Pass Gstar in if we are doing the full ODE system
    '''
    c = x[0] #Cytosolic Ca2+ concentration
    c_tot = x[1] #Total free Ca2+ concentration in cell
    h = x[2] #Deactivation parameter for IP3R
    p = 0 #IP3 concentration
    
    
    #We can choose to either use ODE for IP3 or manually input it
    #ODE
    if(ip3_input_type == None):
        p = x[3] #IP3 concentration
    
    #Explicit input
    else:
        p = get_input_value(ip3_input_type, t)
    
    #Compute ER Ca2+ concentration based on total Ca2+ and Cyt Ca2+
    c_er = (c_tot - c) * cfg.gamma
    
    #First compute some derivative values that will be needed
    #for h and IP3R dynamics
    m_inf = p / (p + cfg.d1)
    n_inf = c / (c + cfg.d5)
    Q2 = cfg.d2 * (p + cfg.d1) / (p + cfg.d3)
    h_inf = Q2 / (Q2 + c)
    tau_h = 1 / (cfg.a2 * (Q2 + c))
    
    #Compute the fluxes through each channel
    J_ip3r = cfg.v_ip3r * (m_inf**3) * (n_inf**3) * (h**3) * (c_er - c)
    J_serca = cfg.v_serca * (c**1.75) / (c**1.75 + cfg.k_serca**1.75)
    J_pmca = cfg.v_pmca * (c**2) / (c**2 + cfg.k_pmca**2)
    J_soc = cfg.v_soc * (cfg.k_soc**4) / (cfg.k_soc**4 + c_er**4)
#     J_soc = cfg.v_soc * (cfg.k_soc**2) / (cfg.k_soc**2 + c_er**2)
    
    #leak fluxes
    J_er_leak = cfg.v_er_leak * (c_er - c) #ER <-> Cyt leak
    J_ecs_add = cfg.v_in - cfg.k_out * c #Cyt <-> extracellular leak
    
    x_dot = np.zeros(len(x))
    
    x_dot[0] = J_ip3r + J_er_leak - J_serca + cfg.delta*(J_ecs_add - J_pmca + J_soc)
    x_dot[1] = cfg.delta*(J_ecs_add - J_pmca + J_soc)
    x_dot[2] = (h_inf - h) / tau_h
    x_dot[3] = ip3_ode_equation(t, x, Gstar)
    
    return x_dot
    
    



def ip3_ode_equation(t, x, Gstar=None):
    '''
    ODE equations for IP3 production and degradation
    This function will be called by ip3_ca_ode_equations if
    manual_ip3 is set to False and we are simulating dynamics
    '''
    c = x[0] #Cytosolic Ca2+ concentration
    c_tot = x[1] #Total free Ca2+ concentration in cell
    h = x[2] #Deactivation parameter for IP3R
    p = x[3] #IP3 concentration 
    
    if(Gstar == None):
        #use a square wave for Gstar for now
        Gstar = pulse_input(t)
        #Gstar = 0
        
    ip3_production = cfg.v_beta*Gstar + cfg.v_delta*((cfg.k_delta)/(1 + p))*((c**2)/(c**2 + cfg.k_plcdelta**2))
    ip3_degradation = cfg.v_3k*((c**4)/(c**4 + cfg.k_d**4))*(p/(p+cfg.k_3)) + cfg.r_5p*p
#     ip3_degradation = cfg.v_3k*((c**2)/(c**2 + k_d**2))*(p/(p+cfg.k_3)) + cfg.r_5p*p
    
    return ip3_production - ip3_degradation

    
    



'''
--------------------------
Input functions
--------------------------
These functions are used to give manual input for glutamate or IP3
The exact shape of input can be modified by changing global parameters
One of the following input string types should be passed in the 
    run_experiment functions to use the corresponding input to run a simulation
Note that to change the shape of the input, certain parameters in cfg can be changed

E.g., cfg.input_duration = 500
    will make a 'pulse' type input last for 500 seconds

!Note: glutamate concentrations are passed in uM (e.g. cfg.input_max=20 means 20uM)

Input types and parameters that can be changed:
    'pulse': used for constant input value, params are
        cfg.input_start:        when the input starts
        cfg.input_duration:     how long input lasts
        cfg.input_max:          how large the input is when on
        cfg.input_min:          how large the input is when off
    'step': used for input that will step up every 150s uniformly to a max
        cfg.step_max_value:     the max value the step will reach
            e.g. if it is 1, the steps will be [0, 0.25, 0.5, 0.75, 0.1, 0.75, 0.5, 0.25, 0]
        cfg.step_time_scale:    how long each step will last
            e.g. if it is 1, each step lasts 150s, if it is 2, each step lasts 300s
    'oscillation': input will oscillate between on and off
        cfg.input_min, cfg.input_max: min/max values of oscillation
        cfg.num_oscilations:      how many oscillations
        cfg.oscillation_on_duration, cfg.oscillation_off_duration: how long each on/off 
                                                           phase will last
    'curve': this is a double exponential curve
        cfg.A, cfg.d_rise, cfg.d_decay, cfg.r_rise, cfg.r_decay: exponential curve parameters
    'exponential_oscillation': same as 'oscillation' but will rise and fall
        exponentially. Uses parameters of both 'oscillation' and 'curve'
    'custom': this is a custom input transient defined by a set of times for switching
        and values to switch to (given by lists or arrays)
        cfg.custom_input_times:     times that the input value switches
        cfg.custom_input_vals:      values that the input switches to at each time
'''

def get_input_value(input_type, t):
    '''
    Helper function - given a certain input type and time t
    return the value for the input
    '''
    if cfg.noise:
        return cfg.input_func(t)

    if(input_type == 'pulse'):
        return pulse_input(t)
    elif(input_type == 'step'):
        return step_input(t)
    elif(input_type == 'oscillation'):
        return oscillation_input(t)
    elif(input_type == 'noisy_oscillation'):
        return noisy_oscillation_input(t)
    elif(input_type == 'curve'):
        return curve_input(t)
    elif(input_type == 'exponential_oscillation'):
        return exponential_oscillation(t)
    elif(input_type == 'exponential_oscillation2'):
        return exponential_oscillation2(t)
    elif(input_type == 'custom'):
        return custom_input(t)
    elif(input_type == 'exponential_train'):
        return exponential_train(t)
    elif(input_type == 'train'):
        return train(t)
    elif(input_type == 'exponential_pulse'):
        return exponential_input(t)
    else:
        return 0




def get_input_plot(input_type, max_step=0.1, noise=False):
    '''
    This function simply returns t and y values to plot
    based on the type of input curve that we used
    E.g., we could call this to get glut_t, glut for the glutamate
    input plot

    If noise is not False, use the value passed as the scaling of the
        normal noise
    '''

    t = np.arange(cfg.t_0, cfg.t_f, max_step)
    y = []
    for x in t:
        y.append(get_input_value(input_type, x))
    y = np.array(y)
    if noise is not False:
        y = y + np.random.normal(scale=noise, size=y.size)
        y = np.clip(y, 0, np.inf)
    return t, y





def pulse_input(t):
    '''
    This function tells us what value glutamate takes at a given time t
    Currently set to square wave
    '''
    if(t > cfg.input_start and t < (cfg.input_start + cfg.input_duration)):
        glut = cfg.input_max
    else:
        glut = cfg.input_min
    return glut





def oscillation_input(t):
    '''
    This function creates an oscillatory glutamate input
    '''
    principle_t = (t) % (cfg.oscillation_on_duration + cfg.oscillation_off_duration)
    
#     print((cfg.oscillation_off_duration + cfg.oscillation_on_duration) * cfg.num_oscillations)
    if(principle_t < cfg.oscillation_on_duration and
        t < (cfg.oscillation_off_duration + cfg.oscillation_on_duration) * cfg.num_oscillations):
        return cfg.input_max
    else:
        return cfg.input_min
    
    



def noisy_oscillation_input(t):
    '''
    This function creates an oscillatory glutamate input, but multiply the spritz
    strength by a normal distributed noise
    '''
    principle_t = (t) % (cfg.oscillation_on_duration + cfg.oscillation_off_duration)
    
    noise = np.random.normal(loc=1, scale=cfg.oscillation_noise_scale)

#     print((cfg.oscillation_off_duration + cfg.oscillation_on_duration) * cfg.num_oscillations)
    if(principle_t < cfg.oscillation_on_duration and
        t < (cfg.oscillation_off_duration + cfg.oscillation_on_duration) * cfg.num_oscillations):
        return cfg.input_max * noise
    else:
        return cfg.input_min



    

def exponential_oscillation(t):
    '''
    This function creates an oscillatory glutamate input where the input grows
    and falls exponentially
    '''
    #compute our IP3 curve helper functions
    oscillation_on_duration = cfg.d_rise + cfg.d_decay
    principle_t = t % (oscillation_on_duration + cfg.oscillation_off_duration)
    
    if(principle_t < oscillation_on_duration and
      t < (oscillation_on_duration + cfg.oscillation_off_duration) * cfg.num_oscillations):
        s_inf = cfg.A / (1 - np.exp(-cfg.r_rise * cfg.d_rise))
        if(principle_t < cfg.d_rise):
            return s_inf * (1 - np.exp(-cfg.r_rise * (principle_t)))
        elif(principle_t >= cfg.d_rise):
            return cfg.A * np.exp(-cfg.r_decay * (principle_t - cfg.d_rise))
    else:
        return 0





def exponential_oscillation2(t):
    '''
    exponential oscillation using the exponential input ramp up and down
    
    note that exponential_input(t) uses cfg.input_start - cfg.input_smoothing as the start point (ramps up
        over the course of cfg.input_smoothing to reach max at cfg.input_start)
        so to get spikes to work, need to make sure we pass in a principle time that allows for negatives
        exponential_input(t) also uses cfg.input_start and cfg.input_duration, so make sure
            cfg.input_start = 0
            cfg.input_duration = 0
        in order for this to work properly and create spikes that directly ramp up and down
        if wanting actually exponential pulse oscillations with flat tops, make sure
            cfg.oscillation_on_duration = cfg.input_duration + 2*cfg.input_smoothing
            cfg.input_start = 0
            
    '''
    principle_t = (t) % (cfg.oscillation_on_duration + cfg.oscillation_off_duration) - cfg.input_smoothing
    
#     print((cfg.oscillation_off_duration + cfg.oscillation_on_duration) * cfg.num_oscillations)
    if(principle_t < cfg.oscillation_on_duration and
        t < (cfg.oscillation_off_duration + cfg.oscillation_on_duration) * cfg.num_oscillations):
        return exponential_input(principle_t)
    else:
        return cfg.input_min





def exponential_train(t):
    '''
    Get exponential input based on a Poisson spike train
    Will automatically generate a train of spike intervals if one is not saved to cfg.train
        we can also generate a new train by calling generate_train_times
    Due to the way exponential_input is written, recommended config values are
    cfg.input_start = 1
    cfg.input_smoothing = 1
    cfg.input_duration = 2
    cfg.input_max = 2
    '''
    try:
        cfg.train
    except:
        cfg.train = generate_train_spike_times()
        
    if t > np.max(cfg.train):
        idx = len(cfg.train) - 1
    else:
        idx = np.argmax(cfg.train > t) - 1
        if idx < 0:
            return cfg.input_min
    
    t = t - cfg.train[idx]
    return exponential_input(t)


def train(t):
    '''Measure train input as the sum of individual inputs in the past
    input_duration time'''
    num_spikes_active = ((cfg.train > t - cfg.input_duration) & (cfg.train < t)).sum()
    output = num_spikes_active * cfg.input_max
    return output

def exponential_input(t):
    '''
    This function will generate an exponential increase and decay glutamate input
    Used only by exponential_train function
    '''
    input_half = cfg.input_duration / 2
    input_end = cfg.input_start + cfg.input_duration
    if(t >= cfg.input_start - cfg.input_smoothing and t < (cfg.input_start)):
        return cfg.input_max * (1 - np.exp((t - (cfg.input_start - cfg.input_smoothing)) / (t - cfg.input_start)))
    elif(t >= (input_end) and t < (input_end + cfg.input_smoothing)):
        return cfg.input_max * np.exp((t - input_end) / (t - (input_end + cfg.input_smoothing)))
    elif(t >= (cfg.input_start) and t <= (input_end)):
        return cfg.input_max
    else:
        return cfg.input_min



def next_time(rate):
    return -math.log(1 - random.random()) / rate

def generate_train(rate, t_f=300):
    total_time = 0
    spikes = []
    while True:
        nxt = next_time(rate)
        total_time += nxt
        if total_time > t_f:
            return np.array(spikes)
        spikes.append(total_time)


def generate_train_spike_times(seed=0):
    """
    Generate a set of Poisson spike train times 
    """
    # Seed the generator for consistent results.
    np.random.seed(seed)

    # We need to start with some sort of r(t). Presumably the asker of this
    # question has access to this.
    r = np.full(4000, 100)

    # Define our time interval, delta t. Use one millisecond as suggested
    # in the paper.
    dt = 0.001

    # Draw 7 random numbers on the interval [0, 1)
    x = np.random.rand(4000)

    # Initialize output.
    spikes = np.zeros_like(r)

    # If x[i] <= r * dt, generate a spike.
    mask = x <= r * dt

    # Set to 1.
    spikes[mask] = 1
    
    interval = 0
    intervals = []
    for spike in spikes:
        if spike == 1:
            intervals.append(interval)
            interval = 0
        else:
            interval += 1
    spike_times = np.cumsum(intervals)
    return spike_times


def curve_input(t):
    '''
    This function tell us what value IP3 takes at a given time t
    if we want to expicitly use IP3 as an input to compute
    the calcium dynamics for a given IP3 transient
    To set the parameters for this curve, use set_ip3_curve()
    '''
    #compute our IP3 curve helper functions
    if(cfg.A != 0):
        s_inf = cfg.A / (1 - np.exp(-cfg.r_rise * cfg.d_rise))
        if(t < cfg.t_star):
            return 0
        elif(t >= cfg.t_star and t < (cfg.t_star + cfg.d_rise)):
            return s_inf * (1 - np.exp(-cfg.r_rise * (t - cfg.t_star)))
        elif(t >= (cfg.t_star + cfg.d_rise)):
            return cfg.A * np.exp(-cfg.r_decay * (t - cfg.t_star - cfg.d_rise))
    else:
        #return steady state IP3
        return 0.056767761
    
    



def step_input(t):
    '''
    This function will create an IP3 input that increases stepwise
    0-50s: 0
    50-200s: 0.125
    200-350s: 0.250
    350-500s: 0.375
    500-650s: 0.5
    650-800s: 0.375
    800-950s: 0.250
    950-1100s: 0.125
    1100s-: 0
    '''
    #times at which the concentration changes
    time_breaks = np.array([0, 50, 200, 350, 500, 650, 800, 950, 1100]) * cfg.step_time_scale
    #values the concentration changes to at each interval
    input_values = np.array([0, 0.25, 0.5, 0.75, 1, 0.75, 0.5, 0.25, 0]) * cfg.step_max_value

    for i in range(len(time_breaks)):
        #check if t is after last time_break
        if(i == len(time_breaks) - 1):
            if(t >= time_breaks[i]):
                return input_values[i]

        #check which interval t is in
        else:
            if(t >= time_breaks[i] and t < time_breaks[i+1]):
                return input_values[i]

    
    
def custom_input(t):
    '''
    Use a custom step input
    The cfg.custom_input_times tells us at what times the input switches
    The cfg.custom_input_values tells us what values the input switches too
    '''
    boolean = np.array(cfg.custom_input_times) > t
    if not np.any(boolean):
        index = len(cfg.custom_input_times) - 1
    else:
        index = np.argmax(boolean) - 1
    return cfg.custom_input_vals[index]



def set_ip3_curve(curve, t=100):
    '''
    This function will set the parameters for our ip3 curve
    curve: the type of curve desired, options are
        'singlepeak'/'SP': IP3 curve that should produce SP
        'multipeak'/'MP': IP3 curve that should produce MP
        'plateau'/'PL': IP3 curve that should produce PL
        'steadystate'/'SS': IP3 set to flat 0.05676
        
    t_star: when to start the IP3 transient, default is 100ms
    '''

    t_star = t
    curves = ['singlepeak', 'SP', 'multipeak', 'MP', 'plateau', 'PL']
    
    if(curve in curves):
        if(curve == 'singlepeak'):
            curve = 'SP'
        elif(curve == 'multipeak'):
            curve = 'MP'
        elif(curve == 'plateau'):
            curve = 'PL'
        cfg.A = ip3_curves[curve]['A']
        cfg.d_rise = ip3_curves[curve]['d_rise']
        cfg.d_decay = ip3_curves[curve]['d_decay']
        cfg.r_rise = ip3_curves[curve]['r_rise']
        cfg.r_decay = ip3_curves[curve]['r_decay']
        
    elif(curve in ['steadystate', 'SS']):
        cfg.A = 0
    else:
        raise Exception('No proper IP3 curve option was selected, check documentation')
        




def set_gpcr_params(set_num):
    '''
    Set params for GPCR to be either the primary (set_num=1) or secondary (set_num=2)
    '''
    if(set_num == 1):
        params = cfg.param_set_1
    elif(set_num == 2):
        params = cfg.param_set_2
    elif(set_num == 'old'):
        params = cfg.param_set_old
    else:
        raise Exception('No valid value given for set_num. Must pass either 1 or 2 or "old"')
        
    for param in params:
        setattr(cfg, param, params[param])

'''

############################

Experiment running functions

############################
These functions will each run an experiment with the given control parameter

!!NOTE: Because these functions save the rewsults 

run_ip3_controlled_experiment: Run ODE with the given IP3 input
run_Gstar_controlled_experiment: Run ODE with the given Gstar input
run_experiment: Run ODE with the given glutamate input

The results of the experiment will be saved to the global space e.g.
t: time value of each step of ODE solve
p, c, h, c_tot: value of variable at each step of ODE solve
t_input: time values used to plot the input parameter
    For example, if controlling IP3, we would plot
    plt.plot(t_input, p)
    to see what the IP3 input looked like
'''

def run_ip3_controlled_experiment(input_type, t_f=1000, max_step=0.1, classification=False):
    '''
    Run an experiment where ip3 is manually given by the specified input type
    If running a classification experiment, we use slightly different initial conds
    '''
    if classification:
        x_0 = cfg.x_02
    else:
        x_0 = cfg.x_0

    sol = scipy.integrate.solve_ivp(ip3_ca_ode_equations, [cfg.t_0, t_f], x_0, 
                                    args=[input_type], max_step=max_step)
    sol['glutamate_input_type'] = input_type
    
    cfg.t = sol.t
    cfg.c = sol.y[0]
    cfg.c_tot = sol.y[1]
    cfg.h = sol.y[2]
    cfg.c_er = (cfg.c_tot - cfg.c) * cfg.gamma

    
    cfg.t_input, cfg.p = get_input_plot(input_type)
    return sol
    



def run_Gstar_controlled_experiment(input_type, t_f=1000, max_step=0.1):
    '''
    Run an experiment where Gstar is manunally given by the specified input type
    '''
    sol = scipy.integrate.solve_ivp(all_ode_equations, [cfg.t_0, t_f], cfg.all_init, 
                                    args=['pulse', '', input_type], max_step=max_step)
    sol['glutamate_input_type'] = input_type
    
    cfg.t = sol.t
    cfg.c = sol.y[0]
    cfg.c_tot = sol.y[1]
    cfg.h = sol.y[2]
    cfg.p = sol.y[3]
    cfg.c_er = (cfg.c_tot - cfg.c) * cfg.gamma
    
    
    cfg.t_input, cfg.Gstar = get_input_plot(input_type)
    return sol



def continue_Gstar_controlled_experiment(input_type='pulse', t_cont=200, max_step=0.1):
    '''
    Continue experiment for t_f amount of time
    '''
    init = [cfg.c[-1], cfg.c_tot[-1], cfg.h[-1], cfg.p[-1], 0, 0, 0, 0]
    t_0 = cfg.t[-1]
    t_f = cfg.t[-1] + t_cont
    
    sol = scipy.integrate.solve_ivp(all_ode_equations, [t_0, t_f], init, 
                                    args=['pulse', '', input_type], max_step=max_step)
    sol['glutamate_input_type'] = input_type
    
    cfg.t = np.append(cfg.t, sol.t)
    cfg.c = np.append(cfg.c, sol.y[0])
    cfg.c_tot = np.append(cfg.c_tot, sol.y[1])
    cfg.h = np.append(cfg.h, sol.y[2])
    cfg.p = np.append(cfg.p, sol.y[3])
    cfg.c_er = np.append(cfg.c_er, (cfg.c_tot - cfg.c) * cfg.gamma)
    
    t_input, Gstar = get_input_plot(input_type)
    cfg.t_input = np.append(cfg.t_input, t_input)
    cfg.Gstar = np.append(cfg.Gstar, Gstar)
    return sol




def run_experiment(input_type='pulse', t_f=1000, max_step=0.5, old_gpcr=False, fluxes=False,
                    fix_c_er=False, noise=False, fix_c_er_func=None, multipliers=False):
    '''
    Run an experiment where glutamate is manually given by the specified input type
    After the dynamical system is simulated, the resulting time-series values will be
        saved into cfg library. For example, you could run
        plt.plot(cfg.t, cfg.c)
        to plot the time resulting time course of cytosolic calcium
    The plotting function plot_experiment_plots automatically plots using cfg values from
        this function
    
    input_type: what type of input should be used for the glutamate. See "Input Functions" section
        for a list of possible strings
    t_f: how long the experiment should be run for
    max_step: how much time between each step of the numerical solver. An experiment with
        max_step=0.5 or 1 is usually precise enough while running fast, but this can be
        decreased for experiments where stimulation is on a fine scale e.g. spritz experiments
    old_gpcr: if True, use gpcr model without downstream Gd2 activator
        IMPORTANT - when setting old_gpcr=True, we will change the parameter set 
        to old ones and then change them back to primary setting, so if secondary 
        is being used please remember to set secondary again after running an old experiment
    if noise is set to a value, add on a normally distributed noise with the passed value
        as the width of the noise
    fix_c_er: A value can be passed to artificially fix the c_er to be a constant value
    fix_c_er_func: A function can be passed to artificially fix c_er to be a function of time
        This function must take a single parameters t and return the c_er value for that time
    fluxes: if True, also compute the calcium channel fluxes with compute_fluxes() function
    '''    
    if(old_gpcr):
        set_gpcr_params('old')
    cfg.t_0 = 0
    cfg.t_f = t_f
    if noise is not False:
        t, y = get_input_plot(input_type, noise=noise)
        cfg.input_func = scipy.interpolate.interp1d(t, y, fill_value=0, bounds_error=False)
        cfg.t_input, cfg.glut = t, y
        cfg.noise = True
    else:
        cfg.noise = False

    sol = scipy.integrate.solve_ivp(all_ode_equations, [cfg.t_0, t_f], cfg.all_init, 
                                  args=[input_type, old_gpcr, False, fix_c_er, fix_c_er_func], max_step=max_step)
    
    cfg.t = sol.t
    cfg.c = sol.y[0] #cytosolic calcium
    cfg.c_tot = sol.y[1] #total calcium
    cfg.c_er = (cfg.c_tot - cfg.c) * cfg.gamma #ER calcium
    cfg.h = sol.y[2] #h parameter of IP3R channel
    cfg.p = sol.y[3] #cytosolic IP3
    cfg.Gstar = sol.y[4]
    cfg.Gd1 = sol.y[5]
    cfg.Gd2 = sol.y[6]
    cfg.lamb = sol.y[7] #downstream GPCR factor that influences Gd2
    cfg.G = 1 - sol.y[4] - sol.y[5] - sol.y[6]
    cfg.Gd = sol.y[5] + sol.y[6] #Gd1 + Gd2


    if noise is False:
        cfg.t_input, cfg.glut = get_input_plot(input_type, max_step=max_step)
    

    if multipliers:
        cfg.c = cfg.c * 1000
        cfg.p = cfg.p * 1000
        cfg.glut = cfg.glut * 1000

    if fluxes:    
        compute_fluxes()

    if(old_gpcr):
        set_gpcr_params(1)

    cfg.noise = False
    return sol




def continue_experiment(input_type='pulse', t_cont=200, max_step=0.1, fluxes=False):
    '''
    Continue simulating an experiment
    
    t_cont: how much more time to simulate for
    '''    
    cfg.t_0 = cfg.t[-1]
    t_f = cfg.t[-1] + t_cont
    cfg.t_f = t_f
    init = [cfg.c[-1], cfg.c_tot[-1], cfg.h[-1], cfg.p[-1], cfg.Gstar[-1], 
            cfg.Gd1[-1], cfg.Gd2[-1], cfg.lamb[-1]]
    sol = scipy.integrate.solve_ivp(all_ode_equations, [cfg.t_0, t_f], init, 
                                  args=[input_type], max_step=max_step)
    
    cfg.t = np.append(cfg.t, sol.t)
    cfg.c = np.append(cfg.c, sol.y[0])
    cfg.c_tot = np.append(cfg.c_tot, sol.y[1])
    cfg.c_er = (cfg.c_tot - cfg.c) * cfg.gamma
    cfg.h = np.append(cfg.h, sol.y[2])
    cfg.p = np.append(cfg.p, sol.y[3])
    cfg.Gstar = np.append(cfg.Gstar, sol.y[4])
    cfg.Gd1 = np.append(cfg.Gd1, sol.y[5])
    cfg.Gd2 = np.append(cfg.Gd2, sol.y[6])
    cfg.lamb = np.append(cfg.lamb, sol.y[7])
    cfg.G = np.append(cfg.G, 1 - sol.y[4] - sol.y[5] - sol.y[6])

    t_input, glut = get_input_plot(input_type, max_step=max_step)
    cfg.t_input = np.append(cfg.t_input, t_input)
    cfg.glut = np.append(cfg.glut, glut)
    
    if fluxes:
        compute_fluxes()

    return sol




def set_init(type='default'):
    '''
    Set initial conditions
    default: default steady state params
    c_t: change initial c_tot to 0.8 of its steady state val
    poisson: change initial to start from after a poisson spike train
    noise: after running with normal noisy of 0.03 sd
    '''
    if type == 'default':
        # cfg.all_init = [0.0951442, 34.841184, 0.673079, 0.056767761, 0, 0, 0, 0]
        cfg.all_init = [0.09013785, 35.744397, 0.66821744, 0.040422910, 0.0, 0.0, 0.0, 0.0]
    
    if type == 'c_t':
        # cfg.all_init = [0.0951442, 34.841184*0.8, 0.673079, 0.056767761, 0, 0, 0, 0]
        cfg.all_init = [0.09013785, 32, 0.66821744, 0.040422910, 0.0, 0.0, 0.0, 0.0]
        # cfg.all_init = [0.09013785, 32.0951, 0.66821744, 0.040422910, 0.0, 0.0, 0.0, 0.0]
    if type == 'poisson':
        #Note that I lost the actual poisson train experiment I used to use
        #but this one seems to work just fine
        # cfg.all_init = [0.11938227778586594, 26.370394600930247, 0.6167209424320252, 0.1442844225467311, 
        #     0.09114690337250154, 0.13599444638425676, 0.3918773661605453, 0.0034478630682041244]
        # cfg.all_init = [0.0951442, 26.370394600930247, 0.6167209424320252, 0.1442844225467311, 0, 0, 0, 0]
        load_experiment('rate_0.2_conc_1_2', verbose=False)
        cfg.all_init = [cfg.c[-1], cfg.c_tot[-1], cfg.h[-1], cfg.p[-1], cfg.Gstar[-1], cfg.Gd1[-1], cfg.Gd2[-1], cfg.lamb[-1]]
    if type == 'noise':
        cfg.all_init = [0.144095, 30.296470, 0.639210, 0.118324, 0.019922, 0.046037, 0.012016, 0.001115]




def run_gpcr_experiment(input_type='pulse', t_f=1000, max_step=0.1):
    '''
    Run an experiment to just analyze GPCR dynamics alone
    '''
    cfg.t_f = t_f
    sol = scipy.integrate.solve_ivp(gpcr_ode_equations, [cfg.t_0, t_f], cfg.gpcr_init,
            args=[input_type], max_step=max_step)

    cfg.t = sol.t
    cfg.Gstar = sol.y[0]
    cfg.Gd1 = sol.y[1]
    cfg.Gd2 = sol.y[2]
    cfg.lamb = sol.y[3]
    cfg.G = 1 - sol.y[0] - sol.y[1] - sol.y[2]
    cfg.t_input, cfg.glut = get_input_plot(input_type, max_step=max_step)

    return sol


def compute_fluxes():
    '''
    Compute the simulated fluxes from the ODEs based on the current
    data loaded in cfg
    '''
    p = cfg.p
    c = cfg.c
    c_er = cfg.c_er
    Gstar = cfg.Gstar
    h = cfg.h
    
    
    #First compute some derivative values that will be needed
    #for h and IP3R dynamics
    m_inf = p / (p + cfg.d1)
    n_inf = c / (c + cfg.d5)
    Q2 = cfg.d2 * (p + cfg.d1) / (p + cfg.d3)
    h_inf = Q2 / (Q2 + c)
    tau_h = 1 / (cfg.a2 * (Q2 + c))
    
    #Compute the fluxes through each channel
    J_ip3r = cfg.v_ip3r * (m_inf**3) * (n_inf**3) * (h**3) * (c_er - c)
    J_serca = cfg.v_serca * (c**1.75) / (c**1.75 + cfg.k_serca**1.75)
    J_pmca = cfg.v_pmca * (c**2) / (c**2 + cfg.k_pmca**2)
    J_soc = cfg.v_soc * (cfg.k_soc**4) / (cfg.k_soc**4 + c_er**4)
#     J_soc = cfg.v_soc * (cfg.k_soc**2) / (cfg.k_soc**2 + c_er**2)
    
    #leak fluxes
    J_er_leak = cfg.v_er_leak * (c_er - c) #ER <-> Cyt leak
    J_ecs_add = cfg.v_in - cfg.k_out * c #Cyt <-> extracellular leak
    
    ip3_production = cfg.v_beta*Gstar + cfg.v_delta*((cfg.k_delta)/(1 + p))*((c**2)/(c**2 + cfg.k_plcdelta**2))
    ip3_degradation = cfg.v_3k*((c**4)/(c**4 + cfg.k_d**4))*(p/(p+cfg.k_3)) + cfg.r_5p*p
    
    cfg.J_ip3r = J_ip3r
    cfg.J_serca = J_serca
    cfg.J_pmca = J_pmca
    cfg.J_soc = J_soc
    cfg.J_er_leak = J_er_leak
    cfg.J_ecs_add = J_ecs_add
    cfg.ip3_production = ip3_production
    cfg.ip3_degradation = ip3_degradation



def save_experiment(name, verbose=True):
    '''
    Save an experiment currently in cfg to data under the given name
    '''
    labels = ['t', 'c', 'c_tot', 'c_er', 'h', 'p', 'Gstar', 'Gd1', 'Gd2', 'Gd', 'lamb', 'G', 't_input', 'glut']
    save_dict = {}
    for label in labels:
        save_dict[label] = getattr(cfg, label)
    pickle.dump(save_dict, open('data/experiments/' + name, 'wb'))
    if(verbose):
        print('Experiment saved')
    
    

    
def load_experiment(name, verbose=True, flux=False):
    '''
    Load an experiment to cfg for plotting and analysis
    '''
    labels = ['t', 'c', 'c_tot', 'c_er', 'h', 'p', 'Gstar', 'Gd1', 'Gd2', 'Gd', 'lamb', 'G', 't_input', 'glut']
    load_dict = pickle.load(open('data/experiments/' + name, 'rb'))
    for label in labels:
        if label in load_dict:
            setattr(cfg, label, load_dict[label])
    if(flux):
        compute_fluxes()
    if(verbose):
        print('Experiment loaded')




def get_stats(variable):
    '''
    Get statistics for a variable from the last run experiment
    For example, pass 'Gstar' to get stats about Gstar
    returns max, integral
    '''
    y = getattr(cfg, variable)
    t = cfg.t
    high = np.max(y)
    integral = np.trapz(y, t)
    
    return high, integral

'''
###########################

Plotting Functions

###########################
These functions help with quickly plotting results from the numerical solver
or from XPP saved data
'''

def plot_experiment_plots(variables, axs=None, add_ylabels=True, add_xlabel=True, plot_input=True, 
                          ylabel_padding=[-0.4, 0.4], legend_label=None, color=None, linestyle='solid',
                          alpha=1, remove_yticks=False, remove_xticks=True, ret_ax=False, multipliers=True,
                          simple_ylabels=False):
    '''
    Plot the solutions of the numerical solver for multiple variables
    Use the passed axs, iterating through them one by one and plotting the variables in the given order
    
    variables: list of variable names e.g. ['Gstar', 'p', 'h', 'c']
    axs: list of axes to plot on. For example a column of plots from plt.subplots
    add_ylabels: whether to add ylabels to the plot
    add_xlabel: whether to add xlabel of time to the last axis
    plot_input: whether first plot is an input (then we will use t_input instead of t)
    ylabel_padding: padding for each ylabel (first is x direction, second is y direction)
    multipliers: whether to multiply glut, c, p by 1000 for plotting purposes
    This should only be used if more than one variable is being plotted
    
    color, linestyle, alpha, remove_yticks, remove_xticks: function the same as their 
        respective matplotlib plotting parameters
    ret_ax: whether to return the axis that is plotted on (useful if you want to make further
        adjustments to the axis after creating it here) 

    Ex.
    run_Gstar_controlled_experiment('pulse')
    fig, ax = plt.subplots(4, 1, figsize=(10, 5))
    plot_experiment_plots(['Gstar', 'p', 'h', 'c'], ax)
    '''
    
    # These keys represent each possible string that can be passed for plotting
    ylabels = {
        # 'glut': r'$\phi$',
        'glut': 'glut \n ($\mu$M)',
        'Gstar': r'$G^*$',
        'G': r'$G$',
        'Gd1': r'$G_\mathrm{d1}$',
        'Gd2': r'$G_\mathrm{d2}$',
        'Gd': r'$G_\mathrm{d,total}$',
        'lamb': r'$\lambda$',
        'p': r'IP$_3$',
        'h': r'$h$',
        # 'c': r'$c$',
        'c': '[Ca$^{2+}$]$_{cyt}$ \n ($\mu$M)',
        'c_tot': r'$c_\mathrm{tot}$',
        # 'c_er': r'$c_\mathrm{ER}$',
        'c_er': '[Ca$^{2+}$]$_{ER}$ \n (mM)',
        'J_ip3r': r'$J_\mathrm{IP3R}$',
        'J_serca': r'$J_\mathrm{SERCA}$',
        'J_pmca': r'$J_\mathrm{PMCA}$',
        'J_soc': r'$J_\mathrm{SOC}$',
        'J_er_leak': r'$J_\mathrm{ER leak}$',
        'J_ecs_add': r'$J_\mathrm{ECS add}$',
        'ip3_production': r'IP$_{3\mathrm{prod}}$',
        'ip3_degradation': r'IP$_{3\mathrm{deg}}$'
    }
    
    simple_ylabels = {
        # 'glut': r'$\phi$',
        'glut': r'$\phi$',
        'Gstar': r'$G^*$',
        'G': r'$G$',
        'Gd1': r'$G_\mathrm{d1}$',
        'Gd2': r'$G_\mathrm{d2}$',
        'Gd': r'$G_\mathrm{d,total}$',
        'lamb': r'$\lambda$',
        'p': r'IP$_3$',
        'h': r'$h$',
        # 'c': r'$c$',
        'c': '[Ca$^{2+}$]$_{cyt}$',
        'c_tot': r'$c_\mathrm{tot}$',
        # 'c_er': r'$c_\mathrm{ER}$',
        'c_er': '[Ca$^{2+}$]$_{ER}$',
        'J_ip3r': r'$J_\mathrm{IP3R}$',
        'J_serca': r'$J_\mathrm{SERCA}$',
        'J_pmca': r'$J_\mathrm{PMCA}$',
        'J_soc': r'$J_\mathrm{SOC}$',
        'J_er_leak': r'$J_\mathrm{ER leak}$',
        'J_ecs_add': r'$J_\mathrm{ECS add}$',
        'ip3_production': r'IP$_{3\mathrm{prod}}$',
        'ip3_degradation': r'IP$_{3\mathrm{deg}}$'
    }
    
    if type(axs) == type(None):
        fig, axs = plt.subplots(len(variables), 1, figsize=(10, 10))

    for i, variable in enumerate(variables):
        y = getattr(cfg, variable)

        if multipliers and variable in ['c', 'p']:
            y = y * 1000
        
        if(plot_input and i==0):
            axs[i].plot(cfg.t_input, y, label=legend_label, color=color, 
                linestyle=linestyle, alpha=alpha)
        else:
            axs[i].plot(cfg.t, y, label=legend_label, color=color, 
                linestyle=linestyle, alpha=alpha)
    
    if(add_ylabels):
        for i, variable in enumerate(variables):
            if simple_ylabels:
                ylabel = simple_ylabels[variable]
                rot = 90
            else:
                ylabel = ylabels[variable]
            
            
            axs[i].set_ylabel(ylabel, rotation='horizontal', ha='left', va='center')
            axs[i].get_yaxis().set_label_coords(ylabel_padding[0], ylabel_padding[1])
            
    if(add_xlabel):
        axs[len(variables) - 1].set_xlabel(r'$t$')

    if(remove_yticks):
        for i, variable in enumerate(variables):
            axs[i].set_yticks([])

    if(remove_xticks):
        for i, variable in enumerate(variables):
            if(i != len(variables) - 1):
                axs[i].set_xticks([])

    if ret_ax:
        return axs


    

def add_abc_to_subaxes(ax, text='A.', left=0, top=1.05):
    '''
    Add ABC type labels to specific subaxes similar to those
    from proplot
    
    Note that if creating a single ax as in fig, ax = pplt.subplots()
    the ax to pass in is add_abc_to_subaxes(ax[0]) since it works
    with subaxes specifically
    '''
    abc_kw = pplt.rc.fill(
        {
            'size': 'abc.size',
            'weight': 'abc.weight',
            'color': 'abc.color',
            'family': 'font.family',
        },
        context=True
    )
    border_kw = pplt.rc.fill(
        {
            'border': 'abc.border',
            'borderwidth': 'abc.borderwidth',
            'bbox': 'abc.bbox',
            'bboxpad': 'abc.bboxpad',
            'bboxcolor': 'abc.bboxcolor',
            'bboxstyle': 'abc.bboxstyle',
            'bboxalpha': 'abc.bboxalpha',
        },
        context=True,
    )
    kw = {'zorder': 3.5, 'transform': ax.transAxes}
    
    ax.text(left, top, text, **abc_kw, **border_kw, **kw)
    



def plot_bifurcation(filename, ax=None, increasing=True, ret=False, jump_lim=0.01):
    '''
    Plot the bifurcation data based on XPPAUTO file (from AUTO > File > Write Points)    
    This will automatically look in the folder 'data/bifurcations/'

    ax: optionally pass in a matplotlib axis to plot on
    increasing: if True, only plot the x axis values that are increasing
        to avoid repeating plotting oscillatory bifurcation curves that tend to double
        back on themselves
    ret: if True, return the dataframe
    jump_lim: amount by which points must jump before considered a split. If things
        are joining weirdly, decrease this value. If things are not being plotted, 
        increase this value

    Ex.
    plot_bifurcation('ip3_ca.dat')
    '''
    if(type(ax) == type(None)):
        ax = plt
    data = load_bifurcation_data(filename)
    
    for i, bifurc_type in enumerate(colors.keys()):
        d = data[data[3] == i+1]
        if(len(d) == 0):
            continue
        splits = np.where(np.abs(d[0].diff()) > jump_lim)[0] #check where the plot splits up and plot 
                                                  #each individual portion of the data
        if(increasing):
            if(d[0].diff().sum() < 0):
                final_index = np.argmin(d[0])
            else:
                final_index = np.argmax(d[0]) #find where the x-axis value no longer increasing
            if(final_index == 0):
                final_index = len(d[0])
#             print(final_index)
            
        if(len(splits) == 0):
            if(increasing):
                end = final_index
            else:
                end = None

            ax.plot(d.iloc[:end, 0], d.iloc[:end, 1], c=colors[bifurc_type], linestyle=linestyles[bifurc_type])
            ax.plot(d.iloc[:end, 0], d.iloc[:end, 2], c=colors[bifurc_type], linestyle=linestyles[bifurc_type])
        else:
            #plot the branches of the bifurcation type
            for j in range(len(splits) + 1):
                if(j == 0):
                    start = 0 #first split starts at 0
                else:
                    start = splits[j - 1]

                if(j == len(splits)):
                    end = None #last split contains all remaining data points
                else:
                    end = splits[j]
                
                if(increasing and (end == None or end > final_index)):
                    end = final_index
                    
                ax.plot(d.iloc[start:end, 0], d.iloc[start:end, 1], c=colors[bifurc_type], linestyle=linestyles[bifurc_type])
                ax.plot(d.iloc[start:end, 0], d.iloc[start:end, 2], c=colors[bifurc_type], linestyle=linestyles[bifurc_type])
    if(ret):
        return data
    
    



def load_bifurcation_data(filename, folder='data/bifurcations/'):
    '''
    Helper function called by plot_bifurcation
    Load the bifurcation data file given by XPP Auto
    in the data/bifurcations folder with the given filename
    '''
    data = []
    with open(folder + filename, 'rb') as f:
        for line in f:
            line = line.split()
            data_dict = {}
            for i, n in enumerate(line):
                data_dict[i] = float(n)
            data.append(data_dict)

    data = pd.DataFrame(data)
    return data




def plot_nullcline(filename, ax=None, ret=False, plot=True):
    '''
    Plot the nullclines from an XPP file
    Will automatically look in the folder 'data/nullclines/' for .dat file

    ax: optional matplotlib axis to plot the nullclines on
    ret: (True/False) whether to return the nullcline data
        'anim': pass this to return the data in specific form to be animated
    plot: whether to actually plot the nullcline

    Note, we can run this with ret=True, plot=False to more easily get
    the nullcline data to work with manually

    '''
    if(type(ax) == type(None)):
        ax = plt
    data = load_nullcline(filename)
    
    if(plot):
        for nc_type in [1, 2]:
            indexes = data[data['type'] == nc_type].index

            if(nc_type == 1):
                color = 'red'
                label = 'p'
            elif(nc_type == 2):
                color = 'green'
                label = 'c'
            for i in indexes: 
                ax.plot([data.loc[i,'x1'], data.loc[i,'x2']], [data.loc[i,'y1'], data.loc[i,'y2']], c=color,
                        label=label)
                label = None
    
    
    if(ret is True):
        return data
    if(ret == 'anim'):
        plot_data = []
        for nc_type in [1, 2]:
            data_block = data[data['type'] == nc_type]
            plot_data.append([data_block['x1'], data_block['y1']])
        return plot_data




def load_nullcline(filename, folder='data/nullclines/'):
    '''
    Helper function called by plot_nullcline
    Gets the nullcline data from data/nullclines
    '''
    data = []
    with open(folder + filename, 'rb') as f:
        nc_type = 1
        even_odd = 0
        for line in f:
            if(line == b'# X-nullcline\n'):
                nc_type = 1
            elif(line == b'# Y-nullcline\n'):
                nc_type = 2
            else:
                line = line.split()
                if(len(line) > 0):
                    #this is a data line
                    if(even_odd == 0):
                        even_odd = 1
                        x1 = float(line[0])
                        y1 = float(line[1])
                    else:
                        even_odd = 0
                        x2 = float(line[0])
                        y2 = float(line[1])
                        data.append([x1, x2, y1, y2, nc_type])

    data = pd.DataFrame(data, columns=['x1', 'x2', 'y1', 'y2', 'type'])
    return data