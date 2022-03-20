'''
#############################

Initiate parameter values

#############################
'''

#----------------------------
#Glutamate -> GPCR parameters
#----------------------------
 
#G <-> G*
# kp = 0.03 #activation rate
# km = 0.04 #inactivation rate
kp = 0.02
km = 0.2

#Old parameter set before lambda is added
# #G* -> Gd1 -> G
# kd1 = 0.01 #homologous (Gd1) deactivation rate
# kr1 = 0.005 #0.003 #homologous reactivation rate
# #G -> Gd2 -> G
# kd2 = 0.003 #0.0025 #heterologous (Gd2) deactivation rate
# kr2 = 0.0007 #0.0004 #heterologous reactivation rate

#New parameter set with lambda (downstream Gd2 activator)
#G* -> Gd1 -> G
kd1 = 0.02 #homologous (Gd1) deactivation rate
kr1 = 0.01 #homologous reactivation rate
#G -> Gd2 -> G
kd2 = 1.2 #heterologous (Gd2) deactivation rate
kr2 = 0.005 #heterologous reactivation rate
#lambda production and degradation
klp = 0.0002 #production
klm = 0.004 #degradation

#----------
#Two parameter sets for GPCR + set for old GPCR model
#----------

param_set_1 = {
    'kp': kp,
    'km': km,
    'kd1': 0.02,
    'kr1': 0.01,
    'kd2': 0.6,
    'kr2': 0.005,
    'klp': 0.0002,
    'klm': 0.004
}

param_set_2 = {
    'kp': kp,
    'km': km,
    'kd1': 0.02,
    'kr1': 0.03,
    'kd2': 0.6,
    'kr2': 0.005,
    'klp': 0.0006,
    'klm': 0.012
}

param_set_old = {
    'kp': kp,
    'km': km,
    'kd1': 0.01,
    'kr1': 0.005,
    'kd2': 0.003,
    'kr2': 0.0007,
    'klp': 0.0002, #arbitrary, not used
    'klm': 0.0004  #arbitrary
}


#-------------------------------
#IP3 -> Ca2+ dynamics parameters
#-------------------------------

gamma = 5.4054 #(cyt vol) / (ER vol)

#----------
#Internal (Cytosol/ER)
#----------
#IP3R channel (Ca2+ ER -> Ca2+ Cyt)
v_ip3r = 0.222 #max IP3R channel flux 

#Li-Rinzel parameters for IP3R channel
d1 = 0.13
d2 = 1.049
d3 = 0.9434
d5 = 0.08234
a2 = 0.04

#SERCA pump (Ca2+ Cyt -> Ca2+ ER)
v_serca = 0.9 #max SERCA flux
k_serca = 0.1 #half-saturation of Ca2+ for SERCA

#Leak (Ca2+ ER <-> Ca2+ Cyt)
v_er_leak = 0.002 #concentration gradient leak

#----------
#External (Cytosol/Extracellular)
#----------
#Leak (Ca2+ Extracellular <-> Ca2+ Cyt)
v_in = 0.05 #constant inward leak (extra -> cyt)
k_out = 1.2 #concentration based outward leak (cyt -> extra)

#PMCA pump (Ca2+ Cyt -> Ca2+ Extracellular)
v_pmca = 10 #max PMCA pump flux
k_pmca = 2.5 #half-saturation of Ca2+ for PMCA

#SOC channel (Ca2+ Extracellular -> Ca2+ Cyt)
v_soc = 1.57 #max SOC channel flux
k_soc = 90 #half-saturation of Ca2+ for SOC


#Ratio of ER to Extracellular (Internal vs External) transmission rates
delta = 0.2


#----------------
#Input parameters
#----------------
#These parmeters are used for different shapes of input

# ------------------------
# Double exponential curve
# ------------------------

#Initiate these variables
A = 0
d_rise = 0
r_rise = 0
d_decay = 0
r_decay = 0
t_star = 0 #when to start the IP3 transient



#--------------------------------------------
#Step input, pulse input and oscillation input
#--------------------------------------------
step_max_value = 1 #maximum value during the ip3 step input
step_time_scale = 1


#Square wave glutamate pulse parameters
input_start = 10
input_duration = 50
input_max = 0.3
input_min = 0
input_smoothing = 10

#Oscillation parameters
num_oscillations = 10
oscillation_on_duration = 50
oscillation_off_duration = 150

oscillation_noise_scale = 0.05

#Custom input
custom_input_times = [0, 50, 100]
custom_input_vals = [0.05, 0.15, 0]

#
#------------------------------
#IP3 Generation and Degradation
#------------------------------
#
#Production parameters
v_beta = 0.2
v_delta = 0.005
k_delta = 1.5
k_plcdelta = 0.1

#Degradation parameters
# v_3k = 2,
v_3k = 0.1
k_d = 0.7
k_3 = 1
r_5p = 0.08
# r_5p = 0.12

#-----------
#Initial conditions
#-----------
t_0 = 0
gpcr_0 = [0, 0, 0] #starting conditions for glutamate GPCR receptor

#calcium starts at steady state
#ca_0 = [0.0949, 34.8645, 0.6731] #starting conditions for calcium transients
ca_0 = [0.0951442, 34.841184, 0.673079] #slightly more precise based on our numerical solve
x_0 = [0.0951442, 34.841184, 0.673079, 0.056767761] #initial steady state conditions for IP3/Ca2+
x_02 = [0.086541496665789, 36.490839775010841, 0.625512446053023, 0]



#initial conditions for full IP3/Ca/GPCR system
# all_init = [0.0951442, 34.841184, 0.673079, 0.056767761, 0, 0, 0, 0]
all_init = [0.09013785, 35.744397, 0.66821744, 0.040422910, 0.0, 0.0, 0.0, 0.0]
# = [c, c_tot, h, p, Gstar, Gd1, Gd2, lambda],
all_init_no_pos_no_neg = [0.08650867, 36.4904174, 0.62551797, 0, 0, 0, 0, 0]
all_init_no_neg = [0.095338464, 34.8106342748, 0.673122891, 0.05723924177, 0, 0, 0, 0]
all_init_no_pos = [0.0865086834, 36.490417408, 0.625517973, 0, 0, 0, 0, 0]

#initial conditions for GPCR system
gpcr_init = [0, 0, 0, 0]

#Final time
t_f = 1000


noise = False
