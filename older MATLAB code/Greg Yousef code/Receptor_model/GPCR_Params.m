%%
%  This file defines all global paramers used in the 
%  GPCR-IP3-Ca (GIC) model 
%%
clear;

global kp km kd1 kr1 kd2 kr2

global v_ip3r gamma v_leak v_in k_out
global d1 d2 d3 d5  a2 
global v_pmca k_pmca v_soc k_soc v_serca delta k_serca

global KPLCdelta Kdelta Vdelta V3k Kd K3
global r5p VB

global glutmax glutmin glutstart glutdur between_pulses Npuls

global Total_time

%% Glutamate parameters for square pulse 

glutmax = 10;
glutmin =  0;
glutstart = 10;
glutdur=750;
%glutend = glutstart+glutdur;
between_pulses=0;

% prompt='Number of Pulses? \n';
% Npuls=input(prompt);
Npuls=1;

Total_time=glutstart+Npuls*(glutdur+between_pulses)+100;

%% Parameters for the GPCR 
kp=.03;   % activation rate from G to Gstar
km=.04;   % deactivation rate from Gstar to G

kd1 = 0.01; % homologous desensitization rate from Gstar to Gd1
kr1 = .005;%.003;   % recovery rate from Gd1 to G

kd2 = 0.003;%0.0025;   % heterologous desensitization rate from G to Gd2
kr2 = 0.0007;%0.0004; % recovery rate from Gd2 to G



%% Parameters for IP3 Dynamics 

% Agonist-independent IP3 production parameters
KPLCdelta = 0.1;
Kdelta = 1.5;
% Maximal rate of IP3production by PLCdelta
Vdelta = 0.01; % decreased from the original value 0.02

% IP3 degradation parameters
V3k = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kd = 0.7; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K3 = 1;

% rate of IP3 degradation by IP-5P
r5p = 0.08;  % increased from the original value 0.04

% Agonist-dependent IP3 production parameters
VB = 0.2; % Maximal rate of IP3production by PLCbeta

%% Ca^2+ Parameters

% Max IP3 receptor flux
v_ip3r=0.222;
% (Cyt vol) / (ER vol)
gamma=5.4054; 

% Leak for ER
v_leak=0.002; 

% Leak for Extracellular Space
v_in=0.05; 
k_out=1.2;

% Li-Rinzel Parameters
d1=0.13; d2=1.049; d3=943.4e-3; d5=0.08234; 
a2=0.04; %adjusted Li-Rinzel Parameter

% PMCA Terms
v_pmca=10; k_pmca=2.5;

% SOCC Terms
v_soc=1.57;k_soc=90;

% SERCA Terms
v_serca=0.9;
 
% Sneyd Parameter
delta=0.2;

% SERCA Ca2+affinity
k_serca=0.1;

