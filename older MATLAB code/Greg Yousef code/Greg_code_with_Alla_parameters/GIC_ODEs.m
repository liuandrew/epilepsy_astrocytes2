%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a file that contains the Handy, Taheri, Pitta ODE system
% Takes in the input times, glutamate peak and decay constant as parameters
% Different spikes can have different peaks and decays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = GIC_ODEs(t, x0)

CaCyt = x0(1); % calcium concentration in the cytosol
CaT = x0(2);   % total free calcium concentration in the cytosol
h = x0(3);     % deactivating variable for the IP3R
ip = x0(4);    % IP3 dynamics

Gstar = x0(5); % activated GPCRs 
Gd1 = x0(6);   % desensitized GPCR from Gstar (homologous) 
Gd2 = x0(7);   % desensitized GPCR from G (heterologous)
G = 1-Gstar-Gd1-Gd2; % unactivated GPCR

%% Define the global parameters 

global v_ip3r gamma v_leak v_in k_out
global d1 d2 d3 d5  a2 
global v_pmca k_pmca v_soc k_soc v_serca delta k_serca

global KPLCdelta Kdelta Vdelta V3k Kd K3
global r5p VB

global kp km kd1 kr1 kd2 kr2

%% Terms for Calcium Dynamics on ER:
CaER = (CaT-CaCyt)*gamma;
% terms for IP3R
minf = ip/(ip + d1);
ninf = CaCyt/(CaCyt + d5);
j_ip3r = v_ip3r*minf^3*ninf^3*h^3*(CaER - CaCyt);

% terms for the h ODE
q2 = d2*(ip+d1)/(ip+d3);
tauh = 1/(a2*(q2+CaCyt));
hinf = q2/(q2+CaCyt);

% Leak Term
j_leak = v_leak*(CaER-CaCyt);

%SERCA Pump
j_serca = v_serca*CaCyt^1.75/(CaCyt^1.75 + k_serca^1.75);

%% Terms for Calcium Dynamics on Plasma Membrane:

%PMCA pump
j_pmca = v_pmca*CaCyt^2/(k_pmca^2 + CaCyt^2);

%SOCC 
j_soc=v_soc*k_soc^4./(k_soc^4+CaER.^4);

%Leak Terms
j_out = k_out*CaCyt;
j_in = v_in;

%% Glutamate Input
% Realistic glutamate input (Decaying exponential)
glut = 0.1;
% for i = 1:length(spike_times)
%     if t > spike_times(i)
%         glut = spike_peaks(i)*exp(-(t-spike_times(i))/spike_decays(i));
%     end
% end

%% Terms governing IP3 Dynamics

% Production: PLCbeta + PLCdelta
IP3_Production = VB*Gstar + Vdelta/(1+ip/Kdelta)*(CaCyt^2/(CaCyt^2+KPLCdelta^2));
exit

% Degradation: Negative Ca^2+ feedback + natural degradation
IP3_Degradation = V3k*CaCyt^4/(CaCyt^4+Kd^4)*ip/(ip+K3)+r5p*ip;

%% ODEs
xdot(1) = j_ip3r-j_serca+j_leak + (j_in-j_out+j_soc-j_pmca).*delta; %ODE for [Ca]cyt
xdot(2) = (j_in-j_out+j_soc-j_pmca)*delta; %ODE for [Ca]ER
xdot(3) = (hinf - h)/tauh; %ODE for h

xdot(4) = IP3_Production - IP3_Degradation; % ODE for IP3

xdot(5) = kp*glut*G - km*Gstar - kd1*Gstar;     % ODE for Gstar
xdot(6) = kd1*Gstar - kr1*Gd1;      % ODE for Gd1
xdot(7) = kd2*Gstar*G - kr2*Gd2;    % ODE for Gd2

xdot=xdot';
end

