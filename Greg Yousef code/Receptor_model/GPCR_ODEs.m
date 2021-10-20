%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains the odes for the G-protein coupled receptor
% submodel. Allows for both homologous and heterologous desensitization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = GPCR_ODEs(t, x0)

Gstar = x0(1);  % activated GPCR
Gd1 = x0(2);    % desensitized GPCR from Gstar (homologous)
Gd2 = x0(3);    % desensitized GPCR from G (heterologous)
 
G = 1-Gstar-Gd1-Gd2; % unactivated GPCR
%% Global parameters
global kp km
global kd1 kr1
global kd2 kr2
global glutmax glutmin glutstart glutdur between_pulses Npuls
%% Stimulus
% glut = 0;
% for i = 1:length(stim_time)
%     if t > stim_time(i)
%         glut = stim_amp(i)*exp(-(t-stim_time(i))/stim_decay(i));
%     end
% end

period=glutdur+between_pulses;
if t>glutstart && mod(t-glutstart,period)< glutdur && t-glutstart<Npuls*period
    glut=glutmax;
else  
    glut=glutmin;
end
%% ODEs
xdot(1) = kp*glut*G - km*Gstar-kd1*Gstar;     % ODE for Gstar
xdot(2) = kd1*Gstar - kr1*Gd1;      % ODE for Gd1
xdot(3) = kd2*Gstar*G - kr2*Gd2;    % ODE for Gd2

xdot=xdot';
end

