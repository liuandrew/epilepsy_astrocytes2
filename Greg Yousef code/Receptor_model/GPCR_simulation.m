%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script that runs GPCR_ODEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear; clc; close all;


GPCR_Params;  % loads the global parameters

% Can be vectors of times, amplitudes and decays (all must be the same
% length)
%stim_time = glutend-glutstart;
stim_amp = glutmax;
stim_decay = Inf;



%% Numerical parameters
% sets the simulation time 
dt = 0.1;
sim_t = [0:dt:Total_time];

% the simulation begins at steady state (based on default parameter values)
init_cond = [0 0 0];
% runs the simulation using ode15s
options = odeset('AbsTol', 10^-6, 'RelTol', 10^-6, 'MaxStep', 0.1);

%% Run the ODE file 



[sim_t,ca_transient1] = ode15s(@GPCR_ODEs, sim_t, init_cond,options);

Gstar = ca_transient1(:,1);
Gd1 = ca_transient1(:,2);
Gd2 = ca_transient1(:,3);
G = 1-Gstar-Gd1-Gd2;

%%
% close all;
% close all
glut = zeros(length(sim_t),1);
period=glutdur+between_pulses;

for k = 1:length(sim_t)
    if sim_t(k)>glutstart && mod(sim_t(k)-glutstart,period)< glutdur && sim_t(k)-glutstart<Npuls*period
        glut(k)=glutmax;
    else  
        glut(k)=glutmin;
    end
end

figure(2)
subplot(5,1,1)
plot(sim_t,glut,'linewidth',1.5)
axis([0 max(sim_t) 0 max(glut)+.01])
%axis([0 1000 0 max(glut)+.01])
%xlabel('Time (s)')
ylabel('[Glut] (\muM)') 
set(gca,'fontsize',8)
hold on
subplot(5,1,2)
plot(sim_t,Gstar,'linewidth',1.5)
%xlabel('Time (s)')
ylabel('G*') 
%set(gca,'fontsize',16)
axis([0 max(sim_t) 0 1])
grid on
hold on
subplot(5,1,3)
plot(sim_t,G,'linewidth',1.5)
axis([0 max(sim_t) 0 1])
%xlabel('Time (s)')
ylabel('G') 
%set(gca,'fontsize',16)
grid 
hold on
subplot(5,1,4)
plot(sim_t,Gd1,'linewidth',1.5)
axis([0 max(sim_t) 0 max(Gd1)+.001])
%set(gca,'fontsize',16)
ylabel('Gd1')
hold on
subplot(5,1,5)
plot(sim_t,Gd2,'linewidth',1.5)
%set(gca,'fontsize',16)
axis([0 max(sim_t) 0 max(Gd2)+.001])
ylabel('Gd2')
hold on

% subplot(4,1,1)
% plot(sim_t,Gstar,'linewidth',1.5)
% set(gca,'fontsize',16)
% xlabel('Gstar')
% 
% subplot(4,1,2)
% plot(sim_t,G,'linewidth',1.5)
% set(gca,'fontsize',16)
% xlabel('G')
% 
% subplot(4,1,3)
% plot(sim_t,Gd1,'linewidth',1.5)
% set(gca,'fontsize',16)
% xlabel('Gd1')
% 
% subplot(4,1,4)
% plot(sim_t,Gd2,'linewidth',1.5)
% set(gca,'fontsize',16)
% xlabel('Gd2')
% 
