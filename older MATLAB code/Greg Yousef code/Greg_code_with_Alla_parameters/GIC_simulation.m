%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script that runs the Handy, Taheri, Pitta 
% ODEs with varying glutamate input in the form of decaying exponential
% The varied parameters are spike_decay and spike_peak
% Then the calcium responses are classifed as SP, PL, MP, LL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

GIC_Params; % loads the global parameters
addpath('./Supporting_Functions')

%This is what Alla had. Might have come from Yousef.
% Can be vectors of times, amplitudes and decays (all must be the same
% length)
%stim_time = glutend-glutstart;
stim_amp = glutmax;
stim_decay = Inf;

%This is what Greg had. For the exponentially decaying pulse
% % glutamate parameters
% spike_times = [20];
% spike_peaks = [5];
% spike_decays = [2];
 spike_times = [0];
 spike_peaks = [0];
 spike_decays = [0];

%% Numerical parameters
% sets the simulation time 
dt = 0.005;
%sim_t = [0:dt:Total_time];
sim_t = [0:dt:2000];

% the simulation begins at steady state (based on default parameter values)
init_cond = [0.0949 34.8645 0.6731 0.056767761 0 0 0];
% runs the simulation using ode15s
%Reduced MaxStep for short pulses
%options = odeset('AbsTol', 10^-6, 'RelTol', 10^-6, 'MaxStep', 0.005); 
options = odeset('AbsTol', 10^-6, 'RelTol', 10^-6, 'MaxStep', 0.1);

%% Run the ODE file 
%[sim_t,ca_transient] = ode15s(@GIC_ODEs, sim_t, init_cond);

CaCyt = ca_transient(:,1);
ip = ca_transient(:,4);
Gstar = ca_transient(:,5);
Gd1 = ca_transient(:,6);
Gd2 = ca_transient(:,7);

% %% Classify the calcium response
% 
% [TrueTroughVal, TrueTroughLoc, PeakVal,PeakLoc]=plotCa_TH(CaCyt,...
%     sim_t, 'time');
% 
% %Prints the calcium response type
% [Result, CaDur, CaAmount, CaLatency, StartOfResp,EndOfResp]=...
%     FourCaResponseTypes_TH(TrueTroughVal,TrueTroughLoc, PeakVal,...
%     PeakLoc, CaCyt, sim_t, 20);
% 
% % save the response type
% all_results = {Result};
% 
%% Plot the data
glut = zeros(length(sim_t),1);
for j = 1:length(sim_t)
    for i = 1:length(spike_times)
        if sim_t(j) > spike_times(i)
            glut(j) = spike_peaks(i)*exp(-(sim_t(j)-spike_times(i))/spike_decays(i));
        end
    end    
end

figure(2);
subplot(3,1,1)
plot(sim_t,CaCyt,'linewidth',1.5)
axis([0 max(sim_t) 0 max(CaCyt)+0.1])
hold on;
plot(sim_t,glut*0.2,'--','linewidth',1.5)
xlabel('Time (s)')
ylabel('[Ca^{2+}] (\muM)') 
set(gca,'fontsize',16)
legend('Ca^{2+}','Glut')

subplot(3,1,2)
plot(sim_t,ip,'linewidth',1.5)
xlabel('Time (s)')
ylabel('[IP_3] (\muM)') 
set(gca,'fontsize',16)

subplot(3,1,3)
plot(sim_t,Gstar,'linewidth',1.5)
xlabel('Time (s)')
ylabel('Gstar') 
set(gca,'fontsize',16)





