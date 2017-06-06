%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testStepCurrent.m 
% --- for simulating step current response waveforms 
% Created: May 16, 2017 by Go Ashida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% time parameters
DT = 0.002; % [ms] time step
T0 = 20.0; % [ms] time without stimulus
T1 = 30.0; % [ms] time with stimulus
N0 = round(T0/DT); % steps 
N1 = round(T1/DT); % steps 
Ntot = N0 + N1 + N0 + N1 + N0 + N1 + N0; % three stimuli
tv = (0:Ntot)*DT; % time vector [ms]; caution:length=Ntot+1

%% stimulus amplitudes (for each model)
ampPIF = [ 350, 500, 650 ]; % [pA]
ampAIF = [ 350, 500, 650 ]; % [pA]
ampOWC = [ 300, 600, 1000 ];% [pA]
ampAWC = [ 200, 400, 600 ]; % [pA]

%% stimulus vector 
iextPIF = [ zeros(1,N0), ones(1,N1)*ampPIF(1), ...
            zeros(1,N0), ones(1,N1)*ampPIF(2), ...
            zeros(1,N0), ones(1,N1)*ampPIF(3), ...
            zeros(1,N0+1) ]; 
iextAIF = [ zeros(1,N0), ones(1,N1)*ampAIF(1), ...
            zeros(1,N0), ones(1,N1)*ampAIF(2), ...
            zeros(1,N0), ones(1,N1)*ampAIF(3), ...
            zeros(1,N0+1) ]; 
iextOWC = [ zeros(1,N0), ones(1,N1)*ampOWC(1), ...
            zeros(1,N0), ones(1,N1)*ampOWC(2), ...
            zeros(1,N0), ones(1,N1)*ampOWC(3), ...
            zeros(1,N0+1) ]; 
iextAWC = [ zeros(1,N0), ones(1,N1)*ampAWC(1), ...
            zeros(1,N0), ones(1,N1)*ampAWC(2), ...
            zeros(1,N0), ones(1,N1)*ampAWC(3), ...
            zeros(1,N0+1) ]; 

spEx = zeros(1,Ntot+1); % dummy synaptic input
spIn = zeros(1,Ntot+1); % dummy synaptic input

%% calling each model 
[spPIF, vPIF] = LSOmodelPIF(spEx, spIn, DT, iextPIF);
[spAIF, vAIF] = LSOmodelAIF(spEx, spIn, DT, iextAIF);
[spOWC, vOWC] = LSOmodelOWC(spEx, spIn, DT, iextOWC);
[spAWC, vAWC] = LSOmodelAWC(spEx, spIn, DT, iextAWC);

%% plotting
figure(101); 

subplot(4,1,1); 
cla; hold on; 
plot(tv,vPIF); 
ylim([-70,-10]);
xlim([0, 4*T0+3*T1]);
title(sprintf('Passive IF model: %.0f / %.0f / %.0f [pA]',...
      ampPIF(1),ampPIF(2),ampPIF(3)));

subplot(4,1,2); 
cla; hold on; 
plot(tv,vAIF); 
ylim([-70,-10]);
xlim([0, 4*T0+3*T1]);
title(sprintf('Active IF model: %.0f / %.0f / %.0f [pA]',...
      ampAIF(1),ampAIF(2),ampAIF(3)));

subplot(4,1,3); 
cla; hold on; 
plot(tv,vOWC); 
ylim([-80,40]);
xlim([0, 4*T0+3*T1]);
title(sprintf('Original Wang-Colburn model: %.0f / %.0f / %.0f [pA]',...
      ampOWC(1),ampOWC(2),ampOWC(3)));

subplot(4,1,4); 
cla; hold on; 
plot(tv,vAWC); 
ylim([-80,40]);
xlim([0, 4*T0+3*T1]);
title(sprintf('Adjusted Wang-Colburn model: %.0f / %.0f / %.0f [pA]',...
      ampAWC(1),ampAWC(2),ampAWC(3)));
xlabel('time [ms]');
ylabel('potential [mV]');
