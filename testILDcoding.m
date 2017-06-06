%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testILDcoding.m 
% --- for simulating ILD coding of LSO models  
% Created: May 15, 2017 by Go Ashida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caution: Running this script may take several minutes (or even longer 
%          depending on your system). Be patient. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% time parameters
DT = 0.002; % [ms] 
Tinit = 80.0; % [ms]
Tmain = 2000.0; % [ms]
Tlast = 20.0; % [ms]
Ninit = round(Tinit/DT); % steps 
Nmain = round(Tmain/DT); % steps 
Nlast = round(Tlast/DT); % steps 
Ntot = Ninit+Nmain+Nlast; 
tv = (0:Ntot)*DT; % time vector [ms]; caution:length=Ntot+1
lmain = logical( [zeros(1,Ninit),ones(1,Nmain),zeros(1,Nlast+1)] );

%% frequency parameters 
iSPL = +35; % [dB] ipsilateral level
iRT  = 30 + 240 ./(1+exp(-(iSPL -20)/6.0) ); % ipsilateral rate 
cSPLs = -10:10:50; % [dB] contralateral levels 
cRTs = 30 + 240 ./(1+exp(-(cSPLs-20)/6.0) ); % contralateral rate 
ILDs = cSPLs-iSPL; % [dB] ILD array 
FQ = 0; % no amplitude modulation 
VS = 0; % no phase-locking 
P0 = 0; % zero initial phase 
Mex = 20; % number of excitatory inputs 
Min = 8;  % number of inhibitory inputs 

%% data array for output rates
RToutCOC = zeros(1,length(ILDs));
RToutEST = zeros(1,length(ILDs));
RToutAST = zeros(1,length(ILDs));
RToutPIF = zeros(1,length(ILDs));
RToutAIF = zeros(1,length(ILDs));
RToutOWC = zeros(1,length(ILDs));
RToutAWC = zeros(1,length(ILDs));

%% main loop 

for i = 1:length(cSPLs)
 
 % spike input vectors
 cRT = cRTs(i) % contralateral (inhibitory) rate 
 spEx = sum( PhaseLock(Mex,length(tv),FQ,VS,iRT,P0,DT), 1 );
 spIn = sum( PhaseLock(Min,length(tv),FQ,VS,cRT,P0,DT), 1 );

 % calling LSO models 
 [spCOC, vCOC] = LSOmodelCOC(spEx, spIn, DT); 
 [spEST, vEST] = LSOmodelEST(spEx, spIn, DT); 
 [spAST, vAST] = LSOmodelAST(spEx, spIn, DT); 
 [spPIF, vPIF] = LSOmodelPIF(spEx, spIn, DT); 
 [spAIF, vAIF] = LSOmodelAIF(spEx, spIn, DT); 
 [spOWC, vOWC] = LSOmodelOWC(spEx, spIn, DT); 
 [spAWC, vAWC] = LSOmodelAWC(spEx, spIn, DT); 

 % getting the main part of the response
 dataCOC = spCOC(lmain); 
 dataEST = spEST(lmain); 
 dataAST = spAST(lmain); 
 dataPIF = spPIF(lmain); 
 dataAIF = spAIF(lmain); 
 dataOWC = spOWC(lmain); 
 dataAWC = spAWC(lmain); 
 
 % calculating output spike rates
 RToutCOC(i) = sum(dataCOC)*1000/Tmain;
 RToutEST(i) = sum(dataEST)*1000/Tmain;
 RToutAST(i) = sum(dataAST)*1000/Tmain;
 RToutPIF(i) = sum(dataPIF)*1000/Tmain;
 RToutAIF(i) = sum(dataAIF)*1000/Tmain;
 RToutOWC(i) = sum(dataOWC)*1000/Tmain;
 RToutAWC(i) = sum(dataAWC)*1000/Tmain;

end

%% display results
sprintf('COC model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutCOC),min(RToutCOC),max(RToutCOC)-min(RToutCOC))
sprintf('EST model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutEST),min(RToutEST),max(RToutEST)-min(RToutEST))
sprintf('AST model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutAST),min(RToutAST),max(RToutAST)-min(RToutAST))
sprintf('PIF model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutPIF),min(RToutPIF),max(RToutPIF)-min(RToutPIF))
sprintf('AIF model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutAIF),min(RToutAIF),max(RToutAIF)-min(RToutAIF))
sprintf('OWC model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutOWC),min(RToutOWC),max(RToutOWC)-min(RToutOWC))
sprintf('AWC model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutAWC),min(RToutAWC),max(RToutAWC)-min(RToutAWC))

%% plotting 
figure(113); 
cla; 

subplot(2,4,1);
plot(ILDs,RToutCOC,'b-');
title('Coincidence Counting Model');
xlim([-45,15]);
ylim([0,180]);
xlabel('ILD [dB]');
ylabel('rate [spikes/sec]');

subplot(2,4,2);
plot(ILDs,RToutEST,'c-');
title('Exponential Stein Model');
xlim([-45,15]);
ylim([0,180]);

subplot(2,4,6)
plot(ILDs,RToutAST,'m-');
title('Alpha Stein Model');
xlim([-45,15]);
ylim([0,180]);

subplot(2,4,3);
plot(ILDs,RToutPIF,'g-');
title('Passive IF Model');
xlim([-45,15]);
ylim([0,180]);

subplot(2,4,7);
plot(ILDs,RToutAIF,'k-');
title('Active IF Model');
xlim([-45,15]);
ylim([0,180]);

subplot(2,4,4);
plot(ILDs,RToutOWC,'y-');
title('Original Wang-Colburn Model');
xlim([-45,15]);
ylim([0,180]);

subplot(2,4,8);
plot(ILDs,RToutAWC,'r-');
title('Adjusted Wang-Colburn Model');
xlim([-45,15]);
ylim([0,180]);
