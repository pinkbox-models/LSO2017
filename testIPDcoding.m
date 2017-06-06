%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testIPDcoding.m 
% --- for simulating IPD/ITD coding of LSO models  
% Created: Apr 27, 2017 by Go Ashida
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

%% input parameters
FQ = 300; 
RT = 180-0.03*FQ; % frequency dependent input rate 
VS = 0.65 * (1-exp( (FQ-2000)/500) ) ./ (1+exp( (FQ-2000)/500) ); % VS
P0 = 0; % zero initial phase 
Mex = 20; % number of excitatory inputs 
Min = 8;  % number of inhibitory inputs 
phases = -180:10:180; % [deg] phase difference between ipsi and contra inputs 

%% data array for output rates
RToutCOC = zeros(1,length(phases));
RToutEST = zeros(1,length(phases));
RToutAST = zeros(1,length(phases));
RToutPIF = zeros(1,length(phases));
RToutAIF = zeros(1,length(phases));
RToutOWC = zeros(1,length(phases));
RToutAWC = zeros(1,length(phases));

%% spike input vectors
spEx  = sum( PhaseLock(Mex,length(tv),FQ,VS,RT,P0,DT), 1 );
spTmp = sum( PhaseLock(Min,length(tv),FQ,VS,RT,P0,DT), 1 );

%% main loop 

for i = 1:length(phases)
 
 % time lag between excitatory and inhibitory inputs 
 pp = phases(i) % [deg] 
 Tlag = 1000.0/FQ*pp/360; % [ms]
 Nlag = round(Tlag/DT); 

 % assign inhibitory vector according to phase difference  
 if(Nlag<=0); 
   spIn = [ zeros(1,-Nlag), spTmp(1:end+Nlag) ];    
 else  
   spIn = [ spTmp(Nlag:end), zeros(1,Nlag) ]; 
 end

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
figure(112); 
cla; 

subplot(2,4,1);
plot(phases,RToutCOC,'b-');
title('Coincidence Counting Model');
xlim([-180,180]);
ylim([0,180]);
xlabel('phase difference [deg]');
ylabel('rate [spikes/sec]');

subplot(2,4,2);
plot(phases,RToutEST,'c-');
title('Exponential Stein Model');
xlim([-180,180]);
ylim([0,180]);

subplot(2,4,6)
plot(phases,RToutAST,'m-');
title('Alpha Stein Model');
xlim([-180,180]);
ylim([0,180]);

subplot(2,4,3);
plot(phases,RToutPIF,'g-');
title('Passive IF Model');
xlim([-180,180]);
ylim([0,180]);

subplot(2,4,7);
plot(phases,RToutAIF,'k-');
title('Active IF Model');
xlim([-180,180]);
ylim([0,180]);

subplot(2,4,4);
plot(phases,RToutOWC,'y-');
title('Original Wang-Colburn Model');
xlim([-180,180]);
ylim([0,180]);

subplot(2,4,8);
plot(phases,RToutAWC,'r-');
title('Adjusted Wang-Colburn Model');
xlim([-180,180]);
ylim([0,180]);

