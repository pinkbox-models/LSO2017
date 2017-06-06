%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testMonaural.m 
% --- for simulating monaural AM coding of LSO models  
% Created: May 16, 2017 by Go Ashida
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
tvmain = tv(lmain);

%% frequency parameters 
FQs = [50, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, ...
       500, 600, 700, 800, 1000, 1200, 1400];
RTs = 180-0.03*FQs; % frequency dependent ipsilateral rate 
VSs = 0.65 * (1-exp( (FQs-2000)/500) ) ./ (1+exp( (FQs-2000)/500) ); % VS
P0 = 0; % zero initial phase 
spontRT = 30.0; % spontaneous rate for contralateral (inhibitory) input 
spontVS = 0; % no phase-locking for contralateral (inhibitory) input 
Mex = 20; % number of excitatory inputs 
Min = 8;  % number of inhibitory inputs 

%% data array for output rates and vector strengths
RToutCOC = zeros(1,length(FQs));
RToutEST = zeros(1,length(FQs));
RToutAST = zeros(1,length(FQs));
RToutPIF = zeros(1,length(FQs));
RToutAIF = zeros(1,length(FQs));
RToutOWC = zeros(1,length(FQs));
RToutAWC = zeros(1,length(FQs));

VSoutCOC = zeros(1,length(FQs));
VSoutEST = zeros(1,length(FQs));
VSoutAST = zeros(1,length(FQs));
VSoutPIF = zeros(1,length(FQs));
VSoutAIF = zeros(1,length(FQs));
VSoutOWC = zeros(1,length(FQs));
VSoutAWC = zeros(1,length(FQs));

%% main loop

for i = 1:length(FQs)
 FQ = FQs(i) 
 RT = RTs(i); 
 VS = VSs(i); 

 % spike input vectors
 spEx = sum( PhaseLock(Mex,length(tv),FQ,VS,RT,P0,DT), 1 );
 spIn = sum( PhaseLock(Min,length(tv),FQ,spontVS,spontRT,P0,DT), 1 );

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

 % calculating output vector strengths 
 data = dataCOC;
 cc = sum( cos( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 ss = sum( sin( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 VSoutCOC(i) = sqrt(cc*cc + ss*ss)/sum(data);

 data = dataEST;
 cc = sum( cos( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 ss = sum( sin( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 VSoutEST(i) = sqrt(cc*cc + ss*ss)/sum(data);

 data = dataAST;
 cc = sum( cos( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 ss = sum( sin( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 VSoutAST(i) = sqrt(cc*cc + ss*ss)/sum(data);

 data = dataPIF;
 cc = sum( cos( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 ss = sum( sin( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 VSoutPIF(i) = sqrt(cc*cc + ss*ss)/sum(data);
 
 data = dataAIF;
 cc = sum( cos( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 ss = sum( sin( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 VSoutAIF(i) = sqrt(cc*cc + ss*ss)/sum(data);

 data = dataOWC;
 cc = sum( cos( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 ss = sum( sin( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 VSoutOWC(i) = sqrt(cc*cc + ss*ss)/sum(data);

 data = dataAWC;
 cc = sum( cos( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 ss = sum( sin( 2*pi*FQ/1000*tvmain(logical(data)) ) ); 
 VSoutAWC(i) = sqrt(cc*cc + ss*ss)/sum(data);
 
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
figure(110); 
cla; 

subplot(2,4,1);
plot(FQs,RToutCOC,'b-');
title('Coincidence Counting Model');
xlim([0,1200]);
ylim([0,180]);
xlabel('AM frequency [Hz]');
ylabel('rate [spikes/sec]');

subplot(2,4,2);
plot(FQs,RToutEST,'c-');
title('Exponential Stein Model');
xlim([0,1200]);
ylim([0,180]);

subplot(2,4,6);
plot(FQs,RToutAST,'m-');
title('Alpha Stein Model');
xlim([0,1200]);
ylim([0,180]);

subplot(2,4,3);
plot(FQs,RToutPIF,'g-');
title('Passive IF Model');
xlim([0,1200]);
ylim([0,180]);

subplot(2,4,7);
plot(FQs,RToutAIF,'k-');
title('Active IF Model');
xlim([0,1200]);
ylim([0,180]);

subplot(2,4,4);
plot(FQs,RToutOWC,'y-');
title('Original Wang-Colburn Model');
xlim([0,1200]);
ylim([0,180]);

subplot(2,4,8);
plot(FQs,RToutAWC,'r-');
title('Adjusted Wang-Colburn Model');
xlim([0,1200]);
ylim([0,180]);

figure(111); 
cla; 

subplot(2,4,1);
semilogx(FQs,20*log10(2*VSoutCOC),'b-'); 
set(gca,'xscale','log');
xlim([50,1500]);
ylim([-6,6]);
xlabel('AM frequency [Hz]');
ylabel('synch gain [dB]');

subplot(2,4,2);
semilogx(FQs,20*log10(2*VSoutEST),'c-'); 
set(gca,'xscale','log');
xlim([50,1500]);
ylim([-6,6]);

subplot(2,4,6);
semilogx(FQs,20*log10(2*VSoutAST),'m-'); 
set(gca,'xscale','log');
xlim([50,1500]);
ylim([-6,6]);

subplot(2,4,3);
semilogx(FQs,20*log10(2*VSoutPIF),'g-'); 
set(gca,'xscale','log');
xlim([50,1500]);
ylim([-6,6]);

subplot(2,4,7);
semilogx(FQs,20*log10(2*VSoutAIF),'k-'); 
set(gca,'xscale','log');
xlim([50,1500]);
ylim([-6,6]);

subplot(2,4,4);
semilogx(FQs,20*log10(2*VSoutOWC),'y-'); 
set(gca,'xscale','log');
xlim([50,1500]);
ylim([-6,6]);

subplot(2,4,8);
semilogx(FQs,20*log10(2*VSoutAWC),'r-'); 
set(gca,'xscale','log');
xlim([50,1500]);
ylim([-6,6]);


