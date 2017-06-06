%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testSynapticInput.m 
% --- for simulating synaptic potentials 
% Created: May 16, 2017 by Go Ashida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% time parameters
DT = 0.002; % [ms] time step
T0 = 5.0; % [ms] time before stimulus
T1 = 15.0; % [ms] time after stimulus
N0 = round(T0/DT); % steps 
N1 = round(T1/DT); % steps 
Ntot = N0 + N1; 
tv = (0:Ntot)*DT; % time vector [ms]; caution:length=Ntot+1

%% stimulus vector
spEx = zeros(1, Ntot+1);
spIn = zeros(1, Ntot+1);

%% data vector
Nmax = 8; % max number of inputs 
vCOCe = zeros(Nmax+1,Ntot+1); 
vCOCi = zeros(Nmax+1,Ntot+1); 
vESTe = zeros(Nmax+1,Ntot+1); 
vESTi = zeros(Nmax+1,Ntot+1); 
vASTe = zeros(Nmax+1,Ntot+1); 
vASTi = zeros(Nmax+1,Ntot+1); 
vPIFe = zeros(Nmax+1,Ntot+1); 
vPIFi = zeros(Nmax+1,Ntot+1); 
vAIFe = zeros(Nmax+1,Ntot+1); 
vAIFi = zeros(Nmax+1,Ntot+1); 
vOWCe = zeros(Nmax+1,Ntot+1); 
vOWCi = zeros(Nmax+1,Ntot+1); 
vAWCe = zeros(Nmax+1,Ntot+1); 
vAWCi = zeros(Nmax+1,Ntot+1); 
spOut = zeros(1,Ntot+1); % dummy spike vector

%% main loop 

for i=0:Nmax
  % excitatory input
  spEx(N0) = i; 
  spIn(N0) = 0;
  [spOut, vCOCe(i+1,:)] = LSOmodelCOC(spEx, spIn, DT);
  [spOut, vESTe(i+1,:)] = LSOmodelEST(spEx, spIn, DT);
  [spOut, vASTe(i+1,:)] = LSOmodelAST(spEx, spIn, DT);
  [spOut, vPIFe(i+1,:)] = LSOmodelPIF(spEx, spIn, DT);
  [spOut, vAIFe(i+1,:)] = LSOmodelAIF(spEx, spIn, DT);
  [spOut, vOWCe(i+1,:)] = LSOmodelOWC(spEx, spIn, DT);
  [spOut, vAWCe(i+1,:)] = LSOmodelAWC(spEx, spIn, DT);

  % inhibitory input
  spEx(N0) = 0; 
  spIn(N0) = i;
  [spOut, vCOCi(i+1,:)] = LSOmodelCOC(spEx, spIn, DT);
  [spOut, vESTi(i+1,:)] = LSOmodelEST(spEx, spIn, DT);
  [spOut, vASTi(i+1,:)] = LSOmodelAST(spEx, spIn, DT);
  [spOut, vPIFi(i+1,:)] = LSOmodelPIF(spEx, spIn, DT);
  [spOut, vAIFi(i+1,:)] = LSOmodelAIF(spEx, spIn, DT);
  [spOut, vOWCi(i+1,:)] = LSOmodelOWC(spEx, spIn, DT);
  [spOut, vAWCi(i+1,:)] = LSOmodelAWC(spEx, spIn, DT);
  
end


%% plotting
figure(102);
cla; 

for i=0:Nmax
  subplot(4,4,1); hold on; 
  plot(tv,vCOCe(i+1,:),'b'); 
  xlim([4,14]);
  title('Coincidence Counting Model');
  subplot(4,4,5); hold on; 
  plot(tv,vCOCi(i+1,:),'r'); 
  xlim([4,14]);
  xlabel('time [ms]');
  
  subplot(4,4,2); hold on; 
  plot(tv,vESTe(i+1,:),'b'); 
  xlim([4,14]);
  title('Exponential Stein Model');
  subplot(4,4,6); hold on; 
  plot(tv,vESTi(i+1,:),'r'); 
  xlim([4,14]);

  subplot(4,4,10); hold on; 
  plot(tv,vASTe(i+1,:),'b'); 
  xlim([4,14]);
  title('Alpha Stein Model');
  subplot(4,4,14); hold on; 
  plot(tv,vASTi(i+1,:),'r'); 
  xlim([4,14]);
  
  subplot(4,4,3); hold on; 
  plot(tv,vPIFe(i+1,:),'b'); 
  xlim([4,14]);
  title('Passive IF Model');
  subplot(4,4,7); hold on; 
  plot(tv,vPIFi(i+1,:),'r'); 
  xlim([4,14]);
  
  subplot(4,4,11); hold on; 
  plot(tv,vAIFe(i+1,:),'b'); 
  xlim([4,14]);
  title('Active IF Model');
  subplot(4,4,15); hold on; 
  plot(tv,vAIFi(i+1,:),'r'); 
  xlim([4,14]);
  
  subplot(4,4,4); hold on; 
  plot(tv,vOWCe(i+1,:),'b'); 
  xlim([4,14]);
  title('Original Wang-Colburn Model');
  subplot(4,4,8); hold on; 
  plot(tv,vOWCi(i+1,:),'r'); 
  xlim([4,14]);
  
  subplot(4,4,12); hold on; 
  plot(tv,vAWCe(i+1,:),'b'); 
  xlim([4,14]);
  title('Adjusted Wang-Colburn Model');
  subplot(4,4,16); hold on; 
  plot(tv,vAWCi(i+1,:),'r'); 
  xlim([4,14]);
  
end


