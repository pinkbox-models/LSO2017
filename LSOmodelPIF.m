function [spOut, vOut] = LSOmodelPIF(spEx, spIn, DT, Iext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passive Integrate-and-fire Model of Lateral Superior Olive 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%  spEx : excitatory input spike vector (number of inputs at each time step)
%  spIn : inhibitory input spike vector (number of inputs at each time step)
%  DT   : size of time step [ms] 
%  Iext (optional) : external current [pA] (will be converted into uA) 
%
% Output
%  spOut : number of output spikes at each time step 
%  vOut  : membrane potential at each time step 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes on the input/output arguments 
% + Iext is an optional input argument to specify (time-varying) external 
%   input current. If this argument is not given (or given as an empty 
%   vector), then the external current is assumed to be zero. 
% + If the external current Iext is a scalar, then the model is assumed 
%   to receive a constant current of size Iext [pA]. 
% + Total simulation time length is determined as below: 
%  - If Iext is not given, given as an empty vector, or given as a scalar, 
%    then the shorter of spEx and spIn determines the length. 
%  - If Iext is a non-empty vectors, then the shortest of spEx, spIn, and 
%    Iext determines the length. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Units to use in the code 
% (all parameters have to be converted into these units before use)
%  T [ms] 
%  I [uA] 
%  V [mV] 
%  g [mS] 
%  C [uF] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes on the model 
%  The passive integrate-and-fire (IF) model of LSO is a (variation of 
%  the) standard leaky IF model with parameters tuned to reproduce 
%  known physiological response properties of LSO neurons. The synaptic 
%  conductance induced is described as an alpha function and fed into the 
%  model. A spike is counted when the model membrane potential reaches 
%  or exceeds the threshold. After each threshold crossing, the potential 
%  is reset and kept there for the refractory period. For more detailed 
%  descriptions, see Ashida et al. (2017). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References 
% Ashida G, Tollin DJ, Kretzberg J (2017) PLoS Comput Biol 13(12): e1005903. 
% "Physiological models of the lateral superior olive" 
% https://doi.org/10.1371/journal.pcbi.1005903
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revisions
% Created (ver 0.9): May 16, 2017 by Go Ashida
% Updated reference info (ver 1.0): Dec 28, 2017 by GA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you find a bug, please report to GA at go.ashida@uni-oldenburg.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Copyright 2017 Go Ashida (go.ashida@uni-oldenburg.de) %%%%%%%%%%%%%
% Permission is hereby granted under the Apache License, Version 2.0; 
% Users of this file must be in compliance with this license, a copy of 
% which may be obtained at http://www.apache.org/licenses/LICENSE-2.0
% This file is provided on an "AS IS" basis, WITHOUT WARRANTIES OR 
% CONDITIONS OF ANY KIND, either express or implied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% simulation parameters
% If Iext is not given, then create it as an empty vector
if ~exist('Iext','var')
 Iext = []; 
end

% Determine the simulation time length 
if length(Iext)>1
 Nsteps = min([length(spEx), length(spIn), length(Iext)]);
else
 Nsteps = min([length(spEx), length(spIn)]);
end

% Make external input vector depending on the input argument 
if length(Iext)>1 % if the Iext is a vector, convert it from [pA] to [uA]
 jext = Iext * 1e-6; 
elseif length(Iext)==1 % if Iext is a scalar, then set a constant input 
 jext = ones(1,Nsteps) * Iext(1) * 1e-6; 
else % if Iext is an empty vector, then use zero input
 jext = zeros(1,Nsteps); 
end

% Create output vectors 
spOut = zeros(1,Nsteps);
vOut = zeros(1,Nsteps);

%% model parameters

% membrane parameters
Cs  = 24.0e-6; % 24.0*1e-6[uF] = 24.0[pF] 
gLL = 26.4e-6; % 26.4*1e-6[mS] = 26.4[nS] 
EL = -60.0; % [mV] leak reversal potential 
Vth = -45.1; % [mV] threshold 
Vreset = -60.0; % [mV] reset potential

% refractory period 
Tref = 1.6; % [ms] 
Nref = round(Tref/DT); % steps
refCount = 0; % refractory counter

% synaptic input parameters
Aex =  3.5; % [nS] excitatory synaptic input amplitude 
Ain = 12.0; % [nS] inhibitory synaptic input amplitude 
Tex = 0.16; % [ms] excitatory synaptic input time constant 
Tin = 0.32; % [ms] inhibitory synaptic input time constant 
Eex =  0.0; % [mV] excitatory synaptic input reversal potential
Ein =-75.0; % [mV] inhibitory synaptic input reversal potential  

% initialization
vOut(1) = Vreset; 

%% synaptic inputs 
gEx = AlphaSynapse(spEx, DT, Aex, Tex); % [nS] 
gIn = AlphaSynapse(spIn, DT, Ain, Tin); % [nS] 

%% main loop 
for t = 1:Nsteps-1

 % calculating all values at time t
 u = vOut(t); % potential 
 iLL = gLL * (EL-u); % leak current
 iEx = gEx(t) * (Eex-u) * 1e-6; % [uA] excitatory synaptic current 
 iIn = gIn(t) * (Ein-u) * 1e-6; % [uA] inhibitory synaptic current
 itot = iLL + iEx + iIn + jext(t); % total current 
  
 % step forward 
 v = vOut(t) + (itot/Cs) * DT; 

 % check for threshold crossing 
 % (1) if in refractory period, then decrement counter and bound the  
 %     membrane potential at the reset potential 
 if(refCount>0); 
  spOut(t) = 0; v = Vreset; refCount = refCount-1;  
 % (2) if threshold is reached, count a spike and set the refractory counter
 elseif(v>=Vth) 
  spOut(t) = 1; refCount = Nref; 
 % (3) if no threshold crossing happened, then no spike output 
 else 
  spOut(t) = 0; 
 end 

 % save membrane potential 
 vOut(t+1) = v; 

end
