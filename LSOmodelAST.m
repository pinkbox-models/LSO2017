function [spOut, vOut] = LSOmodelAST(spEx, spIn, DT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha Stein Model of Lateral Superior Olive 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%  spEx : excitatory input spike vector (number of inputs at each time step)
%  spIn : inhibitory input spike vector (number of inputs at each time step)
%  DT   : size of time step [ms] 
%
% Output
%  spOut : number of output spikes at each time step 
%  vOut  : internal state (virtual membrane potential) at each time step 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes on the input/output arguments 
%  - spEx and spIn should have the same length. Otherwise, the shorter 
%    of them determines the total simulation time length. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes on the model 
%  The alpha Stein model is a revised version of the exponential Stein 
%  model (LSOmodelEST.m) by replacing the exponential decay function with 
%  an alpha function. Namely, each synaptic input to the model neuron is 
%  converted to an alpha function to form a 'virtual membrane potential'. 
%  An output spike is counted when the potential reaches or exceeds the 
%  threshold. For more detailed descriptions, see Ashida et al. (2017). 
%  Calculation of the alpha function is based on the numerical algorithm 
%  by Rotter & Diesmann (1999, Biol Cybern). See also AlphaSynapse.m. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References 
% Ashida G, Tollin DJ, Kretzberg J (2017) Submitted 
%  "Physiological models of the lateral superior olive" 
% Rotter S, Diesmann M (1999) Biol Cybern 81: 381-402
%  "Exact digital simulation of time-invariant linear systems with 
%   applications to neuronal modeling" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revisions
% Created (ver 0.9): May 16, 2017 by Go Ashida
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

%% pre-defined parameters
Tref = 1.6; % [ms] 
Th = 7.3;   % threshold 
Hex = 1.0;  % excitation strength (fixed to one)
Tex = 0.45; % [ms] excitation time constant 
Hin = 1.7;  % inhibition strength 
Tin = 0.63; % [ms] inhibition time constant 

Nref = round(Tref/DT); % steps
Dex = exp(-DT/Tex); % decay factor (at each step) for excitatory inputs
Din = exp(-DT/Tin); % decay factor (at each step) for inhibitory inputs
Aex = Hex * exp(1.0) / Tex; % amplitude conversion factor 
Ain = Hin * exp(1.0) / Tin; % amplitude conversion factor 

%% data vectors
Nsteps = min([length(spEx), length(spIn)]);
spOut = zeros(1,Nsteps);
vOut = zeros(1,Nsteps); % virtual membrane potential 
vEx = 0; % virtual excitatory synaptic input 
vIn = 0; % virtual inhibitory synaptic input 
wEx = 0; % auxiliary (exponential) function for calculating alpha function  
wIn = 0; % auxiliary (exponential0 function for calculating alpha function 
refCount = 0; % refractory counter

%%%%% tech notes for calculating alpha function %%%%%%%%%%%%%%%%%%%%%%%%%%
% wnew = exp(-DT/Tau) * wold (+ input) ..... (exponential decay)
% vnew = exp(-DT/Tau) * vold + DT * wnew ... (alpha function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% main loop 
for t = 1:Nsteps

  % add synaptic inputs at time t 
  if(refCount>0) % if in refractory period, then ignore inputs
   wEx = 0; vEx = 0; 
   wIn = 0; vIn = 0;  
  else 
   wEx = wEx * Dex + Aex * spEx(t);
   vEx = vEx * Dex + wEx * DT; 
   wIn = wIn * Din + Ain * spIn(t);
   vIn = vIn * Din + wIn * DT; 
  end  

  % get the membrane potential at time t and check for threshold crossing 
  vm = vEx - vIn; 
  % (1) if in refractory period, then decrement counter 
  if(refCount>0); 
   spOut(t)=0; refCount = refCount-1;  
  % (2) if threshold is reached, generate a spike, reset synaptic inputs, 
  %     and set the refractory counter
  elseif(vm>=Th) 
   spOut(t)=1; wEx = 0; vEx = 0; wIn = 0; vIn = 0; refCount = Nref; 
  % (3) if no threshold crossing happened, then no spike output 
  else 
   spOut(t)=0; 
  end 

  % save internal state (virtual membrane potential) at time t 
  vOut(t) = vm;

end
