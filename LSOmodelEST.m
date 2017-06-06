function [spOut, vOut] = LSOmodelEST(spEx, spIn, DT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exponential Stein Model of Lateral Superior Olive 
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
%  The prototype of this model was created by Colburn & Moss (1981), 
%  based on the abstract neuron model of Stein (1965,1967). Each synaptic 
%  input to the model neuron is converted to an exponential function to 
%  form a 'virtual membrane potential'. An output spike is counted when 
%  the virtual potential reaches or exceeds the threshold. For more 
%  detailed descriptions, see Ashida et al. (2017). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References 
% Colburn HS, Moss PJ (1981) "Binaural interaction models and mechanisms" 
%  In: Neuronal mechanisms in hearing. New York: Plenum Press. 
% Stein RB (1965,1967) Biophys J 5:173-194. / Biophys J 7: 37-68.
%  "A theoretical analysis of neuronal variability" 
%  "Some models of neuronal variability" 
% Ashida G, Tollin DJ, Kretzberg J (2017) Submitted 
%  "Physiological models of the lateral superior olive" 
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
Th = 5.5;   % threshold 
Hex = 1.0;  % excitation strength (fixed to one)
Tex = 0.70; % [ms] excitation time constant 
Hin = 1.8;  % inhibition strength 
Tin = 0.98; % [ms] inhibition time constant 

Nref = round(Tref/DT); % steps
Dex = exp(-DT/Tex); % decay factor (at each step) for excitatory inputs
Din = exp(-DT/Tin); % decay factor (at each step) for inhibitory inputs

%% data vectors
Nsteps = min([length(spEx), length(spIn)]);
spOut = zeros(1,Nsteps);
vOut = zeros(1,Nsteps); % virtual membrane potential 
vEx = 0; % virtual excitatory synaptic input 
vIn = 0; % virtual inhibitory synaptic input 
refCount = 0; % refractory counter

%% main loop 
for t = 1:Nsteps

  % add synaptic inputs at time t 
  if(refCount>0) % if in refractory period, then ignore inputs
   vEx = 0; 
   vIn = 0; 
  else 
   vEx = vEx * Dex + Hex * spEx(t);
   vIn = vIn * Din + Hin * spIn(t); 
  end  

  % get the membrane potential at time t and check for threshold crossing 
  vm = vEx - vIn; 
  % (1) if in refractory period, then decrement counter 
  if(refCount>0); 
   spOut(t)=0; refCount = refCount-1;  
  % (2) if threshold is reached, generate a spike, reset synaptic inputs, 
  %     and set the refractory counter
  elseif(vm>=Th) 
   spOut(t)=1; vEx = 0; vIn = 0;  refCount = Nref; 
  % (3) if no threshold crossing happened, then no spike output 
  else 
   spOut(t)=0; 
  end 

  % save internal state (virtual membrane potential) at time t 
  vOut(t) = vm;

end
