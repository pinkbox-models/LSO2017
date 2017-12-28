function [spOut, vOut] = LSOmodelAWC(spEx, spIn, DT, Iext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjusted Wang-Colburn Model of Lateral Superior Olive 
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
%  The model is based on the model presented in Wang & Colburn (2012), 
%  but with modifications of the parameters to better replicate known 
%  physiological results. Modifications include: shifting the voltage-
%  dependence of channel kinetics by +5 mV, and revising the parameter 
%  values (capacitance, conductances, reversal potentials, etc.). 
%  More detailed descriptions and analyses of the model can be found 
%  in Ashida et al. (2017). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References 
% Wang L, Colburn HS (2012) JARO 13: 249-267
%  "A modeling study of the responses of the lateral superior olive 
%   to ipsilateral sinusoidally amplitude-modulated tones" 
% Ashida G, Tollin DJ, Kretzberg J (2017) PLoS Comput Biol 13(12): e1005903. 
% "Physiological models of the lateral superior olive" 
% https://doi.org/10.1371/journal.pcbi.1005903
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revisions
% Created (ver 0.9): May 16, 2017 by Go Ashida
% Updated references info (ver. 1.0): Dec. 28, 2017 by GA 
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
gLL = 24.0e-6; % 24.0*1e-6[mS] = 24.0[nS] 
gKL = 15.0e-6; % 15.0*1e-6[mS] = 15.0[nS]
gKH = 440.0e-6; % 440*1e-6[mS] = 440[nS] 
gNa = 4400.0e-6; % 4400*1e-6[mS] = 4400[nS] 
EL = -60.0; % [mV] leak reversal potential 
EK = -75.0; % [mV] potassium reversal potential 
EN = +50.0; % [mV] sodium reversal potential

% synaptic input parameters
Aex =  3.5; % [nS] excitatory synaptic input amplitude 
Ain = 12.0; % [nS] inhibitory synaptic input amplitude 
Tex = 0.16; % [ms] excitatory synaptic input time constant 
Tin = 0.32; % [ms] inhibitory synaptic input time constant 
Eex =  0.0; % [mV] excitatory synaptic input reversal potential
Ein =-75.0; % [mV] inhibitory synaptic input reversal potential  

% temperature factors
Temp_body = 37.0; 
Temp_KLVA = 22.0; 
Temp_KHVA = 22.0;
Temp_Na   = 22.0;
Q10_KL = 3.0; 
Q10_KH = 3.0; 
Q10_Na = 3.0; 
T_KL = ( Q10_KL ^ ((Temp_body-Temp_KLVA)/10.0) ); 
T_KH = ( Q10_KH ^ ((Temp_body-Temp_KHVA)/10.0) ); 
T_Na = ( Q10_Na ^ ((Temp_body-Temp_Na)  /10.0) ); 

% initialization
vRest = -60.3; 
vOut(1) = vRest; 
w = winf(vRest); 
z = zinf(vRest); 
n = ninf(vRest); 
p = pinf(vRest); 
m = minf(vRest); 
h = hinf(vRest); 

% spike detection thresholds
vTh1 = -30.0; % [mV] spike initiation 
vTh2 = -45.0; % [mV] spike end
spikeflag = 0;  % flag for spiking

%% synaptic inputs 
gEx = AlphaSynapse(spEx, DT, Aex, Tex); % [nS] 
gIn = AlphaSynapse(spIn, DT, Ain, Tin); % [nS] 

%% main loop
for t=1:Nsteps-1

 % calculating all values at time t
 u = vOut(t); % potential 
 iLL = gLL * (EL-u); % leak current
 iKL = gKL * (EK-u) * ( w^4 * z ); % KLVA current
 iKH = gKH * (EK-u) * ( 0.85 * n * n + 0.15 * p ); % KHVA current 
 iNa = gNa * (EN-u) * ( m^3 * h ); % Na current
 iEx = gEx(t) * (Eex-u) * 1e-6; % [uA] excitatory synaptic current 
 iIn = gIn(t) * (Ein-u) * 1e-6; % [uA] inhibitory synaptic current
 itot = iLL + iKL + iKH + iNa + iEx + iIn + jext(t); % total current 

 % step forward 
 v = vOut(t) + (itot/Cs) * DT; 
 w = w + T_KL * ( winf(u) - w )/ tauw(u) * DT;
 z = z + T_KL * ( zinf(u) - z )/ tauz(u) * DT;
 n = n + T_KH * ( ninf(u) - n )/ taun(u) * DT;
 p = p + T_KH * ( pinf(u) - p )/ taup(u) * DT;
 m = m + T_Na * ( minf(u) - m )/ taum(u) * DT;
 h = h + T_Na * ( hinf(u) - h )/ tauh(u) * DT;

 % check for spiking 
 if(spikeflag) % if spiking, then look for end of spike 
   spOut(t+1) = 0; 
   if(v<=vTh2)  
    spikeflag = 0; % reset the flag 
   end 
 else % if not spiking, then look for spike threshold crossing 
   if(v>=vTh1) 
    spOut(t+1) = 1; % report spike generation 
    spikeflag = 1; % set the flag 
   else
    spOut(t+1) = 0; 
   end
 end

 % save membrane potential 
 vOut(t+1) = v; 

end % end of main loop

%% end of main function
end % corresponding to the first line of this function declaration

%% internal functions for channel kinetics 
function x = winf(u) % KLVA activation : steady state
 v = u-5.0;
 x = 1.0 / sqrt(sqrt( 1.0 + exp( -(v+48.0)/6.0 ) )); 
end
function x = tauw(u) % KLVA activation : time constant 
 v = u-5.0;
 x = 1.5 + 100.0 / ( 6.0*exp( (v+60.0)/6.0 ) + 16.0*exp( -(v+60.0)/45.0 ) );
end
function x = zinf(u) % KLVA inactivation : steady state
 v = u-5.0;
 x = 0.5 + 0.5 / ( 1.0 + exp( (v+71.0)/10.0 ) ); 
end
function x = tauz(u) % KLVA inactivation : time constant 
 v = u-5.0;
 x = 50.0 + 1000.0 / ( exp( (v+60.0)/20.0 ) + exp( -(v+60.0)/8.0 ) );
end
function x = ninf(u) % KHVA activation 1 : steady state
 v = u-5.0;
 x = 1.0 / sqrt( 1.0 + exp( -(v+15.0)/5.0 ) ); 
end
function x = taun(u) % KHVA activation 1 : time constant 
 v = u-5.0;
 x = 0.7 + 100.0 / ( 11.0*exp( (v+60.0)/24.0 ) + 21.0*exp( -(v+60.0)/23.0 ) );
end
function x = pinf(u) % KHVA activation 2 : steady state
 v = u-5.0;
 x = 1.0 / ( 1.0 + exp( -(v+23.0)/6.0 ) ); 
end
function x = taup(u) % KHVA activation 2 : time constant 
 v = u-5.0;
 x = 5.0 + 100.0 / ( 4.0*exp( (v+60.0)/32.0 ) + 5.0*exp( -(v+60.0)/22.0 ) );
end
function x = minf(u) % Na activation : steady state
 v = u-5.0;
 x = 1.0 / ( 1.0 + exp( -(v+38.0)/7.0 ) ); 
end
function x = taum(u) % Na activation : time constant 
 v = u-5.0;
 x = 0.04 + 10.0 / ( 5.0*exp( (v+60.0)/18.0 ) + 36.0*exp( -(v+60.0)/25.0) );
end
function x = hinf(u) % Na inactivation : steady state
 v = u-5.0;
 x = 1.0 / ( 1.0 + exp( (v+65.0)/6.0 ) ); 
end 
function x = tauh(u) % Na inactivation : time constant 
 v = u-5.0;
 x = 0.6 + 100.0 / ( 7.0*exp( (v+60.0)/11.0 ) + 10.0*exp( -(v+60.0)/25.0) );
end
