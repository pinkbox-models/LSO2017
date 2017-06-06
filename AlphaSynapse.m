function Aout = AlphaSynapse(Sp,DT,Amp,Tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha function for simulating synaptic input waveforms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%  Sp  : spike vector (to show the number of inputs at each step)
%  DT  : size of time step [ms]
%  Amp : input amplitude [arbitrary unit]
%  Tau : time constant [ms]
%
% Output
%  Aout : simulated synaptic input waveform (same unit as Amp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Technical notes
%  This script calculates the alpha function at discrete time steps. 
%  Numerical algorithm based on Rotter & Diesmann (1999, Biol Cybern). 
%  Since the alpha function here is calculated by numerically solving 
%  a second order differential equation (Eq(1) below), we use two 
%  variables, x and y. The resulting alpha function is stored in x, 
%  while y serves as an auxiliary function. Namely, 
%   x: alpha function that satisfies 
%     (Tau^2)(d^2x/dt^2) + (2Tau)(dx/dt) + x = 0 ---Eq(1)
%   y: auxiliary function defined as y = x/Tau + dx/dt
% With the notation above, Eq(1) is equivalent to the system: 
%  dy/dt = - y/Tau     ( exponential decay function )
%  dx/dt = - x/Tau + y ( alpha function )
% Therefore, at each step, we calculate the following: 
%  ynew = exp(-DT/Tau) * yold (+ input)
%  xnew = exp(-DT/Tau) * xold + DT * ynew  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revisions
% Created (ver 0.9): May 16, 2017 by GA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creating output vector
Aout = zeros(1,length(Sp));

% initial values 
x = 0;
y = 0; 

% factors used for calculation
aa = Amp * exp(1.0) / Tau; % amplitude factor for normalizing the input 
ee = exp(-DT/Tau); % decrease factor used at each time step

% step-by-step calculation 
for t=1:length(Sp)
 y = ee*y + aa*Sp(t); 
 x = ee*x + DT*y;
 Aout(t) = x; % store data
end
