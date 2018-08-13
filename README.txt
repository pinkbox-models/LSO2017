----------------------------------------------------------------------------------
 Physiological Models of Lateral Superior Olive Neurons -- Matlab Implementation 
----------------------------------------------------------------------------------

%%% Versions %%% 

+ Ver. 0.9 (Jun. 06, 2017): Initial release of the code on GitHub. 
+ Ver. 1.0 (Dec. 28, 2017): Updated reference info. 

%%% Author %%% 

Go Ashida (University of Oldenburg) go.ashida@uni-oldenburg.de


%%% Contents %%% 

LSOmodelCOC : Coincidence counting model 
LSOmodelEST : Exponential Stein model 
LSOmodelAST : Alpha Stein model 
LSOmodelPIF : Passive integrate-and-fire model 
LSOmodelAIF : Active integrate-and-fire model 
LSOmodelOWC : Original Wang-Colburn model 
LSOmodelAWC : Adjusted Wang-Colburn model 

PhaseLock    : Code for generating phase-locked input sequences 
AlphaSynapse : Code for calculating synaptic input modeled as an alpha function 

testSynapticInput : Sample code for plotting simulated synaptic inputs of each model 
testStepCurrent   : Sample code for plotting step current response of IF and WC models 
testMonaural      : Sample code for plotting monaural rate-MTF and synch-MTF of each model 
testIPDcoding     : Sample code for plotting binaural phase-tuning curve of each model 
testILDcoding     : Sample code for plotting binaural intensity-tuning curve of each model 

+ Notes: See each program file and the reference below for more detailed descriptions. 

%%% Reference %%% 

Ashida G, Tollin DJ, Kretzberg J (2017) 
"Physiological models of the lateral superior olive" 
PLoS Comput Biol 13(12): e1005903. 
https://doi.org/10.1371/journal.pcbi.1005903

