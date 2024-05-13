function [SNReff,bEff] = DWSSFP_SNREfficiencyOscillating(x,optSample,optScanner,bOsc,tau,~)
%%
% calculates DW-SSFP SNR Efficiency
%
% input:
%	x   = array with difference between TR & gradient duration (incorporating deadtime) and the flip angle (ms)
%	optSample   = T1 (ms), T2 (ms), T2* (ms) & D (um^2/ms)
%   optScanner = Maximum readout duration (ms) & Minimum time between end of gradient & Readout (ms)
%	bOsc   = b-value associated with a single diffusion gradient (used for estimation of relative effective b-value) (ms/um2)
%   tau   = Encoding duration
%
% output:
%	SNReff = Estimated SNR Efficiency
%	bEff   = Estimated b-value
%%
%Define array parameters - As TR has to be greater than tau, optimise for the difference
TR=tau+x(1)+optScanner.DeadTime;
alpha=x(2);
%%
%SSFP Signal Estimation
S=FreedDWSSFPOscillating(TR,alpha,optSample.D,optSample.T1,optSample.T2,bOsc);
S0=FreedDWSSFPOscillating(TR,alpha,0,optSample.T1,optSample.T2,bOsc); 
%SSFP Readout Efficiency
Rho=sqrt(min(optScanner.MaxReadout,TR-tau-optScanner.DeadTime)/TR);
%%
%Define Echo Time
TE = (TR-tau-optScanner.DeadTime)/2+tau+optScanner.DeadTime; 
%%
%Estimate T2' decay
T2p=1./((1/optSample.T2s)-(1/optSample.T2));
ET2p=exp(-(TR-TE)./(T2p));
%%
%SSFP SNR Efficiency   
SNReff=S0.*ET2p.*Rho;      
%%
%Estimate effective b-value
bEff=-1./optSample.D.*log(S./S0);
%%
%For optimisation
if nargin==6
    SNReff=[1-SNReff];
end