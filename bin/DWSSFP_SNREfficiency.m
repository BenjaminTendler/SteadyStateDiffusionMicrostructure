function [SNReff,bEff] = DWSSFP_SNREfficiency(x,optSample,optScanner,b)
%%
% calculates DW-SSFP SNR Efficiency
%
% input:
%	x   = array with diffusion gradient (mT/m), duration of diffusion gradient (ms), difference between TR & gradient duration (incorporating deadtime) and the flip angle 
%	optSample   = T1 (ms), T2 (ms), T2* (ms) & D (um^2/ms)
%   optScanner = Maximum readout duration (ms) & Minimum time between end of gradient & Readout (ms)
%	b   = b-value (ms/um2) - if provided will perform optimisation
%
% output:
%	SNReff = Estimated SNR Efficiency
%	bEff   = Estimated b-value
%%
%Define array parameters - As TR has to be greater than tau, optimise for the difference
G=x(1);
tau=x(2);
TR=tau+x(3)+optScanner.DeadTime;
alpha=x(4);
%%
%SSFP Signal Estimation
S=abs(FreedDWSSFP(G,tau,TR,alpha,optSample.D,optSample.T1,optSample.T2)); 
S0=abs(FreedDWSSFP(G,tau,TR,alpha,0,optSample.T1,optSample.T2)); 
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
if nargin==4
    SNReff=[1-SNReff,abs(bEff-b)./max(b,1)];
end