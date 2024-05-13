function [SNReff] = DWSE_SNREfficiencyOscillating(x,optSample,optScanner,tau,~)
%%
% calculates DW-SE SNR Efficiency
%
% input:
%	x   = array with diffusion gradient (mT/m) and repetition time (ms)
%	optSample   = T1 (ms), T2 (ms) & D (um^2/ms)
%   optScanner = Maximum readout duration (ms) & Minimum time between end of gradient & Readout (ms)
%	tau   = encoding duration
%	~   = If variable passed return optimisation estimates
%
% output:
%	SNReff  = Estimated SNR Efficiency
%%
%Define array parameters & correct scalings
Delta=x(1)*10+tau;
TR=x(2)*1000;
%%
TE=Delta*2;
%%
%STE Signal Estimation
S=(1-exp(-TR/optSample.T1)).*exp(-TE/optSample.T2);
%STE Readout Efficiency
Rho=sqrt(max(0,min(2*(Delta-tau-optScanner.DeadTime),optScanner.MaxReadout))/TR);
%STE SNR Efficiency   
SNReff=S.*Rho;  
%%
%For optimisation
if nargin==5
    SNReff=[1-SNReff];
end