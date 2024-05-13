function [SNReff,bEff,Delta] = DWSE_SNREfficiency(x,optSample,optScanner,b,~)
%%
% calculates DW-SE SNR Efficiency
%
% input:
%	x   = array with diffusion gradient (mT/m), duration of diffusion gradient (ms) and repetition time (ms)
%	optSample   = T1 (ms), T2 (ms) & D (um^2/ms)
%   optScanner = Maximum readout duration (ms) & Minimum time between end of gradient & Readout (ms)
%	b   = b-value (ms/um2)
%	~   = If variable passed return optimisation estimates
%
% output:
%	SNReff = Estimated SNR Efficiency
%	bEff    = Estimated b-value
%	Delta    = Estimated Diffusion Time
%%
%Define array parameters & correct scalings
G=x(1)*10;
tau=x(2)*10;
TR=x(3)*1000;
%%
%SE Signal Estimation
[S,bEff,Delta]=DWSE(G,tau,b,TR,optSample.D,optSample.T1,optSample.T2);
%SE Readout Efficiency
Rho=sqrt(max(0,min(2*(Delta-tau-optScanner.DeadTime),optScanner.MaxReadout))/TR);
%SE SNR Efficiency   
SNReff=S.*Rho;  
%%
%For optimisation - Note that whilst b-value is defined for DW-SE, the additional optimisation term prevents a negative diffusion time. It is important when b-values are low, and has a negligable effect when b-values are high.
if nargin==5
    SNReff=[1-SNReff,abs(bEff-b)./max(b,1)];
end