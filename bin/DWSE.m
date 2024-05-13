function [S,bEff,Delta] = DWSE(G,tau,b,TR,D,T1,T2)
%
% calculates DW-SE Signal
%
% input:
%	G   = strength of diffusion gradient (mT/m)
%	tau = duration of diffusion gradient (ms)
%	b   = b-value (ms/um2)
%	TR  = repetition time (ms)
%	D   = diffusion coefficient (um2/ms)
%	T1  = longitudinal relaxation time (ms)
%	T2  = transverse relaxation time (ms)
%
% output:
%	SNReff  = Estimated SNR Efficiency
%	bEff    = Estimated b-value
%	Delta   = Estimated Diffusion time
%%
%Define parameters
gamma = 4258*2*pi;      % Hz/G
TR = TR*10^-3;          % convert to s
tau = tau*10^-3;        % convert to s
b=b*1000;               % Convert to s/mm^2
G = G*10^-2;            % convert to G/mm
T1 = T1*10^-3;          % convert to s
T2 = T2*10^-3;          % convert to s
D=D*10^-3;              %Convert to mm2/s
%%
%Estimate Diffusion Time
Delta=max(tau,min(0.1,max(b./(gamma.^2.*tau.^2*G^2)+tau/3,tau)));
%%
%Define echo time
TE=2*Delta;
%%
%Calculate b value
bEff=(gamma.^2.*tau.^2*G^2).*(Delta-tau/3)/1000;
%%
%Estimate Signal
%S=(1-exp(-TR/T1)).*exp(-TE/T2).*exp(-b.*D);
S=(1-exp(-TR/T1)).*exp(-TE/T2);
%%
%Scale Diffusion Time
Delta=Delta*1000;