function [S,bEff,Tmix] = DWSTE(G,tau,b,TE,TR,D,T1,T2)
%
% calculates DW-STE Signal
%
% input:
%	G   = strength of diffusion gradient (mT/m)
%	tau = duration of diffusion gradient (ms)
%	b   = b-value (ms/um2)
%   TE  = Echo time (ms)
%	TR  = repetition time (ms)
%	D   = diffusion coefficient (um2/ms)
%	T1  = longitudinal relaxation time (ms)
%	T2  = transverse relaxation time (ms)
%
% output:
%	SNReff = Estimated SNR Efficiency
%	bEff    = Estimated b-value
%	Tmix    = Estimated Mixing Time
%%
%Define parameters
gamma = 4258*2*pi;      % Hz/G
TE = TE*10^-3;          % convert to s
TR = TR*10^-3;          % convert to s
tau = tau*10^-3;        % convert to s
b=b*1000;               % Convert to s/mm^2
G = G*10^-2;            % convert to G/mm
T1 = T1*10^-3;          % convert to s
T2 = T2*10^-3;          % convert to s
D=D*10^-3;              %Convert to mm2/s
%%
%Estimate Mixing Time
Tmix=min(max(b./(gamma.^2.*tau.^2*G^2)+tau/3-TE/2,0),TR);
%%
%Calculate b value
bEff=(gamma.^2.*tau.^2*G^2).*(TE/2+Tmix-tau/3)/1000;
%%
%Estimate Signal
%S=1/2.*(1-exp(-(TR-Tmix)/T1)).*exp(-Tmix/T1).*exp(-TE/T2).*exp(-b.*D);
S=1/2.*(1-exp(-(TR-Tmix)/T1)).*exp(-Tmix/T1).*exp(-TE/T2);
%%
%Scale Mixing Time
Tmix=Tmix*1000;