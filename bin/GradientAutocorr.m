function [GwaveAutoCorr]=GradientAutocorr(Gwave,G,dt)
%%
%Define constant - gamma squared
gamma2=71576597699452888; % rad/T
%%
if G==0
    GwaveAutoCorr=0;
else
    %Perform Autocorrelation
    GwaveAutoCorr=autocorr(Gwave,NumLags=length(Gwave)-1);
    %Obtain normalisation coefficient & normalise
    GwaveAutoCorrNorm=dt.^2.*sum(Gwave.*Gwave*gamma2);
    GwaveAutoCorr=GwaveAutoCorr./GwaveAutoCorr(1)*GwaveAutoCorrNorm;
end