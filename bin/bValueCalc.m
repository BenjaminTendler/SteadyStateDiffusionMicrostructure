function [bValue,Gwave] = bValueCalc(opt)    
%%
%Define Diffusion Gradient
Gwave=WaveformCalculate(1,opt.G,opt.tau,opt.TR,opt.SpectraResolution,opt.nOscillations,opt.Waveform);
%%
%Calculates the b-value for a given gradient waveform profile and input options;
bValue=sum(cumsum((opt.gamma/10)*Gwave*(opt.TR*10^-3)./opt.SpectraResolution).^2)*(opt.TR*10^-3)./opt.SpectraResolution/10^3;