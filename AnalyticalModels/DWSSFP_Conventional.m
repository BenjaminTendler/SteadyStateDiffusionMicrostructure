clear all
%%
%Define default parameters
[opt] = ParameterOptions();
%%
%Obtain Analytical solution
SAnalytical=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,opt.D,opt.T1,opt.T2);
%%
%Display
sprintf('Analytical Signal Amplitude = %f x10^-3', abs(SAnalytical)*10^3)
