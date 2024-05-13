clear all
%%
%Define default parameters
[opt] = ParameterOptions();
%%
%Define Diffusion Tensor Components (Taken from a single corpus callosum voxel of the HCP1065 standard-space DTI template in FSL) 
D=[1.1847,-0.0267,-0.0024;-0.0267,0.6334,-0.0099;-0.0024,-0.0099,0.6544];
%%
%Define b-vector (Arbitrary orientation)
bvecs=[1/3;1/3;1/3].^0.5;
%%
%Estimate ADC
ADC=bvecs'*D*bvecs;
%%
%Obtain Analytical solution
SAnalytical=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,ADC,opt.T1,opt.T2);
%%
%Display
sprintf('Analytical Signal Amplitude = %f x10^-3', abs(SAnalytical)*10^3)
