clear all
%%
%Define default parameters
[opt] = ParameterOptions();
%%
%Create an oscillating gradient waveform (here consistent with definition in Aggarwal et al.)
opt.nOscillations=3;    % Number of oscillations per TR
opt.Waveform='Sine';    % Type of Oscillating Gradient ('Sine' or 'Rect')
opt.G=674;              % Diffusion Gradient Amplitude (mT/m)
opt.tau=20;             % Diffusion Gradient Duration (ms)
opt.TR=25;              % Repetition Time (ms)   
%%
%Estimate b-value & synthesise gradient profile of oscillating gradient waveform (note calculation is performed numerically based on opt structure)
[bOsc,GwaveBlock] = bValueCalc(opt);
%%
%Obtain Analytical solution
SAnalytical=FreedDWSSFPOscillating(opt.TR,opt.alpha,opt.D,opt.T1,opt.T2,bOsc);
%%
%Display
sprintf('Analytical Signal Amplitude = %f x10^-3', abs(SAnalytical)*10^3)
%%
%Plot Gradient waveform in a single TR
figure;plot([opt.TR/opt.SpectraResolution:opt.TR/opt.SpectraResolution:opt.TR],GwaveBlock./max(GwaveBlock),'k','linewidth',6)
set(gca,'ytick',[-1,0,1])
xlabel('Time (ms)')
ylabel('Amplitude (normalised)')
title('Gradient Waveform') 
set(findall(gcf,'-property','FontSize'),'FontSize',16)