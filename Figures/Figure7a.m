clear all
%%
%Define Parameters
[opt] = ParameterOptions();
%%
%Define Parameters for Oscillating gradients - Consistent with Aggarwal et al.
opt.nOscillations=3;    % Number of oscillations per TR
opt.Waveform='Sine';    % Type of Oscillating Gradient ('Sine' or 'Rect')
opt.G=674;              % Diffusion Gradient Amplitude (mT/m)
opt.tau=20;             % Diffusion Gradient Duration (ms)
opt.TR=25;              % Repetition Time (ms)   
%% 
%Define gradient block
GwaveBlock=WaveformCalculate(1,opt.G,opt.tau,opt.TR,opt.SpectraResolution,opt.nOscillations,opt.Waveform);
%%
%Plot
figure;plot([opt.TR/opt.SpectraResolution:opt.TR/opt.SpectraResolution:opt.TR],GwaveBlock./max(GwaveBlock),'k','linewidth',6)
set(gca,'ytick',[-1,0,1])
xlabel('Time (ms)')
ylabel('Amplitude (normalised)')
title('Gradient Waveform') 
set(findall(gcf,'-property','FontSize'),'FontSize',16)