clear all
%%
%Define default parameters
[opt] = ParameterOptions();
%%
%Create an oscillating gradient waveform (here consistent with definition in Aggarwal et al.)
opt.nOscillations=3;    % Number of oscillations per TR
opt.Waveform='Sine';    % Type of Oscillating Gradient. Choose from Rect (conventional rectangular gradient), Rect_Sym (symmetric rectangular gradient) or Sine (Sine gradient)
opt.G=674;              % Diffusion Gradient Amplitude (mT/m)
opt.tau=20;             % Diffusion Gradient Duration (ms)
opt.TR=25;              % Repetition Time (ms)   
%%
%Rescale x-axis of b-value distribution
opt.bLim=10;            % Maximim b-value cut-off (ms/um2)
%%
%Generate Pathway Amplitude & Gradient Waveform Dictionaries
%First Dictionary
[PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
%Second Dictionary
[PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
%%
%Generate b-value distribution - Here bMin is the b-value associated with a pair of diffusion gradients
[bValue,SignalAmplitudes] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
%Estimate signal amplitude (assuming Gaussian Diffusion) 
S=sum(SignalAmplitudes.*exp(-bValue.*opt.D));
%%
%Estimate b-value & synthesise gradient profile of oscillating gradient waveform (note calculation is performed numerically based on opt structure)
[bOsc,GwaveBlock] = bValueCalc(opt);
%%
%Obtain Analytical solution
SAnalytical=FreedDWSSFPOscillating(opt.TR,opt.alpha,opt.D,opt.T1,opt.T2,bOsc);
%%
%Compare the estimated DW-SSFP signal amplitude from the b-value distribution to an analytical solution (Freed et al. + Fixed Gradient Duration - Appendix 1) 
sprintf('Estimated Signal Amplitude = %f x10^-3, Analytical Signal Amplitude = %f x10^-3', abs(S)*10^3, abs(SAnalytical)*10^3)
%%
%Plot bar chart (note plotting real component only - give consideration if changing RF phase angle)
figure;bar(bValue,real(SignalAmplitudes),'BarWidth',200)
%Plot formatting
xlim([0,opt.bLim])
set(gca,'ytick',0)
title('DW-SSFP b-value Distribution')
set(gca,'yticklabel',0)
xlabel('b (ms/$\mathrm{\mu m^2}$)','Interpreter', 'latex')
ylabel('A (a. u.)','Interpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
%%
%Plot Gradient waveform in a single TR
figure;plot([opt.TR/opt.SpectraResolution:opt.TR/opt.SpectraResolution:opt.TR],GwaveBlock./max(GwaveBlock),'k','linewidth',6)
set(gca,'ytick',[-1,0,1])
xlabel('Time (ms)')
ylabel('Amplitude (normalised)')
title('Gradient Waveform') 
set(findall(gcf,'-property','FontSize'),'FontSize',16)  
     


