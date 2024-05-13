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
%Generate Pathway Amplitude & Gradient Waveform basis set
%First Dictionary
[PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
%Second Dictionary
[PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
%%
%Generate distribution
[bValue,SignalAmplitudes,bOsc] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
%Estimate signal amplitude (Gaussian Diffusion) 
S=sum(SignalAmplitudes.*exp(-bValue.*opt.D));
%%
%Obtain Analytical solution
SAnalytical=FreedDWSSFPOscillating(opt.TR,opt.alpha,opt.D,opt.T1,opt.T2,bOsc/2);
%%
%Compare the estimated DW-SSFP signal amplitude from the b-value distribution to an analytical solution (Freed et al. + Oscillating Gradients - Appendix 2) 
sprintf('Estimated Signal Amplitude = %f x10^-3, Analytical Signal Amplitude = %f x10^-3', abs(S)*10^3, abs(SAnalytical)*10^3)
%%
%Plot bar chart (note plotting real component only - give consideration if changing RF phase angle)
figure;bar(bValue,real(SignalAmplitudes),'BarWidth',500)
%Plot formatting
xlim([0,10])
yl=ylim;
ylim([0,yl(2)])
set(gca,'ytick',0)
set(gca,'yticklabel',0)
xlabel('b (ms/$\mathrm{\mu m^2}$)','Interpreter', 'latex')
ylabel('A (a. u.)','Interpreter','latex')
title('DW-SSFP b-value Distribution')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
