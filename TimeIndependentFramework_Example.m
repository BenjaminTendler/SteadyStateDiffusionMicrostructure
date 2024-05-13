clear all
%%
%Define default parameters
[opt] = ParameterOptions();
%%
%Generate Pathway Amplitude & Gradient Waveform Dictionaries
%First Dictionary
[PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
%Second Dictionary
[PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
%%
%Generate b-value distribution
[bValue,SignalAmplitudes] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
%Estimate signal amplitude (assuming Gaussian Diffusion) 
S=sum(SignalAmplitudes.*exp(-bValue.*opt.D));
%%
%Obtain Comparison Analytical solution
SAnalytical=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,opt.D,opt.T1,opt.T2);
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

     


