clear all
%%
%Define Parameters
[opt] = ParameterOptions();
%%
%Define threshold array
PrecisionArray=10.^[-18:-3];
%%
%Full distributionan
for k=1:length(PrecisionArray)
    opt.Precision=PrecisionArray(k);
    %Generate Pathway Amplitude & Gradient Waveform basis set
    %First Dictionary
    [PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
    %Second Dictionary
    [PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
    %Generate Histogram
    [bValue,SignalAmplitudes] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
    %Estimate signal amplitude (Gaussian Diffusion)
    S(k)=sum(SignalAmplitudes.*exp(-bValue.*opt.D));
    S0(k)=sum(SignalAmplitudes);
end
%%
%Analytical Solution
SAnalytical=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,opt.D,opt.T1,opt.T2);
S0Analytical=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,0,opt.T1,opt.T2);
%%
%Plot 
figure;plot(log10(PrecisionArray),abs(S./S0),'LineWidth',6);
yline(SAnalytical/S0Analytical,'--','LineWidth',2);
title('Precision of Time-Independent Framework')
xlabel('Minimum Pathway Amplitude (log_{10} Scale)');
ylabel('Diffusion Attenuation (Normalised)');
set(findall(gcf,'-property','FontSize'),'FontSize',16)