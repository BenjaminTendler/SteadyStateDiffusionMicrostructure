clear all
%%
%Define Parameters
[opt] = ParameterOptions();
%%
%Define Parameters for Oscillating gradients
opt.nOscillations=3;    % Number of oscillations per TR
opt.Waveform='Sine';    % Type of Oscillating Gradient. Choose from Rect (conventional rectangular gradient), Rect_Sym (symmetric rectangular gradient) or Sine (Sine gradient)
opt.G=674;              % Diffusion Gradient Amplitude (mT/m)
opt.tau=20;             % Diffusion Gradient Duration (ms)
opt.TR=25;              % Repetition Time (ms)   
%%
%Define Flip Angle array for line plot
FlipArr=[1:179];
%Define Flip Angle array for b-value distributions
FlipArrDist=[1,16,179];
%%
%Estimate b-value associated with a single diffusion gradient. Calculation is performed numerically, based on the opt structure. 
[bOsc] = bValueCalc(opt);  
%%
%%Estimate attenuation relationship as function of flip angle
for k=1:length(FlipArr)   
    SAtt(k)=FreedDWSSFPOscillating(opt.TR,FlipArr(k),opt.D,opt.T1,opt.T2,bOsc)./FreedDWSSFPOscillating(opt.TR,FlipArr(k),0,opt.T1,opt.T2,bOsc);
end
%%
%Plot
figure;
hold all;plot(1:179,abs(SAtt),LineWidth=4)
xlabel('Flip Angle (Deg)')
ylabel('Diffusion Attenuation (Normalised)')
title('Attenuation - DW-SSFP (Oscillating Gradients)')
xlim([1,179])
xticks([10:20:170])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
%%
%Create b-value distributions
for k=1:length(FlipArrDist)
    opt.alpha=FlipArrDist(k);
    %Generate Pathway Amplitude & Gradient Waveform basis set
    %First Dictionary
    [PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
    %Second Dictionary
    [PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
    %Generate Histogram
    [bValue,SignalAmplitudes,bOsc] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
    %Estimate signal amplitude (Gaussian Diffusion)
    S=sum(SignalAmplitudes.*exp(-bValue.*opt.D));
    %Obtain Analytical solution
    %Estimate b-value of input gradient
    SAnalytical=FreedDWSSFPOscillating(opt.TR,opt.alpha,opt.D,opt.T1,opt.T2,bOsc/2);
    sprintf('Estimated Signal Amplitude = %f x10^-3, Analytical Signal Amplitude = %f x10^-3', abs(S)*10^3, abs(SAnalytical)*10^3)
    %Plot bar chart
    figure;bar(bValue,real(SignalAmplitudes),'BarWidth',500)
    %Plot formatting
    xlim([0,10])
    yl=max(abs(SignalAmplitudes))*1.1;
    ylim([0,yl])
    set(gca,'ytick',[0,max(abs(SignalAmplitudes))])
    set(gca,'yticklabel',[0,1])
    xlabel('b (ms/$\mathrm{\mu m^2}$)','Interpreter', 'latex')
    ylabel('A (Normalised)','Interpreter','latex')
    title('DW-SSFP b-value Distribution')
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    text(0.8,0.8,[num2str(opt.alpha),'$^{\circ}$'],'Interpreter', 'latex', 'fontsize', 30,'Units','normalized');
end

