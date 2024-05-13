clear all
%%
%Define Parameters
[opt] = ParameterOptions();
%%
%Define Parameters for Oscillating gradients
opt.nOscillations=3;        % Number of oscillations per TR
opt.Waveform='Sine';        % Type of Oscillating Gradient ('Sine' or 'Rect')
opt.G=674;                  % Diffusion Gradient Amplitude (mT/m)
opt.tau=20;                 % Diffusion Gradient Duration (ms)
opt.TR=25;                  % Repetition Time (ms) 
%%
%Define flip angle array
alpha=10:10:170;
%%
%Calculate analyical and histogram solution as function of flip angle (may take several minutes to run)
for k=1:length(alpha)
    opt.alpha=alpha(k);
    %Generate Pathway Amplitude & Gradient Waveform basis set
    %First Dictionary
    [PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
    %Second Dictionary
    [PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
    %%
    %Generate Histogram
    [bValue,SignalAmplitudes,bOsc] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
    %Estimate signal amplitude (Gaussian Diffusion)
    S(k)=sum(SignalAmplitudes.*exp(-bValue.*opt.D));
    %Obtain Analytical solution
    SAnalytical(k)=FreedDWSSFPOscillating(opt.TR,opt.alpha,opt.D,opt.T1,opt.T2,bOsc(1)/2);
    S0Analytical(k)=FreedDWSSFPOscillating(opt.TR,opt.alpha,0,opt.T1,opt.T2,bOsc(1)/2);
end
%%
%Load Monte Carlo Data
load('MCOutputs/MC_Gaussian_Oscillating_FlipAngle.mat')
%Estimate mean and standard deviation over last 10 TRs
MCmean=mean(S_MC_alpha(end-1001:100:end,:),1);
MCstd=std(S_MC_alpha(end-1001:100:end,:),[],1);
%%
figure;plot(alpha,abs(S./S0Analytical),'LineWidth',6);
hold all;plot(alpha,abs(SAnalytical./S0Analytical),'--','LineWidth',6)
hold all;errorbar(alpha,abs(MCmean./S0Analytical),MCstd./S0Analytical,"MarkerSize",50,'LineStyle','None','Marker','.');
xlim([0,180])
xlabel('Flip Angle (Deg.)')
ylabel('Diffusion Attenuation (Normalised)')
title('Free Gaussian Diffusion (Oscillating Gradients)')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend('Time-Independent framework (Oscillating)','Analytical (Oscillating)','Monte Carlo (Oscillating)','FontSize',14,'location','southwest')