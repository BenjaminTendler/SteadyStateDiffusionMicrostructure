clear all
%%
%Note that this code does not synthesise the Monte Carlo data directly. The MCOutputs folder provides the MC estimates for the default parameter options
%%
%Define Parameters
[opt] = ParameterOptions();
%%
%Define flip angle array
alpha=10:10:170;
%%
%Calculate analyical and b-value distribution solution as function of flip angle (may take several minutes to run)
for k=1:length(alpha)
    opt.alpha=alpha(k);
    %Generate Pathway Amplitude & Gradient Waveform basis set
    %First Dictionary
    [PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
    %Second Dictionary
    [PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
    %%
    %Generate Histogram
    [bValue,SignalAmplitudes] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
    %Estimate signal amplitude (Gaussian Diffusion)
    S(k)=sum(SignalAmplitudes.*exp(-bValue.*opt.D));
    %Obtain Analytical solution
    SAnalytical(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,opt.D,opt.T1,opt.T2);
    S0Analytical(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,0,opt.T1,opt.T2);
end
%%
%Load Monte Carlo Data
load('MCOutputs/MC_Gaussian_FlipAngle.mat')
%Estimate mean and standard deviation over last 10 TRs
MCmean=mean(MC_Gaussian_FlipAngle(end-1001:100:end,:),1);
MCstd=std(MC_Gaussian_FlipAngle(end-1001:100:end,:),[],1);
%%
%Plot figure
figure;plot(alpha,abs(S./S0Analytical),'LineWidth',6);
hold all;plot(alpha,abs(SAnalytical./S0Analytical),'--','LineWidth',6)
hold all;errorbar(alpha,abs(MCmean./S0Analytical),MCstd./S0Analytical,"MarkerSize",50,'LineStyle','None','Marker','.');
xlim([0,180])
xlabel('Flip Angle (Deg.)')
ylabel('Diffusion Attenuation (Normalised)')
title('Free Gaussian Diffusion')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend('Time-Independent Framework','Analytical','Monte Carlo','FontSize',14,'location','southeast')