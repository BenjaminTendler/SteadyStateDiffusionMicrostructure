clear all
%%
%Define Parameters
[opt] = ParameterOptions();
%%
%Define T2 Array
T2Arr=10:10:300;
%%
%Analytical Solution
for k=1:length(T2Arr)
    opt.T2=T2Arr(k);
    SAnalytical(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,opt.D,opt.T1,opt.T2);
    S0Analytical(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,0,opt.T1,opt.T2);
end
%%
%Full distribution
for k=1:length(T2Arr)
    opt.T2=T2Arr(k);
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
%2TP
%Ensure dictionary entries only persist for two periods in transverse plane
opt.nPeriods=1;
for k=1:length(T2Arr)
    opt.T2=T2Arr(k);
    %Generate Pathway Amplitude & Gradient Waveform basis set
    %First Dictionary
    [PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
    %Second Dictionary
    [PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
    %Remove contribution of higher order pathways (e.g. second dictionary) to preserve only the spin & stimulated echo
    PathwaySignalHigherOrder=PathwaySignalHigherOrder*0;
    %Generate Histogram
    [bValue,SignalAmplitudes] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
    %Estimate signal amplitude (Gaussian Diffusion)
    S_2TP(k)=sum(SignalAmplitudes.*exp(-bValue.*opt.D));
    S0_2TP(k)=sum(SignalAmplitudes);
end
%%
%4TP 
for k=1:length(T2Arr)
    opt.T2=T2Arr(k);
    %Identify pathways that persist for four TRs in transverse plane from first dictionary
    opt.nPeriods=2;
    %Generate Pathway Amplitude & Gradient Waveform basis set
    %First Dictionary
    [PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
    %Second Dictionary
    [PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
    %Remove contribution of higher order pathways (e.g. second dictionary)
    PathwaySignalHigherOrder=PathwaySignalHigherOrder*0;
    %Generate Histogram
    [bValue,SignalAmplitudes] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
    %Estimate signal amplitude (Gaussian Diffusion)
    S_4TP(k)=sum(SignalAmplitudes.*exp(-bValue.*opt.D));
    S0_4TP(k)=sum(SignalAmplitudes);
    %Identify pathways that persist for four TRs in transverse plane from combining first and second dictionary
    %Generate Pathway Amplitude & Gradient Waveform basis set
    opt.nPeriods=1;
    %Define option to preserve only pathways that persist for up to four TRs across both dictionaries (this prevents editing the HistogramGenerate file)
    opt.FourTR=1;
    %Generate Pathway Amplitude & Gradient Waveform basis set
    %First Dictionary
    [PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
    %Second Dictionary
    [PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
    %Generate Histogram
    [bValue,SignalAmplitudes] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
    %Estimate final signal amplitude (Gaussian Diffusion) and prevent double contributin of two transverse period pathways
    S_4TP(k)=S_4TP(k)+sum(SignalAmplitudes.*exp(-bValue.*opt.D))-S_2TP(k);
    S0_4TP(k)=S0_4TP(k)+sum(SignalAmplitudes)-S0_2TP(k);
end
%%
%Plot bar chart
figure;plot([1:30]*10/opt.TR,abs(S./S0Analytical),'LineWidth',6);
hold all;plot([1:30]*10/opt.TR,abs(S_2TP./S0_2TP),'LineWidth',6);
hold all;plot([1:30]*10/opt.TR,abs(S_4TP./S0_4TP),'LineWidth',6);
hold all;plot([1:30]*10/opt.TR,abs(SAnalytical./S0Analytical),'--','LineWidth',6);
title('Precision of Time-Independent Framework')
xlabel('T_2/TR (ratio)');
ylabel('Diffusion Attenuation (Normalised)');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend('Time-Independent framework','Two-Transverse Period approximation','Four-Transverse Period approximation','Analytical','FontSize',12);
%%
%Plot where original T2 value lies on graph
opt=ParameterOptions();
xline(opt.T2/opt.TR,'--','LineWidth',2);
