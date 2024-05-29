clear all
%%
%Identify optimal SNR Efficiency as function of b-value
%%
%Define Parameters
[optScanner,optSample,optSSFP,optSE] = ParameterOptionsSNROscillating();
%%
%Define b-value associated with single Gradient instance (ms/um2) - Used to calculate relative effective b-value with DW-SSFP only
bOsc=1;
%%
%Define Encoding Period (equivalent to gradient duration)
tauArr=1:40;
%%
%Define fitting options
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','FunctionTolerance',1E-10,'OptimalityTolerance',1E-10,'StepTolerance',1E-10,'MaxIterations',1E4);
%%
%SSFP Optimisation
for k=length(tauArr):-1:1
    optSSFP.tau=tauArr(k);
    %Define fitting function
    f=@(x)DWSSFP_SNREfficiencyOscillating(x,optSample,optScanner,bOsc,optSSFP.tau,1);
    %Perform fitting
    [fit_out_SSFP(k,:)]=lsqnonlin(f,[optSSFP.TR-optSSFP.tau-optScanner.DeadTime,optSSFP.alpha],[1,1],[2000,179],options);
    %Reconstruct estimated SNR efficiency & effective b-value
    [SNR_eff_SSFP(k),bEff_SSFP(k)]=DWSSFP_SNREfficiencyOscillating(fit_out_SSFP(k,:),optSample,optScanner,bOsc,optSSFP.tau);
    %Define optimal Parameters (G, TR, alpha, q, TE)
    Optimised.SSFP(k,:)=[optSSFP.tau+fit_out_SSFP(k,1)+optScanner.DeadTime,fit_out_SSFP(k,2),optSSFP.tau+fit_out_SSFP(k,1)/2+optScanner.DeadTime];
    %Print output
    fprintf('SNR Efficiency = %.3g%%. Effective b-value = %.3g.\n', SNR_eff_SSFP(k)*100,bEff_SSFP(k))
    fprintf('Gradient Duration = %.3g ms.  flip angle = %.3g degrees. TR = %.3g ms. TE = %.3g ms.\n\n',optSSFP.tau, Optimised.SSFP(k,2), Optimised.SSFP(k,1), Optimised.SSFP(k,3))
    %Update initialisation values
    optSSFP.TR=Optimised.SSFP(k,1);optSSFP.alpha=fit_out_SSFP(k,2);
end
%%
%SE Optimisation
for k=length(tauArr):-1:1
    optSE.tau=tauArr(k);
    %Define fitting function
    f=@(x)DWSE_SNREfficiencyOscillating(x,optSample,optScanner,optSE.tau,1);
    %Perform fitting (using same initalisation as DW-STE) - Note that the fitting algorithm prefers when the input parameters (e.g. G, tau & TR) are of the same scale.
    [fit_out_SE(k,:)]=lsqnonlin(f,[optSE.delta/10-optSE.tau/10,optSE.TR/1000],[0.1/10,100/1000],[200/10,20000/1000],options);
    %Reconstruct estimated SNR efficiency & estimated b-value.
    [SNR_eff_SE(k)]=DWSE_SNREfficiencyOscillating(fit_out_SE(k,:),optSample,optScanner,optSE.tau);
    %Define optimal Parameters (G, tau, TR, Delta, TE, Readout)
    Optimised.SE(k,:)=[fit_out_SE(k,2)*1000,fit_out_SE(k,1)*10+optSE.tau,2*(fit_out_SE(k,1)*10+optSE.tau),(fit_out_SE(k,1)*10+optSE.tau-optScanner.DeadTime-optSE.tau)*2];
    fprintf('SNR Efficiency = %.3g%%.\n', SNR_eff_SE(k)*100)
    fprintf('Gradient Duration = %.3g ms. Diffusion Time = %.3g ms. TE = %.3g ms. TR = %.3g ms. Readout Duration = %.3g ms.\n\n', optSE.tau,Optimised.SE(k,2), Optimised.SE(k,3), Optimised.SE(k,1), Optimised.SE(k,4))
    %Update initialisation values
    optSE.TR=fit_out_SE(k,2)*1000;optSE.delta=fit_out_SE(k,1)*10+optSE.tau;
    k
end
%%
%Plot
figure;plot(tauArr,(SNR_eff_SSFP*100),'linewidth',4)
hold all;plot(tauArr,(SNR_eff_SE*100),'linewidth',4)
legend('DW-SSFP','DW-SE');
title('SNR Efficiency (Oscillating Gradients)')
xlabel('Encoding Duration (ms)');
ylabel('SNR Efficiency (% of M_0)');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
