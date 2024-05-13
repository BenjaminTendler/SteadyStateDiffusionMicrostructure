clear all
%%
%Identify optimal SNR Efficiency as function of Gradient Amplitude
%%
%Define Parameters
[optScanner,optSample,optSSFP,optSTE,optSE] = ParameterOptionsSNR();
%%
%Define T1 array - Note that fitting algorithm will fix the gradient amplitude to GMax during the fitting (defined in ParameterOptionsSNR)
T1Arr=50:50:2000;
%%
%Define b-value to optimise (ms/um2)
b=10;
%%
%Define fitting options
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','FunctionTolerance',1E-10,'OptimalityTolerance',1E-10,'StepTolerance',1E-10,'MaxIterations',1E4);
%%
%SSFP Optimisation
for k=length(T1Arr):-1:1
    optSample.T1=T1Arr(k);
    %Define fitting function
    f=@(x)DWSSFP_SNREfficiency(x,optSample,optScanner,b);
    %Perform fitting
    [fit_out_SSFP(k,:)]=lsqnonlin(f,[optScanner.GMax,optSSFP.tau,optSSFP.TR-optSSFP.tau-optScanner.DeadTime,optSSFP.alpha],[optScanner.GMax,0.1,1,1],[optScanner.GMax,80,100,179],options);
    %Reconstruct estimated SNR efficiency & effective b-value
    [SNR_eff_SSFP(k),bEff_SSFP(k)]=DWSSFP_SNREfficiency(fit_out_SSFP(k,:),optSample,optScanner);
    %Define optimal Parameters (G, tau, TR, alpha, q, TE)
    Optimised.SSFP(k,:)=[fit_out_SSFP(k,1),fit_out_SSFP(k,2),fit_out_SSFP(k,2)+fit_out_SSFP(k,3)+optScanner.DeadTime,fit_out_SSFP(k,4),optScanner.gamma.*fit_out_SSFP(k,2).*fit_out_SSFP(k,1)/(2*pi)/10^8,fit_out_SSFP(k,2)+fit_out_SSFP(k,3)/2+optScanner.DeadTime];
    %Print output
    fprintf('SNR Efficiency = %.3g%%. Target b-value = %.3g. Effective b-value = %.3g.\n', SNR_eff_SSFP(k)*100,b,bEff_SSFP(k))
    fprintf('Gradient Strength = %.3g mT/m. Gradient Duration = %.3g ms. q = %.3g /cm.  flip angle = %.3g degrees. TR = %.3g ms. TE = %.3g ms.\n\n', Optimised.SSFP(k,1), Optimised.SSFP(k,2), Optimised.SSFP(k,5), Optimised.SSFP(k,4), Optimised.SSFP(k,3), Optimised.SSFP(k,6))
    %Update initialisation values
    optSSFP.G=fit_out_SSFP(k,1);optSSFP.tau=fit_out_SSFP(k,2);optSSFP.TR=Optimised.SSFP(k,3);optSSFP.alpha=fit_out_SSFP(k,4);
end;
%%
%STE Optimisation
for k=length(T1Arr):-1:1
    optSample.T1=T1Arr(k);
    %Define fitting function
    f=@(x)DWSTE_SNREfficiency(x,optSample,optScanner,b,1);
    %Perform fitting - Note that the fitting algorithm prefers when the input parameters (e.g. G, tau, TE & TR) are of the same scale. 
    [fit_out_STE(k,:)]=lsqnonlin(f,[optScanner.GMax/10,optSTE.tau/10,(optSTE.TE-2*optSTE.tau)/10,optSTE.TR/1000],[optScanner.GMax/10,1/10,0,100/1000],[optScanner.GMax/10,80/10,50/10,20000/1000],options);
    %Reconstruct estimated SNR efficiency & estimated b-value.
    [SNR_eff_STE(k),bEff_STE(k),Tmix_STE(k)]=DWSTE_SNREfficiency(fit_out_STE(k,:),optSample,optScanner,b);
    %Define optimal Parameters (G, tau, TE, TR, Mixing Time, Readout Duration)
    Optimised.STE(k,:)=[fit_out_STE(k,1)*10,fit_out_STE(k,2)*10,2*fit_out_STE(k,2)*10+fit_out_STE(k,3)*10,fit_out_STE(k,4)*1000,Tmix_STE(k),((2*fit_out_STE(k,2)*10+fit_out_STE(k,3)*10)-2*(optScanner.DeadTime+fit_out_STE(k,2)*10))];
    fprintf('SNR Efficiency = %.3g%%. Target b-value = %.3g. Effective b-value = %.3g.\n', SNR_eff_STE(k)*100,b,bEff_STE(k))
    fprintf('Gradient Strength = %.3g mT/m. Gradient Duration = %.3g ms. Mixing Time = %.3g ms. TE = %.3g ms. TR = %.3g ms.  Readout Duration = %.3g ms.\n\n', Optimised.STE(k,1), Optimised.STE(k,2), Optimised.STE(k,5), Optimised.STE(k,3), Optimised.STE(k,4), Optimised.STE(k,6))
    %Update initialisation values
    optSTE.G=fit_out_STE(k,1)*10;optSTE.tau=fit_out_STE(k,2)*10;optSTE.TE=2*fit_out_STE(k,2)*10+fit_out_STE(k,3)*10;optSTE.TR=fit_out_STE(k,4)*1000;
end
%%
%SE Optimisation
for k=length(T1Arr):-1:1
    optSample.T1=T1Arr(k);
    %Define fitting function
    f=@(x)DWSE_SNREfficiency(x,optSample,optScanner,b,1);
    %Perform fitting (using same initalisation as DW-STE) - Note that the fitting algorithm prefers when the input parameters (e.g. G, tau & TR) are of the same scale.
    [fit_out_SE(k,:)]=lsqnonlin(f,[optScanner.GMax/10,optSE.tau/10,optSE.TR/1000],[optScanner.GMax/10,0.1/10,100/1000],[optScanner.GMax/10,80/10,20000/1000],options);
    %Reconstruct estimated SNR efficiency & estimated b-value.
    [SNR_eff_SE(k),bEff_SE(k),Delta_SE(k)]=DWSE_SNREfficiency(fit_out_SE(k,:),optSample,optScanner,b);
    %Define optimal Parameters (G, tau, TR, Delta, TE, Readout)
    Optimised.SE(k,:)=[fit_out_SE(k,1)*10,fit_out_SE(k,2)*10,fit_out_SE(k,3)*1000,Delta_SE(k),Delta_SE(k)*2,(Delta_SE(k)-optScanner.DeadTime-fit_out_SE(k,2)*10)*2];
    fprintf('SNR Efficiency = %.3g%%. Target b-value = %.3g. Effective b-value = %.3g.\n', SNR_eff_SE(k)*100,b,bEff_SE(k))
    fprintf('Gradient Strength = %.3g mT/m. Gradient Duration = %.3g ms. Diffusion Time = %.3g ms. TE = %.3g ms. TR = %.3g ms. Readout Duration = %.3g ms.\n\n', Optimised.SE(k,1), Optimised.SE(k,2), Optimised.SE(k,4), Optimised.SE(k,5), Optimised.SE(k,3), Optimised.SE(k,6))
    %Update initialisation values
    optSE.G=fit_out_SE(k,1)*10;optSE.tau=fit_out_SE(k,2)*10;optSE.TR=fit_out_SE(k,3)*1000;
end
%%
%Plot
figure;plot(T1Arr,(SNR_eff_SSFP*100),'linewidth',4)
hold all;plot(T1Arr,(SNR_eff_STE*100),'linewidth',4)
hold all;plot(T1Arr,(SNR_eff_SE*100),'linewidth',4)
legend('DW-SSFP','DW-STE','DW-SE','Location','northeast');
title('SNR Efficiency')
xlabel('T_{1} (ms)');
ylabel('SNR Efficiency (% of M_0)');
set(findall(gcf,'-property','FontSize'),'FontSize',16)


