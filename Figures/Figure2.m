clear all
%%
%Define default Parameters
[opt] = ParameterOptions();
%%
%Estimate diffusion attenuation as function of flip angle (1-180o)
FlipArray=[1:179/99:180];
for k=1:length(FlipArray)
    S0_Flip(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,FlipArray(k),0,opt.T1,opt.T2);
    S_Flip(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,FlipArray(k),opt.D,opt.T1,opt.T2);
end
%%
%Estimate diffusion attenuation as function of TR (5-100 ms)
TRArray=[5:95/99:100];
for k=length(TRArray)
    S0_TR(k)=FreedDWSSFP(opt.G,opt.tau,TRArray(k),opt.alpha,0,opt.T1,opt.T2);
    S_TR(k)=FreedDWSSFP(opt.G,opt.tau,TRArray(k),opt.alpha,opt.D,opt.T1,opt.T2);
end
%%
%Estimate diffusion attenuation as function of G (0-300 mT/m)
GArray=[0:300/99:300];
for k=1:length(GArray)
    S0_G(k)=FreedDWSSFP(GArray(k),opt.tau,opt.TR,opt.alpha,0,opt.T1,opt.T2);
    S_G(k)=FreedDWSSFP(GArray(k),opt.tau,opt.TR,opt.alpha,opt.D,opt.T1,opt.T2);
end
%%
%Estimate diffusion attenuation as function of tau (1:20 ms)
tauArray=[1:19/99:20];
for k=1:length(tauArray)
    S0_tau(k)=FreedDWSSFP(opt.G,tauArray(k),opt.TR,opt.alpha,0,k,opt.T2);
    S_tau(k)=FreedDWSSFP(opt.G,tauArray(k),opt.TR,opt.alpha,opt.D,k,opt.T2);
end
%%
%Estimate diffusion attenuation as function of T1 (10-1000 ms)
T1Array=[10:990/99:1000];
for k=1:length(T1Array)
    S0_T1(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,0,T1Array(k),opt.T2);
    S_T1(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,opt.D,T1Array(k),opt.T2);
end
%%
%Estimate diffusion attenuation as function of T2 (1-100 ms)
T2Array=[1:99/99:100];
for k=1:length(T2Array)
    S0_T2(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,0,opt.T1,T2Array(k));
    S_T2(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,opt.D,opt.T1,T2Array(k));
end
%%
%Estimate diffusion attenuation as function of D (0 - 1 um/ms)
DArray=[0:1/99:1];
for k=1:length(DArray)
    S0_D(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,0,opt.T1,opt.T2);
    S_D(k)=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,DArray(k),opt.T1,opt.T2);
end
%%
%Plot Sequence parameter dependencies
figure;plot(S_G./S0_G,LineWidth=4)
hold all;plot(S_tau./S0_tau,LineWidth=4)
hold all;plot(S_Flip./S0_Flip,'--',LineWidth=4)
hold all;plot(S_TR./S0_TR,'--',LineWidth=4)
xlabel('G / \delta / Flip Angle / TR')
ylabel('Diffusion Attenuation (Normalised)')
title('DW-SSFP dependence on sequence parameters')
xlim([1,100])
ylim([0,1])
xticks([1,100])
xticklabels({'      Low','High      '})
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend('G (0-300 mT/m)','\delta (0-20 ms)','Flip Angle (1^o-180^o)','TR (5-100 ms)','Location','northeast','FontSize',10)
%%
%Plot relaxation dependencies
figure;plot(S_D./S0_D,LineWidth=4)
hold all;plot(S_T1./S0_T1,'--',LineWidth=4)
hold all;plot(S_T2./S0_T2,'--',LineWidth=4)
xlabel('D / T_{1} / T_{2}')
ylabel('Diffusion Attenuation (Normalised)')
title('DW-SSFP dependence on sample properties')
xlim([1,100])
ylim([0,1])
xticks([1,100])
xticklabels({'      Low','High      '})
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend('D (0-1 \mum^2/ms)','T_{1} (10-1000 ms)','T_{2} (1-100 ms)','Location','northeast','FontSize',10)