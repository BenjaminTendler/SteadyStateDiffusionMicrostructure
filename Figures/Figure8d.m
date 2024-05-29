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
%Modify default parameter options
opt.Precision=10^-10;      % Cut off for pathway amplitudes 
opt.nPeriods=2;            % Maximim number of TRs in the transverse plane per dictionary element
opt.SpectraResolution=100; % Temporal resolution of TR - Here matched to the MC Simulations
%%
%Define radius of cylinder (um)
R=[0.5:0.5:10];     
%%
%Generate MSD profile
for l=1:length(R)
    [MSD(:,l),t] = MSDAnalytical(R(l),opt);
end
%%
%DWSSFP
%Generate Pathway Amplitude & Gradient Waveform basis set
%First Dictionary
[PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
%Second Dictionary
[PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
%%
%Perform signal estimation - DW-SSFP
for k=1:length(R)
    %Estimate signal attenuation for different cylindrical radii
    [S_DWSSFP(k)]=SAttenuation(opt,MSD(:,k),PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
end
%Estimate S0
[S0_DWSSFP]=SAttenuation(opt,zeros(size(MSD(:,1))),PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
%%
%Load Monte Carlo Data - DW-SSFP
load('MCOutputs/MC_Cylinder_DWSSFP_Oscillating.mat')
%Estimate mean and standard deviation over last 10 TRs
MCmean=mean(MC_Cylinder_DWSSFP_Oscillating(end-1001:100:end,:),1);
MCstd=std(MC_Cylinder_DWSSFP_Oscillating(end-1001:100:end,:),[],1);
MCmeanS0=mean(MC_Cylinder_DWSSFP_Oscillating_S0(end-1001:100:end,1),1);
%%
%DWSE
%Perform signal estimation - DW-SE
opt.PathwayLength=1;      
%Generate Pathway Amplitude & Gradient Waveform basis set
[PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
[PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
%%
%Perform signal estimation - DW-SE
for k=1:length(R)
    %Estimate signal attenuation
    [S_DWSE(k)]=SAttenuation(opt,MSD(:,k),PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
end
%Estimate S0
[S0_DWSE]=SAttenuation(opt,zeros(size(MSD(:,1))),PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
%%
%Load Monte Carlo Data - DW-SE
load('MCOutputs/MC_Cylinder_DWSE_Oscillating.mat')
%%
%Plot
figure;plot(R,abs(S_DWSSFP./S0_DWSSFP),"color","#0072BD",'LineWidth',6);
hold all;errorbar(R,abs(MCmean)./MCmeanS0,MCstd./MCmeanS0,"MarkerSize",50,'LineStyle','None','Marker','.',"MarkerFaceColor","#0072BD","MarkerEdgeColor","#0072BD",'Color',"#0072BD");
hold all;plot(R,abs(S_DWSE./S0_DWSE),"color","#7E2F8E",'LineWidth',6);
hold all;scatter(R,abs(MC_Cylinder_DWSE_Oscillating(end,:)/8),"SizeData",2000,'Marker','.',"MarkerFaceColor","#7E2F8E","MarkerEdgeColor","#7E2F8E",'Color',"#7E2F8E");
xlim([0,10])
ylim([0.8,1.01])
xlabel('Diameter (\mum)')
ylabel('Diffusion Attenuation (normalised)')
title('Restricted Diffusion (Cylinder) (Oscillating Gradients)')
legend('Time-Dependent framework','Monte Carlo (DW-SSFP)','GPA (DW-SE)','Monte Carlo (DW-SE)','FontSize',14,'Location','NorthEast')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
