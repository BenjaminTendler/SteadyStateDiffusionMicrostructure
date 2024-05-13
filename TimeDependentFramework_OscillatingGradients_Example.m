clear all
%%
%Define default parameters
[opt] = ParameterOptions();
%%
%Create an oscillating gradient waveform (here consistent with definition in Aggarwal et al.)
opt.nOscillations=3;            % Number of oscillations per TR
opt.Waveform='Sine';            % Type of Oscillating Gradient ('Sine' or 'Rect')
opt.G=674;                      % Diffusion Gradient Amplitude (mT/m)
opt.tau=20;                     % Diffusion Gradient Duration (ms)
opt.TR=25;                      % Repetition Time (ms)   
%%
%Modify options for time dependent framework - designed to accelerate the estimation 
opt.Precision=10^-10;           % Cut off for pathway amplitudes 
opt.nPeriods=2;                 % Maximim number of TRs in the transverse plane per dictionary element
opt.SpectraResolution=100;      % Temporal resolution of TR
%%
%Define radius of cylinder (um)
R=5;     
%%
%Generate MSD profile
[MSD] = MSDAnalytical(R,opt);
%%
%DWSSFP
%Generate Pathway Amplitude & Gradient Waveform dictionaries
%First Dictionary
[PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
%Second Dictionary
[PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
%%
%Perform signal estimation - DW-SSFP
S_DWSSFP=SAttenuation(opt,MSD',PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
%Estimate S0
S0_DWSSFP=SAttenuation(opt,zeros(size(MSD')),PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
%%
%Calculate Attenuation
sprintf('Attenuation associated with cylinder of radius %.2f micrometers = %.2f', R, abs(S_DWSSFP./S0_DWSSFP))
%%

%---

%%
%Set plot options for power spectrum plot. Note that if you change the input parameters these may need to be adjusted to visualise the plot
xlims=[-300,300];               % Plot limits (Hz)
ylims=[0,100];                  % Plot limits (Amplitude)
GridSize=300;                   % Sets size of grid
ColorLims=[-0.025,0.025];       % Sets colorbar limits
npaths=1E4;                     % Define number of pathways to plot for DW-SSFP (the n pathways that contribute most to the measurement). Set equal to length(Gwave) to plot all pathways (Estimated in Line 50). Note that this can be extremely slow.
%%
%Generate Pathway Amplitude & Gradient Waveforms explicitly
opt.FullWaveform=1;             % Option to calculate the full waveform of each signal-forming pathway rather than the condensed dictionary. Required to synthesise the Power Spectrum plot
[PathwaySignal,Gwave,~] = PathwayOperatorOrder(opt,'First');
%Plot power spectrum - Note plotting real component only, give consideration if editing RF Phase angle
SpectrumPlotDensity(opt,xlims,ylims,npaths,PathwaySignal,Gwave,GridSize,ColorLims)
