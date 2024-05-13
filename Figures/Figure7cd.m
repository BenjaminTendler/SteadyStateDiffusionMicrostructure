clear all
%%
%Define default parameters
[opt] = ParameterOptions();
%%
%Define Parameters for Oscillating gradients
opt.nOscillations=3;    % Number of oscillations per TR
opt.Waveform='Sine';    % Type of Oscillating Gradient ('Sine' or 'Rect')
opt.G=674;              % Diffusion Gradient Amplitude (mT/m)
opt.tau=20;             % Diffusion Gradient Duration (ms)
opt.TR=25;              % Repetition Time (ms)   
%%
%Edit parameters for power spectrum plot
opt.FullWaveform=1;             % Option to calculate the full waveform of each signal-forming pathway rather than the condensed dictionary. Required to synthesise the Power Spectrum plot
opt.nPeriods=2;                 % Maximim number of TRs in the transverse plane (Here 1 = 2-transverse period approximation, 2 = 4-transverse period approximation etc)
opt.SpectraResolution=1000;     % Temporal resolution of TR - Note this needs to be higher versus the conventional case to accurately characterise the oscillating gradient profile.
opt.Precision=1E-10;            % Cut off for pathway amplitudes
npaths=1E4;                     % Define number of pathways to plot for DW-SSFP (the n pathways that contribute most to the measurement). Set equal to length(Gwave) to plot all pathways (Estimated in Line 21). Note that this can be extremely slow.
%%
%DWSSFP
%Set plot options - DW-SSFP. Note that if you change the input parameters these may need to be adjusted to visualise the plot
xlims=[-300,300];               %Plot limits (Hz)
ylims=[0,100];                  %Plot limits (Amplitude)
GridSize=300;                   %Sets size of grid
ColorLims=[-0.025,0.025];       %Sets colorbar limits
%%
%Generate Pathway Amplitude & Gradient Waveforms explicitly
[PathwaySignal,Gwave,~] = PathwayOperatorOrder(opt,'First');
%Plot - DW-SSFP - Note plotting real component only, give consideration if editing RF Phase angle
SpectrumPlotDensity(opt,xlims,ylims,npaths,PathwaySignal,Gwave,GridSize,ColorLims)
%%
%DWSE
%Set plot options - DW-SE. Note that if you change the input parameters these may need to be adjusted to visualise the plot
xlims=[-300,300];         %Plot limits (Hz)
ylims=[0,20];             %Plot limits (Amplitude)
opt.PathwayLength=1;      %Only plot the spin echo pathway by setting the maximum pathway length to 2 TRs.
%%
%Generate Pathway Amplitude & Gradient Waveforms explicitly
[PathwaySignal,Gwave,Trans] = PathwayOperatorOrder(opt,'First');
%Plot - DW-SE
SpectrumPlotDensity(opt,xlims,ylims,npaths,PathwaySignal,Gwave,GridSize)
%Edit title for DWSE
title('DW-SE Power Spectrum')
