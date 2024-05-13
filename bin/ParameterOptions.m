function [opt] = ParameterOptions()
%%
%Sequence parameters
opt.gamma=2*pi*42.58*10^6;  % Gyromagnetic ratio
opt.G=52;                   % Diffusion Gradient Amplitude (mT/m)
opt.tau=13.56;              % Diffusion Gradient Duration (ms)
opt.TR=28;                  % Repetition Time (ms)
opt.alpha=24;               % Flip angle (o)
opt.phi=-pi/2;              % RF Phase (radians)
opt.nOscillations=0;        % Number of oscillations of the diffusion gradient (assumes spoiling independent of gradient)
opt.Waveform='Sine';        % Type of Oscillating Gradient (if opt.nOscillations>0). Choose from 'Sine' or 'Rect'
%%
%Sample Properties
opt.T1=600;                 % T1 (ms)
opt.T2=40;                  % T2 (ms)
opt.D=0.2;                  % Diffusion coefficient (um2/ms)
%%
%Reconstruction parameters - These parameters impact the precision of the signal estimate
opt.nPeriods=5;             % Maximim number of TR pairs in the transverse plane per dictionary element
opt.Precision=1E-28;        % Cut off for pathway amplitudes
opt.PathwayLength=1000;     % Maximum pathway length (number of TRs)
opt.bLim=150;               % b-value cut-off (ms/um2)
opt.SpectraResolution=1000; % Temporal resolution of TR - Higher value equals greater precision of estimations, but longer computation time
opt.bRound=0.001;           % For b-value distribution, round b-values to this value (m/um2 - minimum 0.001 ms/um2).




