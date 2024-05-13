function [optScanner,optSample,optSSFP,optSE] = ParameterOptionsSNROscillating()
%%
%Scanner properties
optScanner.MaxReadout=inf;                       % Maximum Readout Duration (ms)
optScanner.GMax=1000;                            % Maximum Gradient Stength (mT/m)
optScanner.gamma=2*pi*42.58*10^6;                % Gyromagnetic ratio
optScanner.DeadTime=0;                           % Minimum Time between Diffusion Gradient and Readout (ms)
%%
%Sample Properties
optSample.T1=552;                                % T1 (ms)
optSample.T2=26.8;                               % T2 (ms)
optSample.T2s=optSample.T2;                      % T2 star (ms) - used to estimate signal-loss due to non-centred DW-SSFP echo. Set equal to T2 to ignore contribution
optSample.D=0.14;                                % Diffusion coefficient (um2/ms)
%%
%Initial sequence parameters - SSFP
optSSFP.G=optScanner.GMax;                       % Diffusion Gradient Amplitude (mT/m)
optSSFP.tau=10;                                  % Diffusion Gradient Duration (ms)
optSSFP.TR=35;                                   % Repetition Time (ms)
optSSFP.alpha=20;                                % Flip angle (o)
%%
%Initial sequence parameters - SE
optSE.G=optScanner.GMax;                         % Diffusion Gradient Amplitude (mT/m)
optSE.tau=10;                                    % Diffusion Gradient Duration (ms)
optSE.TR=1000;                                   % Repetition Time (ms)
optSE.delta=optSE.tau+1;                         % Diffusion Time (ms)