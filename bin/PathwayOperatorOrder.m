function [PathwaySignal,Gwave,PathwayTrans] = PathwayOperatorOrder(opt,Scenario,Scale)
%%
%Estimate E1 & E2
E1=exp(-opt.TR/opt.T1);
E2=exp(-opt.TR/opt.T2);
%%
%Redefine precision if Second dictionary
if nargin == 3
    opt.Precision=opt.Precision./max(abs(Scale));
end
%%
%Run initialisation - Scenario indicates if First or Second dictionary
[SignalInitial,StateInitial,GwaveInitial,TransInitial] = InitialisationAmplitude(opt.alpha,opt.phi,E1,E2,Scenario);
%%
%Identify Signal-forming pathways
[PathwaySignal,Gwave,PathwayTrans] = PathwayGenerateOrder(SignalInitial,StateInitial,GwaveInitial,TransInitial,opt);