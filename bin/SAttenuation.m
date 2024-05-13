function [S] = SAttenuation(opt,MSD,PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
%%
%Calulate exponential decay term E1
E1=exp(-opt.TR/opt.T1);
%%
%Calculate temporal spacing
dt=(opt.TR*10^-3/opt.SpectraResolution);
%Convert G to G/cm
opt.G=opt.G/10;
%%
%Truncate pathways beyond maximum nPeriods
PathwaySignalFirstOrder(TransFirstOrder>2*opt.nPeriods)=[];
PathwaySignalHigherOrder(TransHigherOrder>2*opt.nPeriods)=[];
GwaveFirstOrder(TransFirstOrder>2*opt.nPeriods)=[];
GwaveHigherOrder(TransHigherOrder>2*opt.nPeriods)=[];
TransFirstOrder(TransFirstOrder>2*opt.nPeriods)=[];
TransHigherOrder(TransHigherOrder>2*opt.nPeriods)=[];
%%
%Calculate gradient block
GwaveBlock=WaveformCalculate(1,opt.G,opt.tau,opt.TR,opt.SpectraResolution,opt.nOscillations,opt.Waveform);
%%
%Obtain ordered indices of higher order pathways
[~,idxHigherOrder] = sort(abs(PathwaySignalHigherOrder),'descend');
%%
%Initialise pathways
Sinit=zeros([1,size(MSD,2)]);
SHigherOrder=0;
%%
%Estimate diffusion weighted signal
parfor k=1:length(PathwaySignalFirstOrder)
    %Calcuate Gradient
    [Gwaveinit]=GwaveEstimate(single(GwaveFirstOrder{k}),1);
    %Extend Gradient for all longitudinal states
    [Gw,PathwaySignal,PathwaySignalAtt]=GwaveAttenuation(opt,Gwaveinit,double(GwaveFirstOrder{k}),1,PathwaySignalFirstOrder(k),E1,GwaveBlock,dt,MSD,[]);
    %Add to Signal
    SAtt=sum([PathwaySignalAtt{:}]);
    Sinit=Sinit+SAtt;
    %Obtain higher order terms via recursion
    for l=1:length(Gw)
        SHigherOrderIt=SAttenuationHigherOrder(opt,PathwaySignal{l},PathwaySignalHigherOrder,Gw{l},GwaveHigherOrder,TransFirstOrder(k),TransHigherOrder,dt,MSD,idxHigherOrder,GwaveBlock,E1,SAtt);
        SHigherOrder = SHigherOrder + SHigherOrderIt;
    end
end
S=Sinit+SHigherOrder;