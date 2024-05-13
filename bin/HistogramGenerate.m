function [bValue,bAmplitude,bOsc] = HistogramGenerate(opt,PathwaySignalFirstOrder,GwaveFirstOrder,PathwaySignalHigherOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder)
%
%Convert bLim to s/mm2
opt.bLim=opt.bLim*1000;
%Convert G to G/cm
opt.G=opt.G/10;
%%
%Initialise b array
bAmplitude=zeros([1,opt.bLim]);
%%
%Define time interval
dt=(opt.TR*10^-3/opt.SpectraResolution);
%%
%Define pathways to evaluate
PathwaySignalFirstOrder(TransFirstOrder>2*opt.nPeriods)=[];
PathwaySignalHigherOrder(TransHigherOrder>2*opt.nPeriods)=[];
GwaveFirstOrder(TransFirstOrder>2*opt.nPeriods)=[];
GwaveHigherOrder(TransHigherOrder>2*opt.nPeriods)=[];
%%
%Define gradient block
GwaveBlock=WaveformCalculate(1,opt.G,opt.tau,opt.TR,opt.SpectraResolution,opt.nOscillations,opt.Waveform);
%%
%Define longitudinal loss function (define decay over 5xMax TRs)
bTrunc=5*opt.PathwayLength;
Loss=exp(-opt.TR/opt.T1)*cosd(opt.alpha);
LossArr=(Loss).^(1:bTrunc);
%%
%Determine the unique longitudinal trajectories
[bLongPathwayUnique]=bLongUnique(GwaveFirstOrder,GwaveHigherOrder);
%Generate weighting functions for longitudinal states
[bWeights,bWeightsPos,bLong] = bWeightAcc(bLongPathwayUnique,bTrunc,LossArr);
%%
%Estimate b for a given trajectory 
%First Order
[bFirstOrder,bOsc] = bHistogramOrderFull(PathwaySignalFirstOrder,GwaveFirstOrder,bWeights,bWeightsPos,bLong,bAmplitude,Loss,GwaveBlock,dt,opt.gamma,opt.bRound*1000);
%Higher order
[bHigherOrder] = bHistogramOrderFull(PathwaySignalHigherOrder,GwaveHigherOrder,bWeights,bWeightsPos,bLong,bAmplitude,Loss,GwaveBlock,dt,opt.gamma,opt.bRound*1000);
%%
%Generate b-value histogram by combining first- and higer-order arrays
%Truncate arrays at bLim
bInit=bFirstOrder(1:opt.bLim);
bHigherOrder=bHigherOrder(1:opt.bLim);
%Perform b-value convolution operation using first and second dictionaries - break when negligable difference with addition of longer pathways
bAmplitude=bInit;
for l=1:opt.PathwayLength
    bTemp=[0,convnfft(bInit,bHigherOrder)];
    bAmplitude=bAmplitude+bTemp(1:length(bAmplitude));
    bInit=bTemp;
    bSum(l)=sum(bAmplitude);
    if l>1
        if abs((bSum(l)-bSum(l-1))./bSum(l-1))<1E-4
            break;
        end
    end
    if isfield(opt,'FourTR')
        break
    end
end
%%
%Define b-values
bValue=(1:opt.bLim)/1000;
bOsc=bOsc/1000;
