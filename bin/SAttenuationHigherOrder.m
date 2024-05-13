function [S] = SAttenuationHigherOrder(opt,PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder,dt,MSD,idxHigherOrder,GwaveBlock,E1,Sinit)

%%
%Initialise signal amplitude
S=zeros([1,size(MSD,2)]);
%%
%Perform correction to only investigate pathways which correspond to the opt.nPeriods condition
idxHigherOrder_nPeriods=idxHigherOrder((TransFirstOrder+TransHigherOrder(idxHigherOrder))<=2*opt.nPeriods);
%%
%Estimate diffusion weighted signal
for l=1:length(idxHigherOrder_nPeriods)
    %Calculate non-diffusion weighted amplitude of pathway
    PathwayAmplitude=PathwaySignalFirstOrder*PathwaySignalHigherOrder(idxHigherOrder_nPeriods(l));
    nTransPaths=TransFirstOrder+TransHigherOrder(idxHigherOrder_nPeriods(l));
    %Add pathway if greater than opt.Precision
    if abs(PathwayAmplitude) > opt.Precision
        %Define gradient waveform of pathway
        Gwave=GwaveHigherOrder{idxHigherOrder_nPeriods(l)};
        %Calcuate Q-vector
        [Gwaveinit]=GwaveEstimate(double(Gwave),1);
        if abs(Gwave(1))==2
            Gwave=[0,Gwave];
        end
        %Extend Q-vector for all longitudinal states
        [Gw,PathwaySignal,PathwaySignalAtt]=GwaveAttenuation(opt,Gwaveinit,Gwave,1,PathwayAmplitude,E1,GwaveBlock,dt,MSD,GwaveFirstOrder);
        %Add to Signal
        SAtt=sum([PathwaySignalAtt{:}]);
        S=S+SAtt;
        %Obtain higher order terms via recursion
        for k=1:length(Gw)
            [SIt]=SAttenuationHigherOrder(opt,PathwaySignal{k},PathwaySignalHigherOrder,[GwaveFirstOrder,Gw{k}],GwaveHigherOrder,nTransPaths,TransHigherOrder,dt,MSD,idxHigherOrder,GwaveBlock,E1,SAtt);
            S=S+SIt;
        end
    else
        break
    end
end