function [b,bPath1] = bHistogramOrderFull(PathwaySignal,GwaveOrder,bWeights,bWeightsPos,bLong,b,Loss,GwaveBlock,dt,gamma,bRound)
%%
%Estimate ratio of Longitudinal scaling
 Gwave=GwaveEstimate([1,-2],GwaveBlock);
 bPath1=sum(cumsum(gamma*Gwave*dt).^2)*dt;
 Gwave=GwaveEstimate([1,-3],GwaveBlock);
 bPath2=sum(cumsum(gamma*Gwave*dt).^2)*dt;
 Ratio=round((bPath2-bPath1)/bRound)*bRound;
%%
%Get length threshold
thr=length(b);
%Calculate b-signal for a given pathway order
for k=1:length(GwaveOrder)
    Gwave=GwaveEstimate(GwaveOrder{k},GwaveBlock);
    %Account for balanced gradient condition to give more precision
    if Ratio>0
        bPath=round(sum(cumsum(gamma*Gwave*dt).^2)*dt/bRound)*bRound;
    else
        bPath=round(round(sum(cumsum(gamma*Gwave*dt).^2)*dt/round(bPath1))*round(bPath1));
    end
    if bPath<thr
        b(bPath)=b(bPath)+PathwaySignal(k);
        %Consider longitudinal periods
        nLongPaths=sum(diff(abs(GwaveOrder{k}))-1);
        %Consider initial longitudinal condition
        if abs(GwaveOrder{k}(1))==2
            nLongPaths=nLongPaths+1;
        end
        if nLongPaths>0
            GwaveInit=GwaveEstimate(GwaveOrder{k},1);
            bPathInit=cumsum(GwaveInit);
            bLongPathway=sort(bPathInit(GwaveInit==0).^2,'descend');
            idx=find(ismember(bLong{nLongPaths}',bLongPathway,'rows'));
            bPos=round((bPath+[bWeightsPos{nLongPaths}{idx}]*Ratio));
            bPos=bPos(bPos<thr);
            %Account for balanced gradient condition. 
            if Ratio>0
                b(bPos)=b(bPos)+[bWeights{nLongPaths}{idx}(1:length(bPos))].*PathwaySignal(k);
            else
                b(bPos(1))=b(bPos(1))+sum([bWeights{nLongPaths}{idx}(1:length(bPos))]).*PathwaySignal(k);
            end
            if sum(bLongPathway==0)>0
                if sum(abs(bLongPathway))==0;
                    b(bPath)=b(bPath)-PathwaySignal(k);
                else
                    b(bPath)=b(bPath)-PathwaySignal(k)+(1./(1-Loss)).^sum(bLongPathway==0)*PathwaySignal(k);
                end
            end
        end
    else
        ;
    end
end