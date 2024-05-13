function [S0,nTrans,nPathways] = PathwayRecursiveS0Long(PathwaySignalFirstOrder,TransFirstOrder,PathwaySignalHigherOrder,TransHigherOrder,GwaveHigherOrder,idxHigherOrder,Precision,loss,nLongPathsInit,PathwayLength,nPeriods)
%%
%Define number of transverse periods
nTrans=zeros([100,1]);
nPathways=zeros([100,1]);
%%
%Perform correction to only investigate pathways which correspond to the nPeriods condition
idxHigherOrdernPeriods=idxHigherOrder((TransFirstOrder+TransHigherOrder(idxHigherOrder))<=2*nPeriods);
%%
%Peform recursive operation to characterise amplitude of each pathway
S0=0;
for l=1:length(idxHigherOrdernPeriods)
    if length(idxHigherOrdernPeriods)==0
        break
    end
    %Calculate amplitude of pathway
    PathwayAmplitude=PathwaySignalFirstOrder*PathwaySignalHigherOrder(idxHigherOrdernPeriods(l));
    Gwave=GwaveHigherOrder{idxHigherOrdernPeriods(l)};
    nTransPaths=TransFirstOrder+TransHigherOrder(idxHigherOrdernPeriods(l));
    if abs(PathwayAmplitude) > Precision
        %Add pathway if greater than Precision
        nTrans(nTransPaths)=nTrans(nTransPaths)+PathwayAmplitude;    
        nPathways(nTransPaths)=nPathways(nTransPaths)+1;
        nLongPaths=sum(diff(abs(Gwave))-1)+nLongPathsInit;
        if abs(Gwave(1))==2
            nLongPaths=nLongPaths+1;
        end
        S0=S0 + PathwayAmplitude;
        [S0It,nTransIt,nPathwaysIt]=PathwayRecursiveS0Long(PathwayAmplitude,nTransPaths,PathwaySignalHigherOrder,TransHigherOrder,GwaveHigherOrder,idxHigherOrder,Precision,loss,nLongPaths,PathwayLength,nPeriods);
        S0=S0+S0It;
        nTrans=nTrans+nTransIt;
        nPathways=nPathways+nPathwaysIt;
        if  nLongPaths>0
            for k=1:PathwayLength
                if abs(PathwayAmplitude).*loss^k>Precision
                    Pathway=PathwayAmplitude.*loss^k;
                    S0=S0+Pathway*nchoosek(nLongPaths+k-1,k);
                    nTrans(nTransPaths)=nTrans(nTransPaths)+Pathway*nchoosek(nLongPaths+k-1,k);
                    nPathways(nTransPaths)=nPathways(nTransPaths)+nchoosek(nLongPaths+k-1,k);
                else
                    break
                end
            end
        end
    else
        break;
    end
end
