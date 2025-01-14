function [S0,nTrans,nPathways] = S0Estimate(opt,PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder)
%%
%Initialise pathways
S0init=0;
S0HigherOrder=0;
%%
%Obtain indices of arrays for pathways to accellerate looping
[~,idxHigherOrder] = sort(abs(PathwaySignalHigherOrder),'descend');
%%
%Define longitudinal loss function
loss=exp(-opt.TR/opt.T1)*cosd(opt.alpha);
%%
%Define arrays for of transverse periods
nTrans=zeros([opt.nPeriods*2,1]);
nPathways=zeros([opt.nPeriods*2,1]);
%%
%Estimate S0
for j=1:length(PathwaySignalFirstOrder)
    if abs(PathwaySignalFirstOrder(j))>0
        nLongPaths=sum(diff(abs(GwaveFirstOrder{j}))-1);
        nTransPaths=TransFirstOrder(j);
        S0init=S0init+PathwaySignalFirstOrder(j);
        nTrans(TransFirstOrder(j))=nTrans(TransFirstOrder(j))+PathwaySignalFirstOrder(j);
        nPathways(TransFirstOrder(j))=nPathways(TransFirstOrder(j))+1;
        [S0HigherOrderIt,nTransIt,nPathwaysIt]=PathwayRecursiveS0Long(PathwaySignalFirstOrder(j),TransFirstOrder(j),PathwaySignalHigherOrder,TransHigherOrder,GwaveHigherOrder,idxHigherOrder,opt.Precision,loss,nLongPaths,opt.PathwayLength,opt.nPeriods);
        S0HigherOrder=S0HigherOrder+S0HigherOrderIt;
        nTrans=nTrans+nTransIt;
        nPathways=nPathways+nPathwaysIt;
        if nLongPaths>0
            for k=1:opt.PathwayLength
                if abs(PathwaySignalFirstOrder(j)).*loss^(k)>opt.Precision
                    Pathway=PathwaySignalFirstOrder(j).*loss^k;
                    S0init=S0init+Pathway*nchoosek(nLongPaths+k-1,k);
                    nTrans(nTransPaths)=nTrans(nTransPaths)+Pathway*nchoosek(nLongPaths+k-1,k);
                    nPathways(nTransPaths)=nPathways(nTransPaths)+nchoosek(nLongPaths+k-1,k);
                else
                    break
                end
            end
        end
    end
end
    %%
%Average
S0=S0init+S0HigherOrder;
