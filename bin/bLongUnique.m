function [bLongPathwayUnique] = bLongUnique(GwaveFirstOrder,GwaveHigherOrder)
%%
%Concatenate arrays
GwaveOrder=[GwaveFirstOrder,GwaveHigherOrder];
%%
%Obtain b-values corresponding to longitudinal states
for k=1:length(GwaveOrder)
    %Consider longitudinal periods
    nLongPaths=sum(diff(abs(GwaveOrder{k}))-1);
    %Consider initial longitudinal condition
    if abs(GwaveOrder{k}(1))==2
        nLongPaths=nLongPaths+1;
    end
    if nLongPaths>0
        GwaveInit=GwaveEstimate(GwaveOrder{k},1);
        bPathInit=cumsum(GwaveInit);
        bLongPathway{k}=sort(bPathInit(GwaveInit==0).^2,'descend');
    end
end
%%
%Sort data by length 
[Length,Index]=sort(cellfun(@length,bLongPathway));
%Obtain unique elements for each length of array
for k=1:Length(end)
    bLongPathwayUnique{k}=unique(vertcat(bLongPathway{Index(Length==k)}),'rows')';
end