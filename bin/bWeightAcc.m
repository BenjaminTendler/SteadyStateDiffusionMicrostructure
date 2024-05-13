function [bWeights,bWeightsPos,bLong] = bWeightAcc(bLong,bTrunc,LossArr)
%%
%For each array, determine number of repeats of each element
for k=1:length(bLong)
    maxbLong(k)=max(bLong{k}(:));
end
bArr=[0:max(maxbLong).^0.5].^2;
bLong{1}=bArr;
for k=1:length(bLong)
    bRepeats{k}=histc(bLong{k},[bLong{1}]',1);
end
%%
%Define arrays for different conditions
ZeroWeight=1./(1-LossArr(1));
%%
%1 Longitudinal period
for k=1:length(bLong{1})
    if bLong{1}(k)==0
        bWeights{1}{k}=ZeroWeight;
        bWeightsPos{1}{k}=0;
    else
        bWeightsPos{1}{k}=1:bTrunc;
        bWeights{1}{k}=zeros(size(bWeightsPos{1}{k}));
        spacing=bLong{1}(k):bLong{1}(k):bTrunc;
        bWeights{1}{k}(spacing)=LossArr(1:length(bWeightsPos{1}{k}(1:length(spacing))));
    end
end
%%
%>1 Longitudinal period
for n=1:length(bLong)-1
    for l=1:size(bLong{n+1},2)
        %Determine where arrays differ by one
        idx=find(sum(abs(bRepeats{n+1}(:,l)-bRepeats{n}),1)==1,1,'first');
        value=sum((bRepeats{n+1}(:,l)-bRepeats{n}(:,idx)).*bLong{1}');
        if sum([bLong{n}(:,idx);value])==0
            bWeights{n+1}{l}=ZeroWeight.*bWeights{n}{idx};
            bWeightsPos{n+1}{l}=0;
        elseif value==0
            bWeightsPos{n+1}{l}=bWeightsPos{n}{idx};
            bWeights{n+1}{l}=ZeroWeight.*bWeights{n}{idx};
        elseif sum(bLong{n}(:,idx))==0
            bWeightsPos{n+1}{l}=bWeightsPos{1}{bLong{1}==value};
            bWeights{n+1}{l}=bWeights{n}{idx}*bWeights{1}{bLong{1}==value};
        else
            noZeros=sum(bLong{n+1}(:,l)==0);
            if noZeros>0
                Scale=ZeroWeight.^noZeros;
            else
                Scale=1;
            end
            bWeightsPos{n+1}{l}=1:bTrunc;
            bWeights{n+1}{l}=zeros(size(bWeightsPos{n+1}{l}));
            spacing=value:value:bTrunc;
            bWeights{n+1}{l}(spacing)=LossArr(1:length(spacing))*Scale;
            bWeights{n+1}{l}=bWeights{n+1}{l}+bWeights{n}{idx};
            temp=bWeights{n+1}{l};
            convmat=convnfft(bWeights{n}{idx},bWeights{1}{bLong{1}==value});
            bWeights{n+1}{l}(2:end)=bWeights{n+1}{l}(2:end)+convmat(1:length(temp)-1);
        end
    end
end