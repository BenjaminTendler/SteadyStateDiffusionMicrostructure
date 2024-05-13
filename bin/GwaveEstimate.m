function [GwaveFull]=GwaveEstimate(Gwave,GwaveBlock)
%%
%Obtain Gwave
GwaveFull=zeros([1,length(GwaveBlock)*abs(Gwave(end))]);
for k=1:length(Gwave)
    GwaveFull((abs(Gwave(k))-1)*length(GwaveBlock)+1:abs(Gwave(k))*length(GwaveBlock))=GwaveBlock*sign(Gwave(k));
end
