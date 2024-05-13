function [Powerspectra,Q]=PowerSpectraEstimate(Gwave,GwaveBlock,Qsize,dt,gamma)
%%
%Obtain Gwave
GwaveFull=zeros([1,length(GwaveBlock)*abs(Gwave(end))]);
for k=1:length(Gwave)
    GwaveFull((abs(Gwave(k))-1)*length(GwaveBlock)+1:abs(Gwave(k))*length(GwaveBlock))=GwaveBlock*sign(Gwave(k));
end
%Estimate Q-vector
Q=cumsum(gamma*GwaveFull)*dt;
Q(end+1:Qsize)=0;
%Estimate power spectrum
Powerspectra=abs(fftshift(fft(fftshift(Q*dt)))).^2;
