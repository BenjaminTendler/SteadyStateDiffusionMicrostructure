function [GwaveInterp]=WaveformCalculate(Gwave,G,delta,TR,nbins,nOscillations,Waveform)
%%
%Define constants
G=G*10^-5;
%%
%Correct for short pulse approximation
if nbins==1
    G=G*delta/TR;
end
%%
%Generate Full Gradient Waveform
GwaveFull=zeros([1,abs(Gwave(end))]);
GwaveFull(abs(Gwave))=1.*sign(Gwave).*G;
%%
%Interpolate
GwaveInterp=zeros([1,round(length(GwaveFull)*nbins)]);
if nbins==1
    GwaveInterp=GwaveFull;
    if nOscillations>0
        GwaveInterp=repelem(GwaveInterp,2*nOscillations);
        GwaveInterp(2:2:end)=GwaveInterp(2:2:end)*-1;
    end
else
    for k=1:length(GwaveFull)
        GwaveBlock=ones([1,floor(nbins)])*GwaveFull(k);
        GwaveBlock(round(delta/TR*nbins+1):end)=0;
        if nOscillations>0
            GwaveBlock=GwaveBlock*0;
            for l=1:4*nOscillations
                OscRange=(round((delta/TR*nbins)/(4*nOscillations)))*(l-1)+1:(round((delta/TR*nbins)/(4*nOscillations)))*(l);
                if strcmp(Waveform,'Rect')
                    flip=repmat([1,1,-1,-1],[1,nOscillations]);
                    GwaveBlock(OscRange)=GwaveFull(k)*flip(l);
                elseif strcmp(Waveform,'Rect_Sym')
                    flip=repmat([1,-1,-1,1],[1,nOscillations]);
                    GwaveBlock(OscRange)=GwaveFull(k)*flip(l);
                elseif strcmp(Waveform,'Sine')
                    if l==1 || l==4*nOscillations
                        GwaveBlock(OscRange)=sin((OscRange-OscRange(1)+1)/(OscRange(end)-OscRange(1)+1)*pi)*GwaveFull(k);
                    else
                        GwaveBlock(OscRange)=sin((OscRange-OscRange(1)+1)/(OscRange(end)-OscRange(1)+1)*pi/2+l*pi/2)*GwaveFull(k);
                    end                  
                else
                    print('Undefined Waveform for oscillating gradients. Please set opt.Waveform to Rect (conventional rectangular gradient), Rect_Sym (symmetric rectangular gradient) or Sine (Sine gradient)');
                    break;
                end
            end
            GwaveBlock((round((delta/TR*nbins)/(2*nOscillations)))*(l)+1:end)=0;
            %Small correction for imperfect refocusing
            GwaveBlock(end)=GwaveBlock(end)-sum(GwaveBlock);
        end
        GwaveInterp((k-1)*floor(nbins)+1:k*floor(nbins))=GwaveBlock;
    end
end
