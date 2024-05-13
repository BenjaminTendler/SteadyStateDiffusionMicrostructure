function [G,Signal,SignalAtt] = GwaveAttenuation(opt,GwaveInit,Gwave,n,PathwaySignal,E1,GwaveBlock,dt,MSD,GwaveFirstOrder)
%%
%Define spacing of gradient elements
Spacing=floor(length(GwaveInit)/max(abs(Gwave)));
%%
%Determine number of longitudinal periods
nLongPeriods=sum(diff(abs(Gwave))-1);
%Identify starting point
LongPeriodsLocation=abs(Gwave(logical(diff(abs(Gwave))-1)));
%%
if nLongPeriods==0
    G{1}=GwaveInit;
    Signal{1}=PathwaySignal;
    GwFull=repelem([GwaveFirstOrder,G{1}],length(GwaveBlock)).*repmat(GwaveBlock,[1,length([GwaveFirstOrder,G{1}])]);
    GAutocorr=GradientAutocorr(GwFull,opt.G,dt);
    SignalAtt{1}=Signal{1}.*exp(sum(GAutocorr.*MSD(1:length(GAutocorr),:)',2)/2)';
else
    %%
    %Perform recursive gradient operation
    ind=1;
    SignalHigherOrder={};
    GHigherOrder={};
    SignalAttHigherOrder={};
    BreakCond=0;
    for l=1:opt.PathwayLength;
        if BreakCond==1;
            break
        end
        if LongPeriodsLocation(n)==0
            Gpath=[zeros([1,Spacing*(l-1)]),GwaveInit];
        else
            Gpath=[GwaveInit(1:LongPeriodsLocation(n)*Spacing),zeros([1,Spacing*l]),GwaveInit(LongPeriodsLocation(n)*Spacing+Spacing+1:end)];
        end
        G{ind}=Gpath;
        if l>1
            Signal{ind}=Signal{ind-1}*E1*cosd(opt.alpha);
           
        else
            Signal{ind}=PathwaySignal;
        end
        GwFull=repelem([GwaveFirstOrder,Gpath],length(GwaveBlock)).*repmat(GwaveBlock,[1,length([GwaveFirstOrder,Gpath])]);
        GAutocorr=GradientAutocorr(GwFull,opt.G,dt);
        Att(ind)=exp(sum(GAutocorr.*MSD(1:length(GAutocorr),:)',2)/2)';
        SignalAtt{ind}=Signal{ind}.*Att(ind);
        if abs(SignalAtt{ind})<opt.Precision
            break
        end
        ind=ind+1;
        if n<nLongPeriods
            if l>1
                Gwave=abs(Gwave)+1;
            end
            [GHigherOrder{ind},SignalHigherOrder{ind},SignalAttHigherOrder{ind}]=GwaveAttenuation(opt,Gpath,Gwave,n+1,Signal{ind-1},E1,GwaveBlock,dt,MSD,GwaveFirstOrder);
        end

    end
    %Remove first pathway (to prevent repeats)
    if n~=1
        G={G{2:end}};
        Signal={Signal{2:end}};
        SignalAtt={SignalAtt{2:end}};
    end
    G=[G,GHigherOrder{:}];
    Signal=[Signal,SignalHigherOrder{:}];
    SignalAtt=[SignalAtt,SignalAttHigherOrder{:}];
end

