function [PathwaySignal,PathwayGwave,PathwayTrans] = PathwayGenerateOrder(Signal,State,Gwave,Trans,opt)
%%
%Set initialisation Parameters
ind=1;ind2=1;ind3=0;
start=1;finish=1;
PathwaySignal=[];
PathwayTrans=[];
PathwayGwave={};
GOrder=State{1}(1);
%%
%Estimate E1 & E2
E1=exp(-opt.TR/opt.T1);
E2=exp(-opt.TR/opt.T2);
%%
%Define signal loss for a transverse pathway
T2Precision=E2*(cosd(opt.alpha/2)^2);
%%
%Generate pathways
for k=1:opt.PathwayLength
    [rotMat] = EPG_Operator(opt.alpha,opt.phi,E1,E2,~isfield(opt,'FullWaveform'));
    for l=start:finish
        for m=1:3
            T2PrecisionPower=abs(abs([State{l}(m)])-1);
            %Estimate if pathways will fall below precision threshold
            if abs(Signal{l}(m).*T2Precision.^T2PrecisionPower)<opt.Precision || Trans{l}(m)>=2*opt.nPeriods
                ;
            else
                %Update signal amplitude, k-order, number of transverse periods & gradient waveform
                ind3=ind3+1;
                Signal{ind+1}(1)=Signal{l}(m)*rotMat(1,m);
                Signal{ind+1}(2)=Signal{l}(m)*rotMat(2,m);
                Signal{ind+1}(3)=Signal{l}(m)*rotMat(3,m);
                State{ind+1}(1)=State{l}(m)+1;
                State{ind+1}(2)=State{l}(m)-1;
                State{ind+1}(3)=State{l}(m);
                Trans{ind+1}(1)=Trans{l}(m)+1;
                Trans{ind+1}(2)=Trans{l}(m)+1;
                Trans{ind+1}(3)=Trans{l}(m);
                Gwave{ind+1}{1}=[Gwave{l}{m},(k+GOrder)];
                Gwave{ind+1}{2}=[Gwave{l}{m},(-k-GOrder)];
                Gwave{ind+1}{3}=Gwave{l}{m};
                %Store signal-forming pathways (k-order = 0)
                if State{ind+1}(2)==0
                    if abs(Signal{ind+1}(2))>opt.Precision
                        PathwaySignal(ind2)=[Signal{ind+1}(2)];
                        PathwayTrans(ind2)=[Trans{ind+1}(2)];
                        PathwayGwave{ind2}=Gwave{ind+1}{2};
                        ind2=ind2+1;
                    end
                    %If generating dictionaries, break after signal is formed
                    if ~isfield(opt,'FullWaveform');
                        Signal{ind+1}(2)=0;
                    end
                end
                ind=ind+1;
            end
        end
        %Reset arrays
        Signal{l}={};
        State{l}={};
        Trans{l}={};
        Gwave{l}={};
    end
    if isempty(Signal) == 1
        break
    end
    %Remove evaluated pathways
    Signal={Signal{finish+1:end}};
    State={State{finish+1:end}};
    Trans={Trans{finish+1:end}};
    Gwave={Gwave{finish+1:end}};
    ind=length(Signal);
    finish=length(Signal);
end



