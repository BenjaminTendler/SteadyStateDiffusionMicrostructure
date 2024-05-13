function [rotMat] = EPG_Operator(alpha,phi,E1,E2,LongitudinalState)
%%
%Define components
TtoT=E2.*(cosd(alpha/2)^2);
T180toT180=E2.*(cosd(alpha/2)^2);
LtoL=E1.*cosd(alpha);
T180toT=exp(1i*2*phi)*E2.*(sind(alpha/2)^2);
LtoT=-1i*exp(1i*phi).*E2.*sind(alpha);
TtoT180=exp(-1i*2*phi)*E2.*(sind(alpha/2)^2);
LtoT180=1i*exp(-1i*phi).*E2.*sind(alpha);
TtoL=-1i./2.*exp(-1i*phi).*E1.*sind(alpha);
T180toL=1i./2*exp(1i*phi).*E1.*sind(alpha);
%%
%Perform summation of longitudinal states
if LongitudinalState==1;
    TtoL=TtoL;
    T180toL=T180toL;
    LtoL=0;
end
%%
%Define matrix
rotMat=[TtoT,T180toT,LtoT;TtoT180,T180toT180,LtoT180;TtoL,T180toL,LtoL];