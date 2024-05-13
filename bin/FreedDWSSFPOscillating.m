function [S] = FreedDWSSFPOscillating(TR,alpha,D,T1,T2,b)
%
% calculates steady-state magnetization with diffusion weighting
% for oscillating diffusion gradient based on the Freed model. (Appendix 2)
%
% input:
%	G   = strength of diffusion gradient (mT/m)
%	TR  = repetition time (ms)
%	alpha = flip angle (degrees)
%	D   = diffusion coefficient (um2/ms)
%	T1  = longitudinal relaxation time (ms)
%	T2  = transverse relaxation time (ms)
%   b   = input b-value (for oscillating gradient per TR)
%
% output:
%	S = Measured signal
%%
%Define parameters
TR = TR*10^-3;          % Convert to s
T1 = T1*10^-3;          % Convert to s
T2 = T2*10^-3;          % Convert to s
alpha = alpha*pi/180;   % Convert to radians
D=D*10^-3;              % Convert to mm2/s
b=b*10^3;               % Convert to s/mm2
%%
%Implementation based on Freed et al: Steady-state free precession experiments and exact treatment of diffusion in a uniform gradient
%J. Chem. Phys. 115, 4249 (2001); https://doi.org/10.1063/1.1389859
cosa = cos(alpha);
sina = sin(alpha);
E1p = @(p)(exp(-TR/T1));
E2p = @(p)(exp(-TR/T2-D*b));
Ap=@(p)(1/2.*(E1p(p)-1).*(1+cosa));
Bp=@(p)(1/2.*(E1p(p)+1).*(1-cosa));
Cp=@(p)(E1p(p)-cosa);
np=@(p)(-E2p(-p).*E2p(p-1).*(Ap(p).^2).*Bp(p-1)./Bp(p));
dp=@(p)((Ap(p)-Bp(p))+E2p(-p-1).*E2p(p).*Bp(p).*Cp(p+1)./Bp(p+1));
ep=@(p)(-E2p(p).*E2p(-p-1).*Bp(p).*Cp(p+1)./Bp(p+1));
%%
%Perform recursive operation
x1=0;
for k=10000:-1:1
    if k==10000
        x1=np(k)./(dp(k)+ep(k));
    else
        x1=np(k)./(dp(k)+x1);
    end
end
%%
%Estimate r1
r1=x1./(E2p(-1).*Bp(0))+(E2p(0).*Cp(1))./Bp(1);
%%
%Calculate Signal
S=abs(r1.*sina.*(1-E1p(0)).*E2p(-1)./(Ap(0)-Bp(0)+E2p(-1).*Cp(0).*r1));