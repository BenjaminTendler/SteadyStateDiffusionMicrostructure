function [MSD,time] = MSDAnalytical(R,opt)
%%
%Implementation based on Eq. 32 in Mortensen et al., Confined Brownian motion tracked with motion blur: estimating diffusion coefficient and size of confining space.
%%
%Convert units 
R=R*10^-3;
D=opt.D*10^-3;
%%
%Define spacing and duration
dt=(opt.TR*10^-3/opt.SpectraResolution);
dur=opt.PathwayLength*opt.TR*10^-3;
%%
%Define first bessel derivitive function
N=100;
dJ1 = @(x) besselj(0,x) - besselj(1,x)./(x + eps);
alpha1 = @(k) fzero(dJ1, [(k-1)*pi, k*pi]);
%%
%Define simulation time (in s) and coefficients
time=0:dt:dur;
tau=R.^2/(D);
%%
%Estimate MSD
MSD=zeros([1,length(time)]);
for n=1:length(time)
    for l=1:length(N)
        MSD(n)=MSD(n)+R.^2*(1-8*exp(-alpha1(l)^2*time(n)/tau)*1/(alpha1(l).^2*(alpha1(l).^2-1)))/2;
    end
end
%%
%Set first point equal to zero
MSD(1)=0;
