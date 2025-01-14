clear all
%%
%Define Parameters
[opt] = ParameterOptions();
%%
%Generate Pathway Amplitude & Gradient Waveform basis set
%First Dictionary
[PathwaySignalFirstOrder,GwaveFirstOrder,TransFirstOrder] = PathwayOperatorOrder(opt,'First');
%Second Dictionary
[PathwaySignalHigherOrder,GwaveHigherOrder,TransHigherOrder] = PathwayOperatorOrder(opt,'Higher',PathwaySignalFirstOrder);
%%
%Estimate Signal as function of number of TRs and number of pathways
[S0,nTrans,nPathways]=S0Estimate(opt,PathwaySignalFirstOrder,PathwaySignalHigherOrder,GwaveFirstOrder,GwaveHigherOrder,TransFirstOrder,TransHigherOrder);
%%
%Define analytical S0
S0Analytical=FreedDWSSFP(opt.G,opt.tau,opt.TR,opt.alpha,0,opt.T1,opt.T2);
%%
%Plot fraction of signal explained by proposed framework thresholded by a given number of transverse TRs
figure;plot([2:2:opt.nPeriods*2],abs(cumsum(nTrans(2:2:opt.nPeriods*2))./S0Analytical*100),'LineWidth',6)
xlabel('no. Transverse Periods');
ylabel('Percentage of Explained Signal (%)')
title(['Precision (non-Diffusion-Weighted Signal)'])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
xticks([2:2:opt.nPeriods*2])