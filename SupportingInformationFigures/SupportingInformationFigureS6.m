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
%Plot figures
figure;
subplot(1,2,1);bar([2:2:10],real(nTrans(2:2:10))./S0)
xlim([0.5,11.5])
ylim([0,1])
xlabel('no. Transverse Periods');
ylabel('Signal fraction (normalised)')
title(['Signal Fraction Contribution'])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
axis square

hold all;subplot(1,2,2);bar([2:2:10],log10(real((nPathways(2:2:10)))))
xlim([0.5,11.5])
xlabel('no. Transverse periods');
ylabel('no. Pathways (log_{10} scale)')
title(['no. Contributing Pathways'])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
axis square