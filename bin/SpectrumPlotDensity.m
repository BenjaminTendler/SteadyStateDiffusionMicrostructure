function [] = SpectrumPlotDensity(opt,xlims,ylims,npaths,PathwaySignal,Gwave,GridSize,ColorLims)
%%
%Estimate E1 & E2
E1=exp(-opt.TR/opt.T1);
E2=exp(-opt.TR/opt.T2);
%%
%Convert G to G/cm
opt.G=opt.G/10;
%%
%Gradient Block Calculations
%Define temporal gradient spacing
dt=(opt.TR*10^-3/opt.SpectraResolution);
%Perform gradient block calculation
GwaveBlock=WaveformCalculate(1,opt.G,opt.tau,opt.TR,opt.SpectraResolution,opt.nOscillations,opt.Waveform);
%%
%Obtain indices of ordered pathway amplitudes
[~,idx] = sort(abs(PathwaySignal),'descend');
%%
%Create grid for Power Spectrum
xGrid=[xlims(1)-(xlims(2)-xlims(1))/GridSize:(xlims(2)-xlims(1))/GridSize:xlims(end)+(xlims(2)-xlims(1))/GridSize];
yGrid=[ylims(1)-((ylims(2)-ylims(1)))/GridSize:((ylims(2)-ylims(1)))/GridSize:ylims(end)+((ylims(2)-ylims(1)))/GridSize];
PowerSpectrumGrid=zeros([length(xGrid),length(yGrid)]);
%%
r = randi([0,1E6],1);%Prevents overlapping figures
figure(r)
%%
%Plot power spectrum density plot
for k=1:min(npaths,length(Gwave))
    %Obtain pathway index
    idxPlot=idx(k);
    %Obtain pathway power spectrum
    [Power]=PowerSpectraEstimate(single(Gwave{idxPlot}),GwaveBlock,opt.SpectraResolution*10^4,dt,opt.gamma);
    %Define x-limits
    x=(-length(Power)/2:length(Power)/2-1)/length(Power)/dt;
    xLocation=dsearchn(xGrid',x(x>=xGrid(1) & x<=xGrid(end))');
    yLocation=dsearchn(yGrid',Power(x>=xGrid(1) & x<=xGrid(end))');
    idxGrid = sub2ind(size(PowerSpectrumGrid),xLocation,yLocation);
    PowerSpectrumGrid(idxGrid)=PowerSpectrumGrid(idxGrid)+PathwaySignal(idxPlot);
    figure(r);imagesc(rot90(real(PowerSpectrumGrid(2:end-1,2:end-1)))); c = gray; c = flipud(c);
    %Add colorbar limits
    if nargin==8
        caxis(ColorLims);
    end
    %Define colormap
    colormap(centered('GyBu'));
end
%%
%Format plot
title('DW-SSFP Power Spectrum')
set(gca,'FontSize',18)
xlabel('Frequency (Hz)', 'Interpreter', 'latex','FontSize',22)
ylabel('$|\tilde{q}(\omega)|^2$', 'Interpreter', 'latex','FontSize',22)
if min(xGrid)<0
    xticks([1,find(xGrid==min(abs(xGrid))),length(xGrid)-2])
    xticklabels([xlims(1),0,xlims(2)])
else
    xticks([1,length(xGrid)-2])
    xticklabels([xlims(1),xlims(2)])
end
yticks([1, length(yGrid)-2])
if ylims(2)<100
    yticklabels([ylims(2),ylims(1)])
else
    yticklabels([ylims(2)/10^floor(log10(ylims(2))),ylims(1)/10^floor(log10(ylims(2)))])
    annotation('textbox', [0.11,1,0.5,0], 'string', ['\times10^',num2str(floor(log10(ylims(2))))],'Fontsize',16)
end
%Colorbar
if nargin==8
    cbh = colorbar ;
    cbh.Ticks = [ColorLims(1),0,ColorLims(2)] ;
    cbh.TickLabels = (num2cell([-1,0,1]));
    ylabel(cbh,'Amplitude (normalised)','FontSize',16,'Rotation',90)
end