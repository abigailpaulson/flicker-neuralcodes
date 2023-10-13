function cf_plot_rippleExample(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%cf_plot_rippleExample
%
%ALP 4/14/2023

%%% define directories for examples to plot
index = [28 200831 3];
f = index(3);

%%% plot


%%% directories
anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(index(1,1)), '_', ...
    num2str(index(1,2)), '\'];

%%% load spiking data. use old PC data for now but I should update...
[singleunits, clusterMetrics, spikeVect] = cf_getspikes(anprocesseddatadir, dirs, params, ...
    index, params.brainReg);
brainReg = {clusterMetrics.brainReg};

%%% load ephys data
bestchdir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\';
load([bestchdir, 'CA1\bestChannel_', num2str(index(1,2)), '.mat'])
ca1ch = bestRippleChan.all; 
ca3ch = bestRippleChan.all; 

ca1 = load([anprocesseddatadir, 'CA1\', num2str(ca1ch), '\eeg', num2str(index(3)), '.mat']);
load([anprocesseddatadir, 'CA1\', num2str(ca1ch), '\ripples', num2str(index(3)), '.mat'])
LFPData = ca1.eeg{index(1)}{index(2)}{index(3)}.data;
LFPtime = 0:1/2000:length(LFPData(1,:))/2000;
RipTs = ripples{index(1)}{index(2)}{index(3)}.samprate;

rippleMid = ripples{index(1)}{index(2)}{index(3)}.midind;
rippleWindows = [rippleMid-0.125*RipTs rippleMid+0.125*RipTs]; % in s, on either side of mid


%%% load rate maps for sorting the spiking
load(['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\spatialmaps_all_500_full_', num2str(index(1,2)), '.mat'])
peakpos = ratemaps.peakpos; 
isPC = ratemaps.includeCell;
PC_peakpos = peakpos(isPC);
PCID = ratemaps.cellIDs(isPC);
[~, iSort] = sort(PC_peakpos);
PCID_sorted = PCID(iSort);

figure(fh)
axes(ax{1})
hold on
% cla

tTime = 0.5; %inS
for t = 1 %first ripple of this day is the example that I want
    plotTime = rippleWindows(t,1):rippleWindows(t,2);
    plotTime = plotTime./RipTs;
    
    %%% get spikes in this period
    isPerSpike = logical(isExcluded(spikeVect(f).spikeTimes, rippleWindows(t,:)./RipTs));
    perSpikeTimes = spikeVect(f).spikeTimes(isPerSpike);
    perSpikeIDs = spikeVect(f).spikeIDs(isPerSpike);
    isPerPC = ismember(perSpikeIDs, PCID);
    PCperSpikeIDs = perSpikeIDs(isPerPC);
    PCperSpikeTimes = perSpikeTimes(isPerPC);
    
    iPlot = 1;
    for c = 1:length(PCID_sorted)
        isCell = PCperSpikeIDs == PCID_sorted(c);
        nSpikes = sum(isCell);
        if nSpikes == 0
            continue
        end
        yVect = iPlot*ones(1,nSpikes);
        plot_raster(PCperSpikeTimes(isCell), yVect)
        iPlot = iPlot+1;
    end
    LFPscale = 30;
    iPlot = 100;
    %     plot(plotTime, iPlot+LFPFilt(timeAroundMid(t,1):timeAroundMid(t,2))./LFPscale)
    %     %scale bar
    %     plot([plotTime(50) plotTime(50)], [iPlot iPlot+100/LFPscale], 'r')
    %     iPlot = 50;
    plot(plotTime, iPlot+LFPData(rippleWindows(t,1):rippleWindows(t,2))./LFPscale, 'k')
    %     %scale bar
    plot([plotTime(50) plotTime(50)], [iPlot iPlot+100/LFPscale], 'r')
    %     xlim([plotTime(1) plotTime(end)])
    plot([plotTime(100) plotTime(300)], [50 50], 'r')
    xticks([])
    yticks([])
    title(num2str(rippleMid(t)./RipTs))
end



end

