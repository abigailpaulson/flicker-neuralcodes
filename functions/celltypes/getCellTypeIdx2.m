function cellTypeInfo = getCellTypeIdx2(dayindex, WFprops,thresholds, iden, savedfiguresdir, plotdays) 
% this function gets the PYR and IN indices
% uses autocorrelogram centerofmass, spikewidth, and mean FR
% separates based on set thresholds 
%SP 4.27.18

%% get index
%eliminate bad neurons
cellTypeParams = [WFprops.ac_mean; WFprops.spikewidth; WFprops.FR_mean]'; %added invert FR_mean ALP 1/12/2020. removed invert ALP 12/14/22
cellTypeParams(isnan(WFprops.spikewidth),:) = nan;
cellTypeParams(isnan(WFprops.ac_mean),:) = nan;
    
%set parameters
INreq1 = find(WFprops.ac_mean > thresholds.acthresh_IN);
INreq2 = find(WFprops.spikewidth < thresholds.widththresh_IN);
INIdx = INreq1(ismember(INreq1, INreq2));

PYRreq1 = find(WFprops.ac_mean < thresholds.acthresh_PYR); %find(~isnan(WFprops.spikewidth));
PYRreq2 = find(WFprops.spikewidth > thresholds.widththres_PYR);
PYRIdx = PYRreq1(ismember(PYRreq1,PYRreq2));

%% plot resulting separation based on index
%plot the averages)
if plotdays
figure; ax1 = subplot(1,2,1); hold on;
end
alignedWF_IN = nan(length(INIdx),600);
for cellIdx = 1:length(INIdx)
    IN = INIdx(cellIdx);
    negpeak = find(WFprops.wf(IN,:) == min(WFprops.wf(IN,:)));
    alignedWF_IN(cellIdx,200-negpeak:200-negpeak+length(WFprops.wf(IN,:))-1) = WFprops.wf(IN,:);
    if plotdays
    plot(1:600, alignedWF_IN);
    end
end
if plotdays
title(['Interneurons (width < ' num2str(thresholds.widththresh_IN) ' ACmean > ' num2str(thresholds.acthresh_IN) ')'])
ax2 = subplot(1,2,2); hold on
end

alignedWF_PYR = nan(length(PYRIdx),600);
for cellIdx = 1:length(PYRIdx)
    PYR = PYRIdx(cellIdx);
    negpeak = find(WFprops.wf(PYR,:) == min(WFprops.wf(PYR,:)));
    alignedWF_PYR(cellIdx,200-negpeak:200-negpeak+length(WFprops.wf(PYR,:))-1) = WFprops.wf(PYR,:);
    if plotdays
    plot(1:600, alignedWF_PYR);
    end
end
if plotdays
title(['Pyramidal (width > ' num2str(thresholds.widththres_PYR) ' ACmean < ' num2str(thresholds.acthresh_PYR) ')'])

linkaxes([ax1 ax2],'xy')
ylim([-1 1])
%save figures
datadir = [savedfiguresdir 'Classification\'];
savefigALP(datadir, ['IN_PYR_classification_indiv_', iden, num2str(dayindex(1)), '_', num2str(dayindex(2))], 'filetype', 'png');
end

%plot the average traces of each cell type
for cellIdx = 1:length(INIdx)
%     IN = INIdx(cellIdx);
%     INwf(cellIdx,:) = WFprops.wf(IN,:);
%     meanwfIN = mean(INwf);
    meanwfIN = nanmean(alignedWF_IN,1);
end
if isempty(INIdx)
    meanwfIN = nan;
end
    
for cellIdx = 1:length(PYRIdx)
%     PYR = PYRIdx(cellIdx);
%     PYRwf(cellIdx,:) = WFprops.wf(PYR,:);
%     meanwfPYR = mean(PYRwf);
    meanwfPYR = nanmean(alignedWF_PYR,1);
end

if isempty(PYRIdx)
    meanwfPYR = nan;
end

datadir = [savedfiguresdir 'Classification\'];

if plotdays
figure; hold on;
plot(meanwfIN, 'LineWidth', 3)
plot(meanwfPYR, 'LineWidth', 3)
title(['Cell Type Averages (width - ' num2str(thresholds.widththres_PYR) ' ACmean - ' num2str(thresholds.acthresh_IN)])
legend('Interneurons','Pyramidal')
ylim([-1 1])

%% save figures and cell type info
savefigALP(datadir, ['IN_PYR_classification_avg_', iden, num2str(dayindex(1)), '_', num2str(dayindex(2))], 'filetype', 'png');
end

celltype = nan(length(cellTypeParams),1);
celltype(PYRIdx) = 1;
celltype(INIdx) = 2;

%store in output structure
cellTypeInfo.celltype = celltype;
cellTypeInfo.PYRidx = PYRIdx;
cellTypeInfo.INidx = INIdx;
cellTypeInfo.params = cellTypeParams;
cellTypeInfo.alignedWF_IN = alignedWF_IN;
cellTypeInfo.alignedWF_PYR = alignedWF_PYR;

end
