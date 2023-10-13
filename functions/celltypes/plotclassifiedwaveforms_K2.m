function [cellTypeInfo] = plotclassifiedwaveforms_K2(dayindex,spikedatadir,savedfiguresdir,allrecProps,thresholds,iden, plotdays)
%SP 3.14.19
%updated ALP 1/16/2020
% 

%load up parameter info
WFprops.ac_mean = cell2mat(allrecProps.autocorr_mean);
WFprops.spikewidth = allrecProps.spikewidth.peak2troughDiff;
WFprops.FR_mean = allrecProps.meanFR;
WFprops.wf = allrecProps.spikewidth.fullWFDiff;

%get indices and plot separation
cellTypeInfo = getCellTypeIdx2(dayindex,WFprops,thresholds, iden, savedfiguresdir, plotdays);

end