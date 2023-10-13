function plotparametersforcelltypeclassification_K2(dayindex, savedfiguresdir,iden, allrecProps)
%SP 3.14.19
% updated ALP 1/16/2020

ac_mean = cell2mat(allrecProps.autocorr_mean);
spikewidth = allrecProps.spikewidth.peak2troughDiff;
FR_mean = allrecProps.meanFR;

%% plot data
figure; clf; hold on;
plot3(ac_mean, spikewidth, FR_mean, 'o')
xlabel('Mean ac (ms)')
ylabel('Spikewidth - peak2trough diff')
zlabel('Mean FR (Hz)')
title(['Mean AC vs. Peak2Trough Diff  vs. Mean FR - ', iden, num2str(dayindex(1)) ' ' num2str(dayindex(2))])
filename = ['MeanFR_MeanAC_peak2troughDiff_A'  num2str(dayindex(1,1)), '_', num2str(dayindex(1,2))];
savefigALP([savedfiguresdir, 'Spikewidth_MeanFR_MeanAC\'], filename, 'filetype', 'pdf')
savefigALP([savedfiguresdir, 'Spikewidth_MeanFR_MeanAC\'], filename, 'filetype', 'fig')


end