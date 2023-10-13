function [statsoutput, statscounter] = plotparametersforcelltypeclassification(dayindex, spikedatadir, savedfiguresdir,iden, allrecProps, statsoutput, statscounter)
%SP 3.14.19

%% load filename
filename = [spikedatadir 'cellTypeProps_' iden num2str(dayindex(1)) '_' num2str(dayindex(2)) '.mat'];
ac_mean = cell2mat(allrecProps.autocorr_mean);
spikewidth = allrecProps.peak2troughDiff;
FR_mean = allrecProps.meanFR;

%% plot data
figure; clf; hold on;
plot3(ac_mean, spikewidth, FR_mean, 'o')
xlabel('Mean ac (ms)')
ylabel('Spikewidth - peak2trough diff')
zlabel('Mean FR (Hz)')
title(['Mean AC vs. Peak2Trough Diff  vs. Mean FR - F' num2str(dayindex(1)) ' ' num2str(dayindex(2))])
filename = [savedfiguresdir 'Spikewidth_MeanFR_MeanAC\MeanFR_MeanAC_peak2troughDiff' iden num2str(dayindex(1)) '_' num2str(dayindex(2))];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png');
statsoutput{statscounter}.descript = 'cell type classification parameters for all animals';
statsoutput{statscounter}.filename = filename;
statsoutput{statscounter}.testtype = 'nan';
statsoutput{statscounter}.acmean = mean(ac_mean);
statsoutput{statscounter}.spikewidthmean = mean(spikewidth);
statsoutput{statscounter}.FRmean = mean(FR_mean);
statscounter = statscounter + 1;

end