function cf_getcelltypes(dirs, params, allindex, metadata)
%cf_classify_celltypes
%
%ALP 12/14/22
indices.allindex = allindex;
indices.dayindex = unique(allindex(:,1:2), 'rows');

%% cell types - CA1
params.thresholds.widththresh_IN = 0.5; % long AC - 0.6
params.thresholds.acthresh_IN = 4.5; %8;
params.thresholds.widththres_PYR = 0.5; %short and long AC
params.thresholds.acthresh_PYR = 4.5; %9;

params.brainReg = 'CA1';
dirs.processedappend = [params.brainReg, '\'];
ct{1} = cf_celltypeclassifier(dirs, indices, params);

%% cell types - CA3
params.thresholds.widththresh_IN = 0.5; % long AC - 0.6
params.thresholds.acthresh_IN = 4.5; %8;
params.thresholds.widththres_PYR = 0.5; %short and long AC
params.thresholds.acthresh_PYR = 4.5; %9;

params.brainReg = 'CA3';
dirs.processedappend = [params.brainReg, '\'];
ct{2} = cf_celltypeclassifier(dirs, indices, params);

%% plot cell types together! 
sw_edges = 0:0.05:1.5; ac_edges = 0:0.25:10; fr_edges = 0:0.5:30;

figure
ax1 = axes('Position',[0.25 0.25 0.6 0.6]);
ax1.TickDir = 'out';
ax2 = axes('Position',[0.08 0.25 0.1 0.6]);
ax2.TickDir = 'out';
ax2.YAxisLocation = 'right';
ax2.YTickLabels = [];
ax2.XDir = 'reverse';
ax3 = axes('Position',[0.25 0.08 0.6 0.1]);
ax3.TickDir = 'out';
ax3.XTickLabels = [];
ax3.XAxisLocation = 'top';
ax3.YDir = 'reverse';

axes(ax1)
hold on
plot(ct{1}.all.ac, ct{1}.all.spikewidth, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2)
plot(ct{2}.all.ac, ct{2}.all.spikewidth, 'diamond', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2)
plot(ct{1}.in.ac, ct{1}.in.spikewidth, 'o', 'MarkerFaceColor', '#00008F', 'MarkerEdgeColor', '#00008F',  'MarkerSize', 2)
plot(ct{2}.in.ac, ct{2}.in.spikewidth, 'diamond', 'MarkerFaceColor', '#00008F', 'MarkerEdgeColor', '#00008F','MarkerSize', 2)
plot(ct{1}.pyr.ac, ct{1}.pyr.spikewidth, 'o', 'MarkerFaceColor', '#BF0000', 'MarkerEdgeColor', '#BF0000','MarkerSize', 2)
plot(ct{2}.pyr.ac, ct{2}.pyr.spikewidth, 'diamond', 'MarkerFaceColor', '#BF0000', 'MarkerEdgeColor', '#BF0000','MarkerSize', 2)
xlim([1 8])
ylim([0 1.5])

allPYRsw = [];
allINsw = [];
allPYRac = [];
allINac = [];
for b = 1:2
    tmpdat = ct{b};
    allPYRsw = [allPYRsw tmpdat.pyr.spikewidth];
    allINsw = [allINsw tmpdat.in.spikewidth];
    allPYRac = [allPYRac tmpdat.pyr.ac];
    allINac = [allINac tmpdat.in.ac];
end
axes(ax2)
hold on
h = histogram(allPYRsw, sw_edges);
h.FaceColor = 'r';
h.Orientation = 'horizontal';
h = histogram(allINsw, sw_edges);
h.FaceColor = 'b';
h.Orientation = 'horizontal';
ylim([0 1.5])
title('spike width (ms)')

axes(ax3)
hold on
h = histogram(allPYRac, ac_edges);
h.FaceColor = 'r';
h = histogram(allINac, ac_edges);
h.FaceColor = 'b';
xlim([1 8])
title('mean of the autocorrelogram (ms)')

makefigurepretty(gcf)

figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\celltypes\';
figname = 'celltypeclassification_CA1CA3_postonly';
savefigALP(figdir, figname, 'filetype', 'pdf')

end

