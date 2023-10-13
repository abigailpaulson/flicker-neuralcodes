function cf_plot_behaviorperformance_zones_table(dirs, params, allindex, metadata, fh, ax, statfid, panelL, tablefilename)
%cf_plot_behaviorperformance_zones
%   plot behavior performance for zone 1 vs. zone 2
%   also to do group 1 and group 2 in distance to reward
%ALP 12/14/2022

%% load datastructures
savedatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\']; 
savefilename = 'behavioranalysis_prepost_recordings_table.mat';
savefilename2 = 'behavioranalysis_prepost_recordings.mat';
load([savedatadir, savefilename])
originaldata = load([savedatadir, savefilename2]);

info = originaldata.info;

%%% get group information from the info structure
gnames = info.gnames;
dnames = info.dnames;
fnames = info.fnames; 
pos_edges = info.params.pos_edges;
units = info.units;
halfcolors = [params.colors.behavior.half];
groupcolors = [params.colors];
xvect = [1:2; 4:5];

%%% get data to include
isPost = strcmp(allTMetrics.timepoint, 'post');
isCorrect = allTMetrics.engaged == 1 & allTMetrics.fullTrial == 1 & allTMetrics.rewarded ==1;

behaviorData = allTMetrics(isPost & isCorrect,:); 
dayZoneData = groupsummary(behaviorData(:,1:end-2), {'day', 'zone'}, 'mean');


%% plot
%%% plots for supplement figure 1
%%% plots for figure 1

figure(fh)
axes(ax{1})
hold on
fname = 'vel_h'; 
patch([0 18 18 0], [0 0 5 5], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
patch([0 -18 -18 0], [0 0 5 5], 'c', 'FaceAlpha', 0.1, 'EdgeColor', 'none')

for z = 1:2
    isZone = dayZoneData.zone == z;
    tmpdat = dayZoneData(isZone,:).(['mean_', fname]);
    mn = nanmean(tmpdat,1);
    stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
    shadedErrorBar(pos_edges(1:end-1), mn, stde, {'Color', halfcolors(z,:)}, 1)
    N(z) = size(tmpdat,1);
end
xlabel('distance to reward zone (deg)')
ylabel('speed (deg/s)')
yticks([0 5 10])
ylim([0 12])
xticks([-81 0  99])
xlim([-81 99])
title(num2str(N))

figure(fh)
axes(ax{2})
hold on
fname = 'lick_nh';
patch([0 18 18 0], [0 0 0.1 0.1], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
patch([0 -18 -18 0], [0 0 0.1 0.1], 'c', 'FaceAlpha', 0.1, 'EdgeColor', 'none')

for z = 1:2
    isZone = dayZoneData.zone == z;
    tmpdat = dayZoneData(isZone,:).(['mean_', fname]);
    mn = nanmean(tmpdat,1);
    stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
    shadedErrorBar(pos_edges(1:end-1), mn, stde, {'Color', halfcolors(z,:)}, 1)
    N(z) = size(tmpdat,1);
end
xlabel('distance to reward zone (deg)')
ylabel('fraction of licks')
xticks([-81 0 99])
xlim([-81 99])
title(num2str(N))

%% group data

dayGroupData = groupsummary(behaviorData(:,1:end-1), {'day', 'group'}, 'mean');

figure(fh)
axes(ax{3})
hold on
fname = 'vel_h'; 
patch([0 18 18 0], [0 0 5 5], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
patch([0 -18 -18 0], [0 0 5 5], 'c', 'FaceAlpha', 0.1, 'EdgeColor', 'none')

for d = 2
    for g = 1:length(gnames)
        tmpdat = [];
        isGroup = strcmp(dayGroupData.group, gnames{g});
        tmpdat = dayGroupData(isGroup,:).(['mean_', fname]);
        mn = nanmean(tmpdat,1);
        stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
        shadedErrorBar(pos_edges(1:end-1), mn, stde, {'Color', groupcolors.(gnames{g}).(dnames{d})}, 1)
        N(g) = size(tmpdat,1);
    end
end
xlabel('distance to reward zone (deg)')
ylabel('speed (deg/s)')
yticks([0 5 10])
ylim([0 12])
xticks([-81 0  99])
xlim([-81 99])
title(['speed ' num2str(N)])


figure(fh)
iPanel = 4;
axes(ax{iPanel})
hold on
fname = 'trial_speed'; %trial speed
for d = 2
    scatterdat = []; bdat = []; bstde = []; cmat = [];
    for g = 1:length(gnames)
        tmpdat = [];
        isGroup = strcmp(dayGroupData.group, gnames{g});
        tmpdat = dayGroupData(isGroup,:).(['mean_', fname]);
        scatterdat{g} = tmpdat;
        bdat(g) = nanmean(tmpdat);
        bstde(g) = nanstd(tmpdat)./sqrt(sum(~isnan(tmpdat)));
        cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
    end
    plotprettypoints(fh, xvect(d,:), scatterdat, 'color', cmat)
    b = bar(xvect(d,:), bdat, 'FaceColor', 'flat');
    b.CData = cmat;
    b.FaceAlpha = 0.6;
    errorbar2(xvect(d,:), bdat, bstde, 0.2, 'k-', 'LineWidth', 0.75);
    
end
ylabel('Speed (deg/s)')
xticks([])
title('trial speed ')
cf_stats2txt2(scatterdat, statfid, panelL{iPanel}, 'days', 'deg/s', 'speed', gnames, tablefilename)



figure(fh)
axes(ax{5})
hold on
hold on
fname = 'lick_nh';
patch([0 18 18 0], [0 0 5 5], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
patch([0 -18 -18 0], [0 0 5 5], 'c', 'FaceAlpha', 0.1, 'EdgeColor', 'none')

for d = 2
    for g = 1:length(gnames)
        tmpdat = [];
        isGroup = strcmp(dayGroupData.group, gnames{g});
        tmpdat = dayGroupData(isGroup,:).(['mean_', fname]);
        mn = nanmean(tmpdat,1);
        stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
        shadedErrorBar(pos_edges(1:end-1), mn, stde, {'Color', groupcolors.(gnames{g}).(dnames{d})}, 1)
        N(g) = size(tmpdat,1);
    end
end
xlabel('distance to reward zone (deg)')
ylabel('fraction of licks')
ylim([0 0.25])
xlim([-81 99])
xticks([-81 0  99])
title(['licking ', num2str(N)])

figure(fh)
iPanel = 6;
axes(ax{iPanel})
hold on
fname = 'licklatency_s'; %lick latency
for d = 2
    scatterdat = []; bdat = []; bstde = []; cmat = [];
    for g = 1:length(gnames)
        tmpdat = [];
       isGroup = strcmp(dayGroupData.group, gnames{g});
        tmpdat = dayGroupData(isGroup,:).(['mean_', fname]);
        scatterdat{g} = tmpdat;
        bdat(g) = nanmean(tmpdat);
        bstde(g) = nanstd(tmpdat)./sqrt(sum(~isnan(tmpdat)));
        cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
    end
    plotprettypoints(fh, xvect(d,:), scatterdat, 'color', cmat)
    b = bar(xvect(d,:), bdat, 'FaceColor', 'flat');
    b.CData = cmat;
    b.FaceAlpha = 0.6;
    errorbar2(xvect(d,:), bdat, bstde, 0.2, 'k-', 'LineWidth', 0.75);
end
xticks([])
ylabel('lick latency (s)')
title('time to first lick')
cf_stats2txt2(scatterdat, statfid, panelL{iPanel}, 'days', 'seconds', 'lick latency', gnames, tablefilename)


% figure(fh)
% axes(ax{1})
% hold on
% for f = 8
%     vdat = [];
%     nm = {'zone1', 'zone2'};
%     for z = 1:2
%         vdat.(nm{z}) = dayZD(1,z).(fnames{f}); 
%     end
%     violinplot_half(vdat, [], 'ViolinColorMat', halfcolors);
%     ylabel(units{f})
%     title(fnames{f})
% end
% figure(fh)
% axes(ax{2})
% hold on
% for f = 12
%     vdat = [];
%     nm = {'zone1', 'zone2'};
%     for z = 1:2
%         vdat.(nm{z}) = dayZD(1,z).(fnames{f}); 
%     end
%     violinplot_half(vdat, [], 'ViolinColorMat', halfcolors);
%     ylabel(units{f})
%     title(fnames{f})
% end
% 
% figure(fh)
% axes(ax{3})
% hold on
% for f = 12
%     vdat = [];
%     nm = {'zone1', 'zone2'};
%     for z = 1:2
%         vdat.(nm{z}) = dayZD(1,z).(fnames{f}); 
%     end
%     violinplot_half(vdat, [], 'ViolinColorMat', halfcolors);
%     ylabel(units{f})
%     title(fnames{f})
% end

%% save data for R statistics




end

