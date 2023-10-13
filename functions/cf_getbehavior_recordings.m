function cf_getbehavior_recordings(dirs, params, allindex, metadata)
%cf_getbehavior_recordings
%   behavior analysis in the annular track during recording sessions
%   note - recording sessions only! 
%   make sure to run PreProcessing for TrialData structure before this
%   script
%ALP 12/7/22

%%% PARAMS
params.pos_edges = -81:3:99; 
params.RZ = [0 18];
params.AZ = [-18+9 0];
RZbins = [find(params.pos_edges == params.RZ(1)) find(params.pos_edges == params.RZ(2))-1];
AZbins = [find(params.pos_edges == params.AZ(1)) find(params.pos_edges == params.AZ(2))-1];

%%% helpful vectors for groups, etc
isGamma = strcmp(metadata.Groups, 'gamma'); 
isRandom = strcmp(metadata.Groups, 'random');
isPre = metadata.FlickerDay < 5;
isPost = metadata.FlickerDay > 5; 
dayID = 1:height(metadata);
gID.gamma.pre = dayID(isGamma & isPre);
gID.gamma.post = dayID(isGamma & isPost);
gID.random.pre = dayID(isRandom & isPre);
gID.random.post = dayID(isRandom & isPost);
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'}; 

%%% intitialize, directories, etc. 
positiondir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\fig_behavior\';
dayindex = unique(allindex(:,1:2), 'rows');
allTData = []; allTMetrics = [];

%%%%% ----- calculate per trial behavior metrics ----- %%%%%
for d = 1:size(dayindex,1)
    trialfilename = ['trialInfo_A',  num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    load([positiondir, trialfilename])
    metrics = [];
    
    %%% analyze trials
    for t = 1:length(trialData)
        time = trialData(t).time; 
        licks = logical(trialData(t).licks);
        theta = trialData(t).theta_d2r;
        rewards = logical(trialData(t).rewards); 
        speed = trialData(t).speed; 
        lickpos = theta(licks); 
        rewardpos = theta(rewards); 
        
        lick_h = histcounts(lickpos, params.pos_edges); 
        [~, ~, theta_bin] = histcounts(theta, params.pos_edges); 
        vel_h = arrayfun(@(x) mean(speed(theta_bin == x), 'omitnan'), 1:length(params.pos_edges)-1, 'UniformOutput', false); 
        vel_h = cell2mat(vel_h); 
        
        inRZ = theta >= params.RZ(1) & theta < params.RZ(2);
        inAZ = theta >= params.AZ(1) & theta < params.AZ(2); 
        inCtrl = ~inRZ & ~inAZ; 
        
        RZtime = time(inRZ);
        if isempty(RZtime)
            RZtime = NaN;
        end
        licktime = time(licks); 
        RZlicktime = licktime(lickpos >= params.RZ(1) & lickpos < params.RZ(2));
        relativelicktime = RZlicktime - RZtime(1); 
        lickDI = (sum(licks(inRZ|inAZ))-sum(licks(inCtrl)))/(sum(licks)); 
        
        licklatency_s = min(relativelicktime(relativelicktime > 0));
        licklatency_pos = min(lickpos(lickpos > 0 & lickpos < 18)); 
        
        trial_speed = mean(speed, 'omitnan'); 
        AZ_speed = mean(speed(inAZ), 'omitnan');
        RZ_speed = mean(speed(inRZ), 'omitnan');
        Ctrl_speed = mean(speed(inCtrl), 'omitnan'); 
        
        speed_slope = (vel_h(AZbins(2)) - vel_h(AZbins(1)))/(AZbins(2)-AZbins(1));
        lick_slope = (lick_h(AZbins(2)) - lick_h(AZbins(1)))/(AZbins(2)-AZbins(1));
        
        metrics(t).lick_h = lick_h';
        metrics(t).vel_h = vel_h';
        metrics(t).lick_nh = lick_h'./sum(lick_h);
        metrics(t).vel_nh = vel_h'./max(vel_h);
        metrics(t).lickDI = lickDI;
        metrics(t).licklatency_s = licklatency_s;
        metrics(t).licklatency_pos = licklatency_pos;
        metrics(t).trial_speed = trial_speed;
        metrics(t).AZ_speed = AZ_speed;
        metrics(t).RZ_speed = RZ_speed;
        metrics(t).Ctrl_speed = Ctrl_speed;
        metrics(t).AZ_speed_slope = speed_slope;
        metrics(t).AZ_lick_slope = lick_slope;
        metrics(t).day = d; 
    end
    %%% save behavior metrics structure for each day
    
    %%% append to large structures
    allTMetrics = [allTMetrics metrics];
    allTData = [allTData trialData];
end

units = {'lick count', 'raw speed (deg/s)', 'fraction of licks', 'norm speed (deg/s)', 'DI', 'latency (s)',...
    'latency (deg)', 'avg speed (deg/s)', 'avg speed (deg/s)', 'avg speed (deg/s)',...
    'avg speed (deg/s)', 'slope', 'slope'}; 

%%%%% ----- split data by zones ----- %%%%%
includeTrials = [allTData.fullTrial] & [allTData.engaged] & [allTData.rewarded];
fnames = fieldnames(allTMetrics);

for f = 1:length(fnames)
    for z = 1:2
        isZone = [allTData.zone] == z;
        ZD(1,z).(fnames{f}) = [allTMetrics(includeTrials & isZone).(fnames{f})];
    end
end

for f = 1:length(fnames)-1
    for z = 1:2
        tmpdat = ZD(1,z).(fnames{f});
        day = ZD(1,z).(fnames{end});
        if f > 4
            day_avg = arrayfun(@(x) nanmean(tmpdat(day == x)), 1:max(day));
        else
            day_avg = arrayfun(@(x) nanmean(tmpdat(:,day == x),2), 1:max(day), 'UniformOutput', false);
            day_avg = cell2mat(day_avg);
        end
        dayZD(1,z).(fnames{f}) = day_avg;
    end
end

%%%%% ----- split data by groups ----- %%%%%
for g = 1:numel(gnames)
    for d = 1:numel(dnames)
        isDay = gID.(gnames{g}).(dnames{d});
        allDays = [allTMetrics.(fnames{end})]; 
        includeDays = ismember(allDays, isDay); 
        for f = 1:length(fnames)
            GD.(gnames{g}).(dnames{d}).(fnames{f}) = [allTMetrics(includeTrials & includeDays).(fnames{f})];
        end
    end
end

for g = 1:numel(gnames)
    for d = 1:numel(dnames)
        tmpdays = GD.(gnames{g}).(dnames{d}).(fnames{end});
        tmpID = unique(tmpdays); 
        for f = 1:length(fnames)-1
            day_avg = [];
            tmpdat = GD.(gnames{g}).(dnames{d}).(fnames{f});
            if f > 4
                day_avg = arrayfun(@(x) nanmean(tmpdat(tmpdays == x)), tmpID);
            else
                day_avg = arrayfun(@(x) nanmean(tmpdat(:,tmpdays == x),2), tmpID, 'UniformOutput', false);
                day_avg = cell2mat(day_avg);
            end
            dayGD.(gnames{g}).(dnames{d}).(fnames{f}) = day_avg;
        end
    end
end

%%%%% ----- SAVE ----- %%%%%
savedatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\']; 
info = [];
info = addhelpfulinfotostruct(info);
info.gnames = gnames;
info.dnames = dnames;
info.fnames = fnames;
info.params = params;
info.units = units;
savefilename = 'behavioranalysis_prepost_recordings.mat';
save([savedatadir, savefilename], 'ZD', 'dayZD', 'GD', 'dayGD', 'allTMetrics', 'allTData', 'info')


%%%%% ----- PLOT show no diff in behavior between zones ----- %%%%%
halfcolors = params.colors.behavior.half;

%%% per trial averages
figure
hold on
tiledlayout(1,length(fnames)-2, 'TileSpacing', 'Compact', 'Padding', 'Compact'); 
hold on
for f = 5:length(fnames)-1
    nexttile
    hold on
    vdat = [];
    nm = {'zone1', 'zone2'};
    for z = 1:2
        vdat.(nm{z}) = ZD(1,z).(fnames{f});
    end
    violinplot_half(vdat, [], 'ViolinColorMat', halfcolors, 'ShowData', false);
    ylabel(units{f})
    title(fnames{f})
end
figname = 'behavior_trials_byzone';
savefigALP(figdir, figname, 'info', 'filetype', 'pdf')

%%% per day averages, probably go into supplement
figure
hold on
tiledlayout(1,length(fnames)-2, 'TileSpacing', 'Compact', 'Padding', 'Compact'); 
hold on
for f = 5:length(fnames)-1
    nexttile
    hold on
    vdat = [];
    nm = {'zone1', 'zone2'};
    for z = 1:2
        vdat.(nm{z}) = dayZD(1,z).(fnames{f}); 
    end
    violinplot_half(vdat, [], 'ViolinColorMat', halfcolors);
    ylabel(units{f})
    title(fnames{f})
end

figname = 'behavior_day_byzone';
savefigALP(figdir, figname, 'info', 'filetype', 'pdf')

f1 = figure;
hold on
ploti = 1;
for f = 1:4
    subplot(1,4,f)
    hold on
    for z=1:2
        tmpdat = dayZD(1,z).(fnames{f})';
        mn = nanmean(tmpdat,1);
        stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
        shadedErrorBar(params.pos_edges(1:end-1), mn, stde, {'Color', halfcolors(z,:)},1)
    end
    ylabel(units{f})
    xlabel('position (deg')
    xlim([-81 99])
    title(fnames{f})
end
figname = 'behavior_histograms_day_byzone';
savefigALP(figdir, figname, 'info', 'filetype', 'pdf')

%%%%% ----- PLOT show no change in behavior between groups ----- %%%%%
%show for day 0/1 and day 9/10 (likely 0/1 will go in the supplement but
%I'm not positive. Zones combined bc no differences
groupcolors = params.colors;
xvect = [1:2; 4:5];

figure
hold on
tiledlayout(1,length(fnames)-2, 'TileSpacing', 'Compact', 'Padding', 'Compact'); 
hold on
for f = 5:length(fnames)-1
    nexttile
    hold on
    vdat = []; cmat = [];
    for d = 1:2
        for g = 1:2
            nm = [gnames{g}, dnames{d}];
            vdat.(nm) = GD.(gnames{g}).(dnames{d}).(fnames{f}); 
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
    end
    vdat.xlabels = {'pre', 'post'};
    violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false);
    ylabel(units{f})
    title(fnames{f})
end
figname = 'behavior_trial_bygroup';
savefigALP(figdir, figname, 'info', 'filetype', 'pdf')

f1 = figure;
hold on
tiledlayout(1,length(fnames)-2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
hold on
for f = 5:length(fnames)-1
    nexttile
    hold on
    for d = 1:length(dnames)
        bdat = []; bstde = []; cmat = []; scatterdat = [];
        for g = 1:length(gnames)
            tmpdat = [];
            tmpdat = dayGD.(gnames{g}).(dnames{d}).(fnames{f});
            scatterdat{g} = tmpdat;
            bdat(g) = nanmean(tmpdat);
            bstde(g) = nanstd(tmpdat)./sqrt(sum(~isnan(tmpdat)));
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
        plotprettypoints(f1, xvect(d,:), scatterdat, 'color', cmat)
        b = bar(xvect(d,:), bdat, 'FaceColor', 'flat');
        b.CData = cmat;
        b.FaceAlpha = 0.6;
        errorbar2(xvect(d,:), bdat, bstde, 0.2, 'k-', 'LineWidth', 0.75);
        
    end
    ylabel(units{f})
    title(fnames{f})
end
figname = 'behavior_day_bygroup';
savefigALP(figdir, figname, 'info', 'filetype', 'pdf')

%%% plot avg licking and velocity histograms
f1 = figure('Position', [148 160 1400 608]);
hold on
ploti = 1;
for f = 1:4
    for d = 1:length(dnames)
        subplot(2,4,ploti)
        hold on
        for g = 1:length(gnames)
            tmpdat = dayGD.(gnames{g}).(dnames{d}).(fnames{f})';
            mn = nanmean(tmpdat,1);
            stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
            shadedErrorBar(params.pos_edges(1:end-1), mn, stde, {'Color', groupcolors.(gnames{g}).(dnames{d})},1)
        end
        ploti = ploti+1;
        ylabel(units{f})
        xlabel('position (deg')
        xlim([-81 99])
        xticks([-81 0 99])
        title(fnames{f})
    end
end
figname = 'behavior_histograms_day_bygroup';
makefigurepretty(f1)
savefigALP(figdir, figname, 'info', 'filetype', 'pdf')


% 
% %%% plots for figure 1
% figure(fh)
% axes(ax{1})
% hold on
% for f = 9
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
% axes(ax{2})
% hold on
% for f = 2
%     for d = 2
%         for g = 1:length(gnames)
%             tmpdat = dayGD.(gnames{g}).(dnames{d}).(fnames{f})';
%             mn = nanmean(tmpdat,1);
%             stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
%             shadedErrorBar(params.pos_edges(1:end-1), mn, stde, {'Color', groupcolors.(gnames{g}).(dnames{d})},1)
%         end
%         ylabel(units{f})
%         xlabel('distance from reward zone (deg)')
%         xlim([-81 99])
%         title(fnames{f})
%     end
% end
% 
% figure(fh)
% axes(ax{3})
% hold on
% for f = 3
%     for d = 2
%         for g = 1:length(gnames)
%             tmpdat = dayGD.(gnames{g}).(dnames{d}).(fnames{f})';
%             mn = nanmean(tmpdat,1);
%             stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
%             shadedErrorBar(params.pos_edges(1:end-1), mn, stde, {'Color', groupcolors.(gnames{g}).(dnames{d})},1)
%         end
%         ylabel(units{f})
%         xlabel('distance from reward zone (deg)')
%         xlim([-81 99])
%         title(fnames{f})
%     end
% end


end

