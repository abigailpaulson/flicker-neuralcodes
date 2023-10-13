function cf_plot_behaviorperformance_groups(dirs, params, allindex, metadata, fh, ax)
%cf_plot_behaviorperformance_groups
%   plot behavior performance for 40Hz vs. random
%ALP 12/14/2022

%% load datastructures
savedatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\']; 
savefilename = 'behavioranalysis_prepost_recordings.mat';
load([savedatadir, savefilename])

%%% get group information from the info structure
gnames = info.gnames;
dnames = info.dnames;
fnames = info.fnames; 
pos_edges = info.params.pos_edges;
units = info.units;
groupcolors = params.colors;
xvect = [1:2; 4:5];

%% plot
%%% plots for figure 1
figure(fh)
axes(ax{1})
hold on
for f = 2
    for d = 2
        for g = 1:length(gnames)
            tmpdat = dayGD.(gnames{g}).(dnames{d}).(fnames{f})';
            mn = nanmean(tmpdat,1);
            stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
            shadedErrorBar(pos_edges(1:end-1), mn, stde, {'Color', groupcolors.(gnames{g}).(dnames{d})},1)
        end
        ylabel(units{f})
        xlabel('distance from reward zone (deg)')
        xlim([-81 99])
    end
end

figure(fh)
axes(ax{2})
hold on
cmat = [];
for f = 12
    for d = 2
        for g = 1:length(gnames)
            tmpdat = [];
            tmpdat = dayGD.(gnames{g}).(dnames{d}).(fnames{f});
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
    ylabel('anticipatory velocity slope')
    xticks([])
end

figure(fh)
axes(ax{3})
hold on
for f = 3
    for d = 2
        for g = 1:length(gnames)
            tmpdat = dayGD.(gnames{g}).(dnames{d}).(fnames{f})';
            mn = nanmean(tmpdat,1);
            stde = nanstd(tmpdat, 0, 1)./sqrt(sum(~isnan(tmpdat(:,1))));
            shadedErrorBar(pos_edges(1:end-1), mn, stde, {'Color', groupcolors.(gnames{g}).(dnames{d})},1)
        end
        ylabel(units{f})
        xlabel('distance from reward zone (deg)')
        xlim([-81 99])
    end
end

figure(fh)
axes(ax{4})
hold on
cmat = [];
for f = 13
    for d = 2
        for g = 1:length(gnames)
            tmpdat = [];
            tmpdat = dayGD.(gnames{g}).(dnames{d}).(fnames{f});
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
    ylabel('anticipatory lick slope')
    xticks([])
end

end

