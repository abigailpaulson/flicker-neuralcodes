   
function cf_spatialmaps_fieldproperties(dirs, params, allindex, metadata, posType)
%cf_spatialmaps_fieldproperties
%ALP 12/19/2022

%%% directories, filenames, etc
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';
mapfilename = ['spatialmaps_preflicker_0_', posType, '_'];
PCfilename = ['spatialmaps_all_500_full_'];
dayindex = unique(allindex(:,1:2), 'rows');

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
groupcolors = params.colors;

%%% helpful stuff about the track
if strcmp(posType, 'full')
    params.pos_edges = 0:2:360;
% RZ = [54 72; 234 252];
RZ = [50 75; 230 255];
%     params.AZ = [234 252];
%     RZbins = [find(params.pos_edges == params.RZ(1)) find(params.pos_edges == params.RZ(2))-1];
%     AZbins = [find(params.pos_edges == params.AZ(1)) find(params.pos_edges == params.AZ(2))-1]
elseif strcmp(posType, 'd2r')
    params.pos_edges = -81:3:99;
    RZ = [0 18];
%     params.AZ = [-18+9 0];
%     RZbins = [find(params.pos_edges == params.RZ(1)) find(params.pos_edges == params.RZ(2))-1];
%     AZbins = [find(params.pos_edges == params.AZ(1)) find(params.pos_edges == params.AZ(2))-1];
end

%%% get per day metrics from the ratemap files

for d = 1:size(dayindex,1)
    dayfilename = [mapfilename, num2str(dayindex(d,2)), '.mat'];
    PCdayfilename = [PCfilename, num2str(dayindex(d,2)), '.mat'];
    load([ratemapdir, dayfilename])
    PCmaps = load([ratemapdir, PCdayfilename]);
    
    tmpdat = ratemaps.fr; 
    findnans = isnan(tmpdat(:,1));
    omit = findnans; 
    allData(d).ratemaps = ratemaps.fr(~omit,:); 
    allData(d).SI = ratemaps.spatialinfo(~omit)';
    allData(d).sparsity = ratemaps.sparsity(~omit);
    allData(d).peakpos = ratemaps.peakpos(~omit);
    allData(d).maxFR = ratemaps.maxFR(~omit)';
    allData(d).meanFR = ratemaps.meanFR(~omit)';
    allData(d).isPC = PCmaps.ratemaps.includeCell(~omit)';
    tmpreg = ones(1,sum(~omit));
    isCA3 = strcmp(ratemaps.brainReg, 'CA3');
    tmpreg(isCA3) = 2;
    allData(d).region = tmpreg;
    allData(d).animal = dayindex(d,1).*ones(1,sum(~omit));
    allData(d).day = d.*ones(1,sum(~omit)); 

end

%%% concatenate all fields by group
fnames = fieldnames(allData);
alldays = 1:length(allData);
for g = 1:numel(gnames)
    for d = 1:numel(dnames)
        isDay = ismember(alldays, gID.(gnames{g}).(dnames{d}));
        for f = 1:numel(fnames)
            if f == 1 %rate maps
                tmp = {allData(isDay).(fnames{f})};
                groupData.(gnames{g}).(dnames{d}).(fnames{f}) = cell2mat(tmp');
            else
                groupData.(gnames{g}).(dnames{d}).(fnames{f}) = [allData(isDay).(fnames{f})];
            end
        end
    end
end

save([ratemapdir, 'ALLDAYS_', mapfilename(1:end), '.mat'], 'allData', 'groupData', 'params', 'RZ')


%%% get averages of some metrics per day
large_edges = 0:18:360;
for g = 1:numel(gnames)
    for d = 1:numel(dnames)
        tmpdays = groupData.(gnames{g}).(dnames{d}).(fnames{end-1});
        tmpID = unique(tmpdays); 
        for f = 2:length(fnames)-3
            day_avg = [];
            tmpdat = groupData.(gnames{g}).(dnames{d}).(fnames{f});
            if f ~= 4
                day_avg = arrayfun(@(x) nanmean(tmpdat(tmpdays == x)), tmpID);
            else
                day_avg = arrayfun(@(x) histcounts(tmpdat(tmpdays == x), large_edges), tmpID, 'UniformOutput', false);
                day_avg = cell2mat(day_avg');
                day_avg = day_avg./(sum(day_avg,2));
%             else
%                 day_avg = arrayfun(@(x) nanmean(tmpdat(:,tmpdays == x),2), tmpID, 'UniformOutput', false);
%                 day_avg = cell2mat(day_avg);
            end
            dayData.(gnames{g}).(dnames{d}).(fnames{f}) = day_avg;
        end
    end
end

%%% plotting
figure(1)
hold on
f = 1;
iPlot = 1;
for d = 1:2
    for g = 1:2
        figure(1)
        hold on
        subplot(2,2,iPlot)
        hold on
        tmp = groupData.(gnames{g}).(dnames{d}).peakpos;
        fr = groupData.(gnames{g}).(dnames{d}).(fnames{f}); 
        [~, iSort] = sort(tmp);
        incl = groupData.(gnames{g}).(dnames{d}).isPC;
        dat = groupData.(gnames{g}).(dnames{d}).(fnames{f})./mean(groupData.(gnames{g}).(dnames{d}).(fnames{f}),2);
        iSort = iSort(incl);
        plotDat = dat(iSort,:);
%         dat = dat(incl,:);
%         iSort = iSort(incl);
        imagesc([0:2:358]+1, 1:size(plotDat), plotDat, [0.5 1.5])
        xlim([0 360])
        xticks([0 360])
        yticks([0 round(size(plotDat,1))])
        ylim([0 size(plotDat,1)])
        title([gnames{g} ' ' dnames{d}])
        colorbar
        xlabel('position (deg)')
        ylabel('place cell count')    
        iPlot = iPlot+1; 
    end
end
figname = ['GROUP_placecellmaps_', posType];
makefigurepretty(gcf)
savefigALP(ratemapdir, figname)

%%% all cells
figure
hold on
tiledlayout(1,length(fnames), 'TileSpacing', 'Compact', 'Padding', 'Compact'); 
hold on
for f = 2:length(fnames)-4
    nexttile
    hold on
    vdat = []; cmat = [];
    for d = 1:2
        for g = 1:2
            nm = [gnames{g}, dnames{d}];
            vdat.(nm) = groupData.(gnames{g}).(dnames{d}).(fnames{f}); 
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
    end
    vdat.xlabels = {'pre', 'post'};
    violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', true);
%     ylabel(units{f})
    title(fnames{f})
end

%%% low fr cells
figure
hold on
tiledlayout(1,length(fnames), 'TileSpacing', 'Compact', 'Padding', 'Compact'); 
hold on
for f = 2:length(fnames)-4
    nexttile
    hold on
    vdat = []; cmat = [];
    for d = 1:2
        datforstats = [];
        for g = 1:2
            nm = [gnames{g}, dnames{d}];
            tmpdat = groupData.(gnames{g}).(dnames{d}).(fnames{f}); 
            incldat = groupData.(gnames{g}).(dnames{d}).meanFR < 5; 
            vdat.(nm) = tmpdat(incldat);
            datforstats{g} = tmpdat(incldat);
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
        stat(d) = ranksum(datforstats{1}, datforstats{2});
    end
    vdat.xlabels = {'pre', 'post'};
    violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false);
%     ylabel(units{f})
    title([fnames{f} num2str(stat)])
end

%%% high fr cells
figure
hold on
tiledlayout(1,length(fnames), 'TileSpacing', 'Compact', 'Padding', 'Compact'); 
hold on
for f = 2:length(fnames)-4
    nexttile
    hold on
    vdat = []; cmat = [];
    for d = 1:2
        for g = 1:2
            nm = [gnames{g}, dnames{d}];
            tmpdat = groupData.(gnames{g}).(dnames{d}).(fnames{f}); 
            incldat = groupData.(gnames{g}).(dnames{d}).meanFR > 5; 
            vdat.(nm) = tmpdat(incldat);
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
    end
    vdat.xlabels = {'pre', 'post'};
    violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false);
%     ylabel(units{f})
    title(fnames{f})
end

%%% RZ cells
figure
hold on
tiledlayout(1,length(fnames), 'TileSpacing', 'Compact', 'Padding', 'Compact'); 
hold on
for f = 2:length(fnames)-4
    nexttile
    hold on
    vdat = []; cmat = [];
    for d = 1:2
        for g = 1:2
            nm = [gnames{g}, dnames{d}];
            tmpdat = groupData.(gnames{g}).(dnames{d}).(fnames{f});
            peakdat = groupData.(gnames{g}).(dnames{d}).peakpos;
            regdat = groupData.(gnames{g}).(dnames{d}).region;
            if size(RZ,1) > 1
                incldat =  (peakdat >= RZ(1,1) & peakdat < RZ(1,2)) | (peakdat >= RZ(2,1) & peakdat < RZ(2,2));
            else
                incldat =  (peakdat >= RZ(1,1) & peakdat < RZ(1,2));
            end
            vdat.(nm) = tmpdat(incldat);
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
    end
    vdat.xlabels = {'pre', 'post'};
    violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', true);
%     ylabel(units{f})
    title(fnames{f})
end

%%% place cells 
yy = {'spatial information (bits/spike)', 'sparsity', 'peak position (deg)', 'peak firing rate (hz)', 'mean firing rate (hz)'};
xxl = {[0 3], [0 1], [-81 99], [0 10], [0 10]};
figure('Position', [1516 241 1685 371])
hold on
for f = 2:length(fnames)-4
    subplot(1,5,f-1)
    hold on
    vdat = []; cmat = []; stat = [];
    for d = 2
        datforstats = [];
        for g = 1:2
            nm = [gnames{g}, dnames{d}];
            tmpdat = groupData.(gnames{g}).(dnames{d}).(fnames{f}); 
            incldat =  groupData.(gnames{g}).(dnames{d}).isPC;
            vdat.(nm) = tmpdat(incldat);
            datforstats{g} = tmpdat(incldat);
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
        stat(d) = ranksum(datforstats{1}, datforstats{2});
    end
    vdat.xlabels = {'pre', 'post'};
    violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false);
%     ylabel(units{f})
    ylabel(yy{f-1})
    ylim(xxl{f-1})
    title([fnames{f}, num2str(stat)])
end
figname = ['GROUP_placecellproperties_', posType];
makefigurepretty(gcf)
savefigALP(ratemapdir, figname)


end

