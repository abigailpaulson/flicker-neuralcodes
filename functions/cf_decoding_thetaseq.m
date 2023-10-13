function cf_decoding_thetaseq(dirs, params, allindex, metadata, posType, cellTypes)
%cf_decoding_thetaseq
%
%ALP 1/5/2023

%%% rewriting files
rewrite.decodingfiles = 1; 

%%% params
params.decodingBinsize = 0.02; %in s
params.timeAroundTrough = 0.2; %in s
params.speedThreshold = 1; %deg/s
params.minTrials = 5;
if strcmp(posType, 'full')
    params.posEdges = 0:2:360;
    params.posBins = 2;
    dType = '360';
    params.nDeg = 360;
    posName = 'theta';
    params.RZ = [54 72; 234 252];
    params.dec_edges = 0:6:360;
else
    params.posEdges = -81:3:99;
    params.posBins = 3;
    dType = '180';
    params.nDeg = 180;
    posName = 'theta_d2r';
    params.RZ = [0 18];
    params.dec_edges = -81:9:99; %could also try 9
end
params.adjustedEdges = -params.nDeg/2:params.posBins:params.nDeg/2;
params.timeEdges = -params.timeAroundTrough:params.decodingBinsize:params.timeAroundTrough;
params.quadrantWindow = 0.16; %in s

%%% define areas of interest, where to get the quadrant ratio, etc
%%% set edges for histograms, areas for quantification, etc
params.PCR_positions{1} = [0-18 0+9];
params.PCR_positions{2} = [9 9+27];

%set up window for calculating the quadrant value
midI = round(length(params.timeEdges)/2); 
binRange = params.quadrantWindow/(2*params.decodingBinsize);
LBinRange = [midI-binRange midI-1];
RBinRange = [midI midI+binRange-1]; %this one should include 0
posMidI = round(length(params.adjustedEdges)/2);
binRange = (length(params.adjustedEdges)-1)/4;
FBinRange = [posMidI posMidI+binRange-1];
BBinRange = [posMidI-binRange posMidI-1];

% quadrant indices, following the convention of Farooq and Dragoi 2019
params.qX = [RBinRange; LBinRange; LBinRange; RBinRange]; %time
params.qY = [FBinRange; FBinRange; BBinRange; BBinRange]; %position

%%% set things as needed
dayindex = unique(allindex(:,1:2), 'rows');

%%% directories
savedatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
trialdatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\';
ripplechandir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\CA1\';
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';

%%% filenames
trialdatafilename = ['trialInfo_', dType, '_A'];
trialspikesfilename = ['trialSpikes_', dType, '_A'];
savefilename = ['decoding_thetaSeq_', dType, '_', cellTypes, '_alltrials'];

%% loop over all days
AllData = table; 
if exist([savedatadir, 'ALL_', savefilename, '.mat'], 'file') && ~rewrite.decodingfiles
    load([savedatadir, 'ALL_', savefilename, '.mat'])
else
    for d = 1:size(dayindex,1)
        index = allindex(allindex(:,2) == dayindex(d,2),:);
        files = index(:,3);
        anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
        dayfilename = [savedatadir, savefilename, '_', num2str(dayindex(d,2)), '.mat'];
        
        if ~exist(dayfilename, 'file') || rewrite.decodingfiles
            disp(['decoding day ', num2str(dayindex(d,2))])
            %%% load trial data
            load([trialdatadir, trialdatafilename, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
            load([trialdatadir, trialspikesfilename, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
            
            %%% load spatial rate maps - training data
            load([ratemapdir, 'spatialmaps_preflicker_0_', posType, '_', num2str(dayindex(d,2)), '.mat'])
            trainingData.ratemaps = ratemaps.fr; %in Hz
            
            %%% incorp cell types options here eventually
            allIDs = trialSpikes.info.clusterID;
            nCells = length(allIDs);
            cellI = 1:length(allIDs);
            inclCells = cellI;
            
            if strcmp(cellTypes, 'PYR')
                ct = trialSpikes.info.cellType;
                isPYR = strcmp(ct, 'PYR');
                avgFR = nanmean(ratemaps.fr,2);
                inclFR = avgFR > 0.2;
                inclCells = cellI(isPYR&inclFR');
                nCells = length(inclCells);
            elseif strcmp(cellTypes, 'PC')
                load([ratemapdir, 'spatialmaps_all_500_full_', num2str(dayindex(d,2)), '.mat'])
                inclFR = ratemaps.includeCell;
                inclCT = ratemaps.isPYR;
                inclCells = cellI(inclFR & inclCT); 
                nCells = length(inclCells);
            elseif strcmp(cellTypes, 'oldPC')
                load(['\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\placefields\CA1\data\placefields_PYR_95threshold_randomshuffle_A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
                tmpIncl = [placefields{dayindex(d,1)}{dayindex(d,2)}.include];
                clear placefields
                if exist(['\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\placefields\CA3\data\placefields_PYR_95threshold_randomshuffle_A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
                    load(['\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\placefields\CA3\data\placefields_PYR_95threshold_randomshuffle_A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
                    tmpIncl = [tmpIncl placefields{dayindex(d,1)}{dayindex(d,2)}.include];
                end
                inclCells = cellI(tmpIncl);
                nCells = length(inclCells); 
            end
            
            %%% get trials to include
            inclTrials = [trialData.engaged] & [trialData.fullTrial] & [trialData.rewarded];
            goodTrials = trialData(inclTrials);
            [goodTrials.spikeTimes] = deal(trialSpikes.data(inclTrials).spikeTimes);
            [goodTrials.spikePosInds] = deal(trialSpikes.data(inclTrials).spikePosInds);
            [goodTrials.spikeIDs] = deal(trialSpikes.data(inclTrials).IDs);
            
            nTrials = sum(inclTrials);
%             if sum(inclTrials) < params.minTrials
%                 continue
%             end

            if nTrials < 1
                continue
            end
            
            if nCells == 0
                continue
            end
            
            %%% get the theta data
            load([ripplechandir, 'bestChannel_', num2str(dayindex(d,2)), '.mat'])
            tmpRChan = bestRippleChan.all;
            theta = loaddatastruct2([anprocesseddatadir, 'CA1\', num2str(tmpRChan), '\EEG\'], dayindex(d,:), 'theta', files);
            
            %%% get testing data for each trial
            [testingData, avgTheta] = cf_gettestingdata_thetaseq(index, goodTrials, theta, params.decodingBinsize, params.timeAroundTrough, params.speedThreshold, posName, nCells, inclCells);
            
            %%% decode position during theta trough
            decodedData = cf_getdecodedposition_thetaseq(trainingData, testingData, nCells, inclCells, ...
                params);
            
            %%% do some calculations on the decoded sequences
            codingRatio = cf_getcodingratio_thetaseq(decodedData, params.PCR_positions, params); 
            
            if isempty(decodedData)
                continue
            end
            decodedData.nTrials = nTrials;
            decodedData.nCells = nCells;
            
            save([savedatadir, savefilename, '_', num2str(dayindex(d,2)), '.mat'], 'testingData', 'decodedData', 'codingRatio', 'avgTheta')
            savefigALP(savedatadir, ['FIG_', savefilename, '_', num2str(dayindex(d,2))], 'filetype', 'pdf')
        else
            load(dayfilename)
        end
        if ~isempty(decodedData)
            %new group data ALP 1/18/23 only include stuff of interest
            appendDat.allDecodedSeq = decodedData.allDecodedSeq;
            appendDat.allTroughPos = decodedData.troughPos;
            appendDat.allTrialID = decodedData.trialID; 
            appendDat.mnDecodedSeq = decodedData.mnDecodedSeq;
            appendDat.mnQR = decodedData.mnQuadrantRatio;
            appendDat.shuffQR = decodedData.shuff_QR_thresh; 
            appendDat.sigDay = decodedData.significantSeq; 
            appendDat.PCR_overposition_trial = codingRatio.PCR_overposition_trial;
            appendDat.PCR_overposition_day = codingRatio.PCR_overposition_day;
            appendDat.PCR_prerew_trial = codingRatio.PCR_prerew_trial;
            appendDat.PCR_prerew_day = codingRatio.PCR_prerew_day;
            appendDat.trough_counts_day = codingRatio.trough_counts_day;
            appednDat.trough_counts_trial = codingRatio.trough_counts_trial; 
            appendDat.mnThetaTrace = avgTheta;
            appendDat.day = d;
            
            if metadata.FlickerDay(d) < 5
                tmpTimepoint = {'pre'};
            else
                tmpTimepoint = {'post'};
            end
            
            allData(d) = appendDat; 
        else
            continue
        end        
        clear appendDat decodedData testingData trainingData avgTheta tmp*
    end
    
    %%% get group data and save
    save([savedatadir, 'ALL_', savefilename], 'allData', '-v7.3')
end

%% concatenate group data by condition and group
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

%%% only include days that meet the significance criteria
for d = 1:size(dayindex,1)
    if isempty(allData(d).sigDay)
        allData(d).sigDay = 0; 
        allData(d).day = d;
    end
end
sigDay = [allData.sigDay]; 
sigDay = logical(sigDay)';

fnames = fieldnames(allData);

%%%%% ----- split data by groups ----- %%%%%
GD = [];
for g = 1:numel(gnames)
    for d = 1:numel(dnames)
        isDay = gID.(gnames{g}).(dnames{d});
        allDays = [allData.(fnames{end})]; 
        includeDays = ismember(allDays, isDay)'; 
        for f = 4:length(fnames)
            tmp = {allData(includeDays&sigDay).(fnames{f})};
           if strcmp(fnames{f}, 'PCR_prerew_trial')
                 GD.(gnames{g}).(dnames{d}).(fnames{f}) = cell2mat(tmp);
           elseif strcmp(fnames{f}, 'mnDecodedSeq')
               GD.(gnames{g}).(dnames{d}).(fnames{f}) = cat(3, tmp{:});
           else
                GD.(gnames{g}).(dnames{d}).(fnames{f}) = cell2mat(tmp');
           end
            clear tmp
        end
        clear isDay allDays includeDays
    end
end

%%% save grouped data for plotting
save([savedatadir, 'GROUPEDDATA_', savefilename, '.mat'], 'GD', 'params', 'groupcolors', 'gnames', 'dnames', 'fnames')

xvect_bar = [1 2];
figure
hold on
f = 11;
for d = 2
    for g = 1:length(gnames)
        scatterdat{g} = GD.(gnames{g}).(dnames{d}).(fnames{f});
        plotdat_mn(g) = mean(scatterdat{g}, 'omitnan');
        plotdat_stde(g) = std(scatterdat{g}, 'omitnan')./sqrt(sum(~isnan(scatterdat{g})));
        colororder(g,:) = groupcolors.(gnames{g}).(dnames{d});
    end
end
plotprettypoints(f, xvect_bar, scatterdat, 'color', colororder)
b = bar(xvect_bar, plotdat_mn, 'FaceColor', 'flat');
b.CData = colororder;
b.FaceAlpha = 0.6;
errorbar2(xvect_bar, plotdat_mn, plotdat_stde, 0.2, 'k-', 'LineWidth', 0.75);
[~, p] = ttest2(scatterdat{1},scatterdat{2});
ylabel('prospective coding ratio')
xticks([])
title(['PCR per sess p = ', num2str(p)])

figure
hold on
f = 9;
plotdat_mn = [];
plotdat_sem = [];
for d = 2
    for g = 1:length(gnames)
        plotdat_mn(g,:) = mean(GD.(gnames{g}).(dnames{d}).(fnames{f}),1,'omitnan');
        plotdat_sem(g,:) = std(GD.(gnames{g}).(dnames{d}).(fnames{f}), 0,1,'omitnan')./sqrt(sum(~isnan(GD.(gnames{g}).(dnames{d}).(fnames{f})(:,1))));
        shadedErrorBar(params.dec_edges(1:end-1)+3, plotdat_mn(g,:), plotdat_sem(g,:), {'Color', groupcolors.(gnames{g}).(dnames{d}), 'LineWidth', 2}, 1);
    end
end

figure
hold on
f = 4;
ploti = 1;
for d = 1:2
    for g = 1:length(gnames)
        subplot(2,2,ploti)
        hold on
        imagesc(mean(GD.(gnames{g}).(dnames{d}).(fnames{f}),3,'omitnan'), [0.015 0.035])
        colorbar
        ploti = ploti+1;
    end
end


% f = figure('Position', [-847 -143 774 995]);
% hold on
% clear plotdat_mn plotdat_std GData RData
% subplot(2,3,3)
% hold on
% GData = [PCRData(isGamma&isPost&inclDay).pcr_pos_avg];
% RData = [PCRData(isRandom&isPost&inclDay).pcr_pos_avg];
% scatterdat{2} = GData(1,:);
% scatterdat{1} = RData(1,:);
% plotdat_mn(2) = mean(scatterdat{2},'omitnan');
% plotdat_mn(1) = mean(scatterdat{1}, 'omitnan');
% plotdat_std(2) = std(scatterdat{2},'omitnan')./sqrt(sum(~isnan(scatterdat{2})));
% plotdat_std(1) = std(scatterdat{1},'omitnan')./sqrt(sum(~isnan(scatterdat{1})));
% colororder(1,:) = groupcolors.random.post;
% colororder(2,:) = groupcolors.gamma.post;
% 
% plotprettypoints(f, xvect_bar, scatterdat, 'color', colororder)
% b = bar(xvect_bar, plotdat_mn, 'FaceColor', 'flat');
% b.CData = colororder;
% b.FaceAlpha = 0.6;
% errorbar2(xvect_bar, plotdat_mn, plotdat_std, 0.2, 'k-', 'LineWidth', 0.75);
% [~, p] = ttest2(scatterdat{1},scatterdat{2});
% ylabel('prospective coding ratio')
% xticks([])
% title(['PCR per sess p = ', num2str(p)])
% 
% 
% clear plotdat_mn plotdat_std GData RData
% % f = figure;
% subplot(2,3,[1 2])
% hold on
% GData = [PCRData(isGamma&isPost&inclDay).pcr_pos]';
% RData = [PCRData(isRandom&isPost&inclDay).pcr_pos]';
% plotdat_mn{2} = mean(GData,1,'omitnan');
% plotdat_mn{1} = mean(RData,1, 'omitnan');
% plotdat_std{2} = std(GData,0,1,'omitnan')./sqrt(sum(~isnan(GData(:,1))));
% plotdat_std{1} = std(RData,0,1,'omitnan')./sqrt(sum(~isnan(RData(:,1))));
% shadedErrorBar(dec_edges(1:end-1)+3, plotdat_mn{1}, plotdat_std{1}, {'Color', groupcolors.random.post, 'LineWidth', 2}, 1);
% shadedErrorBar(dec_edges(1:end-1)+3, plotdat_mn{2}, plotdat_std{2}, {'Color', groupcolors.gamma.post, 'LineWidth', 2}, 1);
% xlabel('position (deg)')
% xlim([-81 99])
% xticks([-81 0 99])
% ylabel('prospective coding ratio')
% title('prospective coding ratio over position')
% figname = ['FIG_prospectivecoding_overPosition_', num2str(nDeg), '_', cellTypes];
% % savefigALP(savedatadir, figname, 'fileType', 'pdf')
% 
% % f = figure;
% subplot(2,3,4)
% hold on
% vdat.randompost = [PCRData(isRandom&isPost&inclDay).pcr_loc_trial];
% vdat.gammapost = [PCRData(isGamma&isPost&inclDay).pcr_loc_trial];
% cmat = [groupcolors.random.post; groupcolors.gamma.post];
% violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false);
% ylabel('RZ prospective coding ratio')
% p = ranksum(vdat.gammapost, vdat.randompost);
% yticks([-1 -0.5 0 0.5 1])
% title(['ranksum p = ', num2str(p)])
% 
% % f = figure;
% % subplot(2,3,5)
% hold on
% cdat.randompost = {PCRData(isRandom&isPost&inclDay).pcr_pos_trial};
% cdat.randompost = cell2mat(cdat.randompost');
% cdat.gammapost = {PCRData(isGamma&isPost&inclDay).pcr_pos_trial};
% cdat.gammapost = cell2mat(cdat.gammapost');
% 
% [~,iSort] = sort(vdat.randompost);
% cdat.randompost_sort = cdat.randompost(iSort,:);
% [~, iSort] = sort(vdat.gammapost); 
% cdat.gammapost_sort = cdat.gammapost(iSort,:); 
% 
% subplot(2,3,5)
% hold on
% imagesc(dec_edges, 1:size(cdat.randompost_sort,1), cdat.randompost_sort, [-2 2])
% xlim([dec_edges(1) dec_edges(end)])
% ylim([1 size(cdat.randompost_sort,1)])
% xlabel('track position (deg)')
% ylabel('trial sorted by RZ PCR')
% colormap(bluewhitered)
% colorbar
% % cmocean('balance', 'pivot', 0)
% title('random')
% subplot(2,3,6)
% hold on
% imagesc(dec_edges, 1:size(cdat.gammapost_sort,1), cdat.gammapost_sort, [-2 2])
% xlim([dec_edges(1) dec_edges(end)])
% ylim([1 size(cdat.gammapost_sort,1)])
% xlabel('track position (deg)')
% ylabel('trial sorted by RZ PCR')
% colormap(bluewhitered)
% colorbar
% % cmocean('balance', 'pivot', 0)
% title('gamma')
% makefigurepretty(gcf)
% savefigALP(savedatadir, ['GROUP_RZ_prospectivecodingratio_', posType, '_', cellTypes])
% 




%% plot 

%%% plot group average over space

%%% get prospective quadrant ratio for the front of the reward zone for each day, plot
%%% the 

end

