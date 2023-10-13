function cf_decoding_thetaseq_table_ctrlspeed(dirs, params, allindex, metadata, posType, cellTypes)
%cf_decoding_thetaseq
%
%ALP 1/5/2023

%%% rewriting files
rewrite.decodingfiles = 0; 

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
    params.AZ = [-36 0];
    params.dec_edges = -81:9:99; %could also try 9
end
params.adjustedEdges = -params.nDeg/2:params.posBins:params.nDeg/2;
params.timeEdges = -params.timeAroundTrough:params.decodingBinsize:params.timeAroundTrough;
params.quadrantWindow = 0.16; %in s

%%% define areas of interest, where to get the quadrant ratio, etc
%%% set edges for histograms, areas for quantification, etc
params.PCR_positions{1} = [0-18 0+9];
params.PCR_positions{2} = [9 9+27];

%%% definie a different area of interes
params.PCR_secondhalf{1} = [18 72];

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
savefilename = ['thetaseqdecoding_alldays_alltrials_230718.mat'];

if ~exist([savedatadir, savefilename]) || rewrite.decodingfiles
%% loop over all days
AllData = table; 
    for d = 1:size(dayindex,1)
        index = allindex(allindex(:,2) == dayindex(d,2),:);
        files = index(:,3);
        anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
        dayfilename = [savedatadir, savefilename, '_', num2str(dayindex(d,2)), '.mat'];
        day = d;
        
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
            %             inclTrials = [trialData.engaged] & [trialData.fullTrial] & [trialData.rewarded];
            %             inclTrials = [trialData.fullTrial]; %want to try getting all data in order to look at unengaged trials
            inclTrials = true(1,length(trialData));
            goodTrials = trialData(inclTrials);
            
            [goodTrials.spikeTimes] = deal(trialSpikes.data(inclTrials).spikeTimes);
            [goodTrials.spikePosInds] = deal(trialSpikes.data(inclTrials).spikePosInds);
            [goodTrials.spikeIDs] = deal(trialSpikes.data(inclTrials).IDs);
            inclTrials = ones(1,length(trialData));
            
            nTrials = sum(inclTrials);
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
            thetatrace = table(day, avgTheta);
            
            %%% decode position during theta trough
            decodedData = cf_getdecodedposition_thetaseq_table(trainingData, testingData, nCells, inclCells, ...
                params);
            
            %%% do some calculations on the decoded sequences
            codingRatio = cf_getcodingratio_thetaseq_table(decodedData, params.PCR_positions, params); 
            
            %%% do some calculations on the decoded sequences
            codingRatio2 = cf_getcodingratio_thetaseq_table(decodedData, params.PCR_secondhalf, params); 
            newVarNames = append("PostR_", fieldnames(codingRatio2)); 
            newVarNames = newVarNames(1:5); %changed from 1:4 ALP 7/18
            codingRatio2.Properties.VariableNames = newVarNames; 
            
            %%% load behavioral metrics
            behaviordatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\'];
            behaviorfilename = ['behavioranalysis_prepost_recordings_table_', num2str(dayindex(d,2)), '.mat'];
            load([behaviordatadir, behaviorfilename], 'metrics')

            %%% figure out if this day has significant sequences or not,
            %%% compare full quadrant ratio on mean day sequence to shuffled data
            sig = cf_getsignificant_thetaseq_table(decodedData, params, metrics); 

            %%% compile data of interest
            tmpData = [metrics(:, end-4:end) metrics(:,25:27) metrics(:, 19) metrics(:,15:18) metrics(:,3:14) decodedData codingRatio codingRatio2];
            tmpData = outerjoin(tmpData, sig, 'MergeKeys', true);
            tmpData = outerjoin(tmpData, thetatrace, 'MergeKeys', true);
            
            if isempty(decodedData)
                continue
            end
            
%             save([savedatadir, savefilename, '_', num2str(dayindex(d,2)), '.mat'], 'testingData', 'decodedData', 'codingRatio', 'avgTheta')
%             savefigALP(savedatadir, ['FIG_', savefilename, '_', num2str(dayindex(d,2))], 'filetype', 'pdf')
        else
%             load(dayfilename)
        end
        
            AllData = [AllData; tmpData];
            
        clear appendDat decodedData testingData trainingData avgTheta tmp* sig codingRatio
    end
    
    save([savedatadir, savefilename], 'AllData', '-v7.3')
    
%% write table for statistics in R
statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_ThetaSequences_perTrial_230615.txt';
inclTrials = AllData.engaged == 1 & AllData.fullTrial == 1 & AllData.rewarded == 1 & AllData.significantSeq == 1;
tmpData = AllData(inclTrials, [1:9,31, 35, 39]); %alp 6/20/23 I think I may have messed up what variables save here.... check with the old version of this structure to be sure. 
% tmpData.date = dayindex(tmpData.day,2);
tmpData.date = tmpData.day; 
writetable(tmpData, fullfile(statsdir, filename))

% dayData = groupsummary(tmpData(inclTrials,:), "day", "mean");
dayData = groupsummary(tmpData(:,[1:2,5:12]), "day", "mean");

% dayData2 = groupsummary(tmpData(inclTrials,[1:2,5:10]), "day", "unique");
A = unique(tmpData.animal);
% newData = join(dayData, tmpData);
[~, iA, ~] = unique(tmpData.day);

dayInfo.timepoint = tmpData.timepoint(iA);
dayInfo.animal = tmpData.animal(iA);
dayInfo.group = tmpData.group(iA);
dayInfo.day = tmpData.day(iA);
dayInfo.date = dayindex(dayInfo.day,2);
dayInfo = struct2table(dayInfo);
newData = join(dayData, dayInfo);

dayData = [];
dayData = newData;
save([savedatadir, 'thetaseqdecoding_includeddays_group_230615.mat'], 'dayData')
filename = 'TableData_ThetaSequences_perDay_230615.txt';
writetable(newData, fullfile(statsdir, filename))
    
else
    load([savedatadir, savefilename])
end

%% set color options
groupcolors = cbrewer('qual', 'Paired', 6); 
gcolors = groupcolors([2,1,4,3],:); %rearranging it cuz of the way that the gramm does the plotting
scolors = groupcolors;
trialscolors = groupcolors(end-1:end,:);

%% try plotting over behavior quartiles
otherdecedges = -81:3:99;
fnames = {'lickDI', 'nLicks', 'licklatency_s', 'duration', 'trial_speed', ...
    'AZ_speed', 'RZ_speed', 'Ctrl_speed', 'AZ_speed_slope', 'AZ_lick_slope'};
figdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 

for f = 1:length(fnames)
    clear g
    tmpData = AllData.(fnames{f});
    qs = quantile(tmpData, 3);
    qs = [min(tmpData) qs max(tmpData)];
    iQ = discretize(tmpData, qs);
    
    g(1,1) = gramm('x', repmat(params.dec_edges(1:end-1), [height(AllData),1]), 'y', AllData.PCR_pos_trial, 'lightness', iQ, 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1);
    g(1,1).facet_grid([], cellstr(AllData.timepoint));
    g(1,1).stat_summary();
    g(1,1).set_names('x','distance to reward (deg)', 'y', 'prospective coding ratio', 'lightness', 'quartile');
    g(1,1).set_title(fnames{f});
    g(1,1).axe_property('YLim', [-0.25 0.25]);
    
    g(2,1) = gramm('x', iQ, 'y', AllData.PCR_loc_trial, 'color', cellstr(AllData.group), 'lightness', iQ, 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & ~isnan(AllData.PCR_loc_trial));
    g(2,1).facet_grid([], cellstr(AllData.timepoint));
%         g(2,1).stat_bin('geom', 'stairs', 'fill', 'transparent');
%     g(2,1).stat_bin('normalization','cdf','geom','stairs');
%     g(2,1).stat_violin();
    g(2,1).stat_boxplot();
    g(2,1).axe_property('YLim', [-0.25 0.25]);

    figure
    g.draw();
%     savefigALP(figdir, ['PCR_overposition_bybehaviorquarts_timepoint_', fnames{f}], 'filetype', 'pdf')
end

%% try plotting control speed separated out by group, and try only plotting speed halves
otherdecedges = -81:3:99;
fnames = {'licklatency_s', 'duration', 'trial_speed', ...
    'AZ_speed', 'RZ_speed', 'Ctrl_speed', 'AZ_speed_slope', 'AZ_lick_slope'};

for f = 1:length(fnames)
    clear g
    inclData = AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & AllData.significantSeq == 1;
    tmpData = AllData.(fnames{f})(inclData);
%     qs = quantile(tmpData, 3);
    qs = median(tmpData, 'omitnan');
    qs = [min(tmpData) qs max(tmpData)];
    iQ = discretize(tmpData, qs);
    PlotData = AllData(inclData,:);
    
    g(1,1) = gramm('x', repmat(params.dec_edges(1:end-1), [height(PlotData),1]), 'y', PlotData.PCR_pos_trial, 'color', cellstr(PlotData.group), 'lightness', iQ, 'subset', PlotData.fullTrial==1 & PlotData.engaged ==1 & PlotData.rewarded == 1 & strcmp(PlotData.timepoint, 'post'), 'column', PlotData.animal);
    g(1,1).facet_grid([], cellstr(PlotData.group));
    g(1,1).stat_summary();
    g(1,1).set_names('x','distance to reward (deg)', 'y', 'prospective coding ratio', 'lightness', 'Data Half');
    g(1,1).set_title(fnames{f});
    g(1,1).axe_property('YLim', [-0.25 0.25]);
    
%     g(2,1) = gramm('x', iQ, 'y', AllData.PCR_loc_trial, 'color', cellstr(AllData.group), 'lightness', iQ, 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & ~isnan(AllData.PCR_loc_trial));
%     g(2,1).facet_grid([], cellstr(AllData.timepoint));
% %         g(2,1).stat_bin('geom', 'stairs', 'fill', 'transparent');
% %     g(2,1).stat_bin('normalization','cdf','geom','stairs');
% %     g(2,1).stat_violin();
%     g(2,1).stat_boxplot();
%     g(2,1).axe_property('YLim', [-0.25 0.25]);

    figure
    g.draw();
%     savefigALP(figdir, ['PCR_overposition_bybehaviorquarts_timepoint_', fnames{f}], 'filetype', 'pdf')
end

%% get animals in the top and bottom speed thresholds
otherdecedges = -81:3:99;
fnames = {'Ctrl_speed'};

for f = 1:length(fnames)
    clear g
    inclData = AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & AllData.significantSeq == 1;
    tmpData = AllData.(fnames{f})(inclData);
%     qs = quantile(tmpData, 3);
    qs = median(tmpData, 'omitnan');
    qs = [min(tmpData) qs max(tmpData)];
    iQ = discretize(tmpData, qs);
    PlotData = AllData(inclData,:);
    
    g(1,1) = gramm('x', PlotData.(fnames{f}), 'color', cellstr(PlotData.group), 'lightness', iQ, 'subset', PlotData.fullTrial==1 & PlotData.engaged ==1 & PlotData.rewarded == 1 & strcmp(PlotData.timepoint, 'post'),'column', PlotData.animal);
    g(1,1).facet_wrap(PlotData.animal);
    g(1,1).stat_bin('edges', qs, 'geom', 'bar')
%     g(1,1).stat_summary();
%     g(1,1).set_names('x','distance to reward (deg)', 'y', 'prospective coding ratio', 'lightness', 'Data Half');
    g(1,1).set_title(fnames{f});
%     g(1,1).axe_property('YLim', [-0.25 0.25]);
    
%     g(2,1) = gramm('x', iQ, 'y', AllData.PCR_loc_trial, 'color', cellstr(AllData.group), 'lightness', iQ, 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & ~isnan(AllData.PCR_loc_trial));
%     g(2,1).facet_grid([], cellstr(AllData.timepoint));
% %         g(2,1).stat_bin('geom', 'stairs', 'fill', 'transparent');
% %     g(2,1).stat_bin('normalization','cdf','geom','stairs');
% %     g(2,1).stat_violin();
%     g(2,1).stat_boxplot();
%     g(2,1).axe_property('YLim', [-0.25 0.25]);

    figure
    g.draw();
%     savefigALP(figdir, ['PCR_overposition_bybehaviorquarts_timepoint_', fnames{f}], 'filetype', 'pdf')
end

%% Manuscript Figure 5
%get a boxplot of prospective coding after the reward zone, per trial, by top and bottom half speed threshold
%6/20/23
figdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 

% want to get everything that I'm interested in into one structure so that
% it is easier to work with
behaviordir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\']; 
behaviorfilename = 'behavioranalysis_prepost_recordings_table_v2.mat';
load(behaviorfilename)

%%% what behavior metric to use to split the data
splitname = 'Ctrl_speed';

%%% get trials to include
isNaNTrial = isnan(AllData.trialNum);
AllData = AllData(~isNaNTrial,:);
inclData = AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & AllData.significantSeq == 1;
tmpData = AllData(inclData,:); 
tmpBhvData = allTMetrics(inclData,:); 

%%% split trials by speed
qs = median(tmpData.(splitname), 'omitnan');
qs = [min(tmpData.(splitname)) qs max(tmpData.(splitname))];
iQ = discretize(tmpData.(splitname), qs);
plotData = tmpData.PostR_PCR_loc_trial; 
velData = tmpBhvData.vel_h; 

clear g
g(1,1) = gramm('x', iQ, 'y', plotData, 'lightness', iQ, 'color', cellstr(tmpData.group)); 
g(1,1).facet_grid([], cellstr(tmpData.group));
g(1,1).stat_violin('width', 2.5);
%g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g(1,1).set_names('x','control speed subset', 'y', 'post reward zone prospective coding ratio', 'lightness', 'Data Half');
g(1,1).set_title({'prospective coding after RZ by speed subsets', 'per trial - post only'});
figure
g.draw(); 

%%% make the stats table for this
postRewardPCR = [];
postRewardPCR.([splitname, '_subset']) = iQ; 
postRewardPCR.PostR_PCR_loc_trial = plotData; 
postRewardPCR.velocity = velData;
postRewardPCR.Ctrl_speed = tmpData.(splitname);
postRewardPCR.PCR_pos_trial = tmpData.PCR_pos_trial;
postRewardPCR.PCR_loc_trial = tmpData.PCR_loc_trial; 
postRewardPCR.animal = tmpData.animal; 
postRewardPCR.group = tmpData.group;

postRPCR = struct2table(postRewardPCR); 
statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_ProspectiveCoding_AfterReward_perTrial_postDataOnly_CtrlSpeed.txt';
writetable(postRPCR, fullfile(statsdir, filename))

%%% ALP to do 6/21/23
behavioredges = -81:3:99;
postR_group = grpstats(postRPCR, {'group', 'Ctrl_speed_subset'}, {'mean', 'sem'});


%%% average of high and low speeds over trials
figure('Position', [242 487 362 185])
hold on
iPlot = 1;
for g = 1:2
    subplot(1,2,g)
    hold on
    for s = 1:2
        mn = []; sem = [];
        mn = postR_group.mean_velocity(iPlot,:);
        sem = postR_group.sem_velocity(iPlot,:);
        shadedErrorBar(behavioredges(1:end-1), mn, sem, {'Color', scolors(iPlot,:)},1);
        iPlot = iPlot+1;
    end
    xlabel('distance to reward zone (deg)')
    ylabel('speed (deg/s)')
    ylim([2 12])
    yticks([2 7 12])
    xticks([-81 0 99])
    xlim([-81 99])
end
makefigurepretty(gcf,1)
filename = 'draftfigure5_velocity_velsplit';
savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

% prospective coding over position split by velocity - per trial
figure('Position', [242 487 362 185])
hold on
iPlot = 1;
for g = 1:2
    subplot(1,2,g)
    hold on
    for s = 1:2
        mn = []; sem = [];
        mn = postR_group.mean_PCR_pos_trial(iPlot,:);
        sem = postR_group.sem_PCR_pos_trial(iPlot,:);
        shadedErrorBar(params.dec_edges(1:end-1), mn, sem, {'Color', scolors(iPlot,:)},1);
        iPlot = iPlot+1;
    end
    xlabel('distance to reward zone (deg)')
    ylabel('prospective coding ratio')
    ylim([-0.25 0.25])
    yticks([-0.25 0 0.25])
    xticks([-81 0 99])
    xlim([-81 99])
end
makefigurepretty(gcf,1)
filename = 'draftfigure5_PCR_overpos_velsplit';
savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

% make this a violin plot for the main figure
gnames = {'gamma', 'random'};
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = [];
for g = 1:2
    isGroup = strcmp(postRPCR.group, gnames{g});
    for s = 1:2
        isPlotTrial = isGroup & (postRPCR.Ctrl_speed_subset == s);
        vdat.([gnames{g}, '_', num2str(s)]) = postRPCR.PostR_PCR_loc_trial(isPlotTrial);
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40);
% violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', true, 'DataPointSize', 2);
xticklabels({'low', 'high', 'low', 'high'})
ylim([-0.5 0.75])
yticks([-0.5 0 0.5])
ylabel('prospective coding ratio')
title('post reward PCR')
makefigurepretty(gcf,1)
filename = 'draftfigure5_postR_PCR_velsplit';
savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

% violin plot of control speeds across groups
gnames = {'gamma', 'random'};
dnames = {'pre', 'post'}; %here post is used for speed set 1 colors
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for g = 1:2
    isGroup = strcmp(postRPCR.group, gnames{g});
    for s = 1:2
        isPlotTrial = isGroup & (postRPCR.Ctrl_speed_subset == s);
        vdat.([gnames{g}, '_', num2str(s)]) = postRPCR.Ctrl_speed(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{s})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40);
% violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', true, 'DataPointSize', 2);
xticklabels({'low', 'high', 'low', 'high'})
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
ylabel('velocity (deg/s)')
title('control speed')
makefigurepretty(gcf,1)
filename = 'draftfigure5_Ctrl_speed_velsplit';
savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

% supplement...
%distribution of control speeds across trials between groups
%other analysis (time below threshold?) per trial 


%supplement 6/28/23 
%   prospective coding in the reward zone comparison between groups split
%   by speed
gnames = {'gamma', 'random'};
dnames = {'pre', 'post'}; %here post is used for speed set 1 colors
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for s = 1:2
    for g = 1:2
        isGroup = strcmp(postRPCR.group, gnames{g});
        isPlotTrial = isGroup & (postRPCR.Ctrl_speed_subset == s);
        vdat.([gnames{g}, '_', num2str(s)]) = postRPCR.PCR_loc_trial(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{s})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', true, 'DataPointSize', 2);
xticklabels({'low', 'high', 'low', 'high'})
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
ylabel('prospective coding ratio')
title('reward related zone prospective coding')
makefigurepretty(gcf,1)
filename = 'draftfigure5_RRZ_PCR_velsplit';
savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

%%% playing around with example ideas for figure 5A
%maybe a heatmap of per trial? fast and slow? Idk
timeEdges = -0.2:0.02:0.2; %copied from params section of cf_decoding_thetaseq. ideally will resave structure with params files
posEdges = -81:3:99; %same as above
adjustedEdges = -180/2:3:180/2; %180 for now
exampledir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\SpeedSubsets_exampleSeq_perTrial\';
for t = 1:height(tmpData)
    figure
    hold on
    tmpSeq = tmpData.PostR_Var5{t};
    imagesc(timeEdges(1:end-1)+0.01, adjustedEdges(1:end-1), tmpData.PostR_Var5{t})
    plot(timeEdges(1:end-1), zeros(1,size(tmpSeq,2)), 'w--')
    plot(zeros(1,size(tmpSeq,1)), adjustedEdges(1:end-1), 'w--')
%     xlim([-0.08 0.08])
    ylim([-45 45])
    yticks([-45 0 45])
    xticks([-0.08 0 0.08])
    
    tmpAn = tmpData.animal(t);
    tmpDay = tmpData.day(t);
    tmpGrp = tmpData.group(t);
    
    if postRPCR.Ctrl_speed_subset(t) == 1
        speedappend = 'SLOW';
    else
        speedappend = 'FAST';
    end
    
    title([num2str(tmpAn), ' ', num2str(tmpDay), ' ', tmpGrp, ' - ', speedappend])
    filename = ['postR_thetaSeq_', num2str(tmpAn), '_', num2str(tmpDay), '_', speedappend, '_', num2str(t)];
    savefigALP(exampledir, filename, 'filetype', 'png', 'date', 1)
    close
end


%%% save data for plotting later
filename = 'figure5_data_ctrlspeed.mat';
save([figdir, filename], 'postRPCR', 'tmpData', '-v7.3')




%% get a boxplot of the post reward zone prospective coding per trial
%by speed threshold

otherdecedges = -81:3:99;
fnames = {'Ctrl_speed'};

for f = 1:length(fnames)
    clear g
    inclData = AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & AllData.significantSeq == 1;
    tmpData = AllData.(fnames{f})(inclData);
%     qs = quantile(tmpData, 3);
    qs = median(tmpData, 'omitnan');
    qs = [min(tmpData) qs max(tmpData)];
    iQ = discretize(tmpData, qs);
    PlotData = AllData(inclData,:);
    
    g(1,1) = gramm('x', PlotData.(fnames{f}),  'color', cellstr(PlotData.group), 'lightness', iQ, 'subset', PlotData.fullTrial==1 & PlotData.engaged ==1 & PlotData.rewarded == 1 & strcmp(PlotData.timepoint, 'post'),'column', PlotData.animal);
    g(1,1).facet_wrap(PlotData.animal);
    g(1,1).stat_bin('edges', qs, 'geom', 'bar')
%     g(1,1).stat_summary();
%     g(1,1).set_names('x','distance to reward (deg)', 'y', 'prospective coding ratio', 'lightness', 'Data Half');
    g(1,1).set_title(fnames{f});
%     g(1,1).axe_property('YLim', [-0.25 0.25]);
    
%     g(2,1) = gramm('x', iQ, 'y', AllData.PCR_loc_trial, 'color', cellstr(AllData.group), 'lightness', iQ, 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & ~isnan(AllData.PCR_loc_trial));
%     g(2,1).facet_grid([], cellstr(AllData.timepoint));
% %         g(2,1).stat_bin('geom', 'stairs', 'fill', 'transparent');
% %     g(2,1).stat_bin('normalization','cdf','geom','stairs');
% %     g(2,1).stat_violin();
%     g(2,1).stat_boxplot();
%     g(2,1).axe_property('YLim', [-0.25 0.25]);

    figure
    g.draw();
%     savefigALP(figdir, ['PCR_overposition_bybehaviorquarts_timepoint_', fnames{f}], 'filetype', 'pdf')
end


%% try plotting behavior by RZq quartile
%%% load behavior data for all trials
behaviordir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\']; 
behaviorfilename = 'behavioranalysis_prepost_recordings_table.mat';
load([behaviordir, behaviorfilename])

%%% get quartiles by PCR value
tmpData = AllData.PCR_loc_trial;
qs = quantile(tmpData, 3);
qs = [min(tmpData) qs max(tmpData)];
iQ = discretize(tmpData, qs);

%%% custom PCR groupings
PCRbins = [-1 -0.1 0.1 1];
[~,~, binID] = histcounts(tmpData, PCRbins);
iQ = binID; 

%%% plot
% ALP 5/18/23 seems like after flicker the extreme values (1st and 4th
% quartiles) have similar velocity behavior. Other potential things of
% interest are that before flicker more prospective coding seems to
% correspond to earlier licking via the histogram - maybe this corresponds
% to getting reward earlier? 
%
%After flicker it seems like the licking bump is more exaggerated after
%flicker, though the initial increase in licking is at the same position.
%Could this indicate moving though the reward zone more quickly? Not sure
clear g
g(1,1) = gramm('x', repmat(params.posEdges(1:end-1), [height(AllData),1]), 'y', allTMetrics.vel_h, 'lightness', iQ, 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1);
g(1,1).facet_grid([], cellstr(AllData.timepoint));
g(1,1).stat_summary();
g(1,1).set_names('x','distance to reward (deg)', 'y', 'velocity', 'lightness', 'quartile');

g(2,1) = gramm('x', repmat(params.posEdges(1:end-1), [height(AllData),1]), 'y', allTMetrics.lick_nh, 'lightness', iQ, 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1);
g(2,1).facet_grid([], cellstr(AllData.timepoint));
g(2,1).stat_summary();
g(2,1).set_names('x','distance to reward (deg)', 'y', 'licking', 'lightness', 'quartile');

figure
g.draw();

%% plot box plots of behavior metrics by pcr quartiles
otherdecedges = -81:3:99;
fnames = {'lickDI', 'nLicks', 'licklatency_s', 'duration', 'trial_speed', ...
    'AZ_speed', 'RZ_speed', 'Ctrl_speed', 'AZ_speed_slope', 'AZ_lick_slope'};
figdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 

for f = 1:length(fnames)
    tmpData = AllData.(fnames{f});
    
    clear g
    g(1,1) = gramm('x', iQ, 'y', tmpData, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & ~isnan(AllData.PCR_loc_trial));
    g(1,1).facet_grid([], cellstr(AllData.timepoint));
    g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
    g(1,1).stat_boxplot();
    g(1,1).set_names('x','RRZ QPCR', 'y', fnames{f});
    g(1,1).set_title(fnames{f});
    
    figure('Position', [376 362 984 547])
    g.draw();
    %     savefigALP(figdir, ['PCR_overposition_bybehaviorquarts_timepoint_', fnames{f}], 'filetype', 'pdf')
end


%% scatter plots
% I'm going to try some scatter plots of PCR vs. behavior metrics
% from looking at plots, it would be interesting to see. Jut going to do
% the ones I have now to start with . would also be nice to try and look at
% where the slowing starts? It seems like there is a fair number of trials
% where the slowing starts later. I wonder if that is essentially just more
% speed in the reward zone. If I wanted to quantify where it starts, I
% could try smoothing the velocity (too noisy w/out) and taking the
% derivative. get the first - slope within a region of interest and find
% the position. from looking at these plots I would think more extreme PCR
% means higher AZ speed and more PCR mean more RZ speed? and more overall
% trial speed
fnames = {'lickDI', 'nLicks', 'licklatency_s', 'duration', 'trial_speed', ...
    'AZ_speed', 'RZ_speed', 'Ctrl_speed', 'AZ_speed_slope', 'AZ_lick_slope'};

for f = 1:length(fnames)
    tmpData = AllData.(fnames{f});
    
    clear g
    g(1,1) = gramm('x', AllData.PCR_loc_trial, 'y', tmpData, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1);
    g(1,1).facet_grid([], cellstr(AllData.timepoint));
    g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
    g(1,1).geom_point('alpha', 0);
    g(1,1).set_names('x','RRZ QPCR', 'y', fnames{f});
     g(1,1).set_title(fnames{f});

    figure('Position', [410.0000 169 1029 577])
    g.draw();
%     savefigALP(figdir, ['Behavior_scatter_byPCRquarts_timepoint_', fnames{f}], 'filetype', 'pdf')
end

%% scatter plots - as 2d histograms - POST ONLY - 40 vs. Random

fnames = {'lickDI', 'nLicks', 'licklatency_s', 'duration', 'trial_speed', ...
    'AZ_speed', 'RZ_speed', 'Ctrl_speed', 'AZ_speed_slope', 'AZ_lick_slope'};

yedgesarray = {[-1:0.1:1], [0:20:1000], [0:0.05:1], [0:10:200], [0:0.5:10], ...
    [0:0.5:10], [0:0.5:10], [0:0.5:10], [-1:0.1:1], [-1:0.1:1]};

for f = 1:length(fnames)
    tmpData = AllData.(fnames{f});
    
    clear g
    g(1,1) = gramm('x', AllData.PCR_loc_trial, 'y', tmpData, 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post'));
    g(1,1).facet_grid([], cellstr(AllData.group));
    g(1,1).stat_bin2d('geom', 'image', 'edges', {-1:0.05:1, yedgesarray{f}});
    g(1,1).set_names('x','RRZ QPCR', 'y', fnames{f});
    g(1,1).set_title([fnames{f}, ' - Post flicker']);

    figure('Position', [410.0000 169 1029 577])
    g.draw();
%     savefigALP(figdir, ['Behavior_byPCR_2dhist_postonly_', fnames{f}], 'filetype', 'pdf')

end

%% calculate spearmans rho for a variety of possible behavior metrics

fnames = {'lickDI', 'nLicks', 'licklatency_s', 'duration', 'trial_speed', ...
    'AZ_speed', 'RZ_speed', 'Ctrl_speed', 'AZ_speed_slope', 'AZ_lick_slope'};

inclData = AllData.fullTrial == 1 & AllData.engaged == 1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & ~isnan(AllData.PCR_loc_trial); 
corrData = AllData(inclData,:);

x = corrData.PCR_loc_trial; 
for f = 1:length(fnames)
    y = [];
    y = corrData.(fnames{f}); 
    [rho, pval] = corr(x,y, 'Type', 'Spearman');

    disp(['spearman correleation for ', fnames{f}, ' rho = ', num2str(rho), '; p = ', num2str(pval)])
end


%% time below threshold in the reward related zone
%this is for the metrics that are only contained in the AllTMetrics table
fnames = {'RRZ_time_below_thresh', 'time_2_nextRZ', 'low_speed_pos', 'high_speed_pos', 'RRZ_time_below_thresh_prop'};

yedgesarray = {[0:5:150], [0:5:150], [-18:2:18], [-18:2:50], [0:0.1:1]};

%%% 2 D plot
for f = 1:length(fnames)
    tmpData = allTMetrics.(fnames{f});
    clear g
    
    g(1,1) = gramm('x', AllData.PCR_loc_trial, 'y', tmpData, 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1);
    g(1,1).facet_grid([], cellstr(AllData.group));
    g(1,1).stat_bin2d('geom', 'image',  'edges', {-1:0.05:1, yedgesarray{f}});
    g(1,1).set_names('x','RRZ QPCR', 'y', fnames{f});
    g(1,1).set_title(fnames{f});
    figure('Position', [410.0000 169 1029 577])
    g.draw();
%     savefigALP(figdir, ['Behavior_byPCR_2dhist_postonly_', fnames{f}], 'filetype', 'pdf')
end

%%% scatter plots
for f = 1:length(fnames)
    tmpData = allTMetrics.(fnames{f});
    
    clear g
    g(1,1) = gramm('x', AllData.PCR_loc_trial, 'y', tmpData, 'color', cellstr(AllData.group), 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1);
    g(1,1).facet_grid([], cellstr(AllData.timepoint));
    g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
    g(1,1).geom_point();
    g(1,1).set_names('x','RRZ QPCR', 'y', fnames{f});
    g(1,1).set_title(fnames{f});
    
    figure('Position', [410.0000 169 1029 577])
    g.draw();
%     savefigALP(figdir, ['Behavior_scatter_byPCRquarts_timepoint_', fnames{f}], 'filetype', 'pdf')
end

%%% box plots
for f = 1:length(fnames)
    tmpData = allTMetrics.(fnames{f});
    
    clear g

    g(1,1) = gramm('x', iQ, 'y', tmpData, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & ~isnan(AllData.PCR_loc_trial));
    g(1,1).facet_grid([], cellstr(AllData.timepoint));
    g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
    g(1,1).stat_boxplot();
    g(1,1).set_names('x','RRZ QPCR quartiles', 'y', fnames{f});
    g(1,1).set_title(fnames{f});    
    
    figure('Position', [410.0000 169 1029 577])
    g.draw();
%     savefigALP(figdir, ['Behavior_boxplot_byPCRbins_', fnames{f}], 'filetype', 'pdf')
end

%% just 40 vs random box plots
fnames = {'RRZ_time_below_thresh', 'postRZ_speed'};
ynames = {'time below threshold (s)', 'speed (deg/s)'};
tnames = {'time below velocity threshold in reward areas', 'speed after the reward zone, 40-55'};
for f = 1:length(fnames)
    tmpData = allTMetrics.(fnames{f});
    
    clear g

    g(1,1) = gramm('x', iQ, 'y', tmpData, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & ~isnan(AllData.PCR_loc_trial) & strcmp(AllData.timepoint, 'post'));
    g(1,1).facet_grid([], cellstr(AllData.group));
    g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
    g(1,1).stat_boxplot('width', 1.5);
    g(1,1).set_names('x',{'reward related zone', 'prospective coding ratio bins'}, 'y', ynames{f});
    g(1,1).set_title(tnames{f});    
    
    figure('Position', [410.0000 169 1029 577])
    g.draw();
%     savefigALP(figdir, ['Behavior_boxplot_byPCRbins_', fnames{f}], 'filetype', 'pdf')
end

%% create table data for testing in R
statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';

PCR_bin = iQ;
PCR_bin_data = table(PCR_bin);
behaviorSeqData = table();
behaviorSeqData = [PCR_bin_data(:,1) AllData(:,31) allTMetrics(:,15:23) AllData(:,1:4) AllData(:,7:9) allTMetrics(:,2:3) AllData(:,36)]; 

filename = 'TableData_ThetaSequences_withBehavior_perTrial.txt';
writetable(behaviorSeqData, fullfile(statsdir, filename))


%% get stats for RRZ time below threshold
%just ranksum for a quick test

inclData = AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & strcmp(AllData.group, 'random');
p = ranksum(allTMetrics.RRZ_time_below_thresh(inclData & PCR_bin ==1), allTMetrics.RRZ_time_below_thresh(inclData & PCR_bin ==2));
p = ranksum(allTMetrics.RRZ_time_below_thresh(inclData & PCR_bin ==2), allTMetrics.RRZ_time_below_thresh(inclData & PCR_bin ==3));



%% get per day per PCR subset to plot
inclData = AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & PCR_bin ~= 0;

goodData = behaviorSeqData(inclData,:);
goodDataDay = grpstats(goodData, {'PCR_bin', 'day', 'group'}, "mean",...
    "DataVars", ["vel_nh", 'RRZ_time_below_thresh', 'time_2_nextRZ', 'low_speed_pos', 'high_speed_pos', 'RRZ_time_below_thresh_prop', 'postRZ_speed']);

clear g
g = gramm('x', repmat(params.posEdges(1:end-1), [height(goodDataDay),1]), 'y', goodDataDay.mean_vel_nh, 'color', cellstr(goodDataDay.group), 'lightness', goodDataDay.PCR_bin);
g(1,1).facet_grid([], cellstr(goodDataDay.group));
g(1,1).stat_summary('type', 'sem');
g(1,1).set_names('x','distance to reward (deg)', 'y', 'velocity(deg/s)', 'lightness', 'PCR_bin');
g(1,1).set_title('velocity by prospective coding ratio bin');
% g(1,1).axe_property('YLim', [-0.25 0.25]);

figure('Position', [410.0000 169 1029 577])
g.draw();

%% get figure of per animal averages across PCR bins
inclData = AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & PCR_bin ~= 0;
goodData = behaviorSeqData(inclData,:);
goodDataAn = grpstats(goodData, {'PCR_bin', 'animal', 'group'}, "mean",...
    "DataVars", ["vel_h", 'RRZ_time_below_thresh', 'time_2_nextRZ', 'low_speed_pos', 'high_speed_pos', 'RRZ_time_below_thresh_prop', 'postRZ_speed']);

%%% get difference between bins 2 and 1 and bins 2 and 3


clear g
g = gramm('x', repmat(params.posEdges(1:end-1), [height(goodDataAn),1]), 'y', goodDataAn.mean_vel_h, 'color', cellstr(goodDataAn.group), 'lightness', goodDataAn.PCR_bin, 'column', goodDataAn.animal);
g(1,1).facet_wrap(goodDataAn.animal);
g(1,1).geom_line();
g(1,1).set_names('x','distance to reward (deg)', 'y', 'velocity(deg/s)', 'lightness', 'PCR_bin');
g(1,1).set_title('velocity by prospective coding ratio bin');
% g(1,1).axe_property('YLim', [-0.25 0.25]);

figure('Position', [410.0000 169 1029 577])
g.draw();

%% get differences 

diffbin = []; diffan = []; diffgroup = []; diffvelocity = [];
animals = unique(goodDataAn.animal);
for an = 1:length(animals)
    isAnimal = goodDataAn.animal == animals(an);
    anData = goodDataAn(isAnimal,:);
    for b = [1,3]
        isBin2 = anData.PCR_bin == 2;
        isBinB = anData.PCR_bin == b;
        
        if sum(isBin2) < 1 || sum(isBinB) < 1
            continue
        end
        
        tmpdiff = anData.mean_vel_h(isBin2,:) - anData.mean_vel_h(isBinB,:);
        
        diffbin = [diffbin; b];
        diffan = [diffan; animals(an)];
        diffgroup = [diffgroup; anData.group(1)];
        diffvelocity = [diffvelocity; tmpdiff];
    end
end

diffData = table(diffbin, diffan, diffgroup, diffvelocity);

clear g
g = gramm('x', repmat(params.posEdges(1:end-1), [height(diffData),1]), 'y', diffData.diffvelocity, ...
    'color', cellstr(diffData.diffgroup), 'lightness', diffData.diffbin, 'column', diffData.diffan);
g(1,1).facet_wrap(diffData.diffan);
g(1,1).geom_hline('yintercept', 0, 'style', 'k-');
g(1,1).geom_line();
g(1,1).set_names('x','distance to reward (deg)', 'y', 'velocity differences(deg/s)', 'lightness', 'PCR bin diff');
g(1,1).set_title('velocity by prospective coding ratio bin');
g(1,1).axe_property('YLim', [-5 5])
% g(1,1).axe_property('YLim', [-0.25 0.25]);
%     g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);

figure('Position', [410.0000 169 1029 577])
g.draw();


%% plot paired difference per session
days = unique(goodDataDay.day); 
binpairs = [2 1; 2 3];

colors{1} = [gcolors(1,:); gcolors(1,:)];
colors{2} = [gcolors(3,:); gcolors(3,:)];

fnames = {'RRZ_time_below_thresh', 'time_2_nextRZ', 'low_speed_pos', 'high_speed_pos', 'postRZ_speed'};
ynames = {'time (seconds)', 'time(seconds)', 'position (deg)', 'position (deg)', 'speed (deg/s)'};

for f = 1:length(fnames)
    
    for g = 1:2
        for b = 1:2
            binnedDay{g}{b} = [];
        end
    end
    
    for d = 1:length(days)
        isDay = goodDataDay.day == days(d);
        dayDat = goodDataDay(isDay,:);
        tmpbins = goodDataDay.PCR_bin(isDay);
        tmpVals = dayDat.(['mean_', fnames{f}]);
        tmpGroup = dayDat(1,:).group;
        if strcmp(tmpGroup, 'gamma')
            g = 1;
        else
            g = 2;
        end
        
        if length(dayDat.PCR_bin) < 3
            continue
        end
        for b = 1:size(binpairs,1)
            tmpDiff = tmpVals(binpairs(b,1)) - tmpVals(binpairs(b,2));
            binnedDay{g}{b} = [binnedDay{g}{b}; tmpDiff];
        end
    end
    
    xvect = [1:2; 4:5]; 
    fh = figure('Position', [558 272 410 382]);
    hold on
    for g = 1:2
        binMn = []; binSem = []; scatterDat = [];
        for b = 1:size(binpairs,1)
            binMn(b) = mean(binnedDay{g}{b}, 'omitnan');
            binSem(b) = nanstd(binnedDay{g}{b})./sqrt(sum(~isnan(binnedDay{g}{b})));
            scatterdat{b} = binnedDay{g}{b};
        end
        
        plotprettypoints(fh, xvect(g,:), scatterdat, 'color', colors{g})
        b = bar(xvect(g,:), binMn, 'FaceColor', 'flat');
        b.CData = colors{g};
        b.FaceAlpha = 0.6;
        errorbar2(xvect(g,:), binMn, binSem, 0.2, 'k-', 'LineWidth', 0.75);
        
    end
    title(fnames{f})
    xticks([1,2,4,5])
    xticklabels({'2-1', '2-3', '2-1', '2-3'})
    ylabel(ynames{f})
    
    clear tmp*
    
end


%% try plotting unengaged trials

%how many are there?
g = gramm('x', AllData.licklatency_s, 'y', AllData.PCR_loc_trial, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1);
% g = gramm('x', AllData.PCR_loc_trial, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1);
g.facet_grid([], cellstr(AllData.timepoint));
g.geom_point();
% g.stat_bin('geom', 'stairs', 'fill', 'transparent');
g.set_names('x','RRZ prospective coding ratio', 'y', 'trial count', 'color','Group');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
figure
g.draw();

% I tested a few different AZ speed thresholds to see if the speed in the
% anticipatory zone (which does seem to be different between groups) is
% driving the differences in prospective coding. when I speed match in this
% way the differences in PCR over the track remain (tested April 2023
% sometime)
g = gramm('x', repmat(params.dec_edges(1:end-1), [height(AllData),1]), 'y', AllData.PCR_pos_trial, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & AllData.AZ_speed > 5);
g.facet_grid([], cellstr(AllData.timepoint));
g.stat_summary();
g.set_names('x','distance to reward (deg)', 'y', 'prospective coding ratio', 'color','Group');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g(1,1).axe_property('YLim', [-0.25 0.25]);
figure
g.draw();

 g = gramm('x', cellstr(AllData.timepoint), 'y', AllData.PCR_pos_trial, 'color', AllData.engaged, 'subset', AllData.fullTrial == 1);
% g.stat_violin('normalization','area');
g.stat_boxplot();
figure
g.draw();
% %% playing around with group stats
% GoodTrials = AllData(AllData.fullTrial&AllData.rewarded&AllData.engaged,:);
% A = grpstats(GoodTrials, ["group", "timepoint"], ["mean", "sem"], "DataVars", 29); 
% B = grpstats(GoodTrials, ["day"], ["mean"], "DataVars", 28:29); 
% 
% tp = metadata.FlickerDay; 
% timepoint(tp <5) = "pre";
% timepoint(tp>5) = "post"; 
% timepoint = convertCharsToStrings(timepoint)';
% 
% DayData = table;
% DayData = [DayData metadata(:, 5)];
% DayData = addvars(DayData, timepoint);
% DayData = [DayData B]; 

%% try plotting unengaged trials 
g(1,1) = gramm('x', AllData.PCR_loc_trial, 'color', AllData.engaged, 'subset', AllData.fullTrial==1);
g(2,1) = copy(g(1));
g(1,1).facet_grid([], cellstr(AllData.timepoint));
% g.geom_point();
g(1,1).stat_bin('geom', 'stairs', 'normalization', 'probability', 'fill', 'transparent');
g(1,1).set_names('x','RRZ prospective coding ratio','color','Engaged');
% g.set_color_options('map', trialscolors, 'n_color', 2, 'n_lightness', 1);

g(2,1).facet_grid([], cellstr(AllData.timepoint));
g(2,1).stat_bin('normalization','cdf','geom','stairs');
g(2,1).set_title('''normalization'',''cdf''','FontSize',10);
g(1,1).set_names('x','RRZ prospective coding ratio','color','Engaged');
figure
g.draw();

g = gramm('x', AllData.AZ_speed_slope, 'color', AllData.engaged, 'subset', AllData.fullTrial==1);
g.facet_grid([], cellstr(AllData.timepoint));
% g.geom_point();
g.stat_bin('geom', 'stairs', 'normalization', 'probability', 'fill', 'transparent');
g.set_names('x','lick discrimination','color','Group');
% g.set_color_options('map', trialscolors, 'n_color', 2, 'n_lightness', 1);
figure
g.draw();

g = gramm('x', AllData.trial_speed, 'color', AllData.engaged, 'subset', AllData.fullTrial==1);
g.facet_grid([], cellstr(AllData.timepoint));
% g.geom_point();
g.stat_bin('geom', 'stairs', 'normalization', 'probability', 'fill', 'transparent');
g.set_names('x','Trial speed','color','Group');
% g.set_color_options('map', trialscolors, 'n_color', 2, 'n_lightness', 1);
figure
g.draw();

%custom engagement grouping
g = gramm('x', repmat(params.dec_edges(1:end-1), [height(AllData),1]), 'y', AllData.PCR_pos_trial, 'color', AllData.engaged, 'subset', AllData.fullTrial==1);
g.facet_grid([], cellstr(AllData.timepoint));
g.stat_summary();
g(1,1).set_names('x','distance to reward (deg)', 'y', 'prospective coding ratio', 'color', 'engaged');
g(1,1).set_title('prospective coding, engaged vs. unengaged trials');
g(1,1).axe_property('YLim', [-0.25 0.25]);
figure
g.draw();

%% figure out how to merge the day averages from the table structure

mnGoodData = cat(3, GoodTrials.mnSeq{:});
mnUnengagedData = cat(3, AllData.mnSeq{AllData.engaged ==0 & AllData.fullTrial == 1});

figure
hold on
subplot(1,2,1)
hold on
imagesc(mean(mnGoodData,3, 'omitnan'), [0.015 0.03])
title('engaged and correct')
colorbar
subplot(1,2,2)
hold on
imagesc(mean(mnUnengagedData,3, 'omitnan'), [0.015 0.03])
colorbar
title('unengaged trials')



end

