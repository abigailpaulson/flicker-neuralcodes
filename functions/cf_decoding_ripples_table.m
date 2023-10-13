function cf_decoding_ripples_table(dirs, params, allindex, metadata, posType, cellTypes)
%cf_decoding_ripples
%   get bayesian decoding during ripples
%ALP 1/9/2023

rewrite.decodingfiles = 0;

%%% set any necessary parameters
params.minRippleDur = 0.015; % in s
params.timeAroundMid = 0.125; % in s, on either side of mid
params.decodingBins = 0.025; %in s
if strcmp(posType, 'full')
    params.posEdges = 0:2:360;
    params.posBins = 2;
    dType = '360';
    nDeg = 360;
    posName = 'theta';
    params.RZ = [54 72; 234 252];
    dec_edges = 0:6:360;
else
    params.posEdges = -81:3:99;
    params.posBins = 3;
    dType = '180';
    params.nDeg = 180;
    posName = 'theta_d2r';
    params.RZ = [0 18];
    dec_edges = -81:9:99; %could also try 9
    params.adjustedEdges = -90:3:90;
end

%%% directories
ripplechdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\'];
positioninfodir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\'];
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';
savedir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_ripples\'];
dayindex = unique(allindex(:,1:2), 'rows');

savefilename = ['rippledecoding_halfratiozones_', dType, '_', cellTypes, '_AllData', '.mat'];
statsdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\'];

if ~exist([savedir, savefilename], 'file') || rewrite.decodingfiles
    
    %%% loop over all days
    AllData = table;
    for d = 1:size(dayindex,1)
        dayData = [];
        disp(['getting ripple analyses for ', num2str(dayindex(d,1)), ' ', num2str(dayindex(d,2))])
        files = allindex(allindex(:,2) == dayindex(d,2),3); %should be only VR files pre flicker or all VR on day 9/10
        sessindex = dayindex(d,:);
        index = repmat(sessindex, length(files), 1);
        index = [index files];
        anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
        
        dayfilename = ['rippledecoding_halfratiozones_', dType, '_', cellTypes, '_' num2str(dayindex(d,2)), '.mat'];
        
        %%% load files
        %best ripple ch
        load([ripplechdir, 'CA1\bestChannel_', num2str(dayindex(d,2)), '.mat'])
        tmpRChan = bestRippleChan.all;
        ripples = loaddatastruct2([anprocesseddatadir, 'CA1\', num2str(tmpRChan), '\'], dayindex(d,:), 'ripples', files);
        nonthetas = loaddatastruct2([anprocesseddatadir, 'CA1\', num2str(tmpRChan), '\'], dayindex(d,:), 'nonthetas', files);
        
        %behaviorInfo
        load([positioninfodir, 'behaviorInfo_A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
        
        %spikes
        [spikes, unitInfo, spikevect] = cf_getspikes(anprocesseddatadir, dirs, params, index, params.brainReg);
        nCells = length(unitInfo);
        cellIDs = 1:1:nCells;
        
        %PC or other cell type information
        if strcmp(cellTypes, 'PC')
            load([ratemapdir, 'spatialmaps_all_500_full_', num2str(dayindex(d,2)), '.mat'])
            inclFR = ratemaps.includeCell;
            inclCT = ratemaps.isPYR;
            inclCells = cellIDs(inclFR & inclCT);
            nCells = length(inclCells);
        elseif strcmp(cellTypes, 'all')
            inclCells = cellIDs;
        end
        
        %%% training data
        %%% load spatial rate maps - training data
        load([ratemapdir, 'spatialmaps_preflicker_0_', posType, '_', num2str(dayindex(d,2)), '.mat'])
        trainingData.ratemaps = ratemaps.fr; %in Hz
        
        for i = 1:size(index,1)
            %%% choose ripples of interest
            tmpRDat= ripples{index(i,1)}{index(i,2)}{index(i,3)};
            lfpsamprate = tmpRDat.samprate;
            tmpRippleI = [tmpRDat.startind tmpRDat.midind tmpRDat.endind];
            if isempty(tmpRippleI)
                continue
            end
            tmpMidI = tmpRippleI(:,2);
            NTInds = [nonthetas{index(i,1)}{index(i,2)}{index(i,3)}.startind nonthetas{index(i,1)}{index(i,2)}{index(i,3)}.endind];
            isNTRipple = isExcluded(tmpMidI, NTInds);
            isNTRipple = logical(isNTRipple);
            tmpSize = tmpRDat.maxthresh;
            
            tmpBehavior = behaviorInfo(index(i,3));
            
            %%% continue to the next loop if no ripples or no behavior data on
            %%% this day
            if sum(isNTRipple) == 0 || isempty(tmpBehavior.(posName))
                continue
            end
            
            tmpDuration = tmpRippleI(:,3) - tmpRippleI(:,1);
            tmpDurationS = tmpDuration./lfpsamprate;
            isLongRipple = tmpDurationS > params.minRippleDur; %must be longer than 15 ms
            
            inclRipple = isNTRipple&isLongRipple;
            rippleInds = tmpRippleI(inclRipple,:); %startind mid ind
            
            rippleMids = [rippleInds(:,2)./lfpsamprate];
            rippleWindows = [rippleInds(:,2)./lfpsamprate-params.timeAroundMid rippleInds(:,2)./lfpsamprate+params.timeAroundMid];
            
            %%% get gamma power for each ripple
            gammapower_CA1_z = cf_getgammadurperiod(anprocesseddatadir, index, rippleInds, lfpsamprate, 'CA1');
            gammapower_CA3_z = cf_getgammadurperiod(anprocesseddatadir, index, rippleInds, lfpsamprate, 'CA3');
            
            %%% get spikes only of included cells
            tmpSpikeIDs = spikevect(index(i,3)).spikeIDs;
            inclSpikes = ismember(tmpSpikeIDs, inclCells);
            tmpSpikeIDs = tmpSpikeIDs(inclSpikes);
            tmpSpikeTimes = spikevect(index(i,3)).spikeTimes(inclSpikes);
            
            tmpSpikeStruct.spiketimes = tmpSpikeTimes;
            tmpSpikeStruct.spikeIDs = tmpSpikeIDs;
            
            decoding_bins = cf_rippledecoding_timebins(index(i,:), trainingData, tmpBehavior.time, tmpBehavior.(posName), tmpSpikeStruct, rippleWindows, rippleMids, params, inclCells, nCells);
            decoding_window = cf_rippledecoding_window(index(i,:), trainingData, tmpBehavior.time, tmpBehavior.(posName), tmpSpikeStruct, rippleWindows, rippleMids, params, inclCells, nCells);
            
            rippleposition = decoding_bins.ripplePos;
            rippleRatio = decoding_bins.rippleRatio;
            rippleDecPos = decoding_window.rippleDecPos;
            rippleRelPos = decoding_window.rippleRelPos;
            rippleDurationS = tmpDurationS(inclRipple);
            rippleSize = tmpSize(inclRipple);
            dayData(i).decoding_bins = decoding_bins;
            dayData(i).decoding_window = decoding_window;
            
            animal = index(i,1).*ones(length(rippleposition),1);
            day = index(i,2).*ones(length(rippleposition),1);
            if metadata.FlickerDay(d) < 5
                timepoint = repmat({'pre'}, length(rippleposition), 1);
            else
                timepoint = repmat({'post'}, length(rippleposition), 1);
            end
            group = repmat(metadata.Groups(d), length(rippleposition),1);
            
            tmpData = table(animal, day, timepoint, group, rippleposition, ...
                rippleRatio, rippleDecPos, rippleRelPos, gammapower_CA1_z, ...
                gammapower_CA3_z, rippleDurationS, rippleSize);
            AllData = [AllData; tmpData];
            
            
            clear tmp* decoding_bins decoding_window rippleMids rippleWindows tmpData
            close all
        end
        
    end
    
    %%% save group data
    %load curr pos data to exclude bad decoding days
    if strcmp(cellTypes, 'all')
        currPos = load('\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_currentposition\Group_decoding_currpos_180_all_test20_210527.mat');
    else
        currPos = load('\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_currentposition\Group_decoding_currpos_180_PC_test20.mat');
    end
    inclDays = [currPos.groupData.include]';
    
    animal = dayindex(:,1);
    day  = dayindex(:,2);
    group = metadata.Groups;
    timepoint = repmat({[]}, length(day),1);
    isPre = metadata.FlickerDay < 5;
    timepoint(isPre) = repmat({'pre'}, sum(isPre),1);
    timepoint(~isPre) = repmat({'post'}, sum(~isPre),1);
    includeDay = inclDays;
    
    
    dayInfo = table(animal, day, group, timepoint, includeDay);
    AllData = outerjoin(AllData, dayInfo, 'MergeKeys', true);
    
    AllData.includeRipple = AllData.includeDay == 1 & ~isnan(AllData.rippleposition) & ~isnan(AllData.rippleRatio) & ~isnan(AllData.rippleDecPos) & ~isnan(AllData.rippleRelPos);
    AllData.RZRipple = AllData.rippleposition >= 0 & AllData.rippleposition<18;
    
    save([savedir, savefilename], 'AllData', 'dayInfo', '-v7.3')
    
    statsfilename = 'TableData_RippleDecoding_SWR.txt';
    writetable(AllData, fullfile(statsdir, statsfilename))
else
    load([savedir, savefilename])
end

%%% get day data
% what are the things I am interested in?
% prop ripples PCR > 0 and prop ripples decoded position difference > 0?
DayData = table;
ripEdges = -99:12:99;
for d = 1:size(dayindex,1)
    incldata = AllData.day == dayindex(d,2) & AllData.includeRipple;
    tmpData = AllData(incldata,:);
    
    prospRipples = sum(tmpData.rippleRatio > 0.1)./length(tmpData.rippleRatio);
    prospRipples_RZ = sum(tmpData.rippleRatio(tmpData.RZRipple) > 0.1)./sum(tmpData.RZRipple); %played around with a few values but changing it back to 0.1 ALP 10/10/23
    
    propForwardRipples = sum(tmpData.rippleRelPos > 0)./length(tmpData.rippleRelPos);
    propForwardRipples_RZ = sum(tmpData.rippleRelPos(tmpData.RZRipple) > 0)./sum(tmpData.RZRipple);
    
    propFarForwardRipples = sum(tmpData.rippleRelPos >= 18)./length(tmpData.rippleRelPos);
    propFarForwardRipples_RZ = sum(tmpData.rippleRelPos(tmpData.RZRipple) >= 18)./sum(tmpData.RZRipple);
    
    propMidRipples = sum(tmpData.rippleRelPos >= 18 & tmpData.rippleRelPos < 36)./length(tmpData.rippleRelPos);
    propMidRipples_RZ = sum(tmpData.rippleRelPos(tmpData.RZRipple) >= 18 & tmpData.rippleRelPos(tmpData.RZRipple) < 36)./sum(tmpData.RZRipple);
    
    propNearRipples = sum(tmpData.rippleRelPos < 18 & tmpData.rippleRelPos >= -18)./length(tmpData.rippleRelPos);
    propFrontRipples = sum(tmpData.rippleRelPos >= 18)./length(tmpData.rippleRelPos);
    propBackRipples = sum(tmpData.rippleRelPos < -18)./length(tmpData.rippleRelPos);
    propFarRipples = sum(tmpData.rippleRelPos > 18 | tmpData.rippleRelPos < -18)./length(tmpData.rippleRelPos);
    
    isRZ = tmpData.RZRipple;
    propNearRipples_RZ = sum(tmpData.rippleRelPos(isRZ) < 18 & tmpData.rippleRelPos(isRZ) >= -18)./length(tmpData.rippleRelPos(isRZ));
    propFarRipples_RZ = sum(tmpData.rippleRelPos(isRZ) > 18 | tmpData.rippleRelPos(isRZ) < -18)./length(tmpData.rippleRelPos(isRZ));
    propFrontRipples_RZ = sum(tmpData.rippleRelPos(isRZ) >= 18)./length(tmpData.rippleRelPos(isRZ));
    propBackRipples_RZ = sum(tmpData.rippleRelPos(isRZ) < -18)./length(tmpData.rippleRelPos(isRZ));
    
    
    %new rel pos
    newRelPos = tmpData.rippleRelPos;
    newRelPos(newRelPos > 90) = newRelPos(newRelPos > 90) - 180;
    newRelPos(newRelPos < -90) = newRelPos(newRelPos < -90) + 180;
    avgRelPos = nanmean(newRelPos);
    avgRelPos_RZ = nanmean(newRelPos(isRZ));
    
    tmpRelPos = tmpData.rippleDecPos(isRZ);
    rippleRelPos_nh = histcounts(tmpRelPos, ripEdges)./sum(isRZ);
    
    nRipples = length(tmpData.rippleRatio);
    nRipples_RZ = sum(tmpData.RZRipple);
    
    animal = dayindex(d,1);
    day = dayindex(d,2);
    
    tmpDayData = table(animal, day, propNearRipples, propFarRipples, prospRipples, ...
        prospRipples_RZ, propForwardRipples, propForwardRipples_RZ, propFarForwardRipples, ...
        propFarForwardRipples_RZ, nRipples, nRipples_RZ, propMidRipples, propMidRipples_RZ, ...
        propNearRipples_RZ, propFarRipples_RZ, propBackRipples, propFrontRipples, propBackRipples_RZ, propFrontRipples_RZ, rippleRelPos_nh, avgRelPos, avgRelPos_RZ);
    
    DayData = [DayData; tmpDayData];
end
DayData = outerjoin(DayData, dayInfo, 'MergeKeys', true);
DayData = DayData(DayData.includeDay & DayData.nRipples > 5,:); 

savefilename = ['rippledecoding_halfratiozones_', dType, '_', cellTypes, '_DayData', '.mat'];
save([savedir, savefilename], 'DayData', 'dayInfo', '-v7.3')

statsfilename = 'TableData_RippleDecoding_SWR_Day.txt';
writetable(DayData, fullfile(statsdir, statsfilename))

%% plotting stuff
groupcolors = cbrewer('qual', 'Paired', 6);
gcolors = groupcolors([2,1,4,3],:); %rearranging it cuz of the way that the gramm does the plotting
trialscolors = groupcolors(end-1:end,:);

%% day plots
%%% day plots
dayfieldnames = {'prospRipples', 'prospRipples_RZ', 'propForwardRipples', 'propForwardRipples_RZ', ...
    'propFarForwardRipples', 'propFarForwardRipples_RZ', 'nRipples', 'nRipples_RZ', 'propMidRipples', ...
    'propMidRipples_RZ', 'propNearRipples', 'propFarRipples', 'propFrontRipples', 'propBackRipples',...
     'propNearRipples_RZ', 'propFarRipples_RZ', 'propFrontRipples_RZ', 'propBackRipples_RZ', 'avgRelPos', 'avgRelPos_RZ'};

for f = 1:2:numel(dayfieldnames)
    fdat = DayData.(dayfieldnames{f});
    g(1,1) = gramm('x', cellstr(DayData.timepoint), 'y', fdat, ...
        'color', cellstr(DayData.group), 'lightness', cellstr(DayData.timepoint), ...
        'subset', DayData.includeDay == 1);
    g(2,1) = copy(g(1,1));
    g(1,1).stat_summary('type', 'sem', 'geom', {'bar', 'black_errorbar'}, 'setylim', true);
    g(1,1).set_title([dayfieldnames{f}, ' per day']);
    g(1,1).set_names('x','timepoint', 'y', 'fraction of ripples' , 'color', 'Group', 'lightness', 'timepoint');
    g(2,1).geom_jitter('width',0.4,'height',0);
    
    fdat = DayData.(dayfieldnames{f+1});
    g(1,2) = gramm('x', cellstr(DayData.timepoint), 'y', fdat, ...
        'color', cellstr(DayData.group), 'lightness', cellstr(DayData.timepoint), ...
        'subset', DayData.includeDay == 1);
    g(2,2) = copy(g(1,2));
    g(1,2).stat_summary('type', 'sem', 'geom', {'bar', 'black_errorbar'}, 'setylim', true);
    g(1,2).set_title([dayfieldnames{f+1}, ' per day']);
    g(1,2).set_names('x','timepoint', 'y', 'fraction of ripples' , 'color', 'Group', 'lightness', 'timepoint');
    g(1,2).set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
    g(2,2).geom_jitter('width',0.4,'height',0);
    
    g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
    
    figure('Position', [365 532 1040 474])
    g.draw();
end

%% plot ripple rel pos over the track
isPost = strcmp(DayData.timepoint, 'post');
PostDayData = DayData(isPost,:);
gnames = {'random', 'gamma'};
figure
hold on

for g = 1:2
    isG = strcmp(PostDayData.group, gnames{g});
    tmpDat = PostDayData.rippleRelPos_nh(isG,:);
    
    tmpMn = nanmean(tmpDat,1);
    tmpSEM = nanstd(tmpDat,0,1)./sqrt(sum(~isnan(tmpDat(:,1))));
    shadedErrorBar(ripEdges(1:end-1), tmpMn, tmpSEM, {'Color', params.colors.(gnames{g}).post, 'LineWidth', 2}, 1)

end


%% Per ripple plots go here
%%% replicate plot here
inclRipples = AllData.includeDay == 1 & ~isnan(AllData.rippleposition) & ~isnan(AllData.rippleRatio) & ~isnan(AllData.rippleDecPos) & ~isnan(AllData.rippleRelPos);
plotData = AllData(inclRipples,:);

%%% plot ripple positions
g = gramm('x', AllData.rippleposition, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', inclRipples);
g.facet_grid([], cellstr(AllData.timepoint));
%g.geom_point();
g.stat_bin('geom', 'stairs', 'fill', 'transparent');
g.set_names('x','ripple position - all ripples', 'color','Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
figure('Position', [365 532 1040 324])
g.draw();

g = gramm('x', AllData.rippleDecPos, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', inclRipples);
g.facet_grid([], cellstr(AllData.timepoint));
%g.geom_point();
g.stat_bin('geom', 'stairs', 'fill', 'transparent');
g.set_names('x','ripple relative position - all ripples', 'color','Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
figure('Position', [365 532 1040 324])
g.draw();

g = gramm('x', AllData.rippleRatio, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', inclRipples);
g.facet_grid([], cellstr(AllData.timepoint));
%g.geom_point();
g.stat_bin('geom', 'stairs', 'fill', 'transparent');
g.set_names('x','ripple PCR ratio - all ripples', 'color','Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
figure
g.draw();

g = gramm('x', AllData.rippleRatio, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', inclRipples & AllData.rippleposition >=0 & AllData.rippleposition < 18);
g.facet_grid([], cellstr(AllData.timepoint));
%g.geom_point();
g.stat_bin('geom', 'stairs', 'fill', 'transparent');
g.set_names('x','ripple PCR ratio - RZ ripples', 'color','Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
figure
g.draw();


gfieldnames = {'rippleSize', 'rippleDurationS', 'gammapower_CA1_z', 'gammapower_CA3_z'};
for i = 1:numel(gfieldnames)
    
    g = gramm('x', AllData.(gfieldnames{i}), 'y', AllData.rippleRatio, 'color', cellstr(AllData.group), 'lightness', cellstr(AllData.timepoint), 'subset', AllData.includeRipple == 1 & AllData.rippleposition >=0 & AllData.rippleposition < 18);
    g.facet_grid([], cellstr(AllData.timepoint));
    g.geom_point();
    %g.stat_glm()
    %g.stat_summary('bin_in',10);
    %g.stat_cornerhist('edges',-4:0.1:10,'aspect',0.5);
    g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
    % g.stat_bin('geom', 'stairs', 'fill', 'transparent');
    g.set_names('x',gfieldnames{i}, 'y',  'ripple PCR', 'color','Group', 'lightness', 'timepoint');
    g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
    figure('Position', [365 532 1040 324])
    g.draw();
end



%% get RZ ripple properties
inclRipplesRZ = ~isnan(AllData.rippleposition) & AllData.rippleposition >=0 & AllData.rippleposition < 18 & strcmp(AllData.timepoint, 'post');
plotData = AllData(inclRipplesRZ,:);
groupcolors = params.colors;

dnames = {'pre', 'post'};
gnames = {'random', 'gamma'};
fnames = {'rippleDurationS', 'rippleSize'};

for f = 1:length(fnames)
    figure
    hold on
    clear vdat
    tmpcolors = [];
    d = 2;
    for g = 1:length(gnames)
        inclDat = strcmp(plotData.group, gnames{g});
        dat = plotData.(fnames{f})(inclDat);
        vdat.([gnames{g}, '_post']) = dat;
        tmpcolors = [tmpcolors; groupcolors.(gnames{g}).(dnames{d})];
    end
    
    violinplot_half(vdat, [], 'ViolinColorMat', tmpcolors, 'BoxWidth', 0.02, 'ShowData', false, 'ViolinAlpha', 0.5);
    ylabel(fnames{f})
    title(fnames{f})
end


%
%
% %%% helpful stuff for plotting
% isGamma = strcmp(metadata.Groups, 'gamma');
% isRandom = strcmp(metadata.Groups, 'random');
% isPre = metadata.FlickerDay < 5;
% isPost = metadata.FlickerDay > 5;
% dayID = 1:height(metadata);
% gID.gamma.pre = dayID(isGamma & isPre & inclDays);
% gID.gamma.post = dayID(isGamma & isPost & inclDays);
% gID.random.pre = dayID(isRandom & isPre & inclDays);
% gID.random.post = dayID(isRandom & isPost & inclDays);
% gnames = {'random', 'gamma'};
% dnames = {'pre', 'post'};

%
% fnames = fieldnames(groupData);
% alldays = 1:length(groupData);
% for g = 1:numel(gnames)
%     for d = 1:numel(dnames)
%         isDay = ismember(alldays, gID.(gnames{g}).(dnames{d}));
%         for f = 1:numel(fnames)
%             tmp = {groupData(isDay).(fnames{f})};
% %             if strcmp(fnames{f}, 'rippleRatio')
% %                 plotData.(gnames{g}).(dnames{d}).(fnames{f}) = cell2mat(tmp);
% %             else
%                 plotData.(gnames{g}).(dnames{d}).(fnames{f}) = cell2mat(tmp');
% %             end
%         end
%     end
% end
%
% %%% plot all ripple positions
% dec_edges = -81:9:99;
% err_edges = -99:18:99;
%
% figure
% hold on
% clear vdat
% tmpcolors = [];
% for d = 2
%     nm = ['rippleDecPos'];
%     nm2 = ['rippleposition'];
%     for g = 1:length(gnames)
%         v_name = [gnames{g}, 'RewardRips'];
%         dec =  plotData.(gnames{g}).(dnames{d}).(nm);
%         pos = plotData.(gnames{g}).(dnames{d}).(nm2);
%         incl = pos > -0 & pos <= 18;
%         posdiff = dec(incl)-pos(incl);
%
%         dec_h = histcounts(dec, dec_edges);
%         dec_nh = dec_h./sum(dec_h);
%
%         err_h = histcounts(posdiff, err_edges);
%         err_nh = err_h./sum(err_h);
%
% %         [density, value] = ksdensity(posdiff);
% %         density = density(value >= min(posdiff) & value <= max(posdiff));
% %         value = value(value >= min(posdiff) & value <= max(posdiff));
% %         value(1) = min(posdiff);
% %         value(end) = max(posdiff);
% %
% %         fill([value value(end:-1:1)], [zeros(size(density)) density(end:-1:1)], groupcolors.(gnames{g}).(dnames{d}), 'LineStyle', 'none');
%         nn = [err_nh zeros(1,length(err_nh))];
%         xx = [err_edges(1:end-1) fliplr(err_edges(1:end-1))];
%         fill(xx, nn, groupcolors.(gnames{g}).(dnames{d}), 'EdgeColor', 'none')
%
%         alpha(0.3)
% %
% %         plotprettypoints_hist(gcf, err_edges, err_nh, posdiff, 'color', colors.(gnames{g})(d,:))
%         datforstats{g} = posdiff;
%
%         clear posdiff density value
%     end
% end
%
% xlabel('decoding error (deg)')
% xlim([-99 99])
% % ylim([0 0.015])
% % yticks([0 0.005 0.01 0.015])
% xticks([-99 0 99])
%
%
% figure
% hold on
% clear vdat
% tmpcolors = [];
% for d = 2
%     nm = ['rippleDecPos'];
%     nm2 = ['rippleposition'];
%
%     for g = 1:length(gnames)
%         v_name = [gnames{g}, 'RewardRips'];
%
%         dec =  plotData.(gnames{g}).(dnames{d}).(nm);
%         pos = plotData.(gnames{g}).(dnames{d}).(nm2);
%
%         incl = pos > 0 & pos <= 18;
%         posdiff = dec(incl)-pos(incl);
%
%         vdat.(v_name) = posdiff;
%         tmpcolors = [tmpcolors; groupcolors.(gnames{g}).(dnames{d})];
%     end
%
% end
%
% text(1.2, -100, {'ranksum p =', num2str(ranksum(vdat.randomRewardRips, vdat.gammaRewardRips))})
%
% % vdat.xlabels = {'post'};
% violinplot_half(vdat, [], 'ViolinColorMat', tmpcolors, 'BoxWidth', 0.02, 'ShowData', true, 'ViolinAlpha', 0.5);
% ylabel('decoding error (deg)')
% title('decoding error, reward zone ripples only')
%
%
% figure
% hold on
% clear vdat
% tmpcolors = [];
% for d = 2
%     nm = ['rippleRatio'];
%     nm2 = ['rippleposition'];
%
%     for g = 1:length(gnames)
%         v_name = [gnames{g}, 'RewardRips'];
%
%         rat =  plotData.(gnames{g}).(dnames{d}).(nm);
%         pos = plotData.(gnames{g}).(dnames{d}).(nm2);
%
%         incl = pos > 0 & pos <= 18;
% %         incl = true(1,length(pos));
%
%         vdat.(v_name) = rat(incl);
%         tmpcolors = [tmpcolors; groupcolors.(gnames{g}).(dnames{d})];
%     end
%
% end
%
% text(1.2, -0.5, {'ranksum p =', num2str(ranksum(vdat.randomRewardRips, vdat.gammaRewardRips))})
%
% % vdat.xlabels = {'post'};
% violinplot_half(vdat, [], 'ViolinColorMat', tmpcolors, 'BoxWidth', 0.02, 'ShowData', true, 'ViolinAlpha', 0.5);
% ylabel('prospective coding ratio')
% title('prospective coding ratio, reward zone ripples only')
%

end

