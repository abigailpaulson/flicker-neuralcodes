function cf_decoding_ripples(dirs, params, allindex, metadata, posType, cellTypes)
%cf_decoding_ripples
%   get bayesian decoding during ripples
%ALP 1/9/2023

rewrite.decodingFiles = 1; 

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

%%% loop over all days
for d = 1:size(dayindex,1)
    dayData = [];
    disp(['getting ripple analyses for ', num2str(dayindex(d,1)), ' ', num2str(dayindex(d,2))])
    files = allindex(allindex(:,2) == dayindex(d,2),3); %should be only VR files pre flicker or all VR on day 9/10
    sessindex = dayindex(d,:);
    index = repmat(sessindex, length(files), 1);
    index = [index files];
    anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
    
    dayfilename = ['rippledecoding_halfratiozones_', dType, '_', cellTypes, '_' num2str(dayindex(d,2)), '.mat'];
    if exist([savedir, dayfilename], 'file') && ~rewrite.decodingFiles 
        load([savedir, dayfilename])
    else
        
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
            
            %%% get spikes only of included cells
            tmpSpikeIDs = spikevect(index(i,3)).spikeIDs;
            inclSpikes = ismember(tmpSpikeIDs, inclCells);
            tmpSpikeIDs = tmpSpikeIDs(inclSpikes);
            tmpSpikeTimes = spikevect(index(i,3)).spikeTimes(inclSpikes);
            
            tmpSpikeStruct.spiketimes = tmpSpikeTimes;
            tmpSpikeStruct.spikeIDs = tmpSpikeIDs;
            
            
            decoding_bins = cf_rippledecoding_timebins(index(i,:), trainingData, tmpBehavior.time, tmpBehavior.(posName), tmpSpikeStruct, rippleWindows, rippleMids, params, inclCells, nCells);
            decoding_window = cf_rippledecoding_window(index(i,:), trainingData, tmpBehavior.time, tmpBehavior.(posName), tmpSpikeStruct, rippleWindows, rippleMids, params, inclCells, nCells);
            
            dayData(i).rippleposition = decoding_bins.ripplePos;
            dayData(i).rippleRatio = decoding_bins.rippleRatio;
            dayData(i).rippleDecPos = decoding_window.rippleDecPos;
            dayData(i).rippleRelPos = decoding_window.rippleRelPos;
            dayData(i).decoding_bins = decoding_bins;
            dayData(i).decoding_window = decoding_window;
            
            clear tmp* decoding_bins decoding_window rippleMids rippleWindows
            close all
        end
        
        save([savedir, 'rippledecoding_halfratiozones_', dType, '_', cellTypes, '_' num2str(dayindex(d,2)), '.mat'], 'dayData')
    end
    if ~isempty(dayData)
        groupData(d).rippleposition = cell2mat({dayData.rippleposition}');
        groupData(d).rippleRatio = cell2mat({dayData.rippleRatio}');
        groupData(d).rippleDecPos = cell2mat({dayData.rippleDecPos}');
        groupData(d).rippleRelPos = cell2mat({dayData.rippleRelPos}');
    end
    
    clear dayData
end

%%% save group data
%load curr pos data to exclude bad decoding days
if strcmp(cellTypes, 'all')
currPos = load('\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_currentposition\Group_decoding_currpos_180_all_test20_210527.mat');
else
    currPos = load('\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_currentposition\Group_decoding_currpos_180_PC_test20.mat');

end
inclDays = [currPos.groupData.include]';

%%% helpful stuff for plotting
isGamma = strcmp(metadata.Groups, 'gamma'); 
isRandom = strcmp(metadata.Groups, 'random');
isPre = metadata.FlickerDay < 5;
isPost = metadata.FlickerDay > 5; 
dayID = 1:height(metadata);
gID.gamma.pre = dayID(isGamma & isPre & inclDays);
gID.gamma.post = dayID(isGamma & isPost & inclDays);
gID.random.pre = dayID(isRandom & isPre & inclDays);
gID.random.post = dayID(isRandom & isPost & inclDays);
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'}; 
groupcolors = params.colors;

fnames = fieldnames(groupData);
alldays = 1:length(groupData);
for g = 1:numel(gnames)
    for d = 1:numel(dnames)
        isDay = ismember(alldays, gID.(gnames{g}).(dnames{d}));
        for f = 1:numel(fnames)
            tmp = {groupData(isDay).(fnames{f})};
%             if strcmp(fnames{f}, 'rippleRatio')
%                 plotData.(gnames{g}).(dnames{d}).(fnames{f}) = cell2mat(tmp);
%             else
                plotData.(gnames{g}).(dnames{d}).(fnames{f}) = cell2mat(tmp');
%             end
        end
    end
end

%%% plot all ripple positions
dec_edges = -81:9:99; 
err_edges = -99:3:99; 

figure
hold on
clear vdat
tmpcolors = [];
tmpH = []; tmpDensity = [];
for d = 2
    nm = ['rippleDecPos'];
    nm2 = ['rippleposition'];
    for g = 1:length(gnames)
        v_name = [gnames{g}, 'RewardRips'];
        dec =  plotData.(gnames{g}).(dnames{d}).(nm);
        pos = plotData.(gnames{g}).(dnames{d}).(nm2); 
        incl = pos > 0 & pos <= 18;
        posdiff = dec(incl)-pos(incl);
        
        dec_h = histcounts(dec, dec_edges);
        dec_nh = dec_h./sum(dec_h);
        
        err_h = histcounts(posdiff, err_edges);
        err_nh = err_h./sum(err_h);
        
        [density, value] = ksdensity(posdiff); 
        density = density(value >= min(posdiff) & value <= max(posdiff));
        value = value(value >= min(posdiff) & value <= max(posdiff));
        value(1) = min(posdiff);
        value(end) = max(posdiff);
% 
%         fill([value value(end:-1:1)], [zeros(size(density)) density(end:-1:1)], groupcolors.(gnames{g}).(dnames{d}), 'LineStyle', 'none');  
        nn = [err_nh zeros(1,length(err_nh))];
        xx = [err_edges(1:end-1) fliplr(err_edges(1:end-1))];
        fill(xx, nn, groupcolors.(gnames{g}).(dnames{d}), 'EdgeColor', 'none')
        
        alpha(0.3)
%         
%         plotprettypoints_hist(gcf, err_edges, err_nh, posdiff, 'color', colors.(gnames{g})(d,:))
        datforstats{g} = posdiff; 
        tmpH = [tmpH; err_nh];
        tmpDensity{g} = density;
        densityValues{g} = value;
        
        clear posdiff density value 
    end
end

xlabel('decoding error (deg)')
xlim([-99 99])
% ylim([0 0.015])
% yticks([0 0.005 0.01 0.015])
xticks([-99 0 99])

figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Figure4\';

figure('Position', [680 820 521 158])
hold on
subplot(2,1,1)
imagesc(err_edges(1:end-1),ones(1,length(tmpH(2,:))), tmpH(2,:), [0 0.08])
colormap(flipud(cmocean('gray')));
xlim([-99 99])
colorbar
xticks([-99 0 99])
yticks([])
title(gnames{2})
subplot(2,1,2)
hold on
imagesc(err_edges(1:end-1),ones(1,length(tmpH(1,:))), tmpH(1,:), [0 0.08])
colorbar
xlim([-99 99])
colorbar
xticks([-99 0 99])
yticks([])
title(gnames{1})
makefigurepretty(gcf)
savefigALP(figdir, 'rippleRelPos_RZ_example1')

figure('Position', [680 820 521 158])
hold on
subplot(2,1,1)
imagesc(densityValues{2},ones(1,length(tmpDensity{2})), tmpDensity{2}, [0.002 0.011])
colormap(flipud(cmocean('gray')));
xlim([-99 99])
colorbar
xticks([-99 0 99])
yticks([])
title(gnames{2})
subplot(2,1,2)
hold on
imagesc(densityValues{1},ones(1,length(tmpDensity{1})), tmpDensity{1}, [0.002 0.011])
colorbar
xlim([-99 99])
box off
xticks([-99 0 99])
yticks([])
title(gnames{1})
makefigurepretty(gcf)
savefigALP(figdir, 'rippleRelPos_RZ_example2')



figure
hold on
clear vdat
tmpcolors = [];
for d = 2
    nm = ['rippleDecPos'];
    nm2 = ['rippleposition'];
    
    for g = 1:length(gnames)
        v_name = [gnames{g}, 'RewardRips'];
        
        dec =  plotData.(gnames{g}).(dnames{d}).(nm);
        pos = plotData.(gnames{g}).(dnames{d}).(nm2); 
        
        incl = pos > 0 & pos <= 18;
        posdiff = dec(incl)-pos(incl);
        
        vdat.(v_name) = posdiff;
        tmpcolors = [tmpcolors; groupcolors.(gnames{g}).(dnames{d})];
    end
    
end

text(1.2, -100, {'ranksum p =', num2str(ranksum(vdat.randomRewardRips, vdat.gammaRewardRips))})

% vdat.xlabels = {'post'};
violinplot_half(vdat, [], 'ViolinColorMat', tmpcolors, 'BoxWidth', 0.02, 'ShowData', true, 'ViolinAlpha', 0.5);
ylabel('decoding error (deg)')
title('decoding error, reward zone ripples only')


figure
hold on
clear vdat
tmpcolors = [];
for d = 1
    nm = ['rippleRatio'];
    nm2 = ['rippleposition'];
    
    for g = 1:length(gnames)
        v_name = [gnames{g}, 'RewardRips'];
        
        rat =  plotData.(gnames{g}).(dnames{d}).(nm);
        pos = plotData.(gnames{g}).(dnames{d}).(nm2); 
        
        incl = pos > 0 & pos <= 18;
%         incl = true(1,length(pos));
        
        vdat.(v_name) = rat(incl);
        tmpcolors = [tmpcolors; groupcolors.(gnames{g}).(dnames{d})];
    end
    
end

text(1.2, -0.5, {'ranksum p =', num2str(ranksum(vdat.randomRewardRips, vdat.gammaRewardRips))})

% vdat.xlabels = {'post'};
violinplot_half(vdat, [], 'ViolinColorMat', tmpcolors, 'BoxWidth', 0.02, 'ShowData', true, 'ViolinAlpha', 0.5);
ylabel('prospective coding ratio')
title('prospective coding ratio, reward zone ripples only')


end

