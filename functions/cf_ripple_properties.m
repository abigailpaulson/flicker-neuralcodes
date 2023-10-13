function cf_ripple_properties(dirs, params, allindex, metadata)
%cf_ripple_properties
%   ripple properties for abby's chronic flicker experiment
%   based on script_SWRprops_VRonly_combineddata
%
%   VR periods only!
% ALP 03/15/2023

params.nonThetaRipples = 1;
savedatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\ripple_properties\'];
%% edges
edges.rippleSize = 3:1:16;
edges.rippleDur = 0:50:500; 
edges.ntDur = 0:5:50;

%% get ripple info per ripple and per file
AllSWRData = table; %initialize empty table
AllNTData = table;
AllFileData = table;
AllDayInfo = table;
for iFile = 1:size(allindex,1)
    %%% index - [1 x 3]
    index = allindex(iFile,:);
    animal = index(1);
    day = index(2);
    [~, iM] = ismember(day, metadata.Date);
    group = metadata.Groups(iM);
    if metadata.FlickerDay(iM) < 5
        timepoint = {'pre'};
    else
        timepoint = {'post'}; 
    end
    
    %%% load best ripple channel for the day
    ripplechandir = ['\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\CA1\'];
    if ~exist([ripplechandir, 'bestChannel_', num2str(index(2)), '.mat'])
        continue
    end
    load([ripplechandir, 'bestChannel_', num2str(index(2)), '.mat'])
    rippleChan = bestRippleChan.all;
    
    %%% set directories
    anperiodsdatadir = [dirs.processeddatadir, params.iden, num2str(index(1)), '_', num2str(index(2)), '\CA1\', num2str(rippleChan), '\'];
    
    %%% get ripple properties
    load([anperiodsdatadir, 'nonthetas', num2str(index(3)), '.mat'])
    load([anperiodsdatadir, 'ripples', num2str(index(3)), '.mat'])
    
    %%% give structures easy names for interfacing
    RipData = ripples{index(1)}{index(2)}{index(3)};
    NTData = nonthetas{index(1)}{index(2)}{index(3)};
    
    %%% if only nontheta ripples...
    if params.nonThetaRipples == 1
        nonthetaper = [NTData.startind NTData.endind];
        isNTRip = isExcluded(RipData.midind, nonthetaper);
        isNTRip = logical(isNTRip);
    else
        isNTRip = true(ones(1,length(RipData.midind)));
    end
    
    %%% only ripples longer than 15 ms
    rippleDurS = (RipData.endind - RipData.startind)./RipData.samprate;
    isLongRip = rippleDurS.*1000 > 15; 
    isInclRip = isNTRip & isLongRip; 
    
    %%% get info about the ripples and nontheta periods overall
    rippleSize = RipData.maxthresh(isInclRip);
    rippleCount = length(RipData.startind(isInclRip));
    rippleDurS = (RipData.endind(isInclRip) - RipData.startind(isInclRip))./RipData.samprate;
    nonThetaDurS = diff(nonthetaper,1,2)./NTData.samprate;
    nonThetaCount = size(nonthetaper,1);
    rippleRateAll = rippleCount/sum(nonThetaDurS);
    
    %%% rate per nontheta period
    perDurByRip = []; 
    nPeriodRip = zeros(size(nonthetaper,1),1);
    for iPeriod = 1:size(nonthetaper,1)
        isPeriodRip = isExcluded(RipData.midind, nonthetaper(iPeriod,:));
        nPeriodRip(iPeriod) = sum(isPeriodRip);
        tempPerDurByRip = nonThetaDurS(iPeriod).*ones(sum(isPeriodRip),1);
        perDurByRip = [perDurByRip; tempPerDurByRip];
        clear tempPerDurByRip
    end
    rippleRatePeriod = nPeriodRip./nonThetaDurS;
    
    %%% rate per by SWR size
    rippleRateSize = histcounts(rippleSize, edges.rippleSize)./sum(nonThetaDurS); 
    rippleSizeProp = histcounts(rippleSize, edges.rippleSize)./rippleCount; 
    
    totnonThetaDurS = sum(nonThetaDurS, 'omitnan');
    
    key1 = 1; %temp for joining
    dayInfo = table(animal, day, timepoint, group, key1);
    
    key1 = ones(length(rippleSize),1); %temp for joining
    PerSWRData = table(rippleSize, rippleDurS, perDurByRip, key1);
    PerSWRData = join(PerSWRData, dayInfo);
    PerSWRData = removevars(PerSWRData, {'key1'});
    
    key1 = ones(length(nPeriodRip),1); %temp for joining
    PerNTData = table(nPeriodRip, nonThetaDurS, rippleRatePeriod, key1);
    PerNTData = join(PerNTData, dayInfo);
    PerNTData = removevars(PerNTData, {'key1'});
    
    key1 = 1; %temp for joining
    PerFileData = table(rippleCount, totnonThetaDurS, key1);
    PerFileData = join(PerFileData, dayInfo);
    PerFileData = removevars(PerFileData, {'key1'});
    
    AllSWRData = [AllSWRData; PerSWRData];
    AllNTData = [AllNTData; PerNTData];
    AllFileData = [AllFileData; PerFileData];
    AllDayInfo = [AllDayInfo; dayInfo];
    
clear ripples RipData nonthetas NTData PerSWRData PerNTData dayInfo
clear -regexp ^ripple ^nonTheta
end

%%% get group summary by day...
AllDayInfo = unique(AllDayInfo, 'rows');
tmpTable1 = groupsummary(AllFileData(:,[1:2,4]), 'day', 'sum');
tmpTable2 = groupsummary(AllSWRData(:,[1:2,5]), 'day', 'mean');
tmpTable = outerjoin(tmpTable1(:,[1,3:end]), tmpTable2(:,[1,3:end]), 'MergeKeys', true);
AllDayData = join(tmpTable, AllDayInfo);

%% save data structures
filename = 'rippleproperties.mat';
save([savedatadir, filename], 'AllSWRData', 'AllNTData', 'AllFileData', 'AllDayData', '-v7.3')

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';

filename = 'TableData_RippleProperties_SWR.txt';
writetable(AllSWRData, fullfile(statsdir, filename));   

filename = 'TableData_RippleProperties_NT.txt';
writetable(AllNTData, fullfile(statsdir, filename));  

filename = 'TableData_RippleProperties_Day.txt';
writetable(AllDayData, fullfile(statsdir, filename));  



%% some plotting to check
savefigdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\ripple_properties\';
groupcolors = cbrewer('qual', 'Paired', 6); 
gcolors = groupcolors([2,1,4,3],:); %rearranging it cuz of the way that the gramm does the plotting
trialscolors = groupcolors(end-1:end,:);


%% per animal plots
%%% ripple duration - all (though exclude outliers > 500ms)
g(1,1) = gramm('x', AllSWRData.animal, 'y', AllSWRData.rippleDurS, 'color', cellstr(AllSWRData.group), ...
    'lightness', cellstr(AllSWRData.timepoint), 'subset', AllSWRData.rippleDurS < 0.5);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('half', 'true', 'fill', 'transparent');
g(1,1).facet_grid([], cellstr(AllSWRData.group));
g(1,1).set_names('x', 'flicker day', 'y', 'ripple duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).facet_grid([], cellstr(AllSWRData.group));
g(2,1).stat_boxplot('width',5,'dodge',0.4);
g(2,1).set_names('x', 'flicker day', 'y', 'ripple duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple duration - <0.5s');
figure('Position', [84 317 1218 610])
g.draw();
savefigALP(savefigdir, 'rippleduration', 'pdf')

%%% ripple duration - NT > 20s
g(1,1) = gramm('x', AllSWRData.animal, 'y', AllSWRData.rippleDurS, 'color', cellstr(AllSWRData.group), ...
    'lightness', cellstr(AllSWRData.timepoint), 'subset', AllSWRData.rippleDurS < 0.5 & AllSWRData.perDurByRip > 20);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('half', 'true', 'fill', 'transparent');
g(1,1).facet_grid([], cellstr(AllSWRData.group));
g(1,1).set_names('x', 'flicker day', 'y', 'ripple duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).facet_grid([], cellstr(AllSWRData.group));
g(2,1).stat_boxplot('width',5,'dodge',0.4);
g(2,1).set_names('x', 'flicker day', 'y', 'ripple duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple duration - <0.5s AND NT>20s');
figure('Position', [84 317 1218 610])
g.draw();
savefigALP(savefigdir, 'rippleduration_longNT', 'pdf')


%%% ripple rate
g(1,1) = gramm('x', AllNTData.animal, 'y', AllNTData.rippleRatePeriod, 'color', cellstr(AllNTData.group), ...
    'lightness', cellstr(AllNTData.timepoint));
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('half', 'true', 'fill', 'transparent');
g(1,1).facet_grid([], cellstr(AllNTData.group));
g(1,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).facet_grid([], cellstr(AllNTData.group));
g(2,1).stat_boxplot('width',3,'dodge',1);
g(2,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple rate');
figure('Position', [84 317 1218 610])
g.draw();
savefigALP(savefigdir, 'ripplerate', 'pdf')


%%% ripple rate - more than 1 ripple per period
g(1,1) = gramm('x', AllNTData.animal, 'y', AllNTData.rippleRatePeriod, 'color', cellstr(AllNTData.group), ...
    'lightness', cellstr(AllNTData.timepoint), 'subset', AllNTData.nPeriodRip > 1);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('half', 'true', 'fill', 'transparent');
g(1,1).facet_grid([], cellstr(AllNTData.group));
g(1,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).facet_grid([], cellstr(AllNTData.group));
g(2,1).stat_boxplot('width',3,'dodge',1);
g(2,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple rate - at least 1 ripple per period');
figure('Position', [84 317 1218 610])
g.draw();
savefigALP(savefigdir, 'ripplerate_1plusripple', 'pdf')


%%% ripple rate - long NT
g(1,1) = gramm('x', AllNTData.animal, 'y', AllNTData.rippleRatePeriod, 'color', cellstr(AllNTData.group), ...
    'lightness', cellstr(AllNTData.timepoint), 'subset', AllNTData.nonThetaDurS > 5);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('half', 'true', 'fill', 'transparent');
g(1,1).facet_grid([], cellstr(AllNTData.group));
g(1,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).facet_grid([], cellstr(AllNTData.group));
g(2,1).stat_boxplot('width',3,'dodge',1);
g(2,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple rate - nonThetaDur > 5s');
figure('Position', [84 317 1218 610])
g.draw();
savefigALP(savefigdir, 'ripplerate_5plusNT', 'pdf')


%%% ripple size - all (though exclude outliers > 500ms)
g(1,1) = gramm('x', AllSWRData.animal, 'y', AllSWRData.rippleSize, 'color', cellstr(AllSWRData.group), ...
    'lightness', cellstr(AllSWRData.timepoint), 'subset', AllSWRData.rippleDurS < 0.5);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('half', 'true', 'fill', 'transparent');
g(1,1).facet_grid([], cellstr(AllSWRData.group));
g(1,1).set_names('x', 'flicker day', 'y', 'ripple size(sd)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).facet_grid([], cellstr(AllSWRData.group));
g(2,1).stat_boxplot('width',5,'dodge',0.4);
g(2,1).set_names('x', 'flicker day', 'y', 'ripple size (sd)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple size - duration<0.5s');
figure('Position', [84 317 1218 610])
g.draw();
savefigALP(savefigdir, 'ripplesize', 'pdf')


%%% ripple size - all (though exclude outliers > 500ms)
g(1,1) = gramm('x', AllSWRData.animal, 'y', AllSWRData.rippleSize, 'color', cellstr(AllSWRData.group), ...
    'lightness', cellstr(AllSWRData.timepoint), 'subset', AllSWRData.rippleDurS < 0.5 & AllSWRData.perDurByRip > 10);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('half', 'true', 'fill', 'transparent');
g(1,1).facet_grid([], cellstr(AllSWRData.group));
g(1,1).set_names('x', 'flicker day', 'y', 'ripple size(sd)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).facet_grid([], cellstr(AllSWRData.group));
g(2,1).stat_boxplot('width',5,'dodge',0.4);
g(2,1).set_names('x', 'flicker day', 'y', 'ripple size (sd)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple size - duration<0.5s & NT dur > 10s');
figure('Position', [84 317 1218 610])
g.draw();
savefigALP(savefigdir, 'ripplesize_10sNT', 'pdf')


%%% non theta duration 
g(1,1) = gramm('x', AllNTData.animal, 'y', AllNTData.nonThetaDurS, 'color', cellstr(AllNTData.group), ...
    'lightness', cellstr(AllNTData.timepoint));
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('half', 'true', 'fill', 'transparent');
g(1,1).facet_grid([], cellstr(AllNTData.group));
g(1,1).set_names('x', 'flicker day', 'y', 'NT period duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).facet_grid([], cellstr(AllNTData.group));
g(2,1).stat_boxplot('width',3,'dodge',1);
g(2,1).set_names('x', 'flicker day', 'y', 'NT period duration (s))', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('nontheta duration');
figure('Position', [84 317 1218 610])
g.draw();
savefigALP(savefigdir, 'nonthetaduration', 'pdf')



%% per group plots 
g = gramm('x', cellstr(AllNTData.timepoint), 'y', AllNTData.rippleRatePeriod, ...
    'color', cellstr(AllNTData.group), 'lightness', cellstr(AllNTData.timepoint));
g(1,1).stat_violin('normalization', 'count', 'fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple rate per period');
figure
g.draw();
savefigALP(savefigdir, 'ripplerate_group', 'pdf')


g = gramm('x', cellstr(AllNTData.timepoint), 'y', AllNTData.rippleRatePeriod, ...
    'color', cellstr(AllNTData.group), 'lightness', cellstr(AllNTData.timepoint), 'subset', AllNTData.nPeriodRip > 1);
g(1,1).stat_violin('normalization', 'count', 'fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple rate - periods w 1+ ripple');
figure
g.draw();
savefigALP(savefigdir, 'ripplerate_1plusripple_group', 'pdf')


g = gramm('x', cellstr(AllSWRData.timepoint), 'y', AllSWRData.rippleDurS, 'color', cellstr(AllSWRData.group), ...
    'lightness', cellstr(AllSWRData.timepoint), 'subset', AllSWRData.rippleDurS < 0.5 & AllSWRData.perDurByRip > 20);
g.stat_violin('normalization', 'count', 'fill', 'transparent');
g.stat_boxplot('width', 0.15);
g.set_names('x', 'flicker day', 'y', 'ripple duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple duration - <0.5s, NT > 20s');
figure
g.draw();
savefigALP(savefigdir, 'rippleduration_longNT_group', 'pdf')


g = gramm('x', cellstr(AllSWRData.timepoint), 'y', AllSWRData.rippleDurS, 'color', cellstr(AllSWRData.group), ...
    'lightness', cellstr(AllSWRData.timepoint), 'subset', AllSWRData.rippleDurS < 0.5);
g.stat_violin('normalization', 'count', 'fill', 'transparent');
g.stat_boxplot('width', 0.15);
g.set_names('x', 'flicker day', 'y', 'ripple duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple duration - <0.5s');
figure
g.draw();
savefigALP(savefigdir, 'rippleduration_group', 'pdf')


%%% ripple size 
g = gramm('x', cellstr(AllSWRData.timepoint), 'y', AllSWRData.rippleSize, 'color', cellstr(AllSWRData.group), ...
    'lightness', cellstr(AllSWRData.timepoint), 'subset', AllSWRData.rippleDurS < 0.5);
g.stat_violin('normalization', 'count', 'fill', 'transparent');
g.stat_boxplot('width', 0.15);
g.set_names('x', 'flicker day', 'y', 'ripple size (sd)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple size - <0.5s');
figure
g.draw();
savefigALP(savefigdir, 'ripplesize_group', 'pdf')


%%% NonThetaDuration
g = gramm('x', cellstr(AllNTData.timepoint), 'y', AllNTData.nonThetaDurS, ...
    'color', cellstr(AllNTData.group), 'lightness', cellstr(AllNTData.timepoint));
g(1,1).stat_violin('fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'nontheta duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple rate - periods w 1+ ripple');
figure
g.draw();
savefigALP(savefigdir, 'nothetaduration_group', 'pdf')


%% per day plots 

%%% ripple rate 
g(1,1) = gramm('x', cellstr(AllDayData.timepoint), 'y', AllDayData.sum_rippleCount./AllDayData.sum_totnonThetaDurS, ...
    'color', cellstr(AllDayData.group), 'lightness', cellstr(AllDayData.timepoint), 'subset', AllDayData.sum_rippleCount > 10);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('normalization', 'count', 'fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).stat_summary('geom',{'bar' 'black_errorbar'},'setylim',true, 'type', 'sem');
g(2,1).set_names('x', 'flicker day', 'y', 'ripple rate (hz)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple rate - days > 10 ripples');
figure
g.draw();
savefigALP(savefigdir, 'ripplerate_dayavg', 'pdf')

%%% ripple size
g(1,1) = gramm('x', cellstr(AllDayData.timepoint), 'y', AllDayData.mean_rippleSize, ...
    'color', cellstr(AllDayData.group), 'lightness', cellstr(AllDayData.timepoint), 'subset', AllDayData.sum_rippleCount > 10);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('normalization', 'count', 'fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'ripple size (sd)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).stat_summary('geom',{'bar' 'black_errorbar'},'setylim',true, 'type', 'sem');
g(2,1).set_names('x', 'flicker day', 'y', 'ripple size(sd)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple size - days > 10 ripples');
figure
g.draw();
savefigALP(savefigdir, 'ripplesize_dayavg', 'pdf')

%%% ripple duration
g(1,1) = gramm('x', cellstr(AllDayData.timepoint), 'y', AllDayData.mean_rippleDurS, ...
    'color', cellstr(AllDayData.group), 'lightness', cellstr(AllDayData.timepoint), 'subset', AllDayData.sum_rippleCount > 10);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('normalization', 'count', 'fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'ripple duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).stat_summary('geom',{'bar' 'black_errorbar'},'setylim',true, 'type', 'sem');
g(2,1).set_names('x', 'flicker day', 'y', 'ripple duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('ripple duration - days > 10ripples');
figure
g.draw();
savefigALP(savefigdir, 'rippleduration_dayavg', 'pdf')

%%% ripple duration
g(1,1) = gramm('x', cellstr(AllDayData.timepoint), 'y', AllDayData.sum_totnonThetaDurS, ...
    'color', cellstr(AllDayData.group), 'lightness', cellstr(AllDayData.timepoint), 'subset', AllDayData.sum_rippleCount > 10);
g(2,1) = copy(g(1,1));
g(1,1).stat_violin('normalization', 'count', 'fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'nontheta duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).stat_summary('geom',{'bar' 'black_errorbar'},'setylim',true, 'type', 'sem');
g(2,1).set_names('x', 'flicker day', 'y', 'non theta duration (s)', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('total nonthetaduration');
figure
g.draw();
savefigALP(savefigdir, 'nonthetaduration_dayavg', 'pdf')

end

