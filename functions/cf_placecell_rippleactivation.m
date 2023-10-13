function cf_placecell_rippleactivation(dirs, params, allindex, metadata)
%cf_placecell_rippleactivation
%   activation/coactivation during sharp-wave ripples
%ALP 3/19/2023

rewrite.files = 1;
%% set stuff up as needed
dayindex = unique(allindex(:,1:2), 'rows');
params.nonThetaRipples = 1;
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';

savedatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\ripples_activation_coactivation\';
savefilename = 'ripple_activation_coactivation_placecells.mat';

if ~exist([savedatadir, savefilename]) || rewrite.files
    %% loop over days
    ActivityData = table;
    CoactivityData = table;
    DayInfo = table;
    for d = 1:size(dayindex,1)
        disp(['getting activity/coactivity for ', num2str(dayindex(d,2))])
        files = allindex(allindex(:,2) == dayindex(d,2),3);
        
        %%% load best ripple channel
        ripplechandir = ['\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\CA1\'];
        if ~exist([ripplechandir, 'bestChannel_', num2str(dayindex(d,2)), '.mat'])
            continue
        end
        load([ripplechandir, 'bestChannel_', num2str(dayindex(d,2)), '.mat'])
        rippleChan = bestRippleChan.all;
        
        %%% set directories
        anprocesseddatadir = [dirs.processeddatadir, params.iden, num2str(dayindex(d,1)), ...
            '_', num2str(dayindex(d,2)), '\'];
        anperiodsdatadir = [dirs.processeddatadir, params.iden, num2str(dayindex(d,1)), ...
            '_', num2str(dayindex(d,2)), '\CA1\', num2str(rippleChan), '\'];
        
        %%% load decoding info to get the ripple position
        decdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\ripples_activation_coactivation\'];
        load(decdir, ['rippledecoding_halratiozones_180_all_', num2str(dayindex(d,2)), '.mat']);
        
        nActive = [];
        nCoactive = [];
        nRipples = [];
        ActivityFile = [];
        for f = 1:length(files)
            index = [dayindex(d,:) files(f)];
            
            %%% get ripple properties
            load([anperiodsdatadir, 'nonthetas', num2str(index(3)), '.mat'])
            load([anperiodsdatadir, 'ripples', num2str(index(3)), '.mat'])
            
            %%% give structures easy names for interfacing
            SWRData = ripples{index(1)}{index(2)}{index(3)};
            NTData = nonthetas{index(1)}{index(2)}{index(3)};
            
            %%% only nontheta ripples and ripples longer than 15 ms
            if params.nonThetaRipples == 1
                nonthetaper = [NTData.startind NTData.endind];
                isNTRip = isExcluded(SWRData.midind, nonthetaper);
                isNTRip = logical(isNTRip);
            else
                isNTRip = true(ones(1,length(SWRData.midind)));
            end
            rippleDurS = (SWRData.endind - SWRData.startind)./SWRData.samprate; %in s
            isLongRip = rippleDurS.*1000 > 15; %ms
            isInclRip = isNTRip & isLongRip;
            rippleTimes = [SWRData.starttime(isInclRip) SWRData.endtime(isInclRip)]; %in s
            
            %%% get ripple position
            ripplePosition = dayData(files(f)).rippleposition;
            
            %%% load spike information
            anspikedir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
            [spikes, unitInfo, spikevect] = cf_getspikes(anspikedir, dirs, params, index, params.brainReg);
            
            if isempty(spikes)
                continue
            end
            
            %%% what kind of cells do we care about? probs place cells
            %load place field information
            load([ratemapdir, 'spatialmaps_all_500_full_', num2str(dayindex(d,2)), '.mat'])
            isPC = ratemaps.includeCell;
            if ~any(isPC)
                continue
            end
            PCid = ratemaps.cellIDs(isPC);
            nCells = length(PCid);
            PCspikes = spikes(files(f)).data(isPC);
            
            %%% get brain regions of the different cells
            regions = ratemaps.brainReg;
            PCregions = regions(isPC);
            
            %%% get activation and coactivation information
            ActivityFile = cf_getactivationprob(PCid, PCspikes, rippleTimes, [], PCregions);
            
            %%% make vectors for each file
            nActive = [nActive ActivityFile.nActive];
            nCoactive = [nCoactive ActivityFile.nCoactive];
            nRipples = [nRipples ActivityFile.nAllRipples];
            
            clear ratemaps SWRData NTData
        end
        
        %%% set up info variables for the table
        animal = dayindex(d,1);
        day = dayindex(d,2);
        group = metadata.Groups(d);
        if metadata.FlickerDay(d) < 5
            timepoint = {'pre'};
        else
            timepoint = {'post'};
        end
        tmpKey = 1;
        nDayRipples = sum(nRipples, 'omitnan');
        
        tmpDayInfo = table(animal, day, group, timepoint, nDayRipples, tmpKey);
        
        %%% activity table
        fracActive = sum(nActive,2,'omitnan')./nDayRipples;
        activeRegions = ActivityFile(1).activeRegions';
        tmpKey = ones(size(nActive,1),1);
        
        tmpActivityData = table(fracActive, activeRegions, tmpKey);
        tmpActivityData = outerjoin(tmpActivityData, tmpDayInfo, 'mergeKeys', true);
        tmpActivityData = removevars(tmpActivityData, 'tmpKey');
        
        %%% coactivity table
        fracCoactive = sum(nCoactive,2,'omitnan')./sum(nRipples, 'omitnan');
        coactiveRegions = ActivityFile(1).coactiveRegions';
        tmpKey = ones(size(nCoactive,1),1);
        
        tmpCoactivityData = table(fracCoactive, coactiveRegions, tmpKey);
        tmpCoactivityData = outerjoin(tmpCoactivityData, tmpDayInfo, 'mergeKeys', true);
        tmpCoactivityData = removevars(tmpCoactivityData, 'tmpKey');
        
        %%% append to full tables
        ActivityData = [ActivityData; tmpActivityData];
        CoactivityData = [CoactivityData; tmpCoactivityData];
        DayInfo = [DayInfo; tmpDayInfo];
        
        
        clear tmp*
    end
    
    %% save the table!
    disp('saving...')
    save([savedatadir, savefilename], 'ActivityData', 'CoactivityData');
    
    statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
    statsfilename = 'TableData_ripples_activation_placecells.txt';
    writetable(ActivityData, fullfile(statsdir, statsfilename));
    
    statsfilename = 'TableData_ripples_coactivation_placecells.txt';
    writetable(CoactivityData, fullfile(statsdir, statsfilename));
else
    load([savedatadir, savefilename])
end

%% get day summary 
%get day averages for active and coactive pairs 
dayActivityData = groupsummary(ActivityData(:,[1:4,7]), {'day', 'activeRegions'}, 'mean'); 
dayActivityData = outerjoin(dayActivityData, DayInfo, 'mergeKeys', true);

dayCoactivityData = groupsummary(CoactivityData(:,[1:4,7]), {'day', 'coactiveRegions'}, 'mean'); 
dayCoactivityData = outerjoin(dayCoactivityData, DayInfo, 'mergeKeys', true);


%% plot!
groupcolors = cbrewer('qual', 'Paired', 6);
gcolors = groupcolors([2,1,4,3],:); %rearranging it cuz of the way that the gramm does the plotting
trialscolors = groupcolors(end-1:end,:);
savefigdir = savedatadir;

%%% ripple duration
g(1,1) = gramm('x', cellstr(ActivityData.timepoint), 'y', ActivityData.fracActive, ...
    'color', cellstr(ActivityData.group), 'lightness', cellstr(ActivityData.timepoint), ...
    'subset', ActivityData.nDayRipples > 10);
g(1,1).stat_violin('fill', 'transparent');
%g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'activation probability', ...
    'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('Activity - all');
figure
g.draw();
savefigALP(savefigdir, 'activationprobability')

clear g
g(1,1) = gramm('x', ActivityData.animal, 'y', ActivityData.fracActive, ...
    'color', cellstr(ActivityData.group), 'lightness', cellstr(ActivityData.timepoint), ...
    'subset', ActivityData.nDayRipples > 10);
%g(1,1).stat_violin('fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'activation probability', ...
    'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('Activity - by animal');
figure
g.draw();
savefigALP(savefigdir, 'activationprobability_peranimal')

clear g
g(1,1) = gramm('x', cellstr(ActivityData.timepoint), 'y', ActivityData.fracActive, ...
    'color', cellstr(ActivityData.group), 'lightness', cellstr(ActivityData.timepoint), ...
    'subset', ActivityData.nDayRipples > 10);
g(1,1).facet_grid([], cellstr(ActivityData.activeRegions));
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'activation probability', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('Activity - by region');
figure
g.draw();
savefigALP(savefigdir, 'activationprobability_byregion')

clear g
g(1,1) = gramm('x', ActivityData.animal, 'y', ActivityData.fracActive, ...
    'color', cellstr(ActivityData.group), 'lightness', cellstr(ActivityData.timepoint), ...
    'subset', ActivityData.nDayRipples > 10);
g(1,1).facet_grid(cellstr(ActivityData.activeRegions), []);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'activation probability', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('Activity - by region and animal');
figure
g.draw();
savefigALP(savefigdir, 'activationprobability_byregion_peranimal')

g(1,1) = gramm('x', cellstr(CoactivityData.timepoint), 'y', CoactivityData.fracCoactive, ...
    'color', cellstr(CoactivityData.group), 'lightness', cellstr(CoactivityData.timepoint), ...
    'subset', CoactivityData.nDayRipples > 10);
%g(1,1).stat_violin('fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'coactivation probability', ...
    'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('Coactivity - all');
figure
g.draw();
savefigALP(savefigdir, 'coactivationprobability')

g(1,1) = gramm('x', CoactivityData.animal, 'y', CoactivityData.fracCoactive, ...
    'color', cellstr(CoactivityData.group), 'lightness', cellstr(CoactivityData.timepoint), ...
    'subset', CoactivityData.nDayRipples > 10);
%g(1,1).stat_violin('fill', 'transparent');
g(1,1).stat_boxplot('width', 0.15);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'coactivation probability', ...
    'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('Coactivity - all');
figure
g.draw();
savefigALP(savefigdir, 'coactivationprobability_peranimal')

g(1,1) = gramm('x', cellstr(CoactivityData.timepoint), 'y', CoactivityData.fracCoactive, ...
    'color', cellstr(CoactivityData.group), 'lightness', cellstr(CoactivityData.timepoint), ...
    'subset', CoactivityData.nDayRipples > 10);
g(1,1).facet_grid([], cellstr(CoactivityData.coactiveRegions));
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'coactivation probability', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('Coactivity - all');
figure
g.draw();
savefigALP(savefigdir, 'coactivationprobability_byregion')

g(1,1) = gramm('x', CoactivityData.animal, 'y', CoactivityData.fracCoactive, ...
    'color', cellstr(CoactivityData.group), 'lightness', cellstr(CoactivityData.timepoint), ...
    'subset', CoactivityData.nDayRipples > 10);
g(1,1).facet_grid(cellstr(CoactivityData.coactiveRegions), []);
g(1,1).stat_boxplot();
g(1,1).set_names('x', 'flicker day', 'y', 'coactivation probability', 'color', 'Group', 'lightness', 'timepoint');
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
g.set_title('Coctivity - all');
figure
g.draw();
savefigALP(savefigdir, 'coactivationprobability_byregion_peranimal')

clear g
g(1,1) = gramm('x', cellstr(dayCoactivityData.timepoint), 'y', dayCoactivityData.mean_fracCoactive, ...
    'color', cellstr(dayCoactivityData.group), 'lightness', cellstr(dayCoactivityData.timepoint), ...
    'subset', dayCoactivityData.mean_nDayRipples > 10);
g(1,1).facet_grid([], cellstr(dayCoactivityData.coactiveRegions));
g(2,1) = copy(g(1,1));
g(1,1).stat_summary('type', 'sem', 'geom', {'bar', 'black_errorbar'}, 'setylim', true);
g(1,1).set_title(['fraction coactive pairs per day']);
g(1,1).set_names('x','timepoint', 'y', 'coactivation probability' , 'color', 'Group', 'lightness', 'timepoint');
g(2,1).geom_jitter('width',0.2,'height',0, 'dodge', 1);
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
figure
g.draw();
savefigALP(savefigdir, 'coactivationprobability_byregion_dayavg')

clear g
g(1,1) = gramm('x', cellstr(dayActivityData.timepoint), 'y', dayActivityData.mean_fracActive, ...
    'color', cellstr(dayActivityData.group), 'lightness', cellstr(dayActivityData.timepoint), ...
    'subset', dayActivityData.mean_nDayRipples > 10);
g(1,1).facet_grid([], cellstr(dayActivityData.activeRegions));
g(2,1) = copy(g(1,1));
g(1,1).stat_summary('type', 'sem', 'geom', {'bar', 'black_errorbar'}, 'setylim', true);
g(1,1).set_title(['fraction Active pairs per day']);
g(1,1).set_names('x','timepoint', 'y', 'activation probability', 'color', 'Group', 'lightness', 'timepoint');
g(2,1).geom_jitter('width',0.2,'height',0, 'dodge', 1);
g.set_color_options('map', gcolors, 'n_color', 2, 'n_lightness', 2);
figure
g.draw();
savefigALP(savefigdir, 'activationprobability_byregion_dayavg')


inclDat = strcmp(CoactivityData.timepoint, 'post') & strcmp(CoactivityData.coactiveRegions, 'CA1CA1') & CoactivityData.nDayRipples > 10;
dataSubset = CoactivityData(inclDat,:);

figure
hold on
groupcolors = params.colors;
vdat.random = dataSubset(strcmp(dataSubset.group, 'random'),:).fracCoactive;
vdat.gamma = dataSubset(strcmp(dataSubset.group, 'gamma'),:).fracCoactive;
cmat = [groupcolors.random.post; groupcolors.gamma.post];
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', true);
ylabel('Coactivity - CA1')
p = ranksum(vdat.gamma, vdat.random);
title(['ranksum p = ', num2str(p)])


end

