function cf_rippleproperties_rate_RZ(dirs, params, allindex, metadata)
%cf_rippleproperties_rate_RZ
%
%ALP 8/29/23
%ALP 10/5/23

%% set up stuff
%%% parameters
params.speedThreshold = 1;
params.minTrials = 0;
params.RZ = [0 18];
params.AZ = [-36 0];
params.posEdges = -81:3:99;
params.VR_samprate = 50; 
dType = '180'; 

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

%% do calculation

for d = 1:size(dayindex,1)
    %%% load trial data
    load([trialdatadir, trialdatafilename, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
    load([trialdatadir, trialspikesfilename, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
    
    %%% get good trials
    inclTrials = [trialData.engaged] & [trialData.fullTrial] & [trialData.rewarded];
    goodTrials = trialData(inclTrials);
    
    RZData = []; 
    for t = 1:length(goodTrials)
        RZperiods = [];
        %%% get time in RZ below a speed threshold for each trials
        isRZ = goodTrials(t).theta_d2r >= params.RZ(1) & goodTrials(t).theta_d2r < params.RZ(2);
        isSlow = goodTrials(t).speed < params.speedThreshold; 
        isInclude = isRZ&isSlow; 
        RZtime = sum(isRZ)./params.VR_samprate;

        if sum(isInclude) > 0
            inRZ = contiguous(isInclude, 1); %get runs of 1s
            inRZ = inRZ{1,2}; % access the runs
            RZperiods = [goodTrials(t).time(inRZ(:,1)) goodTrials(t).time(inRZ(:,2))];
        else
            RZperiods = [];
        end
        
        RZData(t).RZ_time_sum = RZtime;
        RZData(t).RZ_time_per = RZperiods;
        RZData(t).file = goodTrials(t).file;
    end
    
    %%% get best ripple channel
    %%% load best ripple channel for the day
    ripplechandir = ['\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\CA1\'];
    if ~exist([ripplechandir, 'bestChannel_', num2str(dayindex(d,2)), '.mat'])
        continue
    end
    load([ripplechandir, 'bestChannel_', num2str(dayindex(d,2)), '.mat'])
    rippleChan = bestRippleChan.all;
    
    files = unique([goodTrials.file]);
    numRipples = [];
    for f = 1:length(files)
        %%% load ripples
        anperiodsdatadir = [dirs.processeddatadir, params.iden, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\CA1\', num2str(rippleChan), '\'];
        load([anperiodsdatadir, 'ripples', num2str(files(f)), '.mat'])
        
        %%% give structures easy names for interfacing
        RipData = ripples{dayindex(d,1)}{dayindex(d,2)}{files(f)};
        
        %%% get RZ periods per file
        isFile = [RZData.file] == files(f);
        fileRZperiods = cell2mat({RZData(isFile).RZ_time_per}');

        %%% get # ripples in each trial
        isRZripple = isExcluded(RipData.midtime, fileRZperiods); 
        isRZripple = logical(isRZripple); 
        
        %%% ripples per file in RZ below speed threshold
        numRipples(f) = sum(isRZripple);
    end
    
    %%% then will want to get the # of ripples/total time below the velocity
    %%% threshold
    if isempty(numRipples)
        DayData(d).nRipples = 0;
    else
        DayData(d).nRipples = sum(numRipples);
    end
    
    if isempty(RZData)
        DayData(d).RZ_time_s = 0;
    else
        DayData(d).RZ_time_s = sum([RZData.RZ_time_sum]);
    end
    
    if isempty(numRipples)
        DayData(d).rippleRateS = 0;
    else
        DayData(d).rippleRateS = sum(numRipples)/sum([RZData.RZ_time_sum]);
    end
    
    %%% add metadata to the structure
    if metadata.Date(d) == dayindex(d,2)
        DayData(d).group = metadata.Groups(d);
        if metadata.FlickerDay(d) < 5
            DayData(d).timepoint = 'pre';
        elseif metadata.FlickerDay(d) > 5
            DayData(d).timepoint = 'post';
        end
        DayData(d).animal = dayindex(d,1);
        DayData(d).date = dayindex(d,2);
    end    
    
end

RateDataRZ = struct2table(DayData);

%%% save data as a table structure for statistics
statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
statsfilename = 'TableData_RippleRate_RZ.txt';
writetable(RateDataRZ, fullfile(statsdir, statsfilename))

savefilename = 'ripplerate_RZ.mat';
savedatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\ripple_properties\';
save([savedatadir, savefilename], 'RateDataRZ')

%% plot the ripple rate in the reward zone
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'}; 
xvect = [1 2];
groupcolors = params.colors; 
isPost = strcmp(RateDataRZ.timepoint, 'post');
isEnough = RateDataRZ.nRipples > 5; 
RateDataRZ = RateDataRZ(isPost&isEnough,:);
PlotData = RateDataRZ; 

fnames = {'rippleRateS'};
analysisname = {'ripple rate'};

nAnGamma = length(unique(PlotData.animal(strcmp(PlotData.group, 'gamma'))));
nAnRandom = length(unique(PlotData.animal(strcmp(PlotData.group, 'random'))));

disp(['There are ', num2str(nAnGamma) ' animals in 40 Hz RZ ripple rate'])
disp(['There are ', num2str(nAnRandom) ' animals in Random RZ ripple rate'])

fh = figure;
hold on
for f = 1
    for d = 2
        scatterdat = []; bdat = []; bstde = []; cmat = [];
        for g = 1:length(gnames)
            tmpdat = [];
            isGroup = strcmp(RateDataRZ.group, gnames{g});
            tmpdat = RateDataRZ.(fnames{f})(isGroup);
            scatterdat{g} = tmpdat;
            bdat(g) = nanmean(tmpdat);
            bstde(g) = nanstd(tmpdat)./sqrt(sum(~isnan(tmpdat)));
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
        plotprettypoints(fh, xvect(f,:), scatterdat, 'color', cmat)
        b = bar(xvect(f,:), bdat, 'FaceColor', 'flat');
        b.CData = cmat;
        b.FaceAlpha = 0.6;
        errorbar2(xvect(f,:), bdat, bstde, 0.2, 'k-', 'LineWidth', 0.75);
    end

end
ylabel('Ripple abundance (Hz)')
title('ripple rate - RZ')






end

