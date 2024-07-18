function cf_decoding_thetaseq_360dec_180trialsplit_seqlength(dirs, params, allindex, metadata, posType, cellTypes)
%cf_decoding_thetaseq_360dec_180trialsplit_seqlength
%   this version does the theta sequence decoding with the 360 degree
%   position variable but splits the trials by 180degrees
%
%ALP 3/5/2024

%%% rewriting files
rewrite.decodingfiles = 0;

%%% params
% position parameters
if strcmp(posType, 'full')
    params.posEdges = 0:2:360;
    params.posBins = 2;
    dType = '360';
    params.nDeg = 360;
    posName = 'theta';
    params.RZ = [54 72; 234 252];
    params.dec_edges = 0:6:360;
    params.RRZ = [36 64; 216 244];
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

% theta sequence parameters
params.decodingBinsize = 0.02; %in s
params.timeAroundTrough = 0.2; %in s
params.speedThreshold = 1; %deg/s
params.minTrials = 5;
params.adjustedEdges = -params.nDeg/2:params.posBins:params.nDeg/2;
params.timeEdges = -params.timeAroundTrough:params.decodingBinsize:params.timeAroundTrough;
params.quadrantWindow = 0.16; %in s

%%% define areas of interest, where to get the quadrant ratio, etc
%%% set edges for histograms, areas for quantification, etc
params.PCR_positions{1} = [0-18 0+9];
params.PCR_positions{2} = [9 9+27];

%%% define a different area of interest
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
%dayindex = dayindex(dayindex(:,2) ~= 210221,:); %having issues with saving this file, excluding for now cuz pre. ALP 3/7/24
%dayindex = dayindex(dayindex(:,2) ~= 210221,:); %having issues with saving this file, excluding for now cuz pre. ALP 3/7/24

%%% directories
savedatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
trialdatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\';
ripplechandir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\CA1\';
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';

%%% filenames
%trialdatafilename = ['trialInfo_', dType, '_A'];
trialdatafilename = ['trialInfo_180_A'];
trialspikesfilename = ['trialSpikes_180_A'];
savefilename = ['thetaseqdecoding_alldays_alltrials_360dec_180trialsplit.mat'];
AllData = table;
AllTrialData = table;


%% loop over all days
for d = 1:size(dayindex,1)
    index = allindex(allindex(:,2) == dayindex(d,2),:);
    files = index(:,3);
    anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
    tmpfilename = [savedatadir, 'thetaseq_decoding_', dType, '_180trialsplit_'];
    dayfilename = [tmpfilename, num2str(dayindex(d,2)), '.mat'];
    day = d;
    
    if ~exist(dayfilename, 'file') || rewrite.decodingfiles
        disp(['decoding day ', num2str(dayindex(d,2))])
        
        % load trial data
        load([trialdatadir, trialdatafilename, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
        load([trialdatadir, trialspikesfilename, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
        
        
        %%%% ------ get cells of interest ------ %%%%
        %place cells, old versions of this code had options for
        %multiple types of cells
        load([ratemapdir, 'spatialmaps_all_500_full_', num2str(dayindex(d,2)), '.mat'])
        allIDs = trialSpikes.info.clusterID;
        nCells = length(allIDs);
        cellI = 1:length(allIDs);
        inclCells = cellI;
        
        inclPC = ratemaps.includeCell;
        inclCells = cellI(inclPC);
        nCells = length(inclCells);
        
        if nCells == 0
            continue
        end
        
        %%% ------ load theta information ----- %%%%
        load([ripplechandir, 'bestChannel_', num2str(dayindex(d,2)), '.mat'])
        tmpRChan = bestRippleChan.all;
        theta = loaddatastruct2([anprocesseddatadir, 'CA1\', num2str(tmpRChan), '\EEG\'], dayindex(d,:), 'theta', files);
        
        
        % get trials to include
        inclTrials = true(1,length(trialData)); %include all trials, exclude later
        goodTrials = trialData(inclTrials);
        nTrials = sum(inclTrials);
        if nTrials < 1
            continue
        end
        
        % set up spike times, spike IDs
        [goodTrials.spikeTimes] = deal(trialSpikes.data(inclTrials).spikeTimes);
        [goodTrials.spikeIDs] = deal(trialSpikes.data(inclTrials).IDs);
        
        % load spatial rate maps - training data
        load([ratemapdir, 'spatialmaps_preflicker_0_', posType, '_', num2str(dayindex(d,2)), '.mat'])
        trainingData.ratemaps = ratemaps.fr; %in Hz
        
        %%% get testing data for each trial
        [testingData, avgTheta] = cf_gettestingdata_thetaseq(index, goodTrials, theta, ...
            params.decodingBinsize, params.timeAroundTrough, params.speedThreshold, posName, nCells, inclCells);
        thetatrace = table(day, avgTheta);
        
        %%% decode position during theta trough
        decodedData = cf_getdecodedposition_thetaseq_table(trainingData, testingData, nCells, inclCells, ...
            params);
        
        %%% do some calculations on the decoded sequences
        codingRatio = cf_getcodingratio_thetaseq_table(decodedData, params.PCR_positions, params);
        
        %%% do some calculations on the decoded sequences
        %         codingRatio2 = cf_getcodingratio_thetaseq_table(decodedData, params.PCR_secondhalf, params);
        %         newVarNames = append("PostR_", fieldnames(codingRatio2));
        %         newVarNames = newVarNames(1:5); %changed from 1:4 ALP 7/18
        %         codingRatio2.Properties.VariableNames = newVarNames;
        
        %%%% ----- use the 180 data to determine sequence significance
        %%% load behavioral metrics
        behaviordatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\'];
        behaviorfilename = ['behavioranalysis_prepost_recordings_table_', num2str(dayindex(d,2)), '.mat'];
        load([behaviordatadir, behaviorfilename], 'metrics')
        
        %%% figure out if this day has significant sequences or not,
        %%% compare full quadrant ratio on mean day sequence to shuffled data
        sig = cf_getsignificant_thetaseq_table(decodedData, params, metrics);
        
        %save the significance information independently
        sigfilename = ['thetasequence_significance_PC_', num2str(dayindex(d,2)), '.mat'];
        
        %%% compile data of interest
        dayData = [metrics(:, end-4:end) metrics(:,25:27) metrics(:, 19) metrics(:,15:18) metrics(:,3:14) decodedData codingRatio];
        dayData = outerjoin(dayData, sig, 'MergeKeys', true);
        
        if isempty(decodedData)
            continue
        end
        
        save([tmpfilename, num2str(dayindex(d,2)), '.mat'], 'testingData', 'decodedData', 'codingRatio', 'avgTheta')
    else
        load([tmpfilename, num2str(dayindex(d,2)), '.mat'])
%         day = table(dayindex(d,2)*ones(height(decodedData),1), 'VariableNames', {'day'});
%         dayData = [decodedData codingRatio day];
        
    end
    
    load([trialdatadir, trialdatafilename, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
    trialNum = 1:length(trialData);
    tmpTrialData = table([trialData.engaged]', [trialData.fullTrial]', [trialData.rewarded]', ...
        [trialData.an]', [trialData.day]', [trialNum]', 'VariableNames', {'engaged', 'fullTrial', 'rewarded', 'animal', 'day', 'trial'});
    
    
    AllData = [AllData; dayData];
    AllTrialData = [AllTrialData; tmpTrialData];
    
    clear appendDat decodedData testingData trainingData avgTheta tmp* sig codingRatio codingRatio2 tmpData tmpTrialData dayData
end

save([savedatadir, savefilename], 'AllData', '-v7.3')


%% get average sequences at the reward zones
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\decoding_360_pertrial_RZ\';

% per trial and per day out of curiosity
% per trial just for examples
dec_edges = -params.nDeg/2:params.posBins:params.nDeg/2;
dec_edges = dec_edges(1:end-1);

for d = 1:size(dayindex,1)
    dayData = AllData(AllData.day == dayindex(d,2),:);
    for t = 1:height(dayData)
        pos = cell2mat(dayData.troughPos(t));
        seq = cell2mat(dayData.trialSeq(t));
        
        if size(params.RZ,1) == 1 %distance to reward decoding
            isRZ = pos >= params.RRZ(1) & pos < params.RRZ(2);
        elseif size(params.RZ,1) == 2 %360 deg decoding
            isRZ = (pos >= params.RZ(1,1) & pos < params.RZ(1,2)) | (pos >= params.RZ(2,1) & pos < params.RZ(2,2));
        end
        
        RZseq = seq(:,:,isRZ);
        mRZseq{t,1} = mean(RZseq, 3, 'omitnan');
        trial_sigOtherDecoding = getZoneAboveThresh(mRZseq{t}, params.nDeg);
        
        %plot the per trial sequences
        figure
        hold on
        imagesc(params.timeEdges, dec_edges, mRZseq{t})
        colorbar
        xlabel('time from theta trough (s)')
        xlim([-0.2 0.2])
        ylim([-180 180])
        ylabel('relative decoded position')
        title({['theta sequence in the RZ, trial ', num2str(t) ' of ', num2str(height(dayData)), ' - ', num2str(dayindex(d,2))], ...
            [' sig other zone ', num2str(trial_sigOtherDecoding)]})
        figname = ['thetaseqdecoding_', dType, '_', num2str(dayindex(d,2)), '_', num2str(t), '_sig', num2str(trial_sigOtherDecoding)];
        saveas(gcf,[figdir, figname], 'png')
        
        clear pos seq RZseq
    end
    
    dRZseq = cat(3,mRZseq{:});
    mdRZseq = mean(dRZseq, 3, 'omitnan');
    
    sigOtherDecoding = getZoneAboveThresh(mdRZseq, params.nDeg);
    
    %plot!
    figure
    hold on
    imagesc(params.timeEdges, dec_edges, mdRZseq)
    colorbar
    xlabel('time from theta trough (s)')
    xlim([-0.2 0.2])
    ylim([-180 180])
    ylabel('relative decoded position')
    title(['theta sequence in the RZ - sig other zone ', num2str(sigOtherDecoding)])
    figname = ['thetaseqdecoding_', dType, '_', num2str(dayindex(d,2)), '_dayavg'];
    saveas(gcf,[figdir, figname], 'png')
    
    close all
    clear mRZseq dayData
end

end

function out = getZoneAboveThresh(seq, nDeg)

%create a threshold to determine if there is "significant" decoding of
%the other reward zone
vectseq = reshape(seq, 1, []); %unwrap to get the mean
mProb = mean(vectseq); % overall mean
sdevProb = std(vectseq); %overall std dev

% does anything in the other reward region cross the threshold?
%not sure exactly where this should be... going to do 1 zone at the
%beginning and end
if nDeg == 360
    otherRZ = seq([1:9, 172:end],:);
elseif nDeg == 180
    otherRZ = seq([1:6, 175:end],:);
end
vectOtherRZ = reshape(otherRZ,1,[]);
out = any(vectOtherRZ > (2*sdevProb+mProb));

end

