function cf_placecellclassifier(dirs, params, dayindex, files)
%cf_placecellclassifier
%
% ALP 12/19/2022

%% set stuff
%%% params
index = [repmat(dayindex,[length(files),1]) files];
edges = params.place_edges;

%%% rewriting 
rewrite.placecellfiles = 1; 

%%% file names, directories
anspikedir = [dirs.processeddatadir, 'A', num2str(dayindex(1,1)), '_', num2str(dayindex(1,2)), '\'];
positiondir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\';
savedir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';
filename = ['spatialmaps_', params.condition '_', num2str(params.nShuff), '_', params.placecells_postype, '_', num2str(dayindex(1,2)), '.mat'];


if exist([savedir, filename], 'file') && ~rewrite.placecellfiles 
    return
end
%% classify place cells for the recording day
%%% load spikes
[spikes, unitInfo, spikevect] = cf_getspikes(anspikedir, dirs, params, index, params.brainReg);
nCells = length(unitInfo);
cellIDs = 1:1:nCells;

%%% load position information
load([positiondir, 'behaviorInfo_A', num2str(dayindex(1,1)), '_', num2str(dayindex(1,2)), '.mat'])

%%% loop over recording files
position_h = zeros(nCells, length(edges)-1);
occ_h = zeros(1,length(edges)-1);
shuff_position_h = repmat({position_h},1,500);
for f = 1:length(files)
    fID = files(f);
    %ALP commented out below bc index into behvior info using fID, vs
    %positionInfo which needs posID 1/8/2023
%     posindstemp = [behaviorInfo.index];
%     posinds = reshape(posindstemp, [3 length(posindstemp)/3])';
%     posFiles = posinds(:,3);
%     posID = find(posFiles == files(f)); 
%     if isempty(posID)
%         continue
%     end
    
    %%% initialize
    pos_count = []; occ_time = [];
    
    %%% temp spiking vector
    spiketimes = spikevect(fID).spikeTimes;
    spikeIDs = spikevect(fID).spikeIDs;
    
    %%% check indexing into behaviorInfo correctly
    if ~(behaviorInfo(fID).index(1,3) == files(f))
        error('check posInfo indexing')
    end
    
    %%% get position, time, speed
    if strcmp(params.placecells_postype, 'full')
        pos = behaviorInfo(fID).theta;
    elseif strcmp(params.placecells_postype, 'd2r')
        pos = behaviorInfo(fID).theta_d2r;
    end
    
    if isempty(pos)
        disp(['skipping file ', num2str(fID), ' ', num2str(dayindex(1,2))])
        continue
    end
    
    speed = behaviorInfo(fID).speed;
    time = behaviorInfo(fID).time;
    
    pos_count = cf_calcpositionhist(cellIDs, spiketimes, spikeIDs, time, pos, speed, edges, params.speedThreshold, 0);
    
    isMoving = speed > params.speedThreshold;
    
    if params.nShuff > 0
        for s = 1:params.nShuff
            shuff_pos_count{s} = cf_calcpositionhist(cellIDs, spiketimes, spikeIDs, time, pos, speed, edges, params.speedThreshold, 1);
            shuff_position_h{s} = shuff_position_h{s}+shuff_pos_count{s};
        end
    end
    
    %%% calculate occupancy
    vr_samprate = 1/0.02; %bins/s
    occ_counts = histcounts(pos(isMoving), edges);
    occ_time = occ_counts./vr_samprate;
    
    %%% add to full matrix
    position_h = position_h + pos_count;
    occ_h = occ_h + occ_time;
end

%%% smooth
for u = 1:size(position_h,1)
    ratemaps.pos(u,:) = gaussSmooth(position_h(u,:), 2);
end
ratemaps.occ = gaussSmooth(occ_h,2);

if any(ratemaps.occ < 0.1*2)
    ratemaps.pos = NaN(nCells, length(edges)-1);
end

%%% rate map, spatial info, sparsity
ratemaps.fr = ratemaps.pos./ratemaps.occ;
ratemaps.spatialinfo = cf_calcspatialinfo(ratemaps.pos, ratemaps.occ);
nBins = length(edges)-1;

%sparsity Ravassard et al. Science 2013
%0 uniform map, 1 is delta function
for u = 1:size(ratemaps.fr,1)
    ratemaps.sparsity(u) = (1-((1/nBins)*((sum(ratemaps.fr(u,:))^2)/(sum(ratemaps.fr(u,:).^2)))))*(nBins/(nBins-1));
end

%%% shuffles
if params.nShuff > 0
    for s = 1:params.nShuff
        for u = 1:size(position_h,1)
            ratemaps.shuffpos{s}(u,:) = gaussSmooth(shuff_position_h{s}(u,:),2);
        end
%         ratemaps.shuff_fr{s} = ratemaps.shuffpos{s}./ratemaps.occ;
        ratemaps.shuff_fr{s} = ratemaps.shuffpos{s};
        ratemaps.shuff_SI(:,s) = cf_calcspatialinfo(ratemaps.shuff_fr{s}, ratemaps.occ);
    end
    ratemaps.SI_thresh = prctile(ratemaps.shuff_SI, 95, 2);
    ratemaps.SI_thresh_param = params.SI_prctile;
    ratemaps.nShuff = params.nShuff;
end

%%% ouptut structure
ratemaps.cellIDs = cellIDs;
ratemaps.brainReg = {unitInfo.brainReg};

allmap = ratemaps.fr;
[maxfr,iMax] = max(allmap, [], 2, 'omitnan');
[~, iSortAll] = sort(iMax);

ratemaps.peakpos = edges(iMax);
ratemaps.maxFR = maxfr;
ratemaps.meanFR = nanmean(ratemaps.fr,2);

%%% field location
ratemaps.fieldbounds = 0.95*ratemaps.maxFR; %Safaryan and Mehta Nat Neuro 2021
for u = 1:size(ratemaps.fr,1)
    isField = ratemaps.fr(u,:) > ratemaps.fieldbounds(u);
    fieldruns = [];
    if sum(isField) > 0
        fieldruns = contiguous(isField,1);
        fieldruns = fieldruns{1,2};
    end
    locsvect = edges(isField);
    ratemaps.fieldlocations{u} = locsvect;
    ratemaps.fieldruns{u} = fieldruns;
end

%%% does it meet PC criteria?
if params.nShuff
    isThresh = ratemaps.spatialinfo > ratemaps.SI_thresh;
    isFR = (ratemaps.maxFR >= params.minpeakFR) & (ratemaps.meanFR < params.maxmeanFR) & (ratemaps.meanFR >= params.minmeanFR);
    isPYR = strcmp({unitInfo.cellType}, 'PYR');
    includeCell = isThresh & isFR & isPYR';
    
    ratemaps.isThresh = isThresh;
    ratemaps.isFR = isFR;
    ratemaps.isPYR = isPYR';
    ratemaps.includeCell = includeCell;
end

ratemaps.dayindex = dayindex;
ratemaps.files = files;

%%% save ratemaps
save([savedir, filename], 'ratemaps')

%%% plot per day
if params.nShuff > 0
    
    %%% plotting stuff
    isCA1 = strcmp(ratemaps.brainReg, 'CA1')';
    ca1.all_sort = cf_sortbymax(ratemaps.fr(strcmp(ratemaps.brainReg, 'CA1'),:), 'norm2max');
    ca1.include_sort = cf_sortbymax(ratemaps.fr(isCA1&includeCell,:), 'norm2max');
    ca3.all_sort = cf_sortbymax(ratemaps.fr(~isCA1,:), 'norm2max');
    ca3.include_sort = cf_sortbymax(ratemaps.fr(~isCA1&includeCell,:), 'norm2max');
    
    frmap = ratemaps.fr(includeCell,:);
    [maxfrIncl,iMax] = max(frmap, [], 2, 'omitnan');
    [~, iSortIncl] = sort(iMax);
    
    figure('Position', [432 94 617 830])
    hold on
    subplot(3,3,1)
    hold on
    imagesc(edges(1:end-1)+1, 1:nCells, allmap(iSortAll,:)./maxfr(iSortAll), [0 0.75])
    % colormap(flipud(gray))
    xlim([0 360])
    ylim([1 nCells+1])
    xlabel('position (deg)')
    ylabel('cell #')
    title('all cells')
    subplot(3,3,2)
    hold on
    imagesc(edges(1:end-1)+1, 1:sum(includeCell), frmap(iSortIncl,:)./maxfrIncl(iSortIncl), [0 0.75])
    ylim([1 nCells+1])
    xlim([0 360])
    xlabel('position (deg)')
    ylabel('cell #')
    title('place cells')
    subplot(3,3,3)
    hold on
    imagesc(edges(1:end-1)+1, 1:nCells, allmap(iSortAll,:), [0 7])
    ylim([1 nCells+1])
    xlim([0 360])
    xlabel('position (deg)')
    ylabel('cell #')
    title('all cells - rawFR 0 to 7')
    sgtitle(num2str(dayindex(1,2)))
    subplot(3,3,4)
    hold on
    imagesc(edges(1:end-1)+1, 1:size(ca1.all_sort,1), ca1.all_sort, [0 0.75])
    xlim([0 360])
    ylim([0 size(ca1.all_sort,1)+1])
    xlabel('position (deg)')
    ylabel('cell #')
    title('CA1 - all norm')
    subplot(3,3,5)
    hold on
    imagesc(edges(1:end-1)+1, 1:size(ca3.all_sort,1), ca3.all_sort, [0 0.75])
    xlim([0 360])
    ylim([0 size(ca3.all_sort,1)+1])
    xlabel('position (deg)')
    ylabel('cell #')
    title('CA3 - all norm')
    subplot(3,3,7)
    hold on
    imagesc(edges(1:end-1)+1, 1:size(ca1.include_sort,1), ca1.include_sort, [0 0.75])
    xlim([0 360])
    ylim([0 size(ca1.all_sort,1)+1])
    xlabel('position (deg)')
    ylabel('cell #')
    title('CA1 - place cells norm')
    subplot(3,3,8)
    hold on
    imagesc(edges(1:end-1)+1, 1:size(ca3.include_sort,1), ca3.include_sort, [0 0.75])
    xlim([0 360])
    ylim([0 size(ca3.all_sort,1)+1])
    xlabel('position (deg)')
    ylabel('cell #')
    title('CA3 - place cells norm')
    
    figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';
    figname = ['ratemaps_', num2str(dayindex(1,2))];
    savefigALP(figdir, figname, 'filetype', 'pdf')
end




end

