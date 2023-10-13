function out = cf_rippledecoding_window(index, trainingData, time, theta, spikes, rippleWindows, rippleMids, params, inclCells, nCells)
%cf_rippledecoding_timebins
%
%ALP 1/11/2023

%%% get stuff from params
edges = params.posEdges;
nDeg = params.nDeg;
binsize = params.timeAroundMid*2;

%%% directories
figdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_ripples\', num2str(index(1,2)), '\'];
figname = ['rippledecoding_win', num2str(binsize*1000), '_' , num2str(index(1,2)), '_', num2str(index(1,3)), '_'];

%%% testing data
%bins around the ripple time
testingData = cf_gettestingdata_ripples(spikes.spiketimes, spikes.spikeIDs, rippleWindows, binsize, inclCells, nCells);

%%% decode ripple information
% get information about the peak location, try defining a prospective
% coding ratio etc
decodedData = cf_getdecodedposition_ripples(trainingData, testingData);

%%% compare decoded ripple to postion in the VR, save figures, etc
vrtimes = [time(1) time(end)];
vrRipples = isExcluded(rippleMids, vrtimes);
isVRripple = logical(vrRipples);

vrrippleinds = lookup2(rippleMids, time);
ripplePos = theta(vrrippleinds);
ripplePos(~isVRripple) = NaN; 

posBinsize = diff(edges);
posBinsize = posBinsize(1);
adjustedEdges = -nDeg/2:posBinsize:nDeg/2;

rippleDecPos = NaN(size(rippleMids,1),1);
rippleRelPos = NaN(size(rippleMids,1),1);
for i = 1:size(rippleMids)
    goodBin = decodedData.isDecodeBin{i};
    eventSpatialProb = decodedData.spatialProb{i};
    
    if isempty(eventSpatialProb)
        continue
    end
    
    if ~isnan(ripplePos(i))
        [~, iMax] = max(decodedData.spatialProb{i}, [], 'omitnan'); %should be 1 dimension
        rippleDecPos(i) = edges(iMax);
        rippleRelPos(i) = edges(iMax) - ripplePos(i);
    end
    
%     figure('Position', [243 154 1496 815])
%     hold on
%     subplot(1,2,1)
%     hold on
%     imagesc(edges(1:end-1), 1:size(eventSpatialProb,2), eventSpatialProb')
%     xlim([min(edges(1:end-1)) max(edges(1:end-1))])
%     ylim([0.5 1.5])
%     xlabel('decoded position (deg)')
%     ylabel('time bin')
%     colorbar
%     subplot(1,2,2)
%     hold on
%     plot(edges(1:end-1), eventSpatialProb, 'k-')
%     plot(ripplePos(i).*ones(1,2), [0 max(eventSpatialProb)+0.1], 'r-')
%     xlabel('trackposition (deg)')
%     xlim([-81 99])
%     ylabel('probability')
%     sgtitle(['ripple ', num2str(i), ' of ', num2str(size(rippleMids)), ' pos ', num2str(ripplePos(i)), ' relPos ', num2str(rippleRelPos(i))])
%     makefigurepretty(gcf)
%     newfigname = [figname, 'ripple', num2str(i), 'of', num2str(size(rippleMids,1))];
%     savefigALP(figdir, newfigname, 'filetype', 'pdf')
%     close
    clear event*
end

out.testingData = testingData;
out.decodedData = decodedData;
out.ripplePos = ripplePos;
out.rippleDecPos = rippleDecPos; 
out.rippleRelPos = rippleRelPos;





end

