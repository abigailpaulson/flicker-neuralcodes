function out = cf_rippledecoding_timebins(index, trainingData, time, theta, spikes, rippleWindows, rippleMids, params, inclCells, nCells)
%cf_rippledecoding_timebins
%
%ALP 1/11/2023

%%% get stuff from params
edges = params.posEdges;
nDeg = params.nDeg;

%%% directories

figdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_ripples\RZRipples_360\'];
figname = ['rippledecoding_', num2str(nDeg), '_bins', num2str(params.decodingBins*1000), '_' , num2str(index(1,2)), '_', num2str(index(1,3)), '_'];

%%% testing data
%bins around the ripple time
testingData = cf_gettestingdata_ripples(spikes.spiketimes, spikes.spikeIDs, rippleWindows, params.decodingBins, inclCells, nCells);

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
if nDeg == 180
    prospBins = [find(adjustedEdges == 9):1:find(adjustedEdges == 18+45)-1];
    retroBins = [find(adjustedEdges == -18-45):1:find(adjustedEdges == -9)-1];
else
    prospBins = [find(adjustedEdges == 8):1:find(adjustedEdges == 18+44)-1];
    retroBins = [find(adjustedEdges == -18-44):1:find(adjustedEdges == -8)-1];
end

%%% get peak pos from the trainingData for sorting the spiking activity
maps = trainingData.ratemaps(inclCells,:);
[~,iMax] = max(maps,[],2);
[~,iSort] = sort(iMax);
cellIDSort = inclCells(iSort);


rippleRatio = NaN(size(rippleMids,1),1);
for i = 1:size(rippleMids,1)
    
    %%%
    goodBin = decodedData.isDecodeBin{i};
    
    if isempty(goodBin)
        continue
    end
    
    tempcolors = jet(sum(decodedData.isDecodeBin{i}));
    eventTrainingData = decodedData.trainingCounts{i}; 
    eventSpatialProb = decodedData.spatialProb{i};
    [~, eventTrainingPeaks] = max(eventTrainingData, [], 2);
    [~, eventCellOrder] = sort(eventTrainingPeaks);
    
    
    if ~isnan(ripplePos(i))
        posdiff = edges(1:end-1) - ripplePos(i);
        posdiff(posdiff<-(nDeg/2)) = posdiff(posdiff <-(nDeg/2)) + nDeg;
        posdiff(posdiff>(nDeg/2)) = posdiff(posdiff > (nDeg/2)) - nDeg;
        [~, iorder] = sort(posdiff);
        adjustedProb = eventSpatialProb(iorder,:);
        
        goodBinsOnly = eventSpatialProb(:, goodBin);
        probvals = [sum(sum(goodBinsOnly(prospBins,:))) sum(sum(goodBinsOnly(retroBins,:)))];
        rippleRatio(i) = (probvals(1) - probvals(2))/(probvals(1)+probvals(2));
    end
    
%%% figure out if the ripple has significant other zone decoding
%create a threshold to determine if there is "significant" decoding of
%the other reward zone
vectseq = reshape(eventSpatialProb, 1, []); %unwrap to get the mean
mProb = mean(vectseq, 'omitnan'); % overall mean
sdevProb = std(vectseq, 'omitnan'); %overall std dev

% does anything in the other reward region cross the threshold?
%not sure exactly where this should be... going to do 1 zone at the
%beginning and end
if ~isnan(ripplePos(i))
    if nDeg == 360
        otherRZ = adjustedProb([1:9, 172:end],:);
    elseif nDeg == 180
        otherRZ = adjustedProb([1:6, 175:end],:);
    end
    vectOtherRZ = reshape(otherRZ,1,[]);
    otherRZThresh = 2*sdevProb+mProb;
    Event_SigOtherZone(i) = any(vectOtherRZ > (2*sdevProb+mProb));
else
    Event_SigOtherZone(i) = false;
end


%% plotting for examples for the paper - raster plot, heatmap, etc.
% if (ripplePos(i) >= 54 && ripplePos(i) < 72) || (ripplePos(i) >= 234 && ripplePos(i) <252) 
% %     figure('Position', [243 154 1496 815])
% %     hold on
% %     subplot(1,4,1)
% %     hold on
% %     imagesc(edges(1:end-1), 1:size(eventSpatialProb,2), eventSpatialProb')
% %     xlim([min(edges(1:end-1)) max(edges(1:end-1))])
% %     ylim([1 size(eventSpatialProb,2)])
% %     xlabel('decoded position (deg)')
% %     ylabel('time bin')
% %     colorbar
% %     subplot(1,4,2)
% %     hold on
% %     icolor = 1;
% %     for b = 1:length(goodBin)
% %         if goodBin(b)
% %             plot(edges(1:end-1), eventSpatialProb(:,b), '-', 'Color', tempcolors(icolor,:))
% %             icolor = icolor+1;
% %         end
% %     end
% %     xlabel('trackposition (deg)')
% %     %xlim([-81 99])
% %     ylabel('probability')
% %     subplot(1,4,3)
% %     hold on
% %     hold on
% %     xfillvals = adjustedEdges(retroBins);
% %     fill([xfillvals(1) xfillvals(1) xfillvals(end) xfillvals(end)], [0 max(max(adjustedProb)) max(max(adjustedProb)), 0], 'm', 'EdgeColor', 'none')
% %     xfillvals = adjustedEdges(prospBins);
% %     fill([xfillvals(1) xfillvals(1) xfillvals(end) xfillvals(end)], [0 max(max(adjustedProb)) max(max(adjustedProb)), 0], 'm', 'EdgeColor', 'none')
% %     alpha(0.1)
% %     icolor = 1;
% %     for b = 1:length(goodBin)
% %         if goodBin(b)
% %             plot(adjustedEdges(1:end-1), adjustedProb(:,b), '-', 'Color', tempcolors(icolor,:))
% %             icolor = icolor+1;
% %         end
% %     end
% %     plot([adjustedEdges(1) adjustedEdges(end-1)], [otherRZThresh otherRZThresh], 'k--')
% %     xlabel('relative position (deg)')
% %     %xlim([-90 90])
% %     ylabel('probability')
% %     subplot(1,4,4)
% %     hold on
% %     %spiking
% %     iPlot = 1;
% %     isEventSpike = logical(isExcluded(spikes.spiketimes, rippleWindows(i,:)));
% %     eventSpikeTimes = spikes.spiketimes(isEventSpike);
% %     eventSpikeIDs = spikes.spikeIDs(isEventSpike);
% %     for c = 1:length(cellIDSort)
% %         isCellSpike = eventSpikeIDs == cellIDSort(c);
% %         cellSpikeTimes = eventSpikeTimes(isCellSpike);
% %         if isempty(cellSpikeTimes)
% %             continue
% %         end
% %         plot_raster(cellSpikeTimes, iPlot*ones(1,length(cellSpikeTimes)))
% %         iPlot = iPlot +1;
% %     end
% %     xlim([rippleWindows(i,1) rippleWindows(i,2)])
% %     
% %     sgtitle(['ripple ', num2str(i), ' of ', num2str(size(rippleMids)), ' position ', num2str(ripplePos(i)), ' ratio ', num2str(rippleRatio(i)), ' sig other zone ', num2str(Event_SigOtherZone)])
% %     makefigurepretty(gcf)
% %     newfigname = [figname, 'ripple', num2str(i), 'of', num2str(size(rippleMids,1))];
% %     savefigALP(figdir, newfigname, 'filetype', 'png')
% 
%     figure
%     hold on
% %     xfillvals = adjustedEdges(retroBins);
% %     fill([xfillvals(1) xfillvals(1) xfillvals(end) xfillvals(end)], [0 max(max(adjustedProb)) max(max(adjustedProb)), 0], 'm', 'EdgeColor', 'none')
% %     xfillvals = adjustedEdges(prospBins);
% %     fill([xfillvals(1) xfillvals(1) xfillvals(end) xfillvals(end)], [0 max(max(adjustedProb)) max(max(adjustedProb)), 0], 'm', 'EdgeColor', 'none')
% %     alpha(0.1)
%     icolor = 1;
%     for b = 1:length(goodBin)
%         if goodBin(b)
%             plot(adjustedEdges(1:end-1), adjustedProb(:,b), '-', 'Color', tempcolors(icolor,:))
%             icolor = icolor+1;
%         end
%     end
%     plot([adjustedEdges(1) adjustedEdges(end-1)], [otherRZThresh otherRZThresh], 'k--')
%     xlabel('relative position (deg)')
%     %xlim([-90 90])
%     ylabel('probability')
%     sgtitle(['ripple ', num2str(i), ' of ', num2str(size(rippleMids)), ' position ', num2str(ripplePos(i)), ' ratio ', num2str(rippleRatio(i)), ' sig other zone ', num2str(Event_SigOtherZone(i))])
%     newfigname = [figname, 'ripple', num2str(i), 'of', num2str(size(rippleMids,1))];
%     savefigALP(figdir, newfigname, 'filetype', 'png')
%     
%     close
% end
    clear event*
end

out.testingData = testingData;
out.decodedData = decodedData;
out.ripplePos = ripplePos;
out.rippleRatio = rippleRatio;
out.sigOtherZone = Event_SigOtherZone;





end

