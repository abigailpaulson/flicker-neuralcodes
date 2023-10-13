function [outputArg1,outputArg2] = cf_analyzedecoding_ripples(time, theta, decodedData, edges, rippleMids,nDeg)
%cf_anayzedecoding_ripples
%
%ALP 1/11/2023

vrtimes = [time(1) time(end)];
vrRipples = isExcluded(rippleMids, vrtimes);
isVRripple = logical(vrRipples);

vrrippleinds = lookup2(rippleMids, time);
ripplePos = theta(vrrippleinds);
ripplePos(~isVRripple) = NaN; 

posBinsize = diff(edges);
posBinsize = posBinsize(1);
adjustedEdges = -nDeg/2:posBinsize:nDeg/2;
prospBins = [find(adjustedEdges == 18):1:find(adjustedEdges == 18+45)-1];
retroBins = [find(adjustedEdges == -18-45):1:find(adjustedEdges == -18)-1];

for i = 1:size(rippleMids)
    goodBin = decodedData.isDecodeBin{i};
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
    
    figure('Position', [243 154 1496 815])
    hold on
    subplot(1,3,1)
    hold on
    imagesc(edges(1:end-1), 1:size(eventSpatialProb,2), eventSpatialProb')
    xlim([min(edges(1:end-1)) max(edges(1:end-1))])
    ylim([1 size(eventSpatialProb,2)])
    xlabel('decoded position (deg)')
    ylabel('time bin')
    colorbar
    subplot(1,3,2)
    hold on
    icolor = 1;
    for b = 1:length(goodBin)
        if goodBin(b)
            plot(edges(1:end-1), eventSpatialProb(:,b), '-', 'Color', tempcolors(icolor,:))
            icolor = icolor+1;
        end
    end
    xlabel('trackposition (deg)')
    xlim([-81 99])
    ylabel('probability')
    subplot(1,3,3)
    hold on
    hold on
    xfillvals = adjustedEdges(retroBins);
    fill([xfillvals(1) xfillvals(1) xfillvals(end) xfillvals(end)], [0 max(max(adjustedProb)) max(max(adjustedProb)), 0], 'm', 'EdgeColor', 'none')
    xfillvals = adjustedEdges(prospBins);
    fill([xfillvals(1) xfillvals(1) xfillvals(end) xfillvals(end)], [0 max(max(adjustedProb)) max(max(adjustedProb)), 0], 'm', 'EdgeColor', 'none')
    alpha(0.1)
    icolor = 1;
    for b = 1:length(goodBin)
        if goodBin(b)
            plot(adjustedEdges(1:end-1), adjustedProb(:,b), '-', 'Color', tempcolors(icolor,:))
            icolor = icolor+1;
        end
    end
    xlabel('relative position (deg)')
    xlim([-90 90])
    ylabel('probability')
    sgtitle(['ripple ', num2str(i), ' of ', num2str(size(rippleMids)), ' pos ', num2str(ripplePos(i)), ' ratio ', num2str(rippleRatio(i))])

end




end

