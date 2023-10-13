function cf_getplacecells(dirs, params, allindex, metadata)
%cf_getplacecells
%   based on get placefieldsgroup.m and getplacefieldsday.m
%ALP 12/16/2022

%%% params
dayindex = unique(allindex(:,1:2), 'rows');
params.speedThreshold = 1; %deg/s
params.SI_prctile = 95;
params.minmeanFR = 0.2;
params.maxmeanFR = 10;
params.minpeakFR = 1;

if strcmp(params.placecells_postype, 'full')
    params.place_edges = 0:2:360;
else
    params.place_edges = -81:3:99;
end

%%% loop over days
for d = 1:size(dayindex,1)
    index = dayindex(d,:);
    files = allindex(allindex(:,2) == dayindex(d,2),3);
    
    if strcmp(params.condition, 'all')
        recType = metadata.RecordingType{d};
        inclFiles = recType ~= 2;
        files = files(inclFiles);
    end
    
    cf_placecellclassifier(dirs, params, index, files)
end



end

