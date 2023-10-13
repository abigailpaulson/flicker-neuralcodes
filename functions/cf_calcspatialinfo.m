function spatialinfo = cf_calcspatialinfo(FRmap, occ)
%cf_calcspatialinfo
%   for map nUnits x nEdges
%ALP 12/19/22

    posPDF = occ;
for u = 1:size(FRmap,1)
    map = FRmap(u,:);
    posPDF=posPDF./nansum(nansum(posPDF));
    meanrate=nansum(nansum(map.*posPDF));
    NormMap=map./meanrate;
    temp=NormMap.*log2(NormMap);
    informationSpike = nansum(nansum(temp.*posPDF));
    spatialinfo(u) = informationSpike;
% 
% Pi = occ./sum(occ);
% meanlambda = sum(Pi.*map);
% 
% for b = 1:length(map)
%     I(b) = Pi(b).*(map(b)./meanlambda).*log2(map(b)./meanlambda);
% end
end

if size(spatialinfo,1) == 1
    spatialinfo = spatialinfo';
end



end

