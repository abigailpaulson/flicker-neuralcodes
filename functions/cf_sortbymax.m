function out = cf_sortbymax(map, varargin)
%cf_sortbymax
%
%ALP 12/19/2022

norm2max = 0;
for v = 1:numel(varargin)
    if strcmp(varargin, 'norm2max')
        norm2max = 1;
    end
end

[maxfr,iMax] = max(map, [], 2, 'omitnan');
[~, iSort] = sort(iMax);
out = map(iSort,:);

if norm2max
    out = map(iSort,:)./maxfr(iSort);
end




end

