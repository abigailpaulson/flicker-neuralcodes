function [ norm_sorted_rate] = getnormandsortedrate(rate)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
        [maxR, iMax] = max(rate, [], 2);
        [~, iSort] = sort(iMax);
        norm_rate = rate./maxR;
        norm_sorted_rate = norm_rate(iSort,:);
end

