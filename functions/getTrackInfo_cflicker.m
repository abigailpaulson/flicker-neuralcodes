function trackInfo = getTrackInfo_cflicker(tracktype)
% This function defines the track information for the annular paradigm
% SP 11.9.17
%updated for abby's chronic flicker annular track (bigger, and less track
%paradigms)
%ALP 8/14/19
if strcmp(tracktype, 'chronicflicker_annulartrack')
trackInfo.rewardZone = [54 72; 234 252];
trackInfo.zoneLength = 18; 
trackInfo.rewardZoneOffset = 180; 

trackInfo.nonRewardZone = [144 162; 324 342]; %90 deg from reward zone
trackInfo.rewardPerZone = 3;
end

end