function positionInfo = getdistance2reward(positionInfo, trackinfo, recs)
%getdistance2reward add 'distanceToReward' field to positionInfo structure
%   ALP 5/7/2020
%   ALP 12/3/2022 adding comment that this only returns positions with the
%   reward zone at the end of the trak

backedges = [trackinfo.rewardZone(1,2); trackinfo.rewardZone(2,2)]; %back edges of the zone
syms x
pw(x) = piecewise(x<backedges(1), x-backedges(1),...
    ((x>=backedges(1))&(x<backedges(2))), x-backedges(2), x>=backedges(2),...
    x-(backedges(1)+360));

for r = 1:length(recs)
    %adjust theta values
    positionInfo(r).distanceToReward = double(pw(positionInfo(r).thetaSmooth));
    
    %find distance breaks
    switchTheta = find(diff(positionInfo(r).distanceToReward) < -100);
    thetaBreaks = [0; switchTheta; length(positionInfo(r).distanceToReward)];
    positionInfo(r).distanceBreaks = thetaBreaks; 
end



%plot this to check - uncomment if you want to check if
%distance2reward conversion is working
% figure
% hold on
% subplot(1,2,1)
% plot(pw([0:360]))
% xlabel('Position')
% ylabel('Distance to reward')
% title('Piecewise function')
% subplot(1,2,2)
% hold on
% plot(positionInfo(1).thetaSmooth(positionInfo(1).thetaBreaksSmooth(1)+1:positionInfo(1).thetaBreaksSmooth(2)+1),'b')
% plot(positionInfo(1).distanceToReward(positionInfo(1).thetaBreaksSmooth(1)+1:positionInfo(1).thetaBreaksSmooth(2)+1), 'r')
% legend('theta', 'd2reward')
% subplot(1,3,3)
% hold on
% plot(positionInfo(1).thetaSmooth(positionInfo(1).thetaBreaksSmooth(1)+1:positionInfo(1).thetaBreaksSmooth(2)+1), ...
%     positionInfo(1).distanceToReward(positionInfo(1).thetaBreaksSmooth(1)+1:positionInfo(1).thetaBreaksSmooth(2)+1))
% xlabel('Position')
% ylabel('Distance to reward')

end

