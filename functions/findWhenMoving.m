function positionInfo = findWhenMoving(recs, positionInfo, speedThreshold)
%this function finds what times in virmen the animal is moving above a
%certain speed threshold
% SP 12.19.17
% updated ALP 4/24/2020 to work with my data structures

for recIdx = 1:length(recs)
    %find spiketimes when animal was above speed threshold
    isMovingIdx = find(positionInfo(recIdx).speedSmooth > speedThreshold);
    timesMoving = positionInfo(recIdx).timeSmooth(isMovingIdx);
    positionInfo(recIdx).whenMoving = timesMoving;
    positionInfo(recIdx).whenMovingIdx = isMovingIdx;

    %find spiketimes when animal was below speed threshold
    isNotMovingIdx = find(positionInfo(recIdx).speedSmooth < speedThreshold);
    timesNotMoving = positionInfo(recIdx).timeSmooth(isNotMovingIdx);
    positionInfo(recIdx).whenNotMoving = timesNotMoving;
    positionInfo(recIdx).whenNotMovingIdx = isNotMovingIdx;
end

end