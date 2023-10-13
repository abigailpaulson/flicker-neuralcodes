function positionInfoSmooth = smoothPosition(positionInfoRaw, recs);
%this function smooths the raw virmen position data 
%LU 10.25.17 modified
%SP 12.19.17 converted to function

%% set variables
speed = []; speedS = [];
step = 0.02; %resampling every 20ms
WinSize = 11; %smooth 5 sampling points of half window size

%% cycle through recordings
for recIdx = 1:length(recs)
    %time = positionInfoRaw(recIdx).timeElapsed; %%%raw sampling time %%
    %this is wrong as of 11.13.18 when updating the rawpos generatino file
    %don't use time elapsed use raw time - SP
    time = positionInfoRaw(recIdx).timeRaw;
    theta = positionInfoRaw(recIdx).thetaRaw;  %%%raw tracing data
    thetaBreaks = positionInfoRaw(recIdx).thetaBreaks; %%%raw time index for running start/end
    
    if ~isempty(time)
        %redefine theta so distances add up over time (add 360 each trial)
        thetaTransfer = theta;
        for i = 2:length(thetaBreaks)
            temp = (thetaBreaks(i)+1):length(theta);
            thetaTransfer(temp)=thetaTransfer(temp)+360;
        end

        time2Num = round((time(end)-time(1))/step); 
        time2=time(1):step:((time2Num)*step+time(1));
        theta2thetaTransfer = interp1(time,thetaTransfer,time2);  %%%%%re-sampling
        theta2thetaTransfer=smooth2005(theta2thetaTransfer,WinSize,'rlowess'); %%%%%%smoothing
        theta2=mod(theta2thetaTransfer,360);
        runID=round((theta2thetaTransfer-theta2)/360);
        runID(runID<0)=0;
        runID=runID+1;
        for i=1:(length(thetaBreaks)-1)
            StartI(i)=min(find(runID==i));   %%%Run starts
            OverI(i)=max(find(runID==i));    %%%Run ends
        end
        thetaBreaks2=[0;OverI(:)];     %%%%%New thetaBreaks

        % calculate the speed of the animal over time
        speed=abs(diff(theta(:)))./diff(time(:));
        speed(thetaBreaks(2:end))=0;
        speed=[0;speed(:)];

        speed2=abs(diff(theta2(:)))./diff(time2(:));
        speed2(thetaBreaks2(2:end))=0;
        speed2 = [speed2(:)]; %speed2=[0;speed2(:)]; %SP 3.13.18 changed so that speed length would match theta and time

        %save the results into the output data structure
        positionInfoSmooth(recIdx).thetaSmooth = theta2(:);
        positionInfoSmooth(recIdx).timeSmooth = time2(:);
        positionInfoSmooth(recIdx).speed = speed(:);
        positionInfoSmooth(recIdx).speedSmooth = speed2(:);
        positionInfoSmooth(recIdx).thetaBreaksSmooth = thetaBreaks2;

        clear time speed time2 speed2 StartI OverI thetaBreaks2 thetaBreaks
    else
        positionInfoSmooth(recIdx).thetaSmooth = [];
        positionInfoSmooth(recIdx).timeSmooth = [];
        positionInfoSmooth(recIdx).speed = [];
        positionInfoSmooth(recIdx).speedSmooth = [];
        positionInfoSmooth(recIdx).thetaBreaksSmooth = [];
    end
end

end
