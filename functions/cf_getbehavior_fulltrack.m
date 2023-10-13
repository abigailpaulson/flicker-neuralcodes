function cf_getbehavior_fulltrack(dirs, params, allindex, metadata, fh, ax)
%cf_getbehavior_fulltrack
%   simply get the velocity and licking across the entire 360 track
%   probably avg and stde across all recording days
%
%ALP 12/13/2022

%%% params and track info
params.pos_edges = 0:6:360;
trackInfo = getTrackInfo_cflicker('chronicflicker_annulartrack');
RZ = trackInfo.rewardZone;

%%% intitialize, directories, etc.
positiondir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\fig_Behavior\';
dayindex = unique(allindex(:,1:2), 'rows');
flickerDay = metadata.FlickerDay;
isPostDay = flickerDay > 5;

dayindex = dayindex(isPostDay,:);

%%%%% ----- calculate per trial behavior metrics ----- %%%%%
dayLick = []; dayVel = []; dayLick_n = []; dayVel_n = [];

for d = 1:size(dayindex,1)
    trialfilename = ['trialInfo_360_A',  num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    load([positiondir, trialfilename])
    
    anLick = []; anVel = []; anLick_n = []; anVel_n = [];
    %%% analyze trials
    for t = 1:length(trialData)
        if ~trialData(t).fullTrial || ~trialData(t).engaged || ~trialData(t).rewarded
            continue
        end
        time = trialData(t).time;
        licks = logical(trialData(t).licks);
        theta = trialData(t).theta;
        rewards = logical(trialData(t).rewards);
        speed = trialData(t).speed;
        lickpos = theta(licks);
        rewardpos = theta(rewards);
        
        lick_h = histcounts(lickpos, params.pos_edges); 
        [~, ~, theta_bin] = histcounts(theta, params.pos_edges); 
        vel_h = arrayfun(@(x) mean(speed(theta_bin == x), 'omitnan'), 1:length(params.pos_edges)-1, 'UniformOutput', false); 
        vel_h = cell2mat(vel_h); 
        
        lick_nh = lick_h./sum(lick_h);
        vel_nh = vel_h./max(vel_h);
        
        anLick = [anLick; lick_h];
        anVel = [anVel; vel_h];
        anLick_n = [anLick_n; lick_nh];
        anVel_n = [anVel_n; vel_nh];
    end
    
    dayLick = [dayLick; nanmean(anLick,1)];
    dayLick_n = [dayLick_n; nanmean(anLick_n,1)];
    dayVel = [dayVel; nanmean(anVel,1)];
    dayVel_n = [dayVel_n; nanmean(anVel_n,1)]; 
end

%%% averages
mnLick = nanmean(dayLick_n,1);
mnVel = nanmean(dayVel,1);
stdeLick = nanstd(dayLick_n,0,1)./sqrt(sum(~isnan(dayLick_n(:,1))));
stdeVel = nanstd(dayVel,0,1)./sqrt(sum(~isnan(dayVel(:,1))));

%%%% plot with patches for the reward zones
figure(fh)
axes(ax)
hold on
% patch([0 0 RZ(1,2) RZ(1,2)], [1.5 1.7 1.7 1.5], params.colors.behavior.half(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none')
% patch([RZ(1,2) RZ(1,2) RZ(1,2)+180 RZ(1,2)+180], [1.5 1.7 1.7 1.5], params.colors.behavior.half(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none')
% patch([RZ(2,2) RZ(2,2) 360 360], [1.5 1.7 1.7 1.5], params.colors.behavior.half(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch([RZ(1,1) RZ(1,2) RZ(1,2) RZ(1,1)], [0 0 1.4 1.4], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
patch([RZ(2,1) RZ(2,2) RZ(2,2) RZ(2,1)], [0 0 1.4 1.4], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
yyaxis left
hold on
ciplot(mnLick-stdeLick, mnLick+stdeLick, params.pos_edges(1:end-1)+3, 'k')
plot(params.pos_edges(1:end-1)+3, mnLick, 'k-')
ylim([0 0.2])
xticks([0 180 360])
ylabel('fraction of licks')
% shadedErrorBar(params.pos_edges(1:end-1)+3, mnLick, stdeLick, 'k-', 1) %scaled for plotting

yyaxis right
hold on
% plot(params.pos_edges(1:end-1)+3, mnVel)
% shadedErrorBar(params.pos_edges(1:end-1)+3, mnVel+0.6, stdeVel, 'm-', 1)
ciplot(mnVel-stdeVel, mnVel+stdeVel, params.pos_edges(1:end-1)+3, 'm')
plot(params.pos_edges(1:end-1)+3, mnVel, 'm-')
% patch([0 18 18 0], [0 0 1.4 1.4], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
xlabel('position (deg)')
xlim([0 360])
ylabel('speed (deg/s)')
% ylim([0 0.25])

% figure
% hold on
% patch([RZ(1,1) RZ(1,2) RZ(1,2) RZ(1,1)], [0 0 0.5 0.5], 'm', 'FaceAlpha', 0.05, 'EdgeColor', 'none')
% patch([RZ(2,1) RZ(2,2) RZ(2,2) RZ(2,1)], [0 0 0.5 0.5], 'm', 'FaceAlpha', 0.05, 'EdgeColor', 'none')
% plot(params.pos_edges(1:end-1)+3, dayLick_n, '-', 'Color', '#d9d9d9')
% plot(params.pos_edges(1:end-1)+3, mnLick, 'k-', 'LineWidth', 1.5)
% xlabel('position (deg)')
% xlim([0 360])

end

