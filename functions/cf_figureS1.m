function cf_figureS1(dirs, params, allindex, metadata)
%cf_figureS1
%   first supplement figure
%       behavior collapse to d2r
%       recording locations, CA1 and CA3. currently have the probe maps but
%           need a better metric
%ALP 12/13/2022

%% helpful formatting stuff
dayindex = unique(allindex(:,1:2), 'rows');
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Supplement 1\'; 

%% create stats file, rewrite it if it already exists
%this is so that everytime the file is run, the stats are regenerated
figname = 'figureS1';
statsfilename = [figname, '_stats_', num2str(datestr(now, 'yymmdd')), '.txt'];
statsFID = fopen([figdir, statsfilename], 'w');
fprintf(statsFID, '%s\r\n', 'Supplement 1 stats');

tablefilename = [figname, '_stats_'];
tablefilename = [figdir, tablefilename];

%% subplot coordinates
f1 = figure('units','inch','position',[0,0,6.5,6]);

%% behavior example
ax_a = makesubplotwithletter(f1,5,4, 1, 'A', 'spanW', 1.85);
cf_getbehavior_fulltrack(dirs, params, allindex, metadata, f1, ax_a)

%% behavior, by zones and by groups
ax_b = makesubplotwithletter(f1,5,4, 3, 'B', 'spanW', 0.95);
ax_c = makesubplotwithletter(f1,5,4, 4, 'C', 'spanW', 0.95);
ax_d = makesubplotwithletter(f1,5,4, 5, 'D');
ax_e = makesubplotwithletter(f1,5,4, 6, 'E');
ax_f = makesubplotwithletter(f1,5,4, 7, 'F');
ax_g = makesubplotwithletter(f1,5,4, 8, 'G');

cf_plot_behaviorperformance_zones_table(dirs, params, allindex, metadata, ...
    f1, {ax_b, ax_c, ax_d, ax_e, ax_f, ax_g}, statsFID, {'B', 'C', 'D', 'E', 'F', 'G'}, tablefilename)

%% recording location per day
%saving this directly to a pdf bcuz it would be too confusing to have to
%write into the main function
% ax_h = makesubplotwithletter(f1,6,4, 9, 'H', 'nosubplot');
% ax_i = makesubplotwithletter(f1,6,4, 17, 'I', 'nosubplot');

cf_recordinglocation_ripplepower(dayindex, metadata, 'CA1', figdir)
cf_recordinglocation_ripplepower(dayindex, metadata, 'CA3', figdir)

%% final touches
fclose(statsFID); %close the stats file
makefigurepretty(f1,1)

%% save PDF
filename = ['figureS1_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')
end

