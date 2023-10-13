function cf_figureS5(dirs, params, allindex, metadata)
%Figure S5
%ALP 6/30/23

%% helpful stuff
dayindex = unique(allindex(:,1:2), 'rows');
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Supplement5\'; 

%% create stats file, rewrite it if it already exists
%this is so that everytime the file is run, the stats are regenerated
figname = 'figureS5';
statsfilename = [figname, '_captionInfo_', num2str(datestr(now, 'yymmdd')), '.txt'];
statsFID = fopen([figdir, statsfilename], 'w');
fprintf(statsFID, '%s\r\n', ['Figure S5 stats -- ' num2str(datestr(now))]);

tablefilename = [figname, '_stats_'];
tablefilename = [figdir, tablefilename];

%% initialize figure
fh = figure('units','inch','position',[0,0,6.5,6]);
hold on

%% Panel B, C, D are from the theta seq decoding figure
ax_a = makesubplotwithletter(fh,3,2, 1, 'A');
ax_b = makesubplotwithletter(fh,3,2, 2, 'B');

cf_plot_prospectivecoding_byspeed_supplement(dirs, params, allindex, metadata, fh, {ax_a, ax_b}, statsFID, {'A', 'B'}, tablefilename)

%% finish stuff up
%close stats file
fclose(statsFID);

%make figure pretty 
makefigurepretty(fh,1)

%% save the figure as a pdf
filename = ['figureS5_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')


end

