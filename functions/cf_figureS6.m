function cf_figureS6(dirs, params, allindex, metadata)
%supplementary figure 6
%ALP 4/16/24

%% helpful stuff
dayindex = unique(allindex(:,1:2), 'rows');
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Supplement6\'; 

%% create stats file, rewrite it if it already exists
%this is so that everytime the file is run, the stats are regenerated
figname = 'figureS6';
statsfilename = [figname, '_captionInfo_', num2str(datestr(now, 'yymmdd')), '.txt'];
statsFID = fopen([figdir, statsfilename], 'w');
fprintf(statsFID, '%s\r\n', ['Figure S6 stats -- ' num2str(datestr(now))]);

tablefilename = [figname, '_stats_'];
tablefilename = [figdir, tablefilename];

%% initialize figure
fh = figure('units','inch','position',[0,0,6.5,6]);
hold on

%% Panel A is a decoding example, I don't have that yet
ax_a = makesubplotwithletter(fh,3,4, 1, 'A', 'spanW', 1.3);
ax_b = makesubplotwithletter(fh,3,4, 3, 'B', 'spanW', 1.3);
ax_c = makesubplotwithletter(fh,3,4,5, 'C', 'spanW', 2);
ax_d = makesubplotwithletter(fh,3,4, 7,'D');


cf_plot_PCR_engaged_unengaged(dirs, params, allindex, metadata, fh, ...
    {ax_a, ax_b, ax_c, ax_d}, statsFID, {'A', 'B', 'C1', 'C2', 'D'}, tablefilename)

%% finish stuff up
%close stats file
fclose(statsFID);

%make figure pretty 
makefigurepretty(fh,1)

%% save the figure as a pdf
filename = ['figureS6_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')


end

