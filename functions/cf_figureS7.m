function cf_figureS7(dirs, params, allindex, metadata)
%supplementary figure 7
%ALP 4/16/24

%% helpful stuff
dayindex = unique(allindex(:,1:2), 'rows');
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Supplement7\'; 

%% create stats file, rewrite it if it already exists
%this is so that everytime the file is run, the stats are regenerated
figname = 'figureS7';
statsfilename = [figname, '_captionInfo_', num2str(datestr(now, 'yymmdd')), '.txt'];
statsFID = fopen([figdir, statsfilename], 'w');
fprintf(statsFID, '%s\r\n', ['Figure S7 stats -- ' num2str(datestr(now))]);

tablefilename = [figname, '_stats_'];
tablefilename = [figdir, tablefilename];

%% initialize figure
fh = figure('units','inch','position',[0,0,6.5,6]);
hold on

%% Panel A is the proportion of occupied vs unoccupied 
ax_a = makesubplotwithletter(fh,3,4, 1, 'A', 'spanW', 1.9); %need to add this
ax_b = makesubplotwithletter(fh,3,4, 5, 'B', 'spanW', 1.9);
ax_c = makesubplotwithletter(fh,3,4, 7, 'C', 'spanW', 1.9);
ax_d = makesubplotwithletter(fh,3,4, 9, 'D', 'spanW', 1.9);
ax_e = makesubplotwithletter(fh,3,4, 11, 'E', 'spanW', 1.9);


cf_plot_360decoding_unoccpiedzone(dirs, params, allindex, metadata, fh, ...
    {ax_a, ax_b, ax_c, ax_d, ax_e,}, statsFID, {'A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'D1', 'D2', 'E1', 'E2'}, tablefilename)

%% finish stuff up
%close stats file
fclose(statsFID);

%make figure pretty 
makefigurepretty(fh,1)

%% save the figure as a pdf
filename = ['figureS7_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')


end

