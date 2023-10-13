function cf_figureS3(dirs, params, allindex, metadata)
%chronic flicker prospective coding manuscript
%   supplement figure 3
%ALP 1/18/2023

%% helpful stuff
dayindex = unique(allindex(:,1:2), 'rows');
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Supplement3\'; 

%% create stats file, rewrite it if it already exists
%this is so that everytime the file is run, the stats are regenerated
figname = 'figureS3';
statsfilename = [figname, '_stats_', num2str(datestr(now, 'yymmdd')), '.txt'];
statsFID = fopen([figdir, statsfilename], 'w');
fprintf(statsFID, '%s\r\n', ['Supplement 3 stats -- ' num2str(datestr(now))]);

tablefilename = [figname, '_stats_'];
tablefilename = [figdir, tablefilename];

%% initialize figure
fh = figure('units','inch','position',[0,0,6.5,7.5]);

%% Panel A and B are the images of place cells throughout the track
% 360 deegree version is in A and 180 version is in B

ax_a1 = makesubplotwithletter(fh,5,4, 1, 'A', 'spanW', 1.2);
ax_a2 = makesubplotwithletter(fh,5,4, 2, 'A', 'spanW', 1.2);

cf_plot_placecellmaps(dirs, params, allindex, metadata, fh, {ax_a1, ax_a2})


%% place field properties
ax_b = makesubplotwithletter(fh,5,5, 6, 'B');
ax_c = makesubplotwithletter(fh,5,5, 7, 'C');
ax_d = makesubplotwithletter(fh,5,5, 8, 'D');
ax_e = makesubplotwithletter(fh,5,5, 9, 'E');
ax_f = makesubplotwithletter(fh,5,5, 10, 'F');


cf_plot_placecellproperties(dirs, params, allindex, metadata, fh, {ax_b, ax_c, ax_d, ax_e, ax_f}, statsFID, {'B', 'C', 'D', 'E', 'F'}, tablefilename)


%% put panel here wiht the overall theta sequences between the days?
ax_g = makesubplotwithletter(fh, 5,4,9, 'G');

cf_plot_thetaseq_strength(dirs, params, allindex, metadata, fh, {ax_g}, statsFID, {'G'}, tablefilename)


%% theta sequences 180 per trial
ax_h = makesubplotwithletter(fh,5,4,10,'H', 'spanW', 1.9); %line plot through space
ax_i = makesubplotwithletter(fh,5,4,12,'I');
cf_plot_thetadecoding_results_trial(dirs, params, allindex, metadata, fh, {ax_h, ax_i}, statsFID, {'H', 'I'}, tablefilename)


%% current position decoding 
%decoding current position just 180 for now

ax_j1 = makesubplotwithletter(fh,5,4,13,'J'); %line plot through space
ax_j2 = makesubplotwithletter(fh,5,4,14,'J'); %line plot through space
ax_k = makesubplotwithletter(fh, 5,4,15, 'K');

cf_plot_currPosDecoding(dirs, params, allindex, metadata, fh, {ax_j1, ax_j2, ax_k}, statsFID, {'J', 'J', 'K'}, tablefilename)


%% finish stuff up
%close stats file
fclose(statsFID);

%make figure pretty 
makefigurepretty(fh,1)

%% save
filename = ['figureS3_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')



end

