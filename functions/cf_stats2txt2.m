function cf_stats2txt2(data, fid, letter, nUnit, yUnit, analysisname, groupnames, tablefilename)
%cf_stats2txt
%   run stats test and append info to a text file with the name of the
%   figure and the panel name. Add any helpful info you may need!
%ALP 1/18/2023
%   stats2txt2 gets relevant information about the distribution and saves
%   it to a txt file that can be copied for the figure caption
%   no stats tests 
%ALP 3/21/2023

%%%
if mean(data{1}, 'omitnan') < 0.01 
    prec = '%10.2e\n';
else
    prec = '%.4f';
end

%%% get n for each
for g = 1:2
    N(g) = sum(~isnan(data{g}));
end

%%% get percentiles for each
% min, 25, median, 75, max
for g = 1:2
    P(1,g) = min(data{g}, [], 'omitnan');
    P(2:4,g) = prctile(data{g}, [25 50 75]);
    P(5,g) = max(data{g}, [], 'omitnan');
    percentilesString{g} = [num2str(P(1,g), prec), ', ', num2str(P(2,g), prec), ', ', ...
        num2str(P(3,g), prec), ', ', num2str(P(4,g), prec), ', ', num2str(P(5,g), prec)];
end

MedianValue = P(3,:)';
Quarts = P([2,4],:)';

%%% get mean and standard error
for g = 1:2
    M(g) = mean(data{g}, 'omitnan');
    SEM(g) = std(data{g}, 'omitnan')./sqrt(N(g)); 
end

%%% create string to write to file
str2write = [letter, ': ', groupnames{1}, ', ' num2str(M(1), prec), ' +/- ', num2str(SEM(1), prec), ...
    ' ' , yUnit, '; ', groupnames{2}, ', ', num2str(M(2), prec), ' +/- ', num2str(SEM(2), prec), ' ', yUnit, '; STATS GO HERE; ', ...
    groupnames{1}, ', n = ', num2str(N(1), prec), ' ', nUnit, ', ', ...
    analysisname ' percentiles = ', percentilesString{1}, '; ', ...
    groupnames{2}, ' n = ', num2str(N(2), prec), ' ', nUnit, ', ', ...
    analysisname, ' percentiles = ', percentilesString{2}, '.'];

fprintf(fid, '%s\r\n', str2write);
str2write = [''];
fprintf(fid, '%s\r\n', str2write);

%%% make table of info for copying later
newanalysisname = repmat({analysisname}, 2,1);
UnitVal = repmat({nUnit}, 2,1);
MeanVal = M';
SEMVal = SEM';
NVal = N';
LetterVal = repmat({letter}, 2,1);
groupnames = groupnames';
StatsTable = table(LetterVal, newanalysisname, groupnames, NVal, UnitVal, MeanVal, SEMVal, MedianValue, Quarts);
tablefilename = [tablefilename, '_', letter, '.csv']; 
writetable(StatsTable, tablefilename);

%%% get some quan


end

