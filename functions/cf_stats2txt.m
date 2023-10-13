function cf_stats2txt(data, fid, letter, unit, groupnames)
%cf_stats2txt
%   run stats test and append info to a text file with the name of the
%   figure and the panel name. Add any helpful info you may need!
%ALP 1/18/2023

%%% assuming two groups to test on
for g = 1:2
    t(g) = lillietest(data{g});
end
% 
% if t(1) ~= t(2)
%     error('check normality!')
% end

%%% if both are 0, this is normal, use ttest2
if sum(t) < 1 %
    [~, p] = ttest2(data{1}, data{2});
    test = 't-test';
end

%%% if both are 1, data is not normal, use ranksum test
%also do this if one is normal and the othe isn't
if sum(t) >= 1
    [p] = ranksum(data{1}, data{2});
    test = 'ranksum test';
end

%%% get n for each
for g = 1:2
    N(g) = sum(~isnan(data{g}));
end

%%% get star names! 
starDefs = {'n.s.', '* p<0.05', '** p<0.01', '*** p<0.001'};

sigLevel = 1;
if p < 0.05
    sigLevel = 2;
elseif p<0.01
    sigLevel = 3;
elseif p<0.001
    sigLevel = 4;
end

%%% create string to write to file
str2write = [letter, ': ', starDefs{sigLevel}, ' ', test, '; ', groupnames{1}, ' n = ', num2str(N(1)), ' ', unit, '; ' groupnames{2}, ' n = ', num2str(N(2)), ' ', unit];
fprintf(fid, '%s\r\n', str2write);

end

