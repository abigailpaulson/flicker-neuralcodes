function newmat = clustermetrics2mat(clustermetrics, cfieldname, csubfieldname)
%clustermetrics2mat take subfields of clustermetrics and output them into a
%matrix, # units x length(thing of interest). if subfield of interest is 
%stable times, output will be a cellarray. 
%
%   INPUTS:
%       clustermetrics - struct
%       cfieldname - field name of thing of interest, ex: 'WF', 'stable'
%       csubfieldname - subfield name of thing of interest, ex: 'mn'
%
%   common options: 'WF', 'mn'; 'WF', 'std'; 'stable', 'meanFR';
%
% ALP 1/15/2020


if ~strcmp('times', csubfieldname)
    subfieldlength = length(clustermetrics(1).(cfieldname).(csubfieldname));
    newmat = NaN(length(clustermetrics), subfieldlength);
    for u = 1:length(clustermetrics)
        newmat(u,:) = clustermetrics(u).(cfieldname).(csubfieldname);
    end
else
    newmat = {NaN(1,length(clustermetrics))};
    for u = 1:length(clustermetrics)
        newmat{u} = clustermetrics(u).(cfieldname).(csubfieldname);
    end
end

end

