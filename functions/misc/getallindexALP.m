function [allindex, data] = getallindexALP(processeddatadir, spreadsheetdir, varargin)
%
% INPUTS
%   processeddatadir - location of processeddata to save file info
%   spreadsheet dir - location of the experiment spreadsheet
%   varagin - optional inputs. I most often use 'stimulation'
%
% OUTPUT
%   allindex - list of included indices [ANIMAL# YYMMDD FILE#]
%   identifier - list of animal identifier letters as strings
%
% OPTIONS
%   'rewritefileinfo' - 1 to rewrite fileinfo, 0 not to.  default is 0.
%   personalized experiment options to reference which sheet number you
%   want the function to go into, or which stimulation type you want
%
% ALP June 2018
% updated ALP 4/15/2020 

allindex=[];
identifier = [];

rewritefileinfo = 0; abovebrain = 0; ca1 = 0;%set variable options

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'rewritefileinfo'
            rewritefileinfo = varargin{option+1};
        case 'stimulation'
            stimulation = varargin{option+1};
        case 'abovebrain'
            abovebrain = varargin{option+1};

        case 'ca1'
            ca1 = varargin{option+1}; 
    end
end

%Set the stuff you want to be specific to your experiment flags...
headers = {};
spreadsheetfile = spreadsheetdir; %location & name of current spreadsheet
sheetnumber = 'CA1';

if ~exist('stimulation','var')
    [allindex, data] = selectindexALP(spreadsheetfile, sheetnumber, ...
        headers, processeddatadir, rewritefileinfo, 'Include', '~=0') ;
elseif exist('stimulation','var')
    %changed from ~isempty and isempty ALP 02/03/17
    stimeqn=[];
    stimeqn = ['== ''', stimulation, ''''];
    [allindex, data] = selectindexALP(spreadsheetfile, sheetnumber, ...
        headers, processeddatadir, rewritefileinfo, 'Include', '~=0', 'Stimulation', stimeqn);
end

end


