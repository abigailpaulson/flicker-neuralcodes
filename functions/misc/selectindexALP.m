function [index, data] = selectindexALP(spreadsheetfile, sheetnumber, headers, processeddatadir, rewritefileinfo, varargin)
% [index] = selectindex(spreadsheetfile, sheetnumber, headers, varargin)
%
%selects indices of files that correspond to selected values in columns of
%spreadsheet.  All options must be satified for each file for the file to be included If no input
%options selected, all indices will be reported.
%
% INPUTS
%   spreadsheetfile:    filename including location, eg
%   sheetnumber:        sheet in spreadsheet file to read in
%   headers:            names of headers for each column.  WIll be determined
%                       automatically if left blank, but if using a Mac or do not have Excel must be specified.
%                       eg '' or {'animal',	'date',	'cell',	'file'}
%   processeddatadir:   directory where processed data is saved, eg
%                       '/Users/asinger/Desktop/AwakeAutopatchData/', to save fileinfo
%   rewritefileinfo:    default 0.  if 1, rewrite fileinfo##.mat even if it
%                       exists.  if file does not exist it will be written
%   cortex:             1 for cortex recordings, 0 for hpc
%   sizeindex:          how many of first few columns are index, eg 3 or 4
%   fad:                1 if fad analysis
%
% SPREADSHEET FORMAT
%   first 4 columns are indices: [animalLetter animal# YYMMDD file#]
%   first row is value name
%   second row is description of value
%
% INPUT Options
%   name of each column in spreadsheet, must match column name exactly
%   'headername', 'equation', eg 'min Ra', '< 40'
%   for column that is being used to select data, blank rows will be
%   excluded
%
% OUTPUT
%   [animal# YYMMDD cell# file#] or [animal# YYMMDD file#] for fad
% updated ALP 4/15/2020 for chronic flicker experiments 

%read in the spreadsheet
filenumcol = 4;
indcols = [2:4,9:10]; %animal, day, file, flicker day, rectype !!! first 3 entries MUST be animal, day, file 

[data, text, rawData]=xlsread(spreadsheetfile,sheetnumber);
if isempty(headers) & ~isempty(text) %ALP 4/15/2020 shouldn't need to specify bc don't run on mac
    headers = rawData(1,:);
elseif isempty(headers) & isempty(text)
    error('you must specify headers because they cannot be read from file')
end

data = rawData(3:end,:); %ALP 4/15/2020 my experiment spreadsheet the first two rows are header information

%find repeated indices and deal asssign them separate file numbers if have
%different start and end times

%--------find indices that match options---------------
if ~isempty(varargin) %if options are specified

    goodrows = zeros(size(data,1), length(varargin)/2); %will indicate which files to include
    col = 1; %column of goodrows
    for option = 1:2:length(varargin)-1 %for each option
        currentoption =  varargin{option};
        optioneqn = varargin{option+1};
        %find column number to match input options
        [tf loc] = ismember(currentoption, headers);
        if tf == 1
            try
                optiondata = cell2mat(data(:, loc)); %relevant data for this option, eg select column of data matrix
                goodrows(:,col) = eval(['optiondata', optioneqn]); %find rows that meet requirements
            catch
                optiondata = data(:, loc); %relevant data for this option, eg select column of data matrix
                for r = 1:size(optiondata,1)
                    goodrows(r,col) = all(eval(['optiondata{r}', optioneqn])); %find rows that meet requirements %only include if all options satisfied
                end
            end
            col = col+1;
        elseif tf == 0
            error([currentoption ' does not correspond to a header'])
        end
    end
    inclrows = all(goodrows,2); %only include if all options satisfied
    index = cell2mat(data(inclrows,indcols));
    identifier = [data{inclrows,1}]'; 

else isempty(varargin) % if no options specified, select all data indices
    index = cell2mat(data(:,indcols));
end


%---------create fileinfo structure------------

%get info for each index
for r = 1:size(data,1)
    anprocesseddatadir = [processeddatadir, num2str(data{r,1}), num2str(data{r,2}), '_', num2str(data{r,3}) , '\'];

    if ~exist(anprocesseddatadir,'dir')
        mkdir(anprocesseddatadir)
    end

    if ~exist([anprocesseddatadir,'fileinfo', num2str(data{r,filenumcol}),'.mat'], 'file') || rewritefileinfo == 1
        fileinfo = [];
        fileinfo{data{r,indcols(1)}}{data{r,indcols(2)}}{data{r,indcols(3)}} = [];
        for h = 1:size(headers,2) %also refers to column number)
            spaces = [];

            value = data{r,h};
            %makes spaces in headernames into _
            fieldname = headers{h};
            spaces = strfind(headers{h}, ' ');
            if ~isempty(spaces)
                fieldname(spaces) = '_';
            end
%             fileinfo = setfield(fileinfo,fieldname,value); %set each header as a field
            fileinfo{data{r,indcols(1)}}{data{r,indcols(2)}}{data{r,indcols(3)}}.(headers{h}) = data{r,h}; 
        end
        fileinfo{data{r,2}}{data{r,3}}{data{r,4}}.index = cell2mat(data(r,indcols));
        

        %save it
        save([anprocesseddatadir,'fileinfo', num2str(data{r,filenumcol}),'.mat'], 'fileinfo')
    end
end
