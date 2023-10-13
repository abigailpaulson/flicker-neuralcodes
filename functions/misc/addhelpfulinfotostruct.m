function S = addhelpfulinfotostruct(S)
%addhelpfulinfotostruct
%   ALP 12/17/2020

if ~isfield(S, 'saveinfo')
    S.saveinfo.date = datestr(now); 
    st = dbstack;
    S.saveinfo.genFile = {st.file}; 
    S.saveinfo.descript = 'genFile is the file that generated this structure.'; 
end
    
end

