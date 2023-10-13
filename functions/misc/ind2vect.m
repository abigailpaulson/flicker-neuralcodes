function out = ind2vect(ind, v_length, p_width)
%ind2vect   makes a vector of 0s and 1s from vector ind which indicates
%beginning of 1 runs
%   ind - vect of indices where the 1's should start
%   length - total length of vector in indices
%   width - width of pulse in indices
%ALP 12/7/2022

ind = [ind ind+p_width];
longind = arrayfun(@(x) ind(x,1):ind(x,2), 1:size(ind,1), 'UniformOutput', false);
longrewInd = cell2mat(longind);

vect = zeros(1, v_length);

if any(longrewInd)
    vect(longrewInd) = 1;
end

out = vect;

end

