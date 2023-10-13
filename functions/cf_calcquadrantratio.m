function quadrant_ratio = cf_calcquadrantratio(XVals, YVals, Matrix, ratioType)
%cf_calcquadrantratio
%
%ALP 1/5/2023

%qY is first because index into the matrix as rows (Y axis) and then
%columns (x axis)
for q = 1:4
    quadrant_vals(q) = sum(sum(Matrix(YVals(q,1):YVals(q,2), XVals(q,1):XVals(q,2)), 'omitnan'), 'omitnan');
end

if strcmp(ratioType, 'full')
    quadrant_ratio = ((quadrant_vals(1)+quadrant_vals(3))-(quadrant_vals(2)+quadrant_vals(4)))/sum(quadrant_vals);
elseif strcmp(ratioType, 'mod')
    quadrant_ratio = (quadrant_vals(1)-quadrant_vals(3))/(quadrant_vals(1)+quadrant_vals(3));
end


end

