% Post processing procedure for Chen simulation; Increase concentration
% field of species A in phase 2 by threshhold to make it distinguishable in plots

threshold = 8 * 10^(-16);
for i = 2:length(timeStepsSave)
    temp = AconcentrationSave.getdata(i-1);
    grid = AconcentrationSave.grids{i};
    markBasis = sum(reshape(saveXi(gridHyPHMBasis.V0T(:, :), i), gridHyPHMBasis.numT, 3) == 2, 2) > 0.5;
    markSolid = abs(temp) < eps;
    Interp = scatteredInterpolant(gridHyPHMBasis.baryT(:, 1), gridHyPHMBasis.baryT(:, 2), double(markBasis), 'nearest');
    markRefine = logical(Interp(grid.baryT(:, 1), grid.baryT(:, 2)));
    temp(markRefine & markSolid) = threshold;
    AconcentrationSave.setdataGrid(i-1, temp, grid);
end
