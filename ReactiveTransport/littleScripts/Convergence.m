% Data related to Figure 9 and Table 1
% Script to compare simulation results obtained on different meshes
% Calciumtransport is a cell array containing 5x5 instances of the
% HyPHM-class transport related to h=1/25, 1/50, 1/100, 1/200, 1/400 and
% H= 1/15, 1/30, 1/60, 1/120, 1/240

% this might take a while

try
    load('ConvergenceStudies.mat');
catch
    error('unable to load mat file')
end

% Convergence in h for fixed H
for i = 5:5
    for j = 1:4
        diff_h(i, j) = calcDist(CalciumTransport{i, j}.grid, CalciumTransport{i, j}.U.getdata(17), CalciumTransport{i, j + 1}.grid, CalciumTransport{i, j + 1}.U.getdata(17), 0, 0.1);
    end
end


% Convergence in H for fixed h
for i = 1:4
    for j = 5:5
        diff_H(i, j) = calcDist(CalciumTransport{i, j}.grid, CalciumTransport{i, j}.U.getdata(17), CalciumTransport{i + 1, j}.grid, CalciumTransport{i + 1, j}.U.getdata(17), 0, 0.1);
    end
end


%relative Error
for i = 1:5
    for j = 1:5
        diff_Rel(i, j) = calcDist(CalciumTransport{i, j}.grid, CalciumTransport{i, j}.U.getdata(17), CalciumTransport{5, 5}.grid, CalciumTransport{5, 5}.U.getdata(17), 0, 0.1) ...
            / (CalciumTransport{5, 5}.U.getdata(17)' * CalciumTransport{5, 5}.grid.areaT);
    end
end