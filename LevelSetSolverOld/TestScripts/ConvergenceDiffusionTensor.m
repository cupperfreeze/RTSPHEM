% Convergence analysis of effective diffusivity.
% Computes the effective diffusion tensor for multiple refinements of the
% microscopic unit cell and a fixed interface and estimates the order of
% convergence.

numPartitions = 2.^(1:5);
initialLevelSetFunc = @(x, y) 0.24 - norm([x, y], 2);
diffusionTensors = cell(numel(numPartitions), 1);

for i = 1:numel(numPartitions)
    microscaleGrid = FoldedCartesianGrid(2, [-0.5, 0.5, -0.5, 0.5], [5, 5]*numPartitions(i));
    levelSet = arrayfun(initialLevelSetFunc, ...
        microscaleGrid.coordinates(:, 1), microscaleGrid.coordinates(:, 2));

    [A, rhs, isDoF, triangleVolumes] = assembleCellProblem(microscaleGrid, levelSet);
    X = solveSystemFE(microscaleGrid, A, rhs, isDoF);

    diffusionTensors{i} = computeDiffusionTensor(microscaleGrid, X, triangleVolumes);

end

l2Error = zeros(numel(numPartitions)-1, 1);
for i = 1:size(l2Error, 1)

    l2Error(i) = norm(diffusionTensors{i}-diffusionTensors{i + 1});

end

EOC_l2 = zeros(size(l2Error, 1)-1, 1);
for i = 1:size(EOC_l2, 1)

    EOC_l2(i) = log2(l2Error(i)/l2Error(i + 1));

end
