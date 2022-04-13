%Test accuracy of XFEM scheme for diffusion tensor computation
%Compare to mixed formulation solver

global EPS

dimension = 2;
numPartitions = 40;
cellGrid = FoldedCartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
helpGridHyPHM = Grid(cellGrid.coordinates, cellGrid.triangles);
coord = cellGrid.coordinates;
compareNEW = nan(50, 4);
compareXFEM = nan(50, 4);


for i = 20
    radius = 0.01 * (i);
    initialLevelSetFunc = @(x) radius - norm(x, 2);
    coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
    levelSet = cellfun(initialLevelSetFunc, coordCell);

    %Mixed solution
    [gridNew, levelSetNew] = localMeshRefinement(helpGridHyPHM, ...
        levelSet);
    tic

    diffusionMixed = calcDiffusion(gridNew, levelSetNew, 1);

    toc

    %XFEM solution
    tic
    [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, ...
        triangleSurfaces] = assembleCellProblem(cellGrid, ...
        levelSet);
    SOL = solveSystemFE(cellGrid, cellProblemSystemMatrix, rhs, isDoF);
    [diffusionXFEM, porosity] ...
        = computeDiffusionTensor(cellGrid, SOL, ...
        triangleVolumes);
    toc
    compareXFEM(i, :) = reshape(diffusionXFEM, 1, 4);
end


diffusionMixed
diffusionXFEM
