% Apply morphological model to preprecessed data-samples

p = genpath('./VoxelData/SampleBENTHEIMER');
addpath(p);


imSize = [100, 100, 100];
numSamples = 5000;
grid = CartesianGrid3D(3, [0, 1, 0, 1, 0, 1], imSize-1);
numV = uint32(grid.nodes);
coords = grid.coordinates;
neighborIndices = uint32([grid.pull{1, 1, 1}, grid.pull{2, 1, 1}, grid.pull{1, 1, 2}, grid.pull{2, 1, 2}, grid.pull{1, 1, 3}, grid.pull{2, 1, 3}]);


parfor sample = 0:(numSamples -1)
    InName = strcat('out_connected', num2str(sample), ['_', num2str(imSize(1)), 'x', num2str(imSize(2)), 'x', num2str(imSize(3))]);
    im = rwd2mat(InName);

    for PoroStep = 1:5
        saturation = 1 - 0.1 * PoroStep;
        imNew = MorphModel(1-im, saturation);
        OutName = strcat('out_connected', num2str(sample + PoroStep * numSamples));
        [imNew, ~, isConnected] = fillholes3D(numV, coords, neighborIndices, 1-imNew(:), uint8(1));
        mat2rawrwd(reshape(1 - imNew, imSize), OutName, 1, [0, 0, 0], imSize./imSize (1), '');
    end
end
