% Prepare Data from multiple input files
% Inserts file name of all inputs to be considered in 'Names' in the
% rwd or tif format

p=genpath('Preprocessing');
addpath(p);
p=genpath('VoxelData');
addpath(p);

MiniCubeSize = [100,100,100];
SlideShare =1;                                                              % Sliding Frame Parameter in share of 'MiniCubeSize'
Names = {'Bentheimer_2d25um_binary_1000x1000x1000.rwd'};
relPath = './VoxelData/SampleBENTHEIMER/';
%Names = {'Berea_2d25um_binary_1000x1000x1000.rwd'};
%relPath = './VoxelData/SampleBEREA/';
%Names = {'CastleGate_2d25um_binary_1000x1000x1000.rwd'};
%relPath = './VoxelData/SampleCASTLE/';
grid = CartesianGrid3D(3,[0,1,0,1,0,1], MiniCubeSize-1);
counter = 0;                                                                % For subsample labeling
locations = cell(numel(Names),1);

for i=1:numel(Names)
    % extract name
    name = Names{i};
    % find char number related to last 'x'
    for j =numel(name):-1:1
        char = name(j);
        if char == 'x'
            break
        end
    end
    % Extract dimension of the cube
    
    
    % read in file according to format
    if name((end-2):end) == 'rwd'
        image = rwd2mat(name);
    else
        error('unknown file format')
    end
    
    % Perform Preprocessing step on data set
    [counter, locationsLocal] = DataPreparation3D(image, counter, size(image), grid, MiniCubeSize, SlideShare, relPath);
    locations{i} = locationsLocal;
end

%save('Locations.mat', 'locations', '-v7.3');
