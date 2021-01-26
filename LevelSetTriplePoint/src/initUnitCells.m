% initialize Unit Cells
% reduceData -option only saves last geometry state
function out = initUnitCells(cellGrid, numberOfSlices, numTimeSlices, reduceData)

coord = cellGrid.coordinates;
epsilon = sqrt(2) * cellGrid.stepSize(1);
out.epsilon = epsilon;
out.reInit = ones(numberOfSlices, 1);
out.cellGrid = cellGrid;
out.reduceData = reduceData;
out.memoryStep = 1;
out.adapTime = cell(numberOfSlices, 1);
[out.adapTime{:}] = deal(0);
out.adapS = cell(numberOfSlices, 1);
[out.adapS{:}] = deal(zeros(3, 3));
out.saved = cell(0, 1);

if reduceData %compressed data format
    numTimeSlices = 1;
end

radius = 0.2839;
% Initial indicator function for different subdomains
Xi = ones(size(coord, 1), numTimeSlices);
Xi((coord(:, 1) > 0) & (sqrt(coord(:, 1).^2 + coord(:, 2).^2) <= radius), 1) = 2;
Xi((coord(:, 1) <= 0) & (sqrt(coord(:, 1).^2 + coord(:, 2).^2) <= radius), 1) = 3;


%Initial configuration Phi
Phi = nan(size(coord, 1), numTimeSlices);
signed = -nan(cellGrid.nodes, numTimeSlices);
initialPhiFunc = @(x) max(min(-(radius - epsilon - norm(x, 2)), radius + epsilon - norm(x, 2)), (epsilon - norm(x(1) - 10 * eps, inf))*(abs(x(2)) <= radius)-100*(abs(x(2)) > radius));
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), 2);
Phi(:, 1) = cellfun(initialPhiFunc, coordCell);
signed(:, 1) = cellfun(@(x) radius-norm(x), coordCell);

%Calculate distance function d^i only up to this value
out.restrictDist = 3 * epsilon;

%Calibrate Xi to Voronoi Interface
[distFunctions, Xi(:, 1)] = step1(Xi(:, 1), Phi(:, 1), cellGrid, out.restrictDist);

out.Xi = cell(numberOfSlices, 1);
out.Phi = cell(numberOfSlices, 1);
out.signed = cell(numberOfSlices, 1);
out.interfaceLength = cell(numberOfSlices, 1);
out.distFunctions = cell(numberOfSlices, 1); %For internal use only

[out.signed{:}] = deal(signed);
[out.Xi{:}] = deal(Xi);
[out.Phi{:}] = deal(Phi);
[out.distFunctions{:}] = deal(distFunctions);

interfaceLength = nan(numTimeSlices, 2);
interfaceLength(1, :) = radius * pi;
[out.interfaceLength{:}] = deal(interfaceLength);

end