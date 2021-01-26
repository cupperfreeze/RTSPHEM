% initialize Unit Cells
% reduceData -option only saves last geometry state
% outX, outY contain side length of underlying rectangles
function out = initUnitCellsDataSet(cellGrid, numberOfSlicesX, numberOfSlicesY, numTimeSlices, reduceData, outX, outY)

numberOfSlices = numberOfSlicesX * numberOfSlicesY;
coord = cellGrid.coordinates;
epsilon = 2 * cellGrid.stepSize(1);
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

out.Xi = cell(numberOfSlices, 1);
out.Phi = cell(numberOfSlices, 1);
out.signed = cell(numberOfSlices, 1);
out.interfaceLength = cell(numberOfSlices, 1);

for i = 1:numberOfSlicesX
    for j = 1:numberOfSlicesY
        % Initial indicator function for different subdomains
        Xi = ones(size(coord, 1), numTimeSlices);
        Xi((coord(:, 1) > 0) & (abs(coord(:, 2)) < outY(i, j)) & (coord(:, 1) < outX(i, j)), 1) = 2;
        Xi((coord(:, 1) <= 0) & (abs(coord(:, 2)) < outY(i, j)) & (coord(:, 1) > -outX(i, j)), 1) = 3;

        %Initial configuration Phi
        Phi = nan(size(coord, 1), numTimeSlices);
        signed = -nan(cellGrid.nodes, numTimeSlices);

        Phi(:, 1) = max(min([epsilon + outX(i, j) - max(abs([coord(:, 1), (outX(i, j) + epsilon) / (outY(i, j) + epsilon) * coord(:, 2)]), [], 2), ...
            -(-epsilon + outX(i, j) - max(abs([coord(:, 1), (outX(i, j) - epsilon) / (outY(i, j) - epsilon) * coord(:, 2)]), [], 2))], [], 2), ...
            (epsilon - abs(coord(:, 1) - 10 * eps)).*(abs(coord(:, 2)) <= outX(i, j))-100*(abs(coord(:, 2)) > outX(i, j)));
        signed(:, 1) = outX(i, j) - max(abs([coord(:, 1), outX(i, j) / outY(i, j) * coord(:, 2)]), [], 2);

        %Calculate distance function d^i only up to this value
        out.restrictDist = 3 * epsilon;


        out.signed{(j-1)*numberOfSlicesX+i} = deal(single(signed));
        out.Xi{(j-1)*numberOfSlicesX+i} = deal(uint8(Xi));
        out.Phi{(j-1)*numberOfSlicesX+i} = deal(single(Phi));

        interfaceLength = nan(numTimeSlices, 2);
        interfaceLength(1, :) = 2 * (outX(i, j) + outY(i, j));
        out.interfaceLength{(j-1)*numberOfSlicesX+i} = interfaceLength;
    end
end

localXi = cell(numberOfSlices, 1);
localdistFunctions = cell(numberOfSlices, 1);

parfor i = 1:numberOfSlices
    [localdistFunctions{i}, localXi{i}] = step1(out.Xi{i}(:, 1), out.Phi{i}(:, 1), cellGrid, out.restrictDist);
end
out.Xi = localXi;
out.distFunctions = localdistFunctions;

end
