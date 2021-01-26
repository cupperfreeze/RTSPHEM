function [partialDerivative] = levelSetPartialDerivative(grid, levelSet, direction)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numNodes = grid.nodes;
forwardIndex = grid.getForward(1:numNodes, direction);
backwardIndex = grid.getBackward(1:numNodes, direction);
notNaN = ~isnan(forwardIndex) & ~isnan(backwardIndex);

partialDerivative = NaN(size(levelSet));
partialDerivative(notNaN) = (levelSet(forwardIndex(notNaN)) ...
    -levelSet(backwardIndex(notNaN))) / (2 * grid.stepSize(direction));


end
