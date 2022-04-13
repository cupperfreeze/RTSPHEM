function [phi] = initializeFMMOld(grid, lsf)
% Initialization for compution of signed distance function of the interface
% LSF=0.

% Find all nodes to be initialized.
isInitial = findInterfaceNodes(grid, lsf);


% Initialize all node values to positive or negative infinity. The sign
% indicates, in which subdomain (which are divided via the zero level set)
% the node lies.
phi = sign(lsf) .* Inf(numel(lsf), 1);
phi(isnan(phi)) = 0;

delta = Inf(2, sum(isInitial));
distMin = Inf(grid.dimension, sum(isInitial));
initialIndex = find(isInitial);

for d = 1:grid.dimension
    forwardIndex = grid.getForward(initialIndex, d);
    nonNaNForward = ~isnan(forwardIndex);

    backwardIndex = grid.getBackward(initialIndex, d);
    nonNaNBackward = ~isnan(backwardIndex);

    % Solve the linearized equation
    %   phi(NEIGHBOUR) = phi(NODE) + DELTA * (phi(NEIGHBOUR) - phi(NODE))
    % to obtain an approximate location of the root of the level set
    % function on the line from NODE to NEIGHBOUR.
    % This solution is only valid if 0 <= DELTA <= 1. Otherwise
    delta(1, nonNaNForward) = lsf(initialIndex(nonNaNForward)) ...
        ./ (lsf(initialIndex(nonNaNForward)) ...
        -lsf(forwardIndex(nonNaNForward)));
    delta(2, nonNaNBackward) = lsf(initialIndex(nonNaNBackward)) ...
        ./ (lsf(initialIndex(nonNaNBackward)) ...
        -lsf(backwardIndex(nonNaNBackward)));
    delta(delta < 0 | delta > 1) = Inf;
    delta(isnan(delta)) = 0;

    distMin(d, :) = min(delta) * grid.stepSize(d);
end

% Calculate the approximate distance of the nodes to the interfaces
% (see Adalsteinsson and Sethian, 'The Fast Construction of Extension
% Velocities in Level Set Methods').
if (grid.dimension == uint8(1))
    distSquaredInv = (1 ./ distMin.^2)';
else
    distSquaredInv = sum(1./(distMin.^2))';
end
phi(initialIndex) = sign(lsf(initialIndex)) .* ...
    sqrt(1./distSquaredInv);

% Check whether a node is assigned a value iff it is marked as initial.
assert(all(isinf(phi(~isInitial))) & all(~isinf(phi(isInitial))));

end
