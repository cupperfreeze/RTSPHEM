function levelSetOld = levelSetEquationTimeStep(t, tOld, levelSetOld, ...
    grid, normalSpeed, varargin)

%Perform a single explicit Euler step of the level-set equation for a given
%normal velocity field. If normalSpeed is a scalar the velocity field is
%considered constant.
%Default order of the upwind discretization of the gradient
%is 2 but can be reduced by varargin.


%     assert( all( levelSet >= levelSetOld - EPS ) );
if numel(normalSpeed) == 1;
    normalSpeed = normalSpeed * ones(size(levelSetOld));
end

order = 2; %default order
if length(varargin) == 1
    order = varargin{1};
    if order > 2
        print('order not supported');
        return
    end
end


dt = t - tOld;

dim = grid.dimension;
h = grid.stepSize;

Notinf = ~isinf(levelSetOld);
index = find(Notinf); %indizes for banded- representation

numNodes = numel(index);
%    numNodes = grid.nodes;

%     coordCell = mat2cell( grid.coordinates, ones(1, numNodes), dim );
%     timeCell = num2cell( t * ones(numNodes,1) );
%     normalSpeed = cellfun( speed, timeCell, coordCell )

deltaPlusSq = zeros(numNodes, 1); %Initialize arrays
deltaMinusSq = zeros(numNodes, 1);
%     advectiveMovement = zeros( numNodes, 1 );
forwardIndex = zeros(numNodes, dim);
forwardDiff = zeros(numNodes, dim);
backwardIndex = zeros(numNodes, dim);
backwardDiff = zeros(numNodes, dim);

forwardIndex2 = zeros(numNodes, dim);
forwardDiff2 = zeros(numNodes, dim);
forwardSwitch = zeros(numNodes, dim);
forwardSwitch2 = zeros(numNodes, dim);
backwardIndex2 = zeros(numNodes, dim);
backwardDiff2 = zeros(numNodes, dim);
backwardSwitch = zeros(numNodes, dim);
backwardSwitch2 = zeros(numNodes, dim);

for d = 1:dim
    % Assemble term for movement in normal direction
    %forwardIndex(:,d) = grid.getForward( 1:numNodes, d );
    %forwardIndex2(:,d) = grid.getForward( forwardIndex(:,d), d );
    forwardIndex(:, d) = grid.pull{1, 1, d}(index); %use indexes from grid property
    notNaN = (~isnan(forwardIndex(:, d)));


    forwardDiff(notNaN, d) = (levelSetOld(forwardIndex(notNaN, d)) ...
        -levelSetOld(index(notNaN))) ./ h(d);

    %calculate second order correction
    if order == 2
        forwardIndex2(:, d) = grid.pull{1, 2, d}(index);
        notNaN2 = ~isnan(forwardIndex2(:, d));
        forwardDiff2(notNaN2, d) = (+levelSetOld(forwardIndex2(notNaN2, d)) ...
            -2 * levelSetOld(forwardIndex(notNaN2, d)) + levelSetOld(index(notNaN2))) ./ h(d) ./ 2;

        forwardSwitch(notNaN2, d) = (order - 1) * (levelSetOld(forwardIndex2(notNaN2, d)) <= levelSetOld(forwardIndex(notNaN2, d)));
        forwardSwitch2(notNaN2, d) = (order - 1) * (levelSetOld(forwardIndex2(notNaN2, d)) >= levelSetOld(forwardIndex(notNaN2, d)));
    end

    %backwardIndex(:,d) = grid.getBackward( 1:numNodes, d );
    %backwardIndex2(:,d) = grid.getBackward( backwardIndex(:,d), d );


    backwardIndex(:, d) = grid.pull{2, 1, d}(index);
    notNaN = (~isnan(backwardIndex(:, d)));
    backwardDiff(notNaN, d) = (levelSetOld(index(notNaN)) ...
        -levelSetOld(backwardIndex(notNaN, d))) ./ h(d);

    %calculate second order correction
    if order == 2
        backwardIndex2(:, d) = grid.pull{2, 2, d}(index);
        notNaN2 = ~isnan(backwardIndex2(:, d));

        backwardDiff2(notNaN2, d) = (levelSetOld(index(notNaN2)) ...
            -2 * levelSetOld(backwardIndex(notNaN2, d)) + levelSetOld(backwardIndex2(notNaN2, d))) ./ h(d) ./ 2;


        backwardSwitch(notNaN2, d) = (order - 1) * (levelSetOld(backwardIndex2(notNaN2, d)) <= levelSetOld(backwardIndex(notNaN2, d)));
        backwardSwitch2(notNaN2, d) = (order - 1) * (levelSetOld(backwardIndex2(notNaN2, d)) >= levelSetOld(backwardIndex(notNaN2, d)));
    end

    %Apply upwind scheme
    deltaPlusSq = deltaPlusSq + ...
        max([backwardDiff(:, d) + backwardSwitch(:, d) .* backwardDiff2(:, d), ...
        -(forwardDiff(:, d) - forwardSwitch(:, d) .* forwardDiff2(:, d)), zeros(numNodes, 1)], ...
        [], 2).^2;

    deltaMinusSq = deltaMinusSq + ...
        max([-backwardDiff(:, d) - backwardSwitch2(:, d) .* backwardDiff2(:, d), ...
        +(forwardDiff(:, d) - forwardSwitch2(:, d) .* forwardDiff2(:, d)), zeros(numNodes, 1)], ...
        [], 2).^2;


end

levelSetOld(index) = levelSetOld(index) - dt * ( ...
    max(normalSpeed(index), zeros(numNodes, 1)) .* sqrt(deltaPlusSq) + ...
    min(normalSpeed(index), zeros(numNodes, 1)) .* sqrt(deltaMinusSq));

end
