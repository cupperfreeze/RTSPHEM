%Full test case for 4 different phases
close all
%parameter
StepSizeTime = 0.002;
numTimeSteps = 50;

%setup grid
dimension = 2;
numPartitions = 100;
cellGrid = CartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
coord = cellGrid.coordinates;
temp = numPartitions + 1; %temporary variable for plotting

% Initial indicator function for different subdomains
Xi = ones(size(coord, 1), 1);
Xi((coord(:, 2) < 0) & (coord(:, 1) <= 0)) = 2;
Xi((coord(:, 2) < 0) & (coord(:, 1) > 0)) = 3;
Xi((coord(:, 2) >= 0) & (coord(:, 1) > 0)) = 4;

% Velocities at interfaces
% interfaceVelocities(i,j) = velocity between Xi=i and Xi=j
interfaceVelocities = [0, 0.5, -2, 1; 0, 0, 2, 1; 0, 0, 0, 1; 0, 0, 0, 0];
interfaceVelocities = interfaceVelocities - interfaceVelocities'; %add lower half

%Initial configuration Phi
epsilon = 2 / numPartitions;
initialPhiFunc = @(x) max((epsilon-norm([x(1) + eps], 2))*(x(2) < 0)-10*(x(2) >= 0), epsilon-norm([x(2) + eps], 2));
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
Phi = cellfun(initialPhiFunc, coordCell);

%Calculate distance function d^i only up to this value
restrictDist = 3 * epsilon;

for i = 1:numTimeSteps
    %Calculate distance function from each subdomain
    [distFunctions, Xi] = step1General(Xi, Phi, cellGrid, restrictDist);


    %Reconstruct Vonoroi Interface and velocity extension
    [dV, S] = step234General(distFunctions, cellGrid, Xi, interfaceVelocities, restrictDist);

    %Update Phi, but only every 10 steps
    if mod(i, 1) == 0
        Phi = epsilon - dV;
    end

    %Evolve Phi
    Phi = levelSetEquationTimeStep(StepSizeTime, 0, Phi, ...
        cellGrid, S', 1);

end

%plot evolved functions
[distFunctions, Xi] = step1General(Xi, Phi, cellGrid, restrictDist);
[interfaceLength, coordTripel] = evaluateInterfaceTest(cellGrid, Xi, distFunctions, false)
[dV, ~] = step234General(distFunctions, cellGrid, Xi, interfaceVelocities, 4*epsilon);

figure
contour(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Phi, temp, temp), linspace(-epsilon, epsilon, 10));
title('Current Phi Contour');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Phi, temp, temp));
title('Current Phi');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Xi, temp, temp));
title('Current Xi');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(dV, temp, temp));
title('Current dV');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(S, temp, temp));
title('Current S');

% Chart=Xi; delta=0.05;
% for i=1:numel(Xi)
%     x=coord(i,1); y=coord(i,2);
%     dist1 =min( sqrt((x-coord(Xi>1,1)).^2+(y-coord(Xi>1,2)).^2));
%     dist2 =min( sqrt((x-coord(Xi==1,1)).^2+(y-coord(Xi==1,2)).^2));
%     Chart(i) = (dist1>delta) * 0 + (dist2>delta) * 1 + (0<dist2)*(dist2<=delta)*(0.5+dist2/(2*delta)) + (0<dist1)*(dist1<=delta)*(0.5-dist1/(2*delta));
% end
% surf(reshape(Chart,151,151),'EdgeColor','None')
% shading interp
% gca('FontSize',16);
% title('Smooth interpolation')