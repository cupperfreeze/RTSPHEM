%Do reconstruction steps 1-4 of VIIM

close all
tic
%setup grid
dimension = 2;
numPartitions = ceil(101);
cellGrid = CartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
coord = cellGrid.coordinates;
temp = numPartitions + 1;

% indicator function for different subdomains
% Xi = ones(size(coord,1),1);
% Xi((coord(:,2)<0) & (coord(:,1)<=0)) = 2;
% Xi((coord(:,2)<0) & (coord(:,1)>0)) = 3;
Xi = ones(size(coord, 1), 1);
Xi((coord(:, 1).^2+coord(:, 2).^2 < 0.16) & (coord(:, 1) <= 0)) = 2;
Xi((coord(:, 1).^2+coord(:, 2).^2 < 0.16) & (coord(:, 1) >= 0)) = 3;

%Initial configuration Phi
epsilon = 0.05;
%initialPhiFunc = @(x) max((epsilon - norm([x(1)],2))*(x(2)<0) -10*(x(2)>=0),epsilon - norm([x(2)],2));
initialPhiFunc = @(x) max(min(-(0.4 - epsilon - norm(x, 2)), 0.4 + epsilon - norm(x, 2)), (epsilon - norm(x(1) - 10 * eps, inf))*(abs(x(2)) <= 0.4)-100*(abs(x(2)) > 0.4));
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
Phi = cellfun(initialPhiFunc, coordCell);

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Phi, temp, temp), 'EdgeColor', 'none');
set(gca, 'FontSize', 18);
title('Phi');

%Calculate distance function from each subdomain
restrictDist = inf;
[distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);

for i = 1:3
    figure
    surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(distFunctions(:, i), temp, temp), 'EdgeColor', 'none');
    set(gca, 'FontSize', 18);
    %  title(['Distances d^' num2str(i)]);
end

%Reconstruct Vonoroi Interface and velocity extension
interfaceVelocities = [0, 2, 3; 0, 0, 0; 0, 0, 0]; % velocity of interface between phases i<j
interfaceVelocities = interfaceVelocities - interfaceVelocities';
[dV, S] = step234(distFunctions, cellGrid, Xi, interfaceVelocities, restrictDist);

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(dV, temp, temp), 'EdgeColor', 'none');
set(gca, 'FontSize', 18);
title('Distance dV');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(S, temp, temp), 'EdgeColor', 'none');
set(gca, 'FontSize', 18);
title('Normal Velocity field S');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Xi, temp, temp), 'EdgeColor', 'none');
set(gca, 'FontSize', 18);
title('Indicator Xi');

toc
