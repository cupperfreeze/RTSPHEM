% Interpolate precomputed permeability data using stepsize Xlength to the
% finer stepsize fine

function [InterpolX, InterpolY] = createLookUp(Xlength, fine)

% numPartitions = 50;
% PermsX = zeros(length(Xlength));
% PermsY = zeros(length(Xlength));
%
% for i=1:length(Xlength)
%     i
%     for j=1:i
%
%
%         dimension = 2;
%         cellGrid = FoldedCartesianGrid( dimension, ...
%         kron( ones(1,dimension), [-0.5, 0.5] ), ...
%         numPartitions*ones(1,dimension) );
%         coord = cellGrid.coordinates;
%         helpGridHyPHM = Grid(cellGrid.coordinates, cellGrid.triangles);
%
%
%         initialLevelSetFunc = @(x) Xlength(i)-norm([x(1),  ...
%                                           Xlength(i)/Xlength(j)*x(2)],inf);
%         coordCell = mat2cell( coord, ones(1,cellGrid.nodes), dimension );
%         levelSet = cellfun( initialLevelSetFunc, coordCell );
%       %  surf(reshape(levelSet,51,51))
%
%             temp = computePermeabilityTensor(helpGridHyPHM,levelSet);
%             PermsX(i,j) = temp(1);
%             PermsY(i,j) = temp(4);
%             PermsX(j,i) = temp(4);
%             PermsY(j,i) = temp(1);
%
%     end
%
% end

% load precomputed data obtained by the precedure above
load('PermLookupRect.mat');
[temp1, temp2] = meshgrid(Xlength, Xlength);
[temp3, temp4] = meshgrid(fine);
InterpolX = interp2(temp1, temp2, PermsX, temp3, temp4);
InterpolY = interp2(temp1, temp2, PermsY, temp3, temp4);
end