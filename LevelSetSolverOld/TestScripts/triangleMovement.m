% Level - set simulation, contracting triangles

dim = 2;
n = 64; % Number of partitions in each direction
endTime = 0.2;
dt = 0.5 / n;

g = FoldedCartesianGrid(dim, kron(ones(1, dim), [-0.5, 0.5]), n*ones(1, dim));
coord = g.coordinates;

lsf = @(x) max( ...
    min(min(x(1) + 0.25, x(2) + 0.25), 0.25 - norm(x + 0.25, 1)), ...
    min(min(-(x(1) - 0.25), -(x(2) - 0.25)), 0.25 - norm(-(x - 0.25), 1)));
coordCell = mat2cell(coord, ones(1, g.nodes), dim);
initialData = cellfun(lsf, coordCell);

X = g.reshape(coord(:, 1));
Y = g.reshape(coord(:, 2));
plotInitial = g.reshape(initialData);

speed = @(t, x) 0.25; %0 * (  norm( x, Inf ) <= 0.3 );
%variableSpeed = @(t,x) speed(t,x);%* ( ( t <= 0.5*endTime ) - ( t > 0.5*endTime ) );

velocity = [];
% for i = 1:g.numberOfNodes
%     norm1 = norm( coord(i,:) - 0.3 );
%     norm2 = norm( coord(i,:) + 0.3 );
%     if ( norm1 < norm2 )
%         velocity(i,:) = [ min( coord(i,1), 0.3 ) - 0.3, ...
%         min( coord(i,2), 0.3 ) - 0.3 ];
%     else
%         velocity(i,:) = [ max( coord(i,1), -0.3 ) + 0.3, ...
%         max( coord(i,2), -0.3 ) + 0.3 ];
%     end
%     %velocity(i,:) = coord(i,:) - 0.3;
% %     velocity(i,:) = [ max( coord(i,1), -0.3 ) + 0.3, ...
% %         max( coord(i,2), -0.3 ) + 0.3 ];
%     %velocity(i,:) = velocity(i,:) / norm( velocity(i,:) );
% end
% velocity( isnan(velocity) ) = 0;

fig = figure;
%plotGrid( g, fig );
[~, defHandle] = contour(X, Y, plotInitial, zeros(1, 2), 'LineColor', 'r');
set(gca, 'DataAspectRatio', [1, 1, 1]);
xlabel('x');
ylabel('y');
hand = copyobj(defHandle, gca);
set(hand, 'Zdata', plotInitial, 'LineColor', 'b');
%print(fig, '-dpng', 'levelSet0.png');

levelSet = solveLevelSetEquationOld(g, initialData, speed, velocity, ...
    endTime, dt, 'SaveTime', 'all', 'ContourHandle', hand);
