function [diffusion, porosity] = computeDiffusionTensor(grid, ...
    cellSolutions, triangleVolumes)

% Compute the homogenized diffusion tensor and porosity.
%
% The homogenized diffusion tensor is computed from solutions \eta_j of cell
% problems via
%       diffusion(i,j) = \int_{Y_l} ( \delta_{i,j} + \partial_i \eta_j ) dy.

assert(isa(grid, 'CartesianGrid'));
assert(abs(size(cellSolutions, 2) - 2) < eps);

diffusion = zeros(2);
%    numTriangles = size( grid.triangles, 1 );
%
%      verticesActive = grid.triangles(  find(triangleVolumes( : )>eps), : );
%     verticesActive = unique(verticesActive(:));
%     assert( all( ~isnan( cellSolutions(verticesActive , 1 ) ) ), 'calculation of diffusion tensor failed' );            %improved SG

%     for tri = 1:numTriangles
%
%         if ( abs( triangleVolumes( tri ) ) < EPS )
%             continue;
%         end
%
% %         if ( abs( triangleVolumes( tri ) - max( triangleVolumes ) ) > EPS )
% %             continue;
% %         end
%
%         diffusion(1,1) = diffusion(1,1) + triangleVolumes( tri );
%         diffusion(2,2) = diffusion(2,2) + triangleVolumes( tri );
%
%         vertices = grid.triangles( tri, : );
%       %  assert( all( ~isnan( cellSolutions( vertices, 1 ) ) ), num2str( tri ) );
%
%         % TODO: Modify edge orientations so that different orientations in one
%         % triangle are allowed (e.g. [1 -1 1]).
%         diffusion(1,1) = diffusion(1,1) + triangleVolumes( tri ) ...
%             * grid.edgeOrientation( tri, 1 ) ...
%             * ( cellSolutions( vertices(2), 1 ) ...
%                 - cellSolutions( vertices(1), 1 ) ) / grid.stepSize(1);
%
%         diffusion(1,2) = diffusion(1,2) + triangleVolumes( tri ) ...
%             * grid.edgeOrientation( tri, 1 ) ...
%             * ( cellSolutions( vertices(2), 2 ) ...
%                 - cellSolutions( vertices(1), 2 ) ) / grid.stepSize(1);
%
%         diffusion(2,1) = diffusion(2,1) + triangleVolumes( tri ) ...
%             * grid.edgeOrientation( tri, 1 ) ...
%             * ( cellSolutions( vertices(3), 1 ) ...
%                 - cellSolutions( vertices(1), 1 ) ) / grid.stepSize(2);
%
%         diffusion(2,2) = diffusion(2,2) + triangleVolumes( tri ) ...
%             * grid.edgeOrientation( tri, 1 ) ...
%             * ( cellSolutions( vertices(3), 2 ) ...
%                 - cellSolutions( vertices(1), 2 ) ) / grid.stepSize(2);
%
%     end
%
% d=diffusion;

indizes = find(abs(triangleVolumes(:)) >= eps);
diffusion(1, 1) = sum(triangleVolumes(indizes));
diffusion(2, 2) = sum(triangleVolumes(indizes));

vertices = grid.triangles(indizes, :);
assert(all(~isnan(cellSolutions(vertices, 1))), 'calculation of diffusion tensor failed');

% TODO: Modify edge orientations so that different orientations in one
% triangle are allowed (e.g. [1 -1 1]).
diffusion(1, 1) = diffusion(1, 1) + sum(triangleVolumes(indizes) ...
    .*grid.edgeOrientation(indizes, 1) ...
    .*(cellSolutions(vertices(:, 2), 1) ...
    -cellSolutions(vertices(:, 1), 1))) / grid.stepSize(1);

diffusion(1, 2) = sum(triangleVolumes(indizes) ...
    .*grid.edgeOrientation(indizes, 1) ...
    .*(cellSolutions(vertices(:, 2), 2) ...
    -cellSolutions(vertices(:, 1), 2))) / grid.stepSize(1);

diffusion(2, 1) = sum(triangleVolumes(indizes) ...
    .*grid.edgeOrientation(indizes, 1) ...
    .*(cellSolutions(vertices(:, 3), 1) ...
    -cellSolutions(vertices(:, 1), 1))) / grid.stepSize(2);

diffusion(2, 2) = diffusion(2, 2) + sum(triangleVolumes(indizes) ...
    .*grid.edgeOrientation(indizes, 1) ...
    .*(cellSolutions(vertices(:, 3), 2) ...
    -cellSolutions(vertices(:, 1), 2))) / grid.stepSize(2);


porosity = sum(triangleVolumes);
%diffusion = diffusion ./ sum( triangleVolumes );

%d-diffusion
end
