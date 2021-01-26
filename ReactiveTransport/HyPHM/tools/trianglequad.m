%> @file trianglequad.m Construct quadrature points and weights.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Construct quadrature points P and weights W on a triangle specified by
%> its vertices verts.  The quadrature rules used here can be found in
%> @ref EG2000, p. 360.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @retval P      Quadrature points              [#points x 2]
%> @retval W      Quadrature Weights             [#points x 2]
%> @param  coordV Vertices of triangle           [3 x 2]
%> @param  type   Type of quadrature (see table) [string]
%> @param  isVisualize (optional) Visualization  [bool]

function [P, W] = trianglequad(coordV, type, isVisualize)

assert(isequal(size(coordV), [3, 2]))
assert(isa(type, 'char'))

V = mat2cell(coordV, [1, 1, 1], 2); % vertices [1 x 2]
S = abs(det([coordV, ones(3, 1)])) / 2; % surface area

% a row of M has the following structure:
% [bc1 bc2 bc3 weight], bck the kth barycentric coordinate
switch type
    case '11'
        M = [1 / 3, 1 / 3, 1 / 3, S];
    case '23'
        M = [0, 1 / 2, 1 / 2, 1 / 3 * S; ...
            1 / 2, 0, 1 / 2, 1 / 3 * S; ...
            1 / 2, 1 / 2, 0, 1 / 3 * S];
    case '34'
        M = [1 / 3, 1 / 3, 1 / 3, -9 / 16 * S; ...
            1 / 5, 1 / 5, 3 / 5, 25 / 48 * S; ...
            1 / 5, 3 / 5, 1 / 5, 25 / 48 * S; ...
            3 / 5, 1 / 5, 1 / 5, 25 / 48 * S];
    case '37'
        M = [1 / 3, 1 / 3, 1 / 3, 9 / 20 * S; ...
            1 / 2, 1 / 2, 0, 2 / 15 * S; ...
            1 / 2, 0, 1 / 2, 2 / 15 * S; ...
            0, 1 / 2, 1 / 2, 2 / 15 * S; ...
            1, 0, 0, 1 / 20 * S; ...
            0, 1, 0, 1 / 20 * S; ...
            0, 0, 1, 1 / 20 * S];
    case '46'
        a1 = 0.445948490915965;
        a2 = 0.091576213509771;
        w1 = 0.223381589678010;
        w2 = 0.109951743655322;
        M = [a1, a1, 1 - 2 * a1, w1 * S; ...
            a1, 1 - 2 * a1, a1, w1 * S; ...
            1 - 2 * a1, a1, a1, w1 * S; ...
            a2, a2, 1 - 2 * a2, w2 * S; ...
            a2, 1 - 2 * a2, a2, w2 * S; ...
            1 - 2 * a2, a2, a2, w2 * S];
    case '57'
        a1 = 0.101286507323456; % (6-sqrt(15))/21
        a2 = 0.470142064105115; % (6+sqrt(15))/21
        w1 = 0.125939180544827; % (155-sqrt(15))/1200
        w2 = 0.132394152788506; % (155+sqrt(15))/1200
        M = [1 / 3, 1 / 3, 1 / 3, 9 / 40 * S; ...
            a1, a1, 1 - 2 * a1, w1 * S; ...
            a1, 1 - 2 * a1, a1, w1 * S; ...
            1 - 2 * a1, a1, a1, w1 * S; ...
            a2, a2, 1 - 2 * a2, w2 * S; ...
            a2, 1 - 2 * a2, a2, w2 * S; ...
            1 - 2 * a2, a2, a2, w2 * S];
    case '612'
        a1 = 0.063089014491502;
        a2 = 0.249286745170910;
        a = 0.310352451033785;
        b = 0.053145049844816;
        c = 1 - a - b;
        w1 = 0.050844906370206;
        w2 = 0.116786275726378;
        w = 0.082851075618374;
        M = [a1, a1, 1 - 2 * a1, w1 * S; ...
            a1, 1 - 2 * a1, a1, w1 * S; ...
            1 - 2 * a1, a1, a1, w1 * S; ...
            a2, a2, 1 - 2 * a2, w2 * S; ...
            a2, 1 - 2 * a2, a2, w2 * S; ...
            1 - 2 * a2, a2, a2, w2 * S; ...
            a, b, c, w * S; ...
            b, a, c, w * S; ...
            b, c, a, w * S; ...
            a, c, b, w * S; ...
            c, a, b, w * S; ...
            c, b, a, w * S];
    otherwise
        error('HyPHM: Quadrature rule not known.  See documentation.')
end

% QUADRATURE WEIGHTS
W = M(:, 4);
% check weights
assert(abs(sum(W) / S - 1) < 1E-14, 'HyPHM: Weights are no partition of unity [sum/area = %f].  Check quadrature formula.', sum(W)/S);

% QUADRATUR POINTS
dim = length(W);
P = zeros(dim, 2);
for k = 1:dim
    P(k, :) = M(k, 1) * V{1} + M(k, 2) * V{2} + M(k, 3) * V{3};
end

% VISUALIZATION (optional)
if exist('isVisualize', 'var') && isVisualize == true
    line([coordV(:, 1); coordV(1, 1)], [coordV(:, 2); coordV(1, 2)], 'Color', 'r', 'LineWidth', 2)
    for k = 1:dim
        text(P(k, 1), P(k, 2), int2str(k))
        line([P(k, 1), P(k, 1)], [P(k, 2), P(k, 2)], [0, W(k)], 'Color', 'k', 'LineWidth', 2)
    end
end


end
