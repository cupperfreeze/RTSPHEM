%> @file RT0toP1P1slice.m Mapping @f$\mathbb{RT}_0(\mathcal{T}) \rightarrow \mathbb{P}_1(\mathcal{T})^2@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Evaluates the flux for each node of the grid
%> where the global @f$\mathbb{RT}_0@f$ solution is given.
%>
%> The cartesian flux on each
%> vertex for each triangle (there's a jump on triangle boundaries) is
%> evaluated via getCartesianCoords().  Then the global flux @f$\vec{q}@f$  on the
%> vertex  @f$V=\vec{x}_V@f$  is interpolated by
%>
%>   @f[\forall V,\quad \vec{q}(\vec{x}_V) := \sum_{T\supset V} |T|\, \vec{q}\big|_{T}(\vec{x}_V) \Big/  \sum_{T'\supset V} |T'| @f]
%>
%> i.e., for each vertex, the flux is the averaged flux on this node wrt. the
%> adjacent triangles weighted by the resp. volume of the triangle.
%>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>   @param g        Grid [AbstractGrid].
%>   @param v_RT0    Vector-valued solution in @f$\mathbb{RT}_0(\mathcal{T})@f$ basis @f$[\#E \times 1]@f$
%>   @retval v_P1P1  Vector-valued solution in @f$\mathbb{P}_1(\mathcal{T})^2@f$ basis @f$[\#V \times 2]@f$

function [v_P1P1] = RT0toP1P1slice(g, v_RT0)
assert(isa(g, 'AbstractGrid'))
RT0.checkdata(g, v_RT0)

v_P1P1 = zeros(g.numV, 2); % cartesian velocity on each vertex (interpolated due to jumps) [#V x 2]
vW = zeros(g.numV, 1); % the weights by which is divided (area of surrounding triangles)
vTV = zeros(g.numT, 3, 2); % cartesian velocity on each vertex for each triangle [#T x 3 x 2]

for kT = 1:g.numT
    edgeNums = g.E0T(kT, :);
    for kV = 1:3
        vTV(kT, kV, 1:2) = RT0.getCartesianCoords(g, v_RT0(edgeNums)', kT, squeeze(g.coordV0T(kT, kV, :)))';
    end
end

% Evaluation of velocity, triangle oriented assembly
for kT = 1:g.numT
    areaT = g.areaT(kT);
    for kV = 1:3
        v_P1P1(g.V0T(kT, kV), 1) = v_P1P1(g.V0T(kT, kV), 1) + vTV(kT, kV, 1) * areaT;
        v_P1P1(g.V0T(kT, kV), 2) = v_P1P1(g.V0T(kT, kV), 2) + vTV(kT, kV, 2) * areaT;
        vW(g.V0T(kT, kV)) = vW(g.V0T(kT, kV)) + areaT;
    end
end

v_P1P1(:, 1) = v_P1P1(:, 1) ./ vW;
v_P1P1(:, 2) = v_P1P1(:, 2) ./ vW;

% % quick visualisation via matlab
% hold on
% for ig.coordV = 1 : size(g.coordV, 1)
%   node = g.coordV(ig.coordV, :)';
%   q = flux(ig.coordV, 1:2)';
%   line([node(1) , q(1)+node(1)], [node(2) , q(2)+node(2)])
% end

return

end
