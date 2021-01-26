%> @file P0P0toP1P1slice.m Mapping @f$\mathbb{P}_0(\mathcal{T})^2 \rightarrow \mathbb{P}_1(\mathcal{T})^2@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> For each vertex, the velocity is computed by taking the mean velocities
%> of surrounding triangles.
%>
%> <h2>Computation</h2>
%> @f[\forall V,\quad \vec{q}(\vec{x}_V) := \sum_{T\supset V} |T|\, \vec{q}\big|_{T} \Big/  \sum_{T'\supset V} |T'| @f]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param   g        A grid [AbstractGrid].
%> @param   v_P0P0   Vector-valued slice in @f$\mathbb{P}_0(\mathcal{T})^2@f$ @f$[\#T\times2]@f$.
%> @retval  v_P1P1   Vector-valued slice in @f$\mathbb{P}_1(\mathcal{T})^2@f$ @f$[\#V\times2]@f$.
%


function [v_P1P1] = P0P0toP1P1slice(g, v_P0P0)

assert(isa(g, 'AbstractGrid'))
P0P0.checkdata(g, v_P0P0)

v_P1P1 = zeros(g.numV, 2); % cartesian velocity on each vertex (interpolated due to jumps) [#V x 2]
vW = zeros(g.numV, 1); % the weights by which is divided (area of surrounding triangles)

% Evaluation of velocity, triangle oriented assembly
for kT = 1:g.numT
    areaT = g.areaT(kT);
    for kV = 1:3
        v_P1P1(g.V0T(kT, kV), 1) = v_P1P1(g.V0T(kT, kV), 1) + v_P0P0(kT, 1) * areaT;
        v_P1P1(g.V0T(kT, kV), 2) = v_P1P1(g.V0T(kT, kV), 2) + v_P0P0(kT, 2) * areaT;
        vW(g.V0T(kT, kV)) = vW(g.V0T(kT, kV)) + areaT;
    end
end

v_P1P1(:, 1) = v_P1P1(:, 1) ./ vW;
v_P1P1(:, 2) = v_P1P1(:, 2) ./ vW;

end
