%> @file P0toP1slice.m Mapping @f$\mathbb{P}_0(\mathcal{T}) \rightarrow \mathbb{P}_1(\mathcal{T})@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> For each vertex, the value is computed by taking the mean values
%> of surrounding triangles.
%>
%> <h2>Computation</h2>
%> See P0P0toP1P1slice.m
%>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param   g      A grid [AbstractGrid].
%> @param   v_P0   Scalar-valued slice in @f$\mathbb{P}_0(\mathcal{T})@f$ @f$[\#T\times1]@f$.
%> @retval  v_P1   Scalar-valued slice in @f$\mathbb{P}_1(\mathcal{T})@f$ @f$[\#V\times1]@f$.

function [v_P1] = P0toP1slice(g, v_P0)

assert(isa(g, 'AbstractGrid'))
P0.checkdata(g, v_P0)

v_P1 = zeros(g.numV, 1); % cartesian velocity on each vertex (interpolated due to jumps) [#V x 1]
vW = zeros(g.numV, 1); % the weights by which is divided (area of surrounding triangles)

% Evaluation of velocity, triangle oriented assembly
for kT = 1:g.numT
    areaT = g.areaT(kT);
    for kV = 1:3
        v_P1(g.V0T(kT, kV)) = v_P1(g.V0T(kT, kV)) + v_P0(kT) * areaT;
        vW(g.V0T(kT, kV)) = vW(g.V0T(kT, kV)) + areaT;
    end
end

v_P1 = v_P1 ./ vW;

end
