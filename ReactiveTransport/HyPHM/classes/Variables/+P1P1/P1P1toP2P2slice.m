%> @file P1P1toP2P2slice.m Mapping @f$\mathbb{P}_1(\mathcal{T})^2 \rightarrow \mathbb{P}_2(\mathcal{T})^2@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> For each edge, the velocity is computed by taking the mean velocities
%> of surrounding vertices.
%>
%> <h2>Computation</h2>
%> @f[\forall E:\quad \vec{q}(\vec{x}_E) := \big(\vec{q}(\vec{x}_{V_1}) + \vec{q}(\vec{x}_{V_2})\big)\big/ 2@f]
%> @f$V_k@f$ denoting the two vertices of @f$E@f$.
%>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param   g        A grid [ AbstractGrid ].
%> @param   v_P1P1   Vector-valued slice in @f$\mathbb{P}_1(\mathcal{T})^2@f$ @f$[\#V\times2]@f$.
%> @retval  v_P2P2   Vector-valued slice in @f$\mathbb{P}_2(\mathcal{T})^2@f$ @f$[\#V+\#E\times2]@f$.
%

function [v_P2P2] = P1P1toP2P2slice(g, v_P1P1)

assert(isa(g, 'AbstractGrid'))
P1P1.checkdata(g, v_P1P1)

v_P2P2 = [v_P1P1; zeros(g.numE, 2)]; % cartesian velocity on each vertex [#V+#E x 2]

v_P2P2(g.numV+1:end, 1) = (v_P1P1(g.V0E(:, 1), 1) + v_P1P1(g.V0E(:, 2), 1)) / 2;
v_P2P2(g.numV+1:end, 2) = (v_P1P1(g.V0E(:, 1), 2) + v_P1P1(g.V0E(:, 2), 2)) / 2;

%% Remark: The above two lines are the fast way to compute the following:
% for kE = 1 : g.numE
%  idxV = g.V0E(kE, :)';
%  v_P2P2(g.numV + kE, :) = (v_P1P1(idxV(1), :) + v_P1P1(idxV(2), :))/2;
% end

end
