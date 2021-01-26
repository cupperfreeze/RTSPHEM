%> @file P2P2toRT0slice.m Mapping @f$\mathbb{P}_2(\mathcal{T})^2 \rightarrow \mathbb{RT}_0(\mathcal{T})@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @todo Increase Approximation order
%>
%> <h2>Computation</h2>
%> For each edge, the global @f$\vec{RT}_0(T)@f$ coordinate is computed by
%> calculation the velocity over the edge, ie
%>
%>
%> @f[\forall E,\quad q_E := \vec{v}(\vec{x}_E) \cdot \vec{\nu}_E @f]
%>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param   g        A grid [ AbstractGrid ].
%> @param   v_P2P2   Vector-valued slice in @f$\mathbb{P}_2(\mathcal{T})^2@f$ @f$[\#V+\#E\times2]@f$.
%> @retval  q_RT0    Representation vector of vector-valued slice in @f$\mathbb{RT}_0(\mathcal{T})@f$ @f$[\#E\times1]@f$.

function [q_RT0] = P2P2toRT0slice(g, v_P2P2)

% printline(-1, 'HyPHM: Be careful -- You lose accuracy when mapping P2P2 to RT0! ')
assert(isa(g, 'AbstractGrid'))
P2P2.checkdata(g, v_P2P2)

% the local RT0-coordinate is the scalar product of the velocity at the
% edge barycenter and the edge unit normal
q_RT0 = v_P2P2(g.numV+1:end, 1) .* g.nuE(:, 1) + v_P2P2(g.numV+1:end, 2) .* g.nuE(:, 2);

end
