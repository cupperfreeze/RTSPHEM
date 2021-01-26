%> @file P0P0toP2P2slice.m Mapping @f$\mathbb{P}_0(\mathcal{T})^2 \rightarrow \mathbb{P}_2(\mathcal{T})^2@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> <h2>Implementation</h2>
%> Sequential call of @c P0P0toP1P1slice and @c P1P1toP2P2slice.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param   g        A grid [AbstractGrid].
%> @param   v_P0P0   Vector-valued slice in @f$\mathbb{P}_0(\mathcal{T})^2@f$ @f$[\#T\times2]@f$.
%> @retval  v_P2P2   Vector-valued slice in @f$\mathbb{P}_2(\mathcal{T})^2@f$ @f$[\#V+\#E\times2]@f$.

function [v_P2P2] = P0P0toP2P2slice(g, v_P0P0)

assert(isa(g, 'AbstractGrid'))
P0P0.checkdata(g, v_P0P0)

v_P2P2 = P1P1.P1P1toP2P2slice(g, P0P0.P0P0toP1P1slice(g, v_P0P0));

end
