%> @file multiplyRT0P0_P2P2.m Mapping @f$\mathbb{RT}_0(\mathcal{T})\times \mathbb{P}_0(\mathcal{T})\ni(\vec{q},c)\mapsto\vec{q}c\in \mathbb{P}_2(\mathcal{T})^2 @f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Requires RT0.RT0toP0P0slice.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param  grid Variable of AbstractGrid.
%> @param  c @f$\mathbb{P}_0@f$-data [#T x 1]
%> @param  q @f$\mathbb{RT}_0@f$-data [#E x 1]
%> @retval qc @f$\mathbb{P}_2^2@f$-data [#V+#E x 2]

function qc = multiplyRT0P0_P2P2(grid, q, c)

g = grid;
assert(isa(g, 'AbstractGrid'))

RT0.checkdata(g, q)
P0.checkdata(g, c)

q_P0P0 = RT0.RT0toP0P0slice(g, q); % [#T x 2]

qc = P0P0.P0P0toP2P2slice(g, ...
    [c .* q_P0P0(:, 1), ...
    c .* q_P0P0(:, 2)]);

P2P2.checkdata(g, qc)

end
