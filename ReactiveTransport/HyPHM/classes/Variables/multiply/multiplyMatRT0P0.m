%> @file multiplyMatRT0P0.m Mapping @f$\mathbb{R}^{2\times 2}\times \mathbb{RT}_0(\mathcal{T})\times \mathbb{P}_0(\mathcal{T})\ni(\mathbf{A},\vec{q},c)\mapsto \mathbf{A}\vec{q}c\in \mathbb{RT}_0(\mathcal{T}) @f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Requires P1P1toP2P2slice.m and RT0toP1P1slice.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Variable of AbstractGrid.
%> @param A Matrix [2x2]
%> @param q @f$\mathbb{RT}_0@f$-data [#E x 1]
%> @retval Aq Product Aq, @f$\mathbb{RT}_0@f$-data [#E x 1]

function Aqc = multiplyMatRT0P0(grid, A, q, c)

g = grid;

assert(isa(g, 'AbstractGrid'))
assert(isequal(size(A), [2, 2]))
RT0.checkdata(g, q)
P0.checkdata(g, c)

Aqc = zeros(g.numE, 1);

% mapping to P1P1
q_P1P1 = RT0.RT0toP1P1slice(g, q);
% multiplication with c which is also mapped to P1
c_P1 = P0.P0toP1slice(g, c);
q_P1P1(:, 1) = q_P1P1(:, 1) .* c_P1;
q_P1P1(:, 2) = q_P1P1(:, 2) .* c_P1;
% mapping to P2P2 to evaluate values on the edges
q_P2P2 = P1P1.P1P1toP2P2slice(g, q_P1P1);
% take the normal flux on the edges to finally compute the RT0-coordinates
q_E = q_P2P2(g.numV+1:end, :);
for kE = 1:g.numE
    Aqc(kE) = dot(A*q_E(kE, :)', g.nuE(kE, :));
end

end
