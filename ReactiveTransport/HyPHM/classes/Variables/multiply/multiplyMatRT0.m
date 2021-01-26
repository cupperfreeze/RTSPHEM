%> @file multiplyMatRT0.m Mapping @f$\mathbb{R}^{2\times 2}\times \mathbb{RT}_0(\mathcal{T})\ni(\mathbf{A},\vec{q})\mapsto \mathbf{A}\vec{q}\in \mathbb{RT}_0(\mathcal{T}) @f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Requires P1P1toP2P2slice.m and RT0toP1P1slice.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Variable of AbstractGrid.
%> @param A Matrix [2x2]
%> @param q @f$\mathbb{RT}_0@f$-data [#E x 1]
%> @retval Aq Product Aq, @f$\mathbb{RT}_0@f$-data [#E x 1]

function Aq = multiplyMatRT0(grid, A, q)

g = grid;

assert(isa(g, 'AbstractGrid'))
assert(isequal(size(A), [2, 2]))
RT0.checkdata(g, q)

Aq = zeros(g.numE, 1);

q_P2P2 = P1P1.P1P1toP2P2slice(g, RT0.RT0toP1P1slice(g, q));
q_E = q_P2P2(g.numV+1:end, :);

switch g.numE > 20000 % go parallel when array reaches certain size
    case true
        parfor kE = 1:g.numE
            Aq(kE) = dot(A*q_E(kE, :)', g.nuE(kE, :));
        end
    case false
        for kE = 1:g.numE
            Aq(kE) = dot(A*q_E(kE, :)', g.nuE(kE, :));
        end
end

end
