%> @file locMatNum.m Assembly of a local assembly matrix under consideration of the global orientation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param  this  problem class [ Transport ]
%> @param  kT    number of element  [ scalar ]
%> @param  Dinv  the inverse diffusion tensor [ 2x2 ]
%>
%> @retval stiff  local assembly matrix @f$[\int_T D^-1 \vec{phi}_k \cdot\vec{phi}_j]_{j,k}@f$ [3 x 3]
%>

function stiff = locMatNum(this, kT, Dinv)

% INITIALIZATION
g = this.grid;
areaT = g.areaT(kT); % area of triangle [scalar]
areaE = g.areaE(g.E0T(kT, :)); % lengths of edges [3 x 1]
coordV = squeeze(g.coordV0T(kT, :, :))'; % coords of vertices [2 x 3]
quadvals = zeros(3); % weights of quadrature [3 x 3]
sigE0T = g.sigE0T(kT, :)'; % edge signs [3 x 1]

% QUADRATURE POINTS (edge midpoints)
%         | P1(1)  P2(1)  P3(1) |
%         | P1(2)  P2(2)  P3(2) |
baryE = squeeze(g.baryE0T(kT, :, :))';

% ASSEMBLY
for i = 1:3 % rows stiff
    for j = 1:3 % columns stiff
        % QUADRATURE (1)
        quadvals(i, j) = 0;
        for k = 1:3 % loop over quadrature points
            quadvals(i, j) = quadvals(i, j) + (Dinv * (baryE(:, k) - coordV(:, i)))' * (baryE(:, k) - coordV(:, j));
        end
    end
end

stiff = (areaE .* sigE0T) * (areaE .* sigE0T)' .* quadvals / (12 * areaT);

end

%% ...and here's the math:
% (1)
% $  \int_T \phi_i * \phi_j dx =  $
% $  |E_i||E_j|/(4 |T|) * \int_T (x - x_i) * (x - x_j) dx $
% where the integral is evaluated with
% \int f(x) dx = 1/3 * |T| * (f(P1) + f(P2) + f(P3))
