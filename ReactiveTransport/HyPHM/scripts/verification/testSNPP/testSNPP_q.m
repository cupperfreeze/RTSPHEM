%> @file testSNPP_q.m Exact mass fluxes @f$\vec{q}^+@f$ and @f$\vec{q}^-@f$ for both species @f$c^+@f$ and @f$c^-@f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param k 1 for @f$\vec{q}^+@f$, 2 for @f$\vec{q}^-@f$
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description @f$[2\times 1]@f$

function ret = testSNPP_q(k, t, x) % k has to be first argument for parallelization
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

sx = sin(pi*x(1));
cx = cos(pi*x(1));
sy = sin(pi*x(2));
cy = cos(pi*x(2));

switch k
    case 1
        ret = t * [pi * sx + t / pi * sx * cx - t * cx * cx * sy; ...
            t / pi * cx * cy + t * sx * cx * cy];
    case 2
        ret = -t * [t / pi * sx * sy + t * cx * sy * sy; ...
            pi * cy + t / pi * sy * cy - t * sx * sy * cy];
    otherwise
        error(msg)
end

end
