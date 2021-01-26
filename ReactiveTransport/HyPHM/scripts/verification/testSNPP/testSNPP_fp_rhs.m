%> @file testSNPP_fp_rhs.m The whole right hand side for the Stokes problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description @f$[2\times 1]@f$

function ret = testSNPP_fp_rhs(t, x) % k has to be first argument for parallelization
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

cx = cos(pi*x(1));
cy = cos(pi*x(2));
sx = sin(pi*x(1));
sy = sin(pi*x(2));

ret = [-(1 + 2 * t * pi^2) * cx * sy + pi / 2 * sin(2 * pi * x(1)); ...
    (1 + 2 * t * pi^2) * sx * cy + pi / 2 * sin(2 * pi * x(2))];

% % stationary
%  ret =  [-(2*t*pi^2)*cx*sy + pi/2*sin(2*pi*x(1));
%          (2*t*pi^2)*sx*cy + pi/2*sin(2*pi*x(2))];
end
