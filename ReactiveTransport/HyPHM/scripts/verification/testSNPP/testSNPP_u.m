%> @file testSNPP_u.m Exact velocity distribution @f$\vec{u}@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description @f$[2\times 1]@f$

function ret = testSNPP_u(t, x) % k has to be first argument for parallelization
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

sx = sin(pi*x(1));
cx = cos(pi*x(1));
sy = sin(pi*x(2));
cy = cos(pi*x(2));

ret = t * [-cx * sy; ...
    sx * cy];

end
