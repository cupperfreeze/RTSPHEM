%> @file testDNPP_phi.m Exact potential distribution @f$\phi@f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description [scalar]

function ret = testDNPP_phi(t, x)
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

ret = t / pi^2 * (cos(pi * x(1)) - sin(pi * x(2)));
end
