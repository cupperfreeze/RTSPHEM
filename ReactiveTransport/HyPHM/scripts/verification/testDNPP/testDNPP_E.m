%> @file testDNPP_E.m Exact electrostatic field @f$\vec{E}@f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description @f$[2 \times 1]@f$

function ret = testDNPP_E(t, x)
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

ret = t / pi * [sin(pi * x(1)); cos(pi * x(2))];

end
