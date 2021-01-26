%> @file testDNPP_s.m Source/sink term  @f$s^+@f$ and @f$s^-@f$ for Nernst-Planck problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param k 1 for @f$c^{+}@f$, 2 for @f$c^{-}@f$
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description [scalar]

function ret = testDNPP_s(k, t, x) % k has to be first argument for parallelization
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

cx = cos(pi*x(1));
cy = cos(pi*x(2));
sx = sin(pi*x(1));
sy = sin(pi*x(2));

switch k
    case 1
        ret = (1 + pi^2 * t) * cx + t^2 * (cx^2 - sx^2 - cx * sy + pi * sx * cx * sy);
    case 2
        ret = (1 + pi^2 * t) * sy + t^2 * (sy^2 - cy^2 - cx * sy + pi * sx * cy * cy);
    otherwise
        error(msg)
end

end
