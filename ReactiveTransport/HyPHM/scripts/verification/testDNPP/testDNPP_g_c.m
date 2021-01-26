%> @file testDNPP_g_c.m Flux boundary data for the Nernst-Planck problem for both species @f$c^+@f$ and @f$c^-@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param k 1 for @f$c^{+}@f$, 2 for @f$c^{-}@f$
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description [scalar]

function ret = testDNPP_g_c(k, t, x) % k has to be first argument for parallelization
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

sx = sin(pi*x(1));
cx = cos(pi*x(1));
sy = sin(pi*x(2));
% cy = cos(pi*x(2));

switch k
    case 1 % c+
        if x(1) == 0 % west
            ret = t^2 * sy;
        elseif x(1) == 1 % east
            ret = -t^2 * sy;
        elseif ismember(x(2), [0, 1]) % north, south
            ret = -t^2 * sx * cx - t^2 / pi * cx;
        else
            ret = NaN;
        end
    case 2 % c-
        if ismember(x(1), [0, 1]) % east, west
            ret = t^2 * sy^2;
        elseif ismember(x(2), [0, 1]) % north, south
            ret = pi * t;
        else
            ret = NaN;
        end
    otherwise % k not 1 or 2
        error(msg)
end

end
