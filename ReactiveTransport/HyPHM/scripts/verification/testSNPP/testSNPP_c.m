%> @file testSNPP_c.m Exact concentration distribution for both species @f$c^+@f$ and @f$c^-@f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param k 1 for @f$c^{+}@f$, 2 for @f$c^{-}@f$
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description [scalar]

function ret = testSNPP_c(k, t, x) % k has to be first argument for parallelization
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

switch k
    case 1
        ret = t * cos(pi * x(1)); % c+ = t cx
    case 2
        ret = t * sin(pi * x(2)); % c- = t sy
    otherwise
        error(msg)
end

end
