%> @file testSNPP_g_phi.m Neumann boundary data for the Poisson problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description [scalar]

function ret = testSNPP_g_phi(t, x) % k has to be first argument for parallelization
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

if ismember(x(1), [0, 1]) % east, west
    ret = 0;
elseif ismember(x(2), [0, 1]) % north, south
    ret = -t / pi;
else
    ret = NaN;
end

end
