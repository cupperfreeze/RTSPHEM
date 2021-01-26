%> @file div.m The symbolic divergence of a vector-valued function @f$\vec{f}@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> For an example see grad.m.
%>
%> @param f vector-valued funtion @f$[\text{dim}\times1] @f$ of type sym
function ret = div(f, v)

assert(isequal(size(f), size(v)), 'HyPHM: f and v must have the same dimensions.');
assert(iscolumn(f), 'HyPHM: f and v must be column vectors.');

dim = length(v);

ret = sym(0);
for k = 1:dim
    ret = ret + diff(f(k), v(k));
end

end
