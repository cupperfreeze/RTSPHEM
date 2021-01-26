%> @file dt.m The symbolic differentation wrt time @f$t@f$ of a scalar of vector-valued function @f$f@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> For an example see grad.m.
%>
%> @param f funtion @f$[1] @f$ or @f$[\text{dim}\times1] @f$of type sym
%> @param t time variable @f$[1] @f$ of type sym
function ret = dt(f, t)

assert(isa(f, 'sym') && isa(t, 'sym'), 'HyPHM: f and t must be of class symbolic.');
assert(isscalar(t), 'HyPHM: t must be a scalar.');
assert(iscolumn(f), 'HyPHM: f must be a column vector.');

ret = diff(f, t);

end
