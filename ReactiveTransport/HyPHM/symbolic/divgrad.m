%> @file divgrad.m The symbolic laplace of a scalar-valued function @f$f@f$ (requires div.m and grad.m).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> For an example see grad.m.
%>
%> @param f scalar of vector-valued funtion @f$[1] @f$ or @f$[\text{dim}\times1]@f$ of type @c sym
%> @param x vector of variables @f$[\text{dim}\times1]@f$ of type @c sym
function ret = divgrad(f, x)

assert(isa(f, 'sym') && isa(x, 'sym'), 'HyPHM: f and x must be of class symbolic.');
assert(iscolumn(x), 'HyPHM: x must be a column vector.');
assert(iscolumn(f), 'HyPHM: f must be a column vector.');

dimf = length(f);
ret = sym(zeros(dimf, 1));
for k = 1:dimf
    ret(k) = div(grad(f(k), x), x); % WHY DOES THIS THROW AN ERROR?!  THIS WORKS IN THE TERMINAL!
end

end
