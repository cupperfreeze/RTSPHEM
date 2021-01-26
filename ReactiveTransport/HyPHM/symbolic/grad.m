%> @file grad.m The symbolic gradient of a scalar-valued function @f$f@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @b Example:
%>
%> @code
%> syms x y real  % real variables x and y
%> xvec = [x; y]; % vector of variables
%> u = sin(2*pi*x)*y^2; % some function u = u(x, y)
%> grad(u, xvec)
%> @endcode
%>
%> @param f scalar-valued funtion @f$[1]@f$ of type @c sym
%> @param x vector of variables @f$[\text{dim}\times1] @f$ of type @c sym
function ret = grad(f, x)

assert(isa(f, 'sym') && isa(x, 'sym'), 'HyPHM: f and x must be of class symbolic.');
assert(iscolumn(x), 'HyPHM: x must be a column vector.');
assert(isscalar(f), 'HyPHM: f must be a scalar.');

dim = length(x);

ret = sym(zeros(dim, 1));
for k = 1:dim
    ret(k) = diff(f, x(k));
end

end
