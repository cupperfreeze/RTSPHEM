%> @file intT.m Integrate a function @f$f\in L(T)@f$ over @f$T@f$ by quadrature.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Integrate a function @f$f\in L(T)@f$ over @f$T@f$ by quadrature, where
%> @f$f@f$ is given as function_handle and @f$T@f$ is a triangle specified
%> by a Grid @c g and the triangle number or a list of vertex coordinates (soon).
%>
%> The quadrature rule has to be specified by @c type, see
%> trianglequad.m.
%>
%>
%>
%> @code
%> f = @(x) sin(2*pi*x(1)) * cos(2*pi*x(2));
%> g = Grid([0,0;1,0;1,1;0,1], [1,2,3; 1,3,4]);
%> g.visualize('numT')
%> % the integral of f over T_1 is equal to zero.
%> kT = 1;
%> ret = intT(g, kT, f, '23')
%> @endcode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param   g      Grid [Grid ]
%> @param   kT     Global triangle number on g [scalar]
%> @param   intgd  To be integrated scalar function f: R^2 -> R [function_handle]
%> @param   type   Quadrature type [char]
%>
%> @retval  ret    Approximated integral [scalar]
%

function ret = intT(g, kT, intgd, type)

assert(nargin == 4, 'HyPHM: intT expects (g, kT, intgd, type).')
assert(isa(g, 'AbstractGrid'))
assert(isa(intgd, 'function_handle'))
assert(isa(type, 'char'))

% call trianglequad() to obtain quadrature points and weights
coordV = squeeze(g.coordV0T(kT, :, :));
[P, weights] = trianglequad(coordV, type, 0);

ret = 0.0;
for kQ = 1:size(P, 1)
    ret = ret + weights(kQ) * intgd(P(kQ, :)');
end

return

end
