%> @file +P1/localbasis.m Barycentric coordinates of point @f$\vec{x}@f$ on triangle @f$T@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Barycentric coordinates of point @f$\vec{x}@f$ on triangle @f$T@f$, which is also the vector of (values of)
%> local basis functions of @f$\mathbb{P}_1(T)@f$ evaluated in @f$\vec{x}@f$. If @f$\vec{x}\notin T@f$ then zero is returned
%> and a warning is thrown.
%>
%> The node order is as follows:
%> @code
%>          N3
%>                   N2
%>
%>             N1
%> @endcode
%>
%> @b Example (Usage):
%> Let the triangle @f$T@f$ have the local node numbers @f$1@f$ to @f$3@f$.
%> Let @f$\lbrace\varphi_k\rbrace_k@f$ be the <i>basis</i> of
%> @f$\mathbb{P}_1(T)@f$
%> and let @f$f_k@f$ denote the <i>local coordinates</i> of some
%> <i>discrete function</i> @f$f_h\in \mathbb{P}_1(\mathcal{T}_h)@f$ (cf. Variable.type). This method returns
%> the values of each basis function at the point @f$\vec{x}@f$
%> @f[  (\varphi_1(\vec{x}), \ldots,\varphi_3(\vec{x}))^T  @f]
%> for a fixed point @f$\vec{x}\in T@f$.  It holds
%> @f[  f_h(\vec{x})|_T = \begin{pmatrix}\varphi_1(\vec{x})\\\vdots\\\varphi_3(\vec{x})\end{pmatrix} \cdot  \begin{pmatrix}f_1\\\vdots\\f_3\end{pmatrix} = \sum_{k=1}^3 \varphi_k(\vec{x})\,f_k @f]
%> with @f$\sum_{k=1}^3 \varphi_k(\vec{x}) = 1@f$ for all @f$\vec{x}@f$ (partition of unity).
%>
%>
%> @b Example (Visualization):
%> @code
%> g = Grid([0,0;1,0;0,1],[1,2,3]); % the reference triangle
%> dx = 25E-3; % sample fineness
%> [X, Y] = meshgrid(0:dx:1);
%> Z = cell(3, 1); % one component for every (local) basis function
%> for j = 1 : 3
%>   Z{j} = zeros(1/dx+1);
%> end
%> for k = 1 : 1/dx+1
%>   for ell = 1 : 1/dx+1
%>     tmp = P1.localbasis(g, 1, [X(k, ell); Y(k, ell)]);
%>     for j = 1 : 3
%>       Z{j}(k, ell) = tmp(j);
%>     end
%>   end
%> end
%> hold on
%> surf(X, Y, Z{1});
%> surf(X, Y, Z{2});
%> surf(X, Y, Z{3});
%> @endcode
%>
%> @param grid  Instance of Grid.
%> @param kT   Global triangle number of @f$T@f$ [ scalar ]  .
%> @param x    Evaluation point [ 2 x 1 ] .
%> @retval ret Barycentric coordinates of  @f$\vec{x}@f$ on @f$T@f$ [ 3 x 1 ] .


function ret = localbasis(grid, kT, x)

assert(nargin == 3)
assert(isa(grid, 'AbstractGrid'))
assert(isequal(size(x), [2, 1]))

g = grid;
ret = zeros(3, 1);

coordV = squeeze(g.coordV0T(kT, :, :));
D = [coordV, ones(3, 1)];
detD = det(D);

for kV = 1:3
    C = D;
    C(kV, 1:2) = x;
    ret(kV) = det(C);
end

ret = ret / detD;

ret(abs(ret) < 1E-14) = 0; % erase arithmetical error!!!

if any(ret < 0) % this is only possible when arithmetical errors are suppressed!
    warning('HyPHM:P1:localbasis', 'HyPHM: Requested local basis valued not on considered triangle.  kT = %d, x = [%.2e;%.2e], basis = [%.2e;%.2e;%.2e].  You are outside the support of the function; returning basis = [0;0;0] instead.', kT, x(1), x(2), ret(1), ret(2), ret(3))
    ret = zeros(3, 1);
end

end
