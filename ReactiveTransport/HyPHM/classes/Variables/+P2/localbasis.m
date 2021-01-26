%> @file +P2/localbasis.m Returns the vector of (values of) local (polynomial quadratic) basis functions of @f$P_2(T)@f$ evaluated in @f$\vec{x}@f$.
%> If @f$\vec{x}\notin T@f$ zero ist returned and a warning thrown. Depends on P1/localbasis.m.
%>
%> The node order is as follows:
%> @code
%>          N3
%>               N4
%>           N6        N2
%>                 N5
%>             N1
%> @endcode
%>
%> @b Example (Usage):
%> See example in P1/localbasis.m.
%>
%>
%> @b Example (Visualization):
%> First, confirm that @f$\phi_k(\vec{x}_j)=\delta_{kj}@f$, @f$\delta@f$ denoting the Kronecker-Delta:
%> @code
%> g = Grid([0,0;1,0;0,1],[1,2,3]); % the reference triangle
%> P2.localbasis(g, 1, [0;0])
%> P2.localbasis(g, 1, [1;0])
%> P2.localbasis(g, 1, [0;1])
%> P2.localbasis(g, 1, [.5;.5])
%> P2.localbasis(g, 1, [0;.5])
%> P2.localbasis(g, 1, [.5;0])
%> @endcode
%> Check the partition of unity by application on an intermediate point.
%>
%> The visualization can be performed by
%> @code
%>  g = Grid([0,0;1,0;0,1],[1,2,3]); % thereference triangle
%>  dx = 25E-3; % sample fineness
%>  [X, Y] = meshgrid(0:dx:1);
%>  Z = cell(6, 1); % one component for every (local) basis function
%>  for j = 1 : 6
%>    Z{j} = zeros(1/dx+1);
%>  end
%>  for k = 1 : 1/dx+1
%>    for ell = 1 : 1/dx+1
%>      tmp = P2.localbasis(g, 1, [X(k, ell); Y(k, ell)]);
%>      for j = 1 : 6
%>        Z{j}(k, ell) = tmp(j);
%>      end
%>    end
%>  end
%>
%>  hold on
%>  for j = 1 : 6
%>  surf(X, Y, Z{j});
%>  end
%> @endcode
%>
%> @param grid  Instance of Grid.
%> @param kT   Global triangle number of @f$T@f$ [ scalar ]  .
%> @param x    Evaluation point [ 2 x 1 ] .
%> @retval ret @f$P_2@f$-coordinates of @f$\vec{x}@f$ on @f$T@f$ [ 6 x 1 ] .


function ret = localbasis(grid, kT, x)

assert(nargin == 3)
assert(isa(grid, 'AbstractGrid'))
assert(isequal(size(x), [2, 1]))

ret = zeros(6, 1);

lam = +P1.localbasis(grid, kT, x); % [3 x 1] barycentric coordinates of x in T or zero if x not in T

for kV = 1:3
    ret(kV) = lam(kV) * (2 * lam(kV) - 1);
end
ret(6) = 4 * lam(1) * lam(2);
ret(4) = 4 * lam(2) * lam(3);
ret(5) = 4 * lam(1) * lam(3);

end
