%> @file circumcircleTest.m Test if a point lies within a passing through three given points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Test if a point @f$\vec{p}@f$ lies within a passing through three given points @f$\vec{a}, \vec{b}, \vec{c}@f$.
%>
%> @f[\mathtt{ret} = \left(\det\begin{bmatrix}1&a_1&a_2&|\vec{a}|^2\\1&b_1&b_2&|\vec{b}|^2\\
%> 1&c_1&c_2&|\vec{c}|^2\\1&p_1&p_2&|\vec{p}|^2\\\end{bmatrix}\,\det\begin{bmatrix}
%> 1&a_1&a_2\\1&b_1&b_2\\1&c_1&c_2\end{bmatrix} < 0\right)@f]
%>
%>
%>
%>
%>
%> @param a first point on the circle's boundary [2x1]
%> @param b second point on the circle's boundary [2x1]
%> @param c third point on the circle's boundary [2x1]
%> @param p this point is either inside or outside the circle [2x1]
%> @retval ret [1x1 logical]

function ret = circumcircleTest(a, b, c, p)

delta = [ones(4, 1), [a'; b'; c'; p'], [norm(a); norm(b); norm(c); norm(p)].^2];
gamma = [ones(3, 1), [a'; b'; c']];
ret = det(gamma) * det(delta) < 0;

end
