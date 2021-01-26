%> @file getRT0coords.m Determination of the local @f$\vec{RT}_0@f$-coordinates of an elementwise constant vector given in cartesian coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param  g     Instance of AbstractGrid.
%> @param  q     2D-vector of @f$P_0(T)^2@f$ [2 x 1]
%> @param  kT    Triangle number.
%> @retval a     Local @f$RT_0(T)@f$-coordinates of q [1 x 3]
function a = getRT0coords(g, q, kT)
assert(isa(g, 'AbstractGrid'))

a = zeros(1, 3);
for k = 1:3
    a(k) = q' * g.nuE(g.E0T(kT, k), :)';
end

end
