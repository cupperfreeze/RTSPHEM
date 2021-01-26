%> @file +C/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x is a dummy which will not be sampled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> There's only one dof for the whole domain per time step.
function ret = setfh2data(time, ~, fun2)

ret = fun2(time, [inf; inf]);

end
