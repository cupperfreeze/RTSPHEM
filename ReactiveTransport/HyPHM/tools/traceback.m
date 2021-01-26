%> @file traceback.m  Computation of the characteristic trace-back within a time-interval @f$[t_1, t_2]@f$ for a fixed point @f$\vec{x}@f$ by successive explicit Euler or explicit Heun stepping.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Characteristic trace back @f$\vec{x}^n = \vec{x}^n(t, \vec{x})@f$ on
%> @f$[t_{n-1}, t_n] \times \Omega@f$
%> (the upper index @f$n@f$ in @f$\vec{x}^n@f$ relates to the time interval @f$[t_{n-1}, t_n]@f$
%> and is not to be confused with the evaluation in @f$t_n@f$!)
%>
%> @f[ \partial_t \vec{x}^n = \mathrm{charVelocity}(t, \vec{x}^n) \quad \text{in}\quad [t_{n-1}, t_n] \times \Omega@f]
%> @f[ \vec{x}^n(t_n, x) = \vec{x}  \quad \text{on} \quad   {t_n} \times \Omega@f]
%>
%> @param    t1       [ 1 x 1, double ]          evaluation time/old time
%> @param    t2       [ 1 x 1, double ]          current time, start of the traceback
%> @param    x        [ 2 x 1, double ]          evaluation point
%> @param    fun      [ 1 x 1, function_handle ] characteristic velocity charVel(t, x)
%> @param    numSteps [ 1 x 1, double ]          number of approximation steps
%> @param    method   [ 1 x n, char ]            approximation method {'euler', 'heun'}
%> @retval   y        [ 2 x 1, double ]          evaluated traceback point

function y = traceback(t1, t2, x, fun, numSteps, method)

if t1 >= t2 || ~isa(fun, 'function_handle') || ~isequal(size(x), [2, 1])
    error('MFEM:wrong_args', 'error in traceback()')
end

y = x;
t = t2;
dt = (t2 - t1) / numSteps;

for i = 1:numSteps
    yold = y;
    told = t;

    switch lower(method)
        case 'euler'
            y = yold - fun(told, yold) * dt;
        case 'heun'
            K1 = fun(told, yold);
            K2 = fun(told-dt, yold-K1*dt);
            y = yold - 0.5 * dt * (K1 + K2);
    end

    t = told - dt;
end

return

end
