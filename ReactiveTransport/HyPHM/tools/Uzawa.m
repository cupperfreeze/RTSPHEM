function [sol, iter] = Uzawa(A, B, f, fp, p, n, scale)
%Apply Uzawa's algorithm to solve the saddle point problem arising from
%Stokes equation in 2D. A and B correspond to the assembled matrices for a
%single dimension (half the dimension of rhs). If available substitution
%operations are performed in the tailored mex method 'backSubs'.

try
    backSubs(sparse(1), sparse(1), 1);
    backSubsFound = true;
catch
    backSubsFound = false;
end
[R, ~, P] = chol(A); % Cholesky decomposition of Laplace part
maxiter = 1000;
tol = 10^(-7);

fp = -scale .* fp;
n = uint32(n/2);
temp = f + B' * (scale .* p);
u = zeros(2*n, 1);
w = zeros(2*n, 1);

if backSubsFound
    u = backSubs(R, P, temp);
else
    u(1:n, 1) = P * (R \ (R' \ (P' * temp(1:n))));
    u((n+1):2*n, 1) = P * (R \ (R' \ (P' * temp((n+1):2 * n))));
end


r = scale .* (B * u) + fp;
d = r;
iter = 0;

while iter < 10 || mod(iter, 2) > 0 || (norm([[A * u(1:n); A * u(n+1:2*n)] - B' * (scale .* p) - f; (-scale .* (B * u) - fp)]) > tol * norm([u; p]) & iter < maxiter)
    temp = B' * (scale .* d);
    if backSubsFound
        w = backSubs(R, P, temp);
    else
        w(1:n, 1) = P * (R \ (R' \ (P' * temp(1:n))));
        w((n+1):2*n, 1) = P * (R \ (R' \ (P' * temp((n+1):2 * n))));
    end

    rTilde = scale .* (B * w);
    rho = (r' * d) / (d' * rTilde);

    p = p - rho * d;
    u = u - rho * w;

    rnew = r - rho * rTilde;

    beta = (rnew' * rnew) / (r' * r);
    d = rnew + beta * d;
    r = rnew;
    iter = iter + 1;

end
%norm([A*u-B'*p-f; -B*u-fp])
assert(iter < maxiter, 'Uzawa did not converge');

sol = [u; scale .* p];
