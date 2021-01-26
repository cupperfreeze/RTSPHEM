% Vector of the basis functions of P2 evaluated in x.  The ordering of the
% nodes is as descibed below.
%
% Input
%   g      [ Grid  ]  grid
%   kT     [ 2 x 1 ]  global triangle number of T
%   x      [ 2 x 1 ]  evaluation point
%   varargin
%     'x'       return dphi/dx instead of lambda
%     'y'       return dphi/dy instead of lambda
% Output
%   ret    [ 6 x 1 ]  value of the kth quadratic basis function at
%                          point x.  The nodes are defined as follows:
%
%                       N3
%                                N4
%                                        N2
%                      N5
%                              N6
%
%                     N1
%
% Hint: You can define a function_handle for a fixed triangle T as:
%       phi_fh = @(x) phi(g, kT, x, varargin);
%
% See also lambda
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%        Department of Mathematics, University of Erlangen-Nuremberg,
%        G E R M A N Y

function ret = phi(g, kT, x, varargin)

assert(isequal(size(x), [2, 1]))
assert(nargin <= 4)
assert(isa(g, 'AbstractGrid'))

lam = lambda(g, kT, x);

ret = zeros(6, 1);

if any(ismember(varargin, 'x')) % derivative d phi/dx
    lamx = lambda(g, kT, x, 'x');
    ret(1:3) = lamx .* (4 * lam);
    ret(6) = 4 * (lam(1) * lamx(2) + lam(2) * lamx(1));
    ret(4) = 4 * (lam(2) * lamx(3) + lam(3) * lamx(2));
    ret(5) = 4 * (lam(1) * lamx(3) + lam(3) * lamx(1));
elseif any(ismember(varargin, 'y')) % derivative d phi/dy
    lamy = lambda(g, kT, x, 'y');
    ret(1:3) = lamy .* (4 * lam);
    ret(6) = 4 * (lam(1) * lamy(2) + lam(2) * lamy(1));
    ret(4) = 4 * (lam(2) * lamy(3) + lam(3) * lamy(2));
    ret(5) = 4 * (lam(1) * lamy(3) + lam(3) * lamy(1));
else % phi
    for kV = 1:3
        ret(kV) = lam(kV) * (2 * lam(kV) - 1);
    end
    ret(6) = 4 * lam(1) * lam(2);
    ret(4) = 4 * lam(2) * lam(3);
    ret(5) = 4 * lam(1) * lam(3);
end

end
