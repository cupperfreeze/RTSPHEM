% This function can be used as function handle (eg for integration) which
% returns the kth barycentric coordinate of a given point x on triangle T.
%
% ret = lambdak_fh(g, kT, x, varargin)
%
% Input
%   g      [ Grid  ]  grid
%   kT     [ 2 x 1 ]  global triangle number of T
%   x      [ 2 x 1 ]  evaluation point
%   k      {1, 2, 3}  component number
%   varargin
%     'x'       return dlambda/dx instead of lambda
%     'y'       return dlambda/dy instead of lambda
%
% Output
%   ret    [ 3 x 1 ]  barycentric coords of x on T
%
%                       N3
%                                N2
%
%                     N1
%
% Hint: You can define a function_handle for a fixed triangle T as:
%       lambda_fh = @(x) lambda(g, kT, x, varargin);
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%        Department of Mathematics, University of Erlangen-Nuremberg,
%        G E R M A N Y

function ret = lambdak_fh(g, kT, x, k, varargin)

assert(isequal(size(x), [2, 1]))
assert(nargin <= 5)
assert(isa(g, 'Grid'))

verts = squeeze(g.coordV0T(kT, :, :));
D = [verts, ones(3, 1)];
detD = det(D);


if any(ismember(varargin, 'x')) % derivative dlambda/dx
    y1 = verts(1, 2);
    y2 = verts(2, 2);
    y3 = verts(3, 2);
    ret = [y2 - y3; y3 - y1; y1 - y2];
    ret = ret(k);
elseif any(ismember(varargin, 'y')) % derivative dlambda/dy
    x1 = verts(1, 1);
    x2 = verts(2, 1);
    x3 = verts(3, 1);
    ret = [x3 - x2; x1 - x3; x2 - x1];
    ret = ret(k);
else % lambda
    D(k, 1:2) = x;
    ret = det(D);
end

ret = ret / detD;

end
