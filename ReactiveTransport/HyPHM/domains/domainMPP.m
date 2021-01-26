%> @file domainMPP.m Unit square bounded domain as use in M++ convergence tests.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Unit square bounded domain as use in M++ convergence tests.
%>
%> @param  a          [ scalar ] x-length of domain
%> @param  b          [ scalar ] y-length of domain
%> @param  h          [ scalar ] mesh fineness (<b>heigth</b> of triangles, <b>not</b> largest edge)
%> @retval g          [ Grid ]   the grid in HyPHM format
%>


function g = domainMPP(a, b, h)

msg = 'HyPHM: Wrong input arguments.';
assert(isscalar(a) && isscalar(b) && isscalar(h), msg)
assert(a > 0 && b > 0 && h > 0, msg)

dim = ceil(1/h); % number of edges per side of the unit square
switch 1 / h == dim
    case true % h is equal to 1/n
        printline(3, 'HyPHM: domainMPP:  Using h = 1/%d', dim);
    otherwise
        printline(3, 'HyPHM: domainMPP:  Using h = 1/%d instead of %.3e', dim, h);
        h = 1 / dim;
end


%%%%%%%%%%%%

%% coordV %%
%%%%%%%%%%%%
[X, Y] = meshgrid(0:h:1);
Xlist = reshape(X, length(X)^2, 1);
Ylist = reshape(Y, length(X)^2, 1);
coordV = [Xlist, Ylist];

%%%%%%%%%

%% V0T %%
%%%%%%%%%

pat1 = [1, dim + 2, 2]; % pattern of "lower-left" triangles
V0T1 = repmat(pat1, dim*(dim + 1), 1) + repmat((0:dim*(dim + 1)-1)', 1, 3);
V0T1(dim+1:dim+1:dim*(dim + 1), :) = [];
pat2 = [dim + 2, dim + 3, 2];
V0T2 = repmat(pat2, dim*(dim + 1), 1) + repmat((0:dim*(dim + 1)-1)', 1, 3);
V0T2(dim+1:dim+1:dim*(dim + 1), :) = [];

%%%%%%%%%%

%% Grid %%
%%%%%%%%%%
g = Grid(coordV, [V0T1; V0T2]);
g.idE(g.baryE(:, 2) == 0) = 1; % south
g.idE(g.baryE(:, 1) == 1) = 2; % east
g.idE(g.baryE(:, 2) == 1) = 3; % north
g.idE(g.baryE(:, 1) == 0) = 4; % west
printline(3, 'Boundary IDs set!')

end
