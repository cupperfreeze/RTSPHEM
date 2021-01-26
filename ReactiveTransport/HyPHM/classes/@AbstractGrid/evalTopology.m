%> @file evalTopology.m Evaluates the topological properties based on @ref BC2005.

%%% EVALUATION OF TOPOLOGY OF THE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following function was based on a code from C. Bahriawati and   %
% C. Carstensen underlying the following copyright:                   %
%                                                                     %
% Copyright 2008 C. Bahriawati, C. Carstensen                         %
% cf. "Three Matlab Implementations of the Lowest-Order Raviart-Thomas%
%      MFEM with a Posteriori Error Control"                          %
% by C.Bahriawati and C. Carstensen                                   %


function [sigE0T, V2T, V2E, V0E, T0E, E0T] ...
    = evalTopology(coordV, V0T)

numT = size(V0T, 1);
numV = size(coordV, 1);

% V2T (deps.: V0T) % THIS IMPLICITELY DEFINES THE SIGNS OF THE EDGES
try
    V2T = fsparse(V0T(:, [1, 2, 3, 1, 2, 3, 1, 2, 3]), ...
        V0T(:, [2, 3, 1, 2, 3, 1, 2, 3, 1]), ...
        [(1:numT)', zeros(numT, 3), (1:numT)', zeros(numT, 3), (1:numT)'], [], [], 1);
catch
    V2T = sparse(V0T(:, [1, 2, 3, 1, 2, 3, 1, 2, 3]), ...
        V0T(:, [2, 3, 1, 2, 3, 1, 2, 3, 1]), ...
        [(1:numT)', zeros(numT, 3), (1:numT)', zeros(numT, 3), (1:numT)'], numV, numV);
end
% V2E (sparse, symmetric, deps.: V2T) % THIS IMPLICITELY DEFINES THE EDGE NUMBERS
[k, ell] = find(triu(V2T + V2T'));
V2E = sparse(k, ell, 1:size(k, 1), numV, numV);
V2E = V2E+V2E';

% V0E, T0E (deps.: V2T, V2E) (see Footnote (1))
%idxE = reshape(full(diag(V2E(V0T(end:-1:1, [1,2,3]), V0T(end:-1:1, [2,3,1])))), numT, 3)';
idxE = reshape(full(V2E(sub2ind(size(V2E), V0T(end:-1:1, [1, 2, 3]), V0T(end:-1:1, [2, 3, 1])))), numT, 3)';

V0E(idxE(:), 1) = reshape(V0T(end:-1:1, [1, 2, 3])', 3*numT, 1);
V0E(idxE(:), 2) = reshape(V0T(end:-1:1, [2, 3, 1])', 3*numT, 1);

%T0E(idxE(:), 1) = reshape(reshape(full(diag(V2T(V0T(end:-1:1,[1,2,3]), V0T(end:-1:1,[2,3,1])))), numT, 3)', 3*numT, 1);
%T0E(idxE(:), 2) = reshape(reshape(full(diag(V2T(V0T(end:-1:1,[2,3,1]), V0T(end:-1:1,[1,2,3])))), numT, 3)', 3*numT, 1);

T0E(idxE(:), 1) = reshape(reshape(full(V2T(sub2ind(size(V2T), V0T(end:-1:1, [1, 2, 3]), V0T(end:-1:1, [2, 3, 1])))), numT, 3)', 3*numT, 1);
T0E(idxE(:), 2) = reshape(reshape(full(V2T(sub2ind(size(V2T), V0T(end:-1:1, [2, 3, 1]), V0T(end:-1:1, [1, 2, 3])))), numT, 3)', 3*numT, 1);


% E0T (deps.: V2E)
%E0T = reshape(full(diag(V2E(V0T(:, [2 3 1]), V0T(:, [3 1 2])))), numT, 3);
E0T = reshape(full(V2E(sub2ind(size(V2E), V0T(:, [2, 3, 1]), V0T(:, [3, 1, 2])))), numT, 3);


% sigE0T (deps.: T0E, E0T)
sigE0T = 1 - 2 * (reshape(T0E(E0T, 2), numT, 3) == repmat((1:numT)', 1, 3));

end % evalTopology


%%%%%%%%%%%%%%%%%%%%%

%% Footnote (1) %%
%%%%%%%%%%%%%%%%%%%%%
% The original idea is:
%
% V0E = zeros(numE, 2);
% T0E = zeros(numE, 2);
% for kT = 1 : numT
%   for k = 1 : 3
%     idxE = V2E(V0T(kT, k), V0T(kT, rem(k,3) + 1));
%     if ~V0E(idxE, 1) % if row idxE of V0E was not yet defined
%       V0E(idxE, :) = [V0T(kT, k),                            V0T(kT, rem(k,3)+1)];
%       T0E(idxE, :) = [V2T(V0T(kT,k), V0T(kT,rem(k,3)+1)),    V2T(V0T(kT,rem(k,3)+1), V0T(kT,k))];
%     end
%   end
% end
%
% Note, that we can omit the if-clause when indexing backwards:
%
% for kT = numT : -1 : 1
%   for k = 1 : 3
%     idxE = V2E(V0T(kT, k), V0T(kT, rem(k,3) + 1));
%     V0E(idxE, :) = [V0T(kT, k),                            V0T(kT, rem(k,3)+1)];
%     T0E(idxE, :) = [V2T(V0T(kT,k), V0T(kT,rem(k,3)+1)),    V2T(V0T(kT,rem(k,3)+1), V0T(kT,k))];
%   end
% end
%
% And this can be vectorized by some indexing magic.
