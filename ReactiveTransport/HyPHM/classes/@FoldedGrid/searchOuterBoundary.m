%> @file searchOuterBoundary.m Initialization of northV, eastV, southV, westV, northE, eastE, southE, westE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Initialization of northV, eastV, southV, westV, northE, eastE, southE, westE.
%>
%> Dependencies
%>
%>   FoldedGrid.rectGrid
%>
%> This function initializes the properties
%>
%>   - FoldedGrid.northV
%>   - FoldedGrid.eastV
%>   - FoldedGrid.southV
%>   - FoldedGrid.westV
%>   - FoldedGrid.northE
%>   - FoldedGrid.eastE
%>   - FoldedGrid.southE
%>   - FoldedGrid.westE
%>
%> where, eg, northV  [1 x #north bdry edges] contains the vertex numbers
%> rel. to north boundary on rectangular grid FoldedGrid.rectGrid

function searchOuterBoundary(this)

assert(isa(this, 'FoldedGrid'))
tol = 1E-14;

%%% fetching information/data
xmin = min(this.rectGrid.coordV(:, 1));
xmax = max(this.rectGrid.coordV(:, 1));
ymin = min(this.rectGrid.coordV(:, 2));
ymax = max(this.rectGrid.coordV(:, 2));

%%% boundary vertex lists, these are enumerated in axes' direction
vertlist = 1:this.rectGrid.numV;
northV = vertlist(abs(this.rectGrid.coordV(:, 2) - ymax) < tol);
eastV = vertlist(abs(this.rectGrid.coordV(:, 1) - xmax) < tol);
southV = vertlist(abs(this.rectGrid.coordV(:, 2) - ymin) < tol);
westV = vertlist(abs(this.rectGrid.coordV(:, 1) - xmin) < tol);

% restore ordering of nodes
[~, idx] = sort(this.rectGrid.coordV(northV, 1));
northV = northV(idx);
[~, idx] = sort(this.rectGrid.coordV(southV, 1));
southV = southV(idx);
[~, idx] = sort(this.rectGrid.coordV(eastV, 2));
eastV = eastV(idx);
[~, idx] = sort(this.rectGrid.coordV(westV, 2));
westV = westV(idx);

%%% check if number of vertices match
assert(length(northV) == length(southV), ...
    'HyPHM: Number of top boundary vertices (%d) does not match the bottom ones (%d)', length(northV), length(southV))
if ~this.isCylinder % only check right and left match, if grid is not folded only vertically
    assert(length(eastV) == length(westV), ...
        'HyPHM: Number of right boundary vertices (%d) does not match the left ones (%d).', length(eastV), length(westV))
end

%%% check if number of vertices is greater than 4 (we need at least 3 edges)
assert(length(northV) > 3 && length(southV) > 3 && length(northV) > 3 && length(southV) > 3, ...
    'HyPHM: Number of edges has to be at least 3 for each boundary.')

%%% check if also coordinates of vertices match


if (this.rectGrid.coordV(northV, 1) - this.rectGrid.coordV(southV, 1)) > 1E-6
    printline(-1, 'Top vertices do not match bottom ones -> increase boundary finness.  Fold the grid anyway, but there may be errors.')
end
if ~this.isCylinder % only check right and left match, if grid is not folded only vertically
    if (this.rectGrid.coordV(eastV, 2) - this.rectGrid.coordV(westV, 2)) > 1E-6
        printline(-1, 'Left vertices do not match right ones -> increase boundary finness.  Fold the grid anyway, but there may be errors.')
    end
end

this.northV = northV;
this.eastV = eastV;
this.southV = southV;
this.westV = westV;


%%% boundary edge lists, these are also enumerated in axes' direction
edgelist = 1:this.rectGrid.numE;
northE = edgelist(this.rectGrid.baryE(:, 2) == ymax);
eastE = edgelist(this.rectGrid.baryE(:, 1) == xmax);
southE = edgelist(this.rectGrid.baryE(:, 2) == ymin);
westE = edgelist(this.rectGrid.baryE(:, 1) == xmin);
% restore ordering of nodes
[~, idx] = sort(this.rectGrid.baryE(northE, 1));
northE = northE(idx);
[~, idx] = sort(this.rectGrid.baryE(eastE, 2));
eastE = eastE(idx);
[~, idx] = sort(this.rectGrid.baryE(southE, 1));
southE = southE(idx);
[~, idx] = sort(this.rectGrid.baryE(westE, 2));
westE = westE(idx);

this.northE = northE;
this.eastE = eastE;
this.southE = southE;
this.westE = westE;

end
