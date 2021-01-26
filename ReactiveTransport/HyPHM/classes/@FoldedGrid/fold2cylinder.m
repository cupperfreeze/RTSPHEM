%> @file fold2cylinder.m Folding a rectangular grid to a vertically
%> periodic grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Folding a rectangular grid to a vertically periodic grid.
%>
%> This function makes a periodic coordV and V0T from an existent
%> rectangular-bounded coordV_R, V0T_R made, eg. with mesh2d.
%> The result can later be user for generation of the topology via
%> makeTopology().
%>
%> @code
%>   PPPPPPPPPPPPPP
%>   B            B
%>   B     BB     B              P: periodic
%>   B    B  B    B
%>   B     BB     B              B: outer and inner boundary
%>   B            B
%>   PPPPPPPPPPPPPP
%> @endcode
%>
%> Dependencies
%>   FoldedGrid.rectGrid
%>   FoldedGrid.fg.northV
%>   FoldedGrid.fg.eastV
%>   FoldedGrid.fg.southV
%>   FoldedGrid.fg.westV
%>
%> Output
%>   FoldedGrid.coordV
%>   FoldedGrid.mapV
%>
%> @warning Apply this method only on a data for a rectangular grid (may
%> contain holes).  Each boundary vertex has to have an opposite vertex.
%> Furthermore, after grid generation, the barycenters of the folded
%> triangles as well as edges have to be redefined/defined.
%>
%> @note Do not change this file, change fold2torus
%> instead and re-create this function.
%>
%> @sa fold2torus


% _R indicates the rectangular grid, no subindex the periodic one.
function coordV = fold2cylinder(fg)

V0T_tmp = fg.rectGrid.V0T; % this is being transformed to V0T

%%% identification of north vertices with related south's
for kV = 1:length(fg.northV)
    V0T_tmp(V0T_tmp == fg.northV(kV)) = fg.southV(kV);
end

% note: V0T_tmp now contains the new triangle information, but with
% references to the old vertices.  A new enumeration has to be made containing
% only used vertices.  Also the V0T_tmp has to be updated wrt this numbering.

%%% BUILDING THE MAPPING MAPVERTS
% oldDeprecatedVertsMark contains the logical indices of the vertices not
% used any more, oldUsedVertsMark the complement
oldDeprecatedVertsMark = ismember(1:fg.rectGrid.numV, fg.northV);
oldUsedVertsMark = ~oldDeprecatedVertsMark;
% going to numerical indexing we obtain
oldUsedVertsIdxs = find(oldUsedVertsMark); % [#new vertices]
% oldDeprecatedVertsIdxs = find(oldDeprecatedVertsMark); % [#old-#new
% vertices]

%%% BUILDING NEW VERTEX-COORDINATE-TABLE
coordV = fg.rectGrid.coordV(oldUsedVertsMark, :);

%%% BUILDUNG UP THE MAPPING MAP_VERTS
% the kth entry of mapV [numV_R x 1] contains the number of the new
% vertex (wrt the new numbering) to which the old vertex k (wrt the old
% numbering) is mapped to.  This mapping is NOT an injection.
mapV = zeros(fg.rectGrid.numV, 1);
mapV(oldUsedVertsIdxs) = find(oldUsedVertsIdxs);
mapV(fg.northV) = mapV(fg.southV);

fg.mapV = mapV;
fg.coordV = coordV;

end
