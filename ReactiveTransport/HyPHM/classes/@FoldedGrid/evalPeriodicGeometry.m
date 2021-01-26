%> @file evalPeriodicGeometry.m This initializes the properties areaE, baryE, vecE, nuE, areaT, baryT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> This initializes the properties
%>
%>   - areaE
%>   - baryE
%>   - vecE
%>   - nuE
%>   - areaT
%>   - baryT
%>
%> of the periodic grid fg via mapping of the geometry data of the
%> underlying rectangular grid fg.rectGrid.


function fg = evalPeriodicGeometry(fg)

assert(isa(fg, 'FoldedGrid'))

% areaE
fg.areaE(fg.mapE) = fg.rectGrid.areaE;
fg.areaE = fg.areaE';

% baryE
baryE1(fg.mapE) = fg.rectGrid.baryE(:, 1);
baryE2(fg.mapE) = fg.rectGrid.baryE(:, 2);
fg.baryE = [baryE1', baryE2'];

% coordV must not be used, since there may be edges with jumping vertices
fg.vecE = zeros(fg.numE, 2);
% for kT = 1 : fg.numT
%   for k = 1 : 3 % locE
%     idxE = fg.E0T(kT, k);
%     v = fg.coordV0T(kT, mod(k, 3) + 1, 1:2) - fg.coordV0T(kT, mod(k+1, 3) + 1, 1:2);
%
%     fg.vecE(idxE, 1) = v(1,1,1);
%     fg.vecE(idxE, 2) = v(1,1,2);
%
%     fg.vecE(idxE, :) = fg.vecE(idxE, :) * fg.sigE0T(kT, k);
%   end
% end
% fg.vecE = -fg.vecE;


for k = 1:3
    idxE = fg.E0T(:, k);
    v = fg.coordV0T(:, mod(k, 3)+1, 1:2) - fg.coordV0T(:, mod(k + 1, 3)+1, 1:2);

    fg.vecE(idxE, 1) = v(:, 1, 1);
    fg.vecE(idxE, 2) = v(:, 1, 2);

    fg.vecE(idxE, :) = fg.vecE(idxE, :) .* fg.sigE0T(:, k);
end

fg.vecE = -fg.vecE;


% nuE (dep: areaE)  ---  ATTENTION, nu does change the sign at fissure
%                        solution: new evaluation of normals
% for iEdge = 1 : fg.numE
%   fg.nuE(iEdge, :)   = fg.vecE(iEdge, :) * [0,-1; 1,0] / fg.areaE(iEdge);
% end

fg.nuE(:, :) = [fg.vecE(:, 2), -fg.vecE(:, 1)] ./ fg.areaE(:);
% areaT, baryT

fg.areaT = fg.rectGrid.areaT;
fg.baryT = fg.rectGrid.baryT;

end
