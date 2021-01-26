%> @file AbstractGrid/visualize.m Visualization of instances of the class Grid or FoldedGrid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @b Options (additional visualization of attributes):
%>
%> <table>
%> <tr>
%>   <th>Visualize</th>
%>   <th>Tag</th>
%> </tr>
%> <tr>
%>   <td>edge normals</td>
%>   <td><code>'edgenormals'</code> or <code>'nuE'</code></td>
%> </tr>
%> <tr>
%>   <td>edge vectors</td>
%>   <td><code>'edgevectors'</code> or <code>'vecE'</code></td>
%> </tr>
%> <tr>
%>   <td>edge numbers</td>
%>   <td><code>'edgenums'</code> or <code>'numE'</code></td>
%> </tr>
%> <tr>
%>   <td>edge IDs</td>
%>   <td><code>'idE'</code></td>
%> </tr>
%> <tr>
%>   <td>vertex numbers</td>
%>   <td><code>'vertnums'</code> or <code>'numV'</code></td>
%> </tr>
%> <tr>
%>   <td>vertex coordinates</td>
%>   <td><code>'coordV'</code></td>
%> </tr>
%> <tr>
%>   <td>triangle numbers</td>
%>   <td><code>'triangnums'</code> or <code>'numT'</code></td>
%> </tr>
%> <tr>
%>   <td>Donald diagram</td>
%>   <td><code>'Donald' or 'donald'</code></td>
%> </tr>
%> </table>
%>
%> <b>Example 1</b>
%> @code
%> g = domainRectangle(0, 1, 0, 1, .4)
%> g.visualize('numV', 'numE', 'Donald', 'numT')
%> @endcode
%> Generates the following figure:
%> @image html images/visualizeGrid1.png
%>
%> <b>Example 2</b>
%> @code
%> g = domainCirce(0, 0, 1, .4)
%> g.visualize('idE', 'vecE', 'nuE')
%> @endcode
%> Generates the following figure:
%> @image html images/visualizeGrid2.png
function visualize(grid, varargin)

%figure('Color',[1 1 1]);
hold on

vertexColor = [.7, .6, .0];
edgeColor = [.0, .5, .0];
donaldColor = [.8, .0, .0];

xmin = min(min(grid.baryE(:, 1)), min(grid.coordV(:, 1)));
xmax = max(max(grid.baryE(:, 1)), max(grid.coordV(:, 1)));
ymin = min(min(grid.baryE(:, 2)), min(grid.coordV(:, 2)));
ymax = max(max(grid.baryE(:, 2)), max(grid.coordV(:, 2)));
xborder = .4 * (xmax - xmin);
yborder = .4 * (ymax - ymin);

axis([xmin - xborder, xmax + xborder, ymin - yborder, ymax + yborder])

if isa(grid, 'FoldedGrid')
    trisurf(grid.rectGrid.V0T, grid.rectGrid.coordV(:, 1), ...
        grid.rectGrid.coordV(:, 2), zeros(size(grid.rectGrid.coordV, 1), 1), ...
        'facecolor', 'none');
else
    trisurf(grid.V0T, grid.coordV(:, 1), grid.coordV(:, 2), zeros(size(grid.coordV, 1), 1), 'facecolor', 'none');
end

axis off
view(0.0, 90.0)
daspect([1, 1, 1])


%%%%%%%%%%%%%%%%%%%%%%%%  Option Edge Normals  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              (have to be plotted before text plots)                     %
if any(ismember(varargin, 'edgenormals')) || any(ismember(varargin, 'nuE'))
    arrow = @(x, y) line(x, y); %&& arrowh(x, y, 'b');                                improved SG
    % scalefactor sc is def. by half of min edge length (all normals have
    % length one)
    sf = .3 * min(grid.areaE);
    for kE = 1:grid.numE
        arrow([grid.baryE(kE, 1), grid.baryE(kE, 1) + sf * grid.nuE(kE, 1)], ...
            [grid.baryE(kE, 2), grid.baryE(kE, 2) + sf * grid.nuE(kE, 2)]);
    end
end
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%  Option Edge Vectors  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              (have to be plotted before text plots)                     %
if any(ismember(varargin, 'edgevectors')) || any(ismember(varargin, 'vecE'))
    arrow = @(x, y) line(x, y, 'Color', 'k') && arrowh(x, y, edgeColor);
    %   for kE = 1 : grid.numE
    %     arrow([grid.baryE(kE, 1), grid.baryE(kE, 1) + sf*grid.nuE(kE, 1)], ...
    %       [grid.baryE(kE, 2), grid.baryE(kE, 2) + sf*grid.nuE(kE, 2)]);
    %   end
    for kE = 1:grid.numE
        arrow([grid.coordV(grid.V0E(kE, 1), 1), ...
            grid.coordV(grid.V0E(kE, 1), 1) + grid.vecE(kE, 1)], ...
            [grid.coordV(grid.V0E(kE, 1), 2), ...
            grid.coordV(grid.V0E(kE, 1), 2) + grid.vecE(kE, 2)]);
    end
end
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Option Donald  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              (have to be plotted before text plots)                     %
if any(ismember(varargin, 'Donald')) || any(ismember(varargin, 'donald'))
    for kT = 1:grid.numT
        weights = {[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]};
        for kE = 1:3
            line([grid.baryT(kT, 1), grid.baryE(grid.E0T(kT, kE), 1)], ...
                [grid.baryT(kT, 2), grid.baryE(grid.E0T(kT, kE), 2)], 'Color', donaldColor);
            dot(grid.coordV(grid.V0T(kT, :), 1), weights{kE});
            text(dot(grid.coordV(grid.V0T(kT, :), 1), weights{kE}), ...
                dot(grid.coordV(grid.V0T(kT, :), 2), weights{kE}), num2str(grid.V0T(kT, kE)), 'Color', donaldColor);
        end

    end
end
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%  Option Vertices Numbers  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
if any(ismember(varargin, 'vertnums')) || any(ismember(varargin, 'numV'))
    for kV = 1:grid.numV
        text(grid.coordV(kV, 1), grid.coordV(kV, 2), int2str(kV), 'Color', vertexColor);
    end
end
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Option Vertices Coordinates  %%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
if any(ismember(varargin, 'coordV'))
    for kV = 1:grid.numV
        text(grid.coordV(kV, 1), grid.coordV(kV, 2), ...
            ['(', num2str(grid.coordV(kV, 1)), ', ', num2str(grid.coordV(kV, 2)), ')'], ...
            'Color', vertexColor);
    end
end
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%  Option Edge IDs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
if any(ismember(varargin, 'idE'))
    for kE = 1:grid.numE
        if grid.idE(kE) ~= 0 % do not print inner edge ids (those are equal to zero)
            text(grid.baryE(kE, 1), grid.baryE(kE, 2), int2str(grid.idE(kE)), 'Color', edgeColor);
        end
    end
end
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Option Edge Numbers  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
if any(ismember(varargin, 'edgenums')) || any(ismember(varargin, 'numE'))
    for kE = 1:grid.numE
        text(grid.baryE(kE, 1), grid.baryE(kE, 2), int2str(kE), 'Color', edgeColor);
    end
end
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%  Option Triangle Numbers  %%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
if any(ismember(varargin, 'triangnums')) || any(ismember(varargin, 'numT'))
    for kT = 1:grid.numT
        text(grid.baryT(kT, 1), grid.baryT(kT, 2), int2str(kT));
    end
end
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % visualize
