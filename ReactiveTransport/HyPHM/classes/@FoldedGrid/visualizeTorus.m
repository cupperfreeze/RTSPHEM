%> @file visualizeTorus.m Similar to visualize() this function displays a FoldedGrid as a torus.

function visualizeTorus(obj, varargin)
assert(isa(obj, 'FoldedGrid'))
assert(isnumeric(obj.coordV))
assert(isnumeric(obj.V0T))

coordV = obj.coordV;
V0T = obj.V0T;
numV = obj.numV;

hold on
% axis equal


% init (we now go 3d)
x = coordV(:, 1);
y = coordV(:, 2);
z = zeros(numV);

xmax = max(x);
ymax = max(y);
xmin = min(x);
ymin = min(y);

% this is the (defined) distance between the last node to the very first
dist2endnodex = (xmax - xmin) / (sum(coordV(:, 2) == ymin) - 1);
dist2endnodey = (ymax - ymin) / (sum(coordV(:, 1) == xmin) - 1);

% if data is not a torus but has to be folded anyway
% dist2endnodex = 0;
% dist2endnodey = 0;

xend = xmax + dist2endnodex;
yend = ymax + dist2endnodey;

%%% folding the grid to a cylinder (y == const)
circumference = xend - xmin;
radius = circumference / 2 / pi;
centerx = xmin;
centerz = radius;
for iVert = 1:numV
    phi = (x(iVert) - xmin) / circumference * 2 * pi;
    x(iVert) = centerx + radius * sin(phi);
    z(iVert) = centerz - radius * cos(phi);
end


%%% folding the grid to a torus (z == const)
centerx = xend;
centery = ymin;

for iVert = 1:numV
    circumference = xend - x(iVert);
    radius = circumference / 2 / pi;
    phi = (y(iVert) - ymin) / (yend - ymin) * 2 * pi;
    x(iVert) = centerx + radius * sin(phi);
    y(iVert) = centery - radius * cos(phi);
end


for iTriang = 1:size(V0T, 1)
    patch(x(V0T(iTriang, :)), y(V0T(iTriang, :)), z(V0T(iTriang, :)), 'r');

    % line(x([V0T(iTriang, 1), V0T(iTriang, 2) ]), ...
    %   y([V0T(iTriang, 1), V0T(iTriang, 2) ]), ...
    %   z([V0T(iTriang, 1), V0T(iTriang, 2) ]));
    % line(x([V0T(iTriang, 2), V0T(iTriang, 3) ]), ...
    %   y([V0T(iTriang, 2), V0T(iTriang, 3) ]), ...
    %   z([V0T(iTriang, 2), V0T(iTriang, 3) ]));
    % line(x([V0T(iTriang, 3), V0T(iTriang, 1) ]), ...
    %   y([V0T(iTriang, 3), V0T(iTriang, 1) ]), ...
    %   z([V0T(iTriang, 3), V0T(iTriang, 1) ]));
end

%%%%%%%%%%%%%%%%%%%%%%  Option Vertices Numbers  %%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
if any(ismember(varargin, 'vertnums'))
    for iVert = 1:numV
        text(x(iVert), y(iVert), z(iVert), int2str(iVert));
    end
end
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
