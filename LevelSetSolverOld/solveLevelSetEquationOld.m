function [levelSet] = solveLevelSetEquation(grid, initialData, speedFun, velocity, ...
    endTime, dt, varargin)
% Solves the level set equation \partial_t \Phi + F | \nabla \Phi | = 0.

global EPS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input

% Possible input parameters:
%   SaveTime: 'all' indicates that all time steps of the evolution should be
%       stored as columns in the output. Otherwise only the level set at the
%       last time step is returned.
%   Reinitialize: Boolean variable. If true, after every time step the level
%       set function is reinitialized to guarantee the property
%       |\nabla L| = 1 almost everywhere.
%   ContourFigure: Figure in which the level set should be plotted.
%   ContourHandle: Handle of previous level set plots in the same figure.
%       This allows to keep previous plots, i.e. for the initial level set
%       contour.
%   Print: Boolean variable. If true and a 'ContourHandle' is specified, the
%       plotted contour will be saved as a PNG image file.
%   PrintFileName: String for file names of PNG files generated if contours
%       are printed. Pattern: [FILENAME][currentTimeStep].png, so i.e.
%       for FILENAME = 'foo' one has file names like 'foo42.png'

parser = inputParser;
saveTimeValidationFun = @(x) validateattributes(x, {'char'}, ...
    {'nonempty', 'scalartext'});
parser.addParameter('SaveTime', 'last', saveTimeValidationFun);
reinitializeValidationFun = @(x) validateattributes(x, {'logical'}, ...
    {'nonempty', 'scalar'});
parser.addParameter('Reinitialize', false, reinitializeValidationFun);
%     parser.addParameter( 'PlotFigure', [] );
%     contourHandleValidationFun = @(x) validateattributes( x, ...
%         {'matlab.graphics.chart.primitive.Contour'}, {'nonempty', 'scalar'} );
parser.addParameter('ContourHandle', []); %, contourHandleValidationFun );
parser.addParameter('ContourFigure', []);
parser.addParameter('Print', false);
parser.addParameter('PrintFileName', '');
normalVelocityModifierValidationFunc = @(x) validateattributes(x, {'function_handle'}, ...
    {'nonempty', 'scalar'});
parser.addParameter('NormalVelocityModifier', [], normalVelocityModifierValidationFunc);
parser.parse(varargin{:});

isSaveTimeAll = strcmpi(parser.Results.SaveTime, 'all');
assert(isSaveTimeAll || strcmpi(parser.Results.SaveTime, 'last'));
toReinitialize = parser.Results.Reinitialize;
contourHandle = parser.Results.ContourHandle;
contourFig = parser.Results.ContourFigure;
toContourPlot = ~isempty(contourHandle);
shouldPrint = parser.Results.Print && toContourPlot;
printFileName = parser.Results.PrintFileName;
normalVelocityModifierFunc = parser.Results.NormalVelocityModifier;

assert(isempty(velocity) || iscolumn(velocity), ...
    'Parameter ''velocity'' must be an empty matrix or a column vector.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

dim = grid.dimension;
numNodes = grid.nodes;
h = grid.stepSize;
coord = grid.coordinates;

t_ = 0;
phi = initialData(:);

if (toReinitialize)
    phi = reinitializeLevelSet(grid, phi);
end
phiOld = phi;

if (isSaveTimeAll)
    levelSet = NaN(numel(initialData), ceil(endTime / dt)+1);
    levelSet(:, 1) = phi;
else
    levelSet = phi;
end

timeStep = 1;
while (t_ <= endTime + EPS)

    assert(all(phi >= phiOld - EPS));

    if (t_ + dt > endTime - EPS)
        dt = endTime - t_;
        if (abs(dt) < EPS)
            break;
        end
    end

    coordCell = mat2cell(grid.coordinates, ones(1, numNodes), dim);
    timeCell = num2cell(t_*ones(numNodes, 1));
    if (isempty(velocity))
        normalSpeed = cellfun(speedFun, timeCell, coordCell);
    else
        levelSetPartialDerivatives = zeros(numel(phi), 1);
        for d = 1:dim
            levelSetPartialDerivatives(:, d) = levelSetPartialDerivative(grid, phi, d);
        end
        normalSpeed = levelSetPartialDerivatives * velocity;
        if (~isempty(normalVelocityModifierFunc))
            normalSpeed = arrayfun(normalVelocityModifierFunc, normalSpeed);
        end
    end

    phi = levelSetEquationTimeStep(t_+dt, t_, phi, grid, ...
        normalSpeed, 1);

    if (toReinitialize)
        phi = reinitializeLevelSet(grid, phi);
    end

    if (toContourPlot)
        plotLS(grid, phi, timeStep, contourHandle, contourFig);
    end

    if (shouldPrint) % && timeStep >= 50 && timeStep <= 70 )
        %             print( contourFig, '-dpng', [ printFileName, ...
        %                 num2str( timeStep ), '.png' ] );
        print(contourFig, '-depsc', [printFileName, ...
            num2str(timeStep), '.eps']);
    end

    t_ = t_ + dt;
    phiOld = phi;
    timeStep = timeStep + 1;

    if (isSaveTimeAll)
        levelSet(:, timeStep) = phi;
    else
        levelSet = phi;
    end

end % while

if (isSaveTimeAll)
    levelSet(:, all(isnan(levelSet))) = [];
end

end

function plotLS(grid, phi, timeStep, handle, fig)
dim = grid.dimension;
if (dim < 2)
    PHI = phi;
else
    PHI = grid.reshape(phi);
end

set(handle, 'Zdata', PHI, 'LineColor', 'b', 'LineWidth', 3);
if (~isempty(fig))
    %         print( fig, '-dpng', [ 'data/levelSet', num2str( timeStep ), '.png' ] );
    %         print( fig, '-dpng', [ 'data/Biofilm/EllipseLevelSet', ...
    %             num2str( timeStep ), '.png' ] );
end
drawnow;
end
