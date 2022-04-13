% Evolve geometry according to two-phase solid with different normal
% interface velocities. Porosity and diffusivity of the arising unit cell
% geometries are calculated

% cf. [5] Section 4.2

dimension = 2;
numPartitions = 64; % Number of partitions in each direction
endTime = 3.5;
timeStepSize = 0.5/numPartitions;

cellGrid = FoldedCartesianGrid( dimension, kron( ones(1,dimension), [-0.5 0.5] ), numPartitions*ones(1,dimension) );
coord = cellGrid.coordinates;

initialLevelSetFunc = @(x) abs( x(1) ) - 0.05;
% lsf = @(x) min( abs( x(1) + 0.25 ) - 0.05, abs( x(1) - 0.25 ) - 0.05 );
coordCell = mat2cell( coord, ones(1, cellGrid.nodes), dimension );
initialData = cellfun( initialLevelSetFunc, coordCell );
clear coordCell;

X = cellGrid.reshape( coord(:,1) );
Y = cellGrid.reshape( coord(:,2) );
clear coord;

plotInitial = cellGrid.reshape( initialData );

% speed = @(t,x) -35 * (x(1)-0.3).^2 * ( x(1) + 0.5 ).^2 * ( x(1) <= 0.3 );
% variableSpeed = @(t,x) speed(t,x) * ( ( t <= 0.5*endTime ) - ( t > 0.5*endTime ) );

fig = figure;
[~, defHandle] = contour( X, Y, plotInitial, zeros(1,2), 'LineColor', 'r', ...
    'LineWidth', 3 );
set(gca,'DataAspectRatio',[1 1 1]);
% xlabel('x');
% ylabel('y');
handle = copyobj( defHandle, gca );
set( handle, 'Zdata', plotInitial, 'LineColor', 'b', 'LineWidth', 3 );
set( gca, 'xtick', [] );
set( gca, 'ytick', [] );
folderName = 'data/IrregularSoil/Periodic/';
% if ( ~exist( folderName, 'dir' ) )
%     mkdir( folderName );
% end
toPrint = false;
printFileName = [ folderName, 'TimeStep_' ];
if( toPrint )
    print( fig, '-depsc', [ printFileName, num2str( 0 ), '.eps' ] );
end
toSave = false; % Flag: save variables to disk
speedFun = @speedFunction;

levelSet = solveLevelSetEquationOld( cellGrid, initialData, @speedFunction, ...
    [], endTime, timeStepSize, 'ContourHandle', handle, 'SaveTime', 'all', ...
    'ContourFigure', fig, 'Print', toPrint, 'PrintFileName', printFileName, ...
    'Reinitialize', false );

clear velocity;
clear fig handle defHandle;

numTimeSlices = size( levelSet, 2 );
diffusionTensors = NaN( 4, numTimeSlices );
porosities = NaN( numTimeSlices, 1 );

for ts = 1:numTimeSlices
    disp( [ 'Time step ', num2str(ts) ] );

    [ stiffnessMatrix, balanceVector, isDoF, triangleVolumes ] = assembleCellProblem( cellGrid, levelSet(:,ts) );
    cellProblemSolutions = solveSystemFE( cellGrid, stiffnessMatrix, balanceVector, isDoF );

    [ diffusion, porosity ] = computeDiffusionTensor( cellGrid, cellProblemSolutions, triangleVolumes );
    diffusionTensors(:,ts) = diffusion(:);
    porosities(ts) = porosity;
end
clear ts diffusion porosity;

folderName = [ 'data/IrregularSoil/Periodic/NumPartitions_', ...
    num2str( numPartitions ), '/' ];
dateString = datestr( now, 'yyyymmdd_HHMM' );
if ( toSave )
    filePath = [ folderName, 'Data_', dateString ];
    save( filePath );
    disp( [ 'Data saved to ', filePath ] );
end
% copyfile( filePath, [ folderName, 'Data' ] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function normalVelocity = speedFunctionTemplate( t, x )

    calciteSpeed = 0.5;
    dolomiteSpeed = 0.05;

    if ( any( abs(x) > 0.5 ) )
        normalVelocity = 0;
        return;
    end
    
    if ( x(1) < -0.15 && x(2) > 0.3 * x(1) + 0.475 )
        normalVelocity = dolomiteSpeed;
        return;
    end
    
    if ( x(1) > 0.1 && x(2) > 0.3 * x(1) + 0.4 )
        normalVelocity = dolomiteSpeed;
        return;
    end
    
    if ( x(1) < 0 && ( ( x(2) < 0.3 * x(1) + 0.4 ) && ( x(2) > 0.2 * x(1) + 0.3 ) ) )
        normalVelocity = dolomiteSpeed;
        return;
    end

    if ( ( x(2) > 0.15 * x(1) + 0.0 ) && ( x(2) < 0.3 * x(1) + 0.25 ) )
        normalVelocity = dolomiteSpeed;
        return;
    end
    
    if ( x(2) < ( 0.3 * ( x(1) < 0 ) + 0.05 * ( x(1) >= 0 ) ) * x(1) - 0.1 )
        if ( x(2) > ( -0.3 * ( x(1) < 0 ) + -0.25 * ( x(1) >= 0 ) ) )
            normalVelocity = dolomiteSpeed;
            return;
        end
    end
    
    if ( x(2) < ( 0.2 * ( x(1) < 0 ) ) * x(1) - 0.35 )
        normalVelocity = dolomiteSpeed;
        return;
    end
    
    normalVelocity = calciteSpeed;

end

function normalVelocity = speedFunction( t, x )

    y = x;
    y( x >= 0 ) = 2 * y( x >= 0 ) - 0.5;
    y( x < 0 ) = -2 * y( x < 0 ) - 0.5;
    
    normalVelocity = speedFunctionTemplate( t, y );

end