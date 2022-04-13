% Perform micro-scale simulation with periodic lateral bounadries
% Similar to [4] Section 5.3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% General setting variables

global EPS Solver;
Solver = 'ilukIterative';
EPS=eps;


lengthXAxis = 0.1 ; % [ dm ]
lengthYAxis = 0.05; % [ dm ]
numPartitionsMicroscale = 100 % Number of partitions 



tic;    % Preprocessing

% Physical parameters
dimension = 2;
spaceScaleFactor = 1; % spaceScaleFactor [length of Y] = 1 [cm]

% Parameters in dm and s
diffusionCoefficient = 0.2 * 1e-4 * 1e-2; % [ dm^2 s^(-1) ]
molarVolume = 64.12 * 1e-3; % [ dm^3 mol^(-1) ]
%rateCoefficientHydrogen = 0.89 * 1e-2; % [ mol dm^(-2) s^(-1) ]
rateCoefficientTST = 4.5e-4 * 1e-2; % [ mol dm^(-2) s^(-1) ]
equilibriumRateConstant = 10^(-16.5);
%inletVelocity = 0.1* 1e-1; % [ dm s^(-1) ]; Pe = 5000
%inletVelocity = 0.001 * 1e-1; % [ dm s^(-1) ]; Pe = 50
inletVelocity = 0.01 * 1e-1; % [ dm s^(-1) ]; Pe = 500
initialHydrogenConcentration =  1e-5; % [ mol dm^(-3) ]

dissolutionReactionRate = @( cHydrogen, cCalcium, cCarbonate, cMagnesium ) ...
    (cHydrogen^0.5 * rateCoefficientTST ) ...
    * ( 1 - cCalcium * cCarbonate^2 * cMagnesium / equilibriumRateConstant );

% Simulation and discretization parameters (comment for convergence tests)

numberOfSlices = 1;
disp( [ 'numSlices = ', num2str( numberOfSlices ) ] );
pecletNumber = ( inletVelocity * lengthXAxis ) / diffusionCoefficient;
fprintf( [ 'Peclet number: ', num2str( pecletNumber ), '\n' ] );

% dt_macro = 1 * 3600; % Pe = 5000
% initialMacroscaleTimeStepSize = 160; % Pe = 50
initialMacroscaleTimeStepSize =60;
% dt_micro = dt_macro / 100; % 0;
numMicroscaleTimeSteps = 10;
endTime = 3600* 600; % [h] Pe = 5000
% endTime = 10 * 60 * 60; % Pe = 50
% endTime = 200 * dt_macro;

%% Computation of time steps in simulation
% stepperType = 'linear';
timeStepperType = 'expmax';

switch ( timeStepperType )
    case 'linear'
        if ( mod( endTime, initialMacroscaleTimeStepSize ) < EPS )
            timeSteps = 0:initialMacroscaleTimeStepSize:endTime;
        else
            timeSteps = [ 0:initialMacroscaleTimeStepSize:endTime, endTime ];
        end
        timeStepSizeFactor = 1;
    case 'exp'
      %  timeStepSizeFactor = sqrt(2);
        timeStepSizeFactor = 2;
        numTimeSteps = floor( log( endTime / initialMacroscaleTimeStepSize ) ...
            / log( timeStepSizeFactor ) );
        timeSteps = [0; timeStepSizeFactor.^( 0:numTimeSteps )' ] ...
            * initialMacroscaleTimeStepSize;
        if ( timeSteps( end ) < endTime )
            timeSteps( end + 1 ) = endTime;
        end
    case 'expmax'
        maximalStep = 0.25*10^5;
        timeStepSizeFactor = 2;
        timeSteps = [0,initialMacroscaleTimeStepSize];
        while timeSteps(end) < endTime
            timeSteps = [timeSteps, min(timeSteps(end)*timeStepSizeFactor, timeSteps(end)+maximalStep)];
        end
        timeSteps = [timeSteps(1:(end-1)), endTime];
        size(timeSteps)
    otherwise
        error( 'Time stepper type not implemented.' );
end
numTimeSlices = numel( timeSteps );
levelSetEvolutionTime = NaN( numTimeSlices, 1 );
cellProblemTime = NaN( numTimeSlices, 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of level set method variables

microscaleGrid = FoldedCartesianGridCylinder( dimension, ...
    [0,lengthXAxis,0,lengthYAxis], ...
    numPartitionsMicroscale * [round(lengthXAxis/lengthYAxis),1] , 2);
coord = microscaleGrid.coordinates;
coordCell = mat2cell( coord, ones( 1, microscaleGrid.nodes ), dimension );


g = @(x) -1;
initialLevelSetFunc = @(x) g(x);
scale = 14;

initialLevelSetDataCells = cell( numberOfSlices, 1 );
initialLevelSetDataCells{1} = cellfun( initialLevelSetFunc, coordCell );

f = @(x,a) 0.0141/scale - norm(x-a,2);
for i=1:2*scale
    i
    for j=1:(scale)
	shift = lengthYAxis/scale/3;
        m = lengthYAxis/scale*(i-0.5);
        n = lengthYAxis/scale*(j-0.5) -shift/2 + mod(i,2)*shift;
        
        indizes = logical(indizesF([m,n], coord, lengthYAxis/scale));
        sub = mat2cell( coord(indizes,:), ones(1,sum(indizes)), dimension );
        initialLevelSetDataCells{1}(indizes,1) = max(cellfun(@(x) f(x, [m, n]), sub), initialLevelSetDataCells{1}(indizes,1));     

	if j==1
 	n1 = lengthYAxis/scale*(scale-0.5) -shift/2 + mod(i,2)*shift;
	n2 = lengthYAxis/scale*(scale+1-0.5) -shift/2 + mod(i,2)*shift;
	indizes = logical(indizesF([m,n1], coord, lengthYAxis/scale));
        sub = mat2cell( coord(indizes,:), ones(1,sum(indizes)), dimension );
        initialLevelSetDataCells{1}(indizes,1) = max(cellfun(@(x) f(x, [m, n2]), sub), initialLevelSetDataCells{1}(indizes,1)); 
	end   

	if j==scale
 	n1 = lengthYAxis/scale*(1-0.5) -shift/2 + mod(i,2)*shift;
	n2 = lengthYAxis/scale*(0-0.5) -shift/2 + mod(i,2)*shift;
	indizes = logical(indizesF([m,n1], coord, lengthYAxis/scale));
        sub = mat2cell( coord(indizes,:), ones(1,sum(indizes)), dimension );
        initialLevelSetDataCells{1}(indizes,1) = max(cellfun(@(x) f(x, [m, n2]), sub), initialLevelSetDataCells{1}(indizes,1)); 
	end   
        
    end
end

[a,b]=meshgrid(0:lengthYAxis/numPartitionsMicroscale:lengthXAxis,0:lengthYAxis/numPartitionsMicroscale:lengthYAxis);
        contour(a,b,reshape(initialLevelSetDataCells{1},2*numPartitionsMicroscale + 1 ,numPartitionsMicroscale+1)',[0,1])
        axis equal

interfaceNormalVelocity = @(cHydrogen, cCalcium, cCarbonate, cMagnesium ) ...
    spaceScaleFactor * molarVolume * dissolutionReactionRate( cHydrogen, ...
    cCalcium, cCarbonate, cMagnesium );

currentTime = 0;

currentLevelSetDataCells = cell( numberOfSlices, 1 );
[ currentLevelSetDataCells{:} ] = deal( initialLevelSetDataCells{1} );
oldLevelSetDataCells = currentLevelSetDataCells;
% TODO Viel zu viele Spalten in levelSetData
% prev: ( ceil( endTime / dt_macro ) + 1 ) columns
levelSetData = NaN( numel( initialLevelSetDataCells{1} ), numTimeSlices );
levelSetData(:,1) = currentLevelSetDataCells{1};
levelSet = cell( numberOfSlices, 1 );
[ levelSet{:} ] = deal( levelSetData );
clear levelSetData;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gridHyPHM =Grid(coord,microscaleGrid.triangles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1 : gridHyPHM.numE
    if (all(gridHyPHM.coordV(gridHyPHM.V0E(i,:),1)<eps) )
        gridHyPHM.idE(i)=4;
    end
    if (all(gridHyPHM.coordV(gridHyPHM.V0E(i,:),1)>lengthXAxis-eps) )
        gridHyPHM.idE(i)=2;
    end
    if (all(gridHyPHM.coordV(gridHyPHM.V0E(i,:),2)<eps) )
        gridHyPHM.idE(i)=1;
    end
    if (all(gridHyPHM.coordV(gridHyPHM.V0E(i,:),2)>lengthYAxis-eps) )
        gridHyPHM.idE(i)=3;
    end
end


gridHyPHM = FoldedGrid(gridHyPHM,'Cylinder');
macroCoordCell = mat2cell( gridHyPHM.baryT, ones( gridHyPHM.numT, 1 ), 2 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of initial effective parameters

diffusionTensors = cell( numberOfSlices, 1 );
[ diffusionTensors{:} ] = deal( NaN( 4, numTimeSlices ) );
porosities = cell( numberOfSlices, 1 );
[ porosities{1:end} ] = deal( NaN( numTimeSlices, 1 ) );
clear numTimeSlices;

surfaceArea = cell( numberOfSlices, 1 );
PorositiyEstimate = sum(levelSet{1}(:,1)<-eps)/numel(levelSet{1}(:,1))

for i = 1:numberOfSlices
    
    [ cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces ] ...
        = assembleCellProblem( microscaleGrid, levelSet{i}(:,1) );
  
    surfaceArea{i}(1) = sum( triangleSurfaces );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of HyPHM variables (cont.)

flowStepper = Stepper(0:1);
transportStepper = Stepper( timeSteps );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stokes
    %Calculate Pressure distribution and Stokes velocity
    disp( ' ' );
    disp( 'Initializing Stokes problem...' );
   
    reducedLevel = levelSet{1}(~microscaleGrid.isFoldedNode,1);
    phead = Variable(gridHyPHM, flowStepper, 'pressure head', 'P1');
    StokesVelocity = Variable(gridHyPHM, flowStepper, 'Stokes Velocity', 'P2P2');
    phead.setdata(0, @(t, x) 0.0);
    StokesVelocity.setdata(0, @(t, x) 0.0);

    StokesL = StokesLEVEL(gridHyPHM, flowStepper, 'Stokes problem');
    StokesL.L.setdata(reducedLevel)
    StokesL.id2D = {4,3,1};
    StokesL.uD.setdata(@(t, x) inletVelocity*(x(1)<EPS)*[1;0] );
    StokesL.F.setdata(@(t, x) 0);
    StokesL.U = StokesVelocity;
    StokesL.P = phead;

    flowStepper.next;
    StokesL.computeLevel('s');
    flowStepper.prev;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
flow = Variable( gridHyPHM, flowStepper, 'Flow', 'RT0' );
helper = StokesL.U.getdata(1);
flow.setdata(helper(gridHyPHM.numV+1:end, 1) .* gridHyPHM.nuE(:, 1) + helper(gridHyPHM.numV+1:end, 2) .* gridHyPHM.nuE(:, 2));
porosityFunc = @(t,x) porosityHelperFun( t, x, porosities, numberOfSlices, 1 );


% Boundary IDs: 1 = down, 2 = right, 3 = up, 4 = left

disp( ' ' );
disp( 'Initializing calcium transport...' );


calciumConcentration = Variable( gridHyPHM, transportStepper, ...
    'Ca^(2+)', 'P0' );
calciumConcentration.setdata( 0, @(t,x) 0 );
% calciumConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
calciumTransport = TransportLEVEL( gridHyPHM, transportStepper, 'Ca^(2+) Transport' );
% calciumTransport.id2D = {2};
calciumTransport.D.setdata( diffusionCoefficient*eye(2));
calciumTransport.id2N = {1,2,3};
calciumTransport.id2F = {4};
calciumTransport.U = calciumConcentration;
% calciumTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% calciumTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
calciumTransport.gF = Variable( gridHyPHM, transportStepper, 'gF', 'P0E' );
calciumTransport.gF.setdata(0, zeros( gridHyPHM.numE, 1 ) );
calciumTransport.A  = Variable( gridHyPHM, transportStepper, 'Porosity', 'C' );
calciumTransport.A.setdata( 1  );
% calciumTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
calciumTransport.C.setdata(0, flow.getdata(1) );
calciumTransport.isUpwind = 'exp';

disp( ' ' );
disp( 'Initializing carbonate transport...' );

carbonateConcentration = Variable( gridHyPHM, transportStepper, ...
    'CO_3^(2-)', 'P0' );
carbonateConcentration.setdata( 0, @(t,x) 0 );
% carbonateConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
carbonateTransport = TransportLEVEL( gridHyPHM, transportStepper, 'CO_3^(2-) Transport' );
% carbonateTransport.id2D = {2};
carbonateTransport.D.setdata( diffusionCoefficient*eye(2));
carbonateTransport.id2N = {1,2,3};
carbonateTransport.id2F = {4};
carbonateTransport.U = carbonateConcentration;
% carbonateTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% carbonateTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
carbonateTransport.gF = Variable( gridHyPHM, transportStepper, 'gF', 'P0E' );
carbonateTransport.gF.setdata(0, zeros( gridHyPHM.numE, 1 ) );
carbonateTransport.A = Variable( gridHyPHM, transportStepper, 'Porosity', 'C' );
carbonateTransport.A.setdata(  1  );

% carbonateTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
carbonateTransport.C.setdata( flow.getdata(1) );
carbonateTransport.isUpwind = 'exp';

disp( ' ' );
disp( 'Initializing magnesium transport...' );

magnesiumConcentration = Variable( gridHyPHM, transportStepper, ...
    'Mg^(2+)', 'P0' );
magnesiumConcentration.setdata( 0, @(t,x) 0 );
% magnesiumConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
magnesiumTransport = TransportLEVEL( gridHyPHM, transportStepper, 'MG^(2+) Transport' );
% magnesiumTransport.id2D = {2};
magnesiumTransport.D.setdata( diffusionCoefficient*eye(2));
magnesiumTransport.id2N = {1,2,3};
magnesiumTransport.id2F = {4};
magnesiumTransport.U = magnesiumConcentration;
% magnesiumTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% magnesiumTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
magnesiumTransport.gF = Variable( gridHyPHM, transportStepper, 'gF', 'P0E' );
magnesiumTransport.gF.setdata(0, zeros( gridHyPHM.numE, 1 ) );
magnesiumTransport.A = Variable( gridHyPHM, transportStepper, 'Porosity', 'C' );
magnesiumTransport.A.setdata( 1  );
% magnesiumTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
magnesiumTransport.C.setdata(0, flow.getdata(1) );
magnesiumTransport.isUpwind = 'exp';
disp( ' ' );
disp( 'Initializing hydrogen transport...' );

hydrogenConcentration = Variable( gridHyPHM, transportStepper, ...
    'H^+_pH', 'P0' );
hydrogenConcentration.setdata( 0, @(t,x) initialHydrogenConcentration );
hydrogenTransport = TransportLEVEL( gridHyPHM, transportStepper, 'H^+ Transport' );
hydrogenTransport.id2N = {1,2,3};
hydrogenTransport.id2F = {4};
% hydrogenTransport.id2D = {4};
hydrogenTransport.U = hydrogenConcentration;
hydrogenTransport.D.setdata( diffusionCoefficient*eye(2));
hydrogenTransport.gF.setdata( ...
    @(t,x) - initialHydrogenConcentration * inletVelocity * ( x(1) < EPS ) );

% hydrogenTransport.gF.setdata( @(t,x) -1e-7 );
% hydrogenTransport.uD.setdata( @(t,x) initialHydrogenConcentration );
hydrogenTransport.A = Variable( gridHyPHM, transportStepper, 'Porosity', 'C' );
hydrogenTransport.A.setdata(  1 );
hydrogenTransport.C.setdata(0, flow.getdata( 1 ) );
hydrogenTransport.isUpwind = 'exp';

speciesCells = { calciumTransport; carbonateTransport; hydrogenTransport; magnesiumTransport };

% calciumDataOld = calciumConcentration.getdata( 0 );
% carbonateDataOld = carbonateConcentration.getdata( 0 );
% hydrogenDataOld = hydrogenConcentration.getdata( 0 );

preprocessingTime = toc;     % Preprocessing

disp( ' ' );
disp( [ 'Preprocessing done in ', num2str( preprocessingTime ), ' seconds.' ] );

hydrogenTransport.L.setdata(0, reducedLevel);
magnesiumTransport.L.setdata(0, reducedLevel);
carbonateTransport.L.setdata(0, reducedLevel);
calciumTransport.L.setdata(0, reducedLevel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time iteration

%while ( t_ <= endTime + EPS )
while transportStepper.next
    timeIterationStep = transportStepper.curstep;
    macroscaleTimeStepSize = transportStepper.curtau;
    disp( '----------------------------------------' );
    disp( [ 'Time step ', num2str( timeIterationStep ) ] );
    disp( [ '  Current time = ', ...
        num2str( currentTime ) ] );
    disp( [ '  Current time step size = ', ...
        num2str( macroscaleTimeStepSize ) ] );

    if ( currentTime + macroscaleTimeStepSize > endTime - EPS )
        macroscaleTimeStepSize = endTime - currentTime;
        if ( abs( macroscaleTimeStepSize ) < EPS )
            break;
        end
    end
    
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   inner corrector
    if timeIterationStep>0
        speciesCells = { calciumTransport; carbonateTransport; hydrogenTransport; magnesiumTransport };
        [out, Corrector] = Continuation(speciesCells, timeIterationStep,levelSet{1}(:,timeIterationStep),20);
        calciumTransport = out{1};
        carbonateTransport = out{2};
        hydrogenTransport = out{3};
        magnesiumTransport = out{4};
    else
       corrector = zeros(gridHyPHM.numT,4); 
    end      
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level set evolution step

    disp( 'Evolution of level set ...' );
    tic;
    coordCell = mat2cell( microscaleGrid.coordinates, ...
        ones( 1, microscaleGrid.nodes ), dimension );
    timeCell = num2cell( currentTime * ones( microscaleGrid.nodes, 1 ) );
    
    
    calciumData = calciumConcentration.getdata( timeIterationStep - 1 );
    carbonateData = carbonateConcentration.getdata( timeIterationStep - 1 );
    magnesiumData = magnesiumConcentration.getdata( timeIterationStep - 1 );
    hydrogenData = hydrogenConcentration.getdata( timeIterationStep - 1 );
 
    
    
     %Convert Concentrations from T --> V
    calciumDataV=zeros(gridHyPHM.numV,1);
    magnesiumDataV=zeros(gridHyPHM.numV,1);
    hydrogenDataV=zeros(gridHyPHM.numV,1);
    carbonateDataV=zeros(gridHyPHM.numV,1);
       
    
     normalSpeed = arrayfun( interfaceNormalVelocity, ...
            hydrogenData  , ...
            calciumData , ...
            carbonateData , ...
            magnesiumData   );
    
    
    for i = 1:numberOfSlices
%         normalSpeed = cellfun( interfaceNormalVelocity, timeCell, coordCell, ...
%             num2cell( calciumSliceAverageCells{i} * ones( microscaleGrid.nodes, 1 ) ), ...
%             num2cell( carbonateSliceAverageCells{i} * ones( microscaleGrid.nodes, 1 ) ) );
        
%         t_micro = currentTime;
%         for j = 1:ceil( dt_macro / dt_micro )
% %             disp( ['    Level set substep ', num2str(j) ] );
%             currentLevelSetDataCells{i} = levelSetEquationTimeStep( t_micro + dt_micro, currentLevelSetDataCells{i}, ...
%                 t_micro, oldLevelSetDataCells{i}, microscaleGrid, normalSpeed );
%             
%             t_micro = t_micro + dt_micro;
%             oldLevelSetDataCells = currentLevelSetDataCells;
%             
%         end


        [~,extensionVelocity] = Phase2Extension(microscaleGrid, gridHyPHM, levelSet{1}(:,timeIterationStep), normalSpeed);
        extensionVelocity = microscaleGrid.synchronizeValues(extensionVelocity)';

	CFL = 1*1/3* 1/numPartitionsMicroscale * lengthYAxis / max(abs(extensionVelocity)) 
         microscaleTimeStepSize = transportStepper.curtau / numMicroscaleTimeSteps;
        microscaleTime = currentTime;
        oldMicroscaleTime = currentTime;
        localIterates = ceil(transportStepper.curtau/CFL)

        for j = 1:localIterates
%             disp( ['    Level set substep ', num2str(j) ] );

            if j<localIterates
		jump = CFL;
	    else
		jump = transportStepper.curtau - (localIterates-1)*CFL;
	    end
            % [] argument is unused in method (needed for implicit methods)
         %   currentLevelSetDataCells{i} = levelSetEquationTimeStep( ...
         %       newMicroscaleTime, ...
         %       oldMicroscaleTime, oldLevelSetDataCells{i}, microscaleGrid, expandNodes(extensionVelocity,microscaleGrid),1 );
            currentLevelSetDataCells{i} = levelSetEquationTimeStep( ...
                jump, ...
                0, oldLevelSetDataCells{i}, microscaleGrid, extensionVelocity,1 );
            oldLevelSetDataCells = currentLevelSetDataCells;
            
        end
        
        levelSet{i}( :, timeIterationStep + 1 ) = currentLevelSetDataCells{i};
    end
    
%     t_old = currentTime;
%     currentTime = currentTime + dt_macro;
%     oldLevelSetDataCells = currentLevelSetDataCells;
    currentTime = currentTime + macroscaleTimeStepSize;
    
    levelSetEvolutionTime( timeIterationStep ) = toc;
    disp( [ '    ... done in ', ...
        num2str( levelSetEvolutionTime( timeIterationStep ) ), ' seconds.' ] );
    
       [ cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces ] ...
        = assembleCellProblem( microscaleGrid, levelSet{1}(:, timeIterationStep + 1) );
    porositiesMICRO(timeIterationStep)=1-sum(triangleVolumes)/lengthXAxis/lengthYAxis;

     reducedLevel = levelSet{1}(~microscaleGrid.isFoldedNode,timeIterationStep+1);	
    hydrogenTransport.L.setdata(timeIterationStep, reducedLevel);
    magnesiumTransport.L.setdata(timeIterationStep, reducedLevel);
    carbonateTransport.L.setdata(timeIterationStep, reducedLevel);
    calciumTransport.L.setdata(timeIterationStep, reducedLevel);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculation of effective parameters
    
    disp( 'Calculation of effective parameters ...' );
    tic;
%     
%     for i = 1:numberOfSlices
%         
%         [ cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, ...
%             triangleSurfaces ] = assembleCellProblem( microscaleGrid, ...
%             levelSet{i}( :, timeIterationStep ) );
%       
%         surfaceArea{i}( timeIterationStep ) = sum( triangleSurfaces );
%                
%     end
    cellProblemTime( timeIterationStep ) = toc;
    disp( [ '    ... done in ', num2str( cellProblemTime( timeIterationStep ) ), ...
        ' seconds.' ] );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% Stokes
    %Calculate Pressure distribution and Stokes velocity
    disp( ' ' );
    disp( 'Initializing Stokes problem...' );
   
    if sum(abs((levelSet{i}( :, timeIterationStep + 1 )>0)-(levelSet{i}( :, timeIterationStep )>0)))>0.5 
    	StokesL.L.setdata(reducedLevel)

    	flowStepper.next;
    	StokesL.computeLevel('s');
    
    	helper = StokesL.U.getdata(1);
    	flow.setdata(helper(gridHyPHM.numV+1:end, 1) .* gridHyPHM.nuE(:, 1) + helper(gridHyPHM.numV+1:end, 2) .* gridHyPHM.nuE(:, 2));
    	flowStepper.prev;
   else
	disp('Stokes skipped')
   end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Macroscopic transport step
    
     
    
%     surfaceAreaFun = @(x) areaHelperFun( -1, x, microscaleGrid, levelSet, numberOfSlices, timeStep );
    surfaceAreaFunc = @(x) areaHelperFun( -1, x, surfaceArea, numberOfSlices, ...
        timeIterationStep );

%    specificSurfaceArea = cellfun( surfaceAreaFunc, macroCoordCell ) * spaceScaleFactor;
    % TODO evaluate dissolutionReactionRate and other related function
    % handles (really needed?)
    Levels=hydrogenTransport.L.getdata(timeIterationStep);
    %%%%%%
    % handle edge orientation
    edgeOrientation = zeros(gridHyPHM.numE,1);
  %  for kt = 1: gridHyPHM.numT
  %      if sum(Levels(gridHyPHM.V0T(kt,:))>-eps)>2.5                        %solid triangle
  %          edgeOrientation(gridHyPHM.E0T(kt,:)) = gridHyPHM.sigE0T(kt,:);
  %      end
  %  end    
        
    
    
    %%%%%%
    
    
    %concentrations P_0(T)--> P_0(E)
    calciumDataE = zeros(gridHyPHM.numE,1);
    magnesiumDataE = zeros(gridHyPHM.numE,1);
    hydrogenDataE = zeros(gridHyPHM.numE,1);
    carbonateDataE = zeros(gridHyPHM.numE,1);
    
        
    %%%%%%
    calciumTransportRhsData = 0 ...
        .* ( sqrt(hydrogenDataE) * rateCoefficientTST ) .* ...
        ( 1 - calciumDataE .* carbonateDataE .*carbonateDataE .*magnesiumDataE ./ equilibriumRateConstant ).* ...
        edgeOrientation;
    
    carbonateTransportRhsData = 0 * calciumTransportRhsData;
    magnesiumTransportRhsData = 0 * calciumTransportRhsData;
    
    hydrogenTransportRhsData = 0 * calciumTransportRhsData;
    
 
    
    calciumTransport.gF.setdata( timeIterationStep, calciumTransportRhsData );
    calciumTransport.C.setdata( timeIterationStep, flow.getdata(1) );
    carbonateTransport.gF.setdata( timeIterationStep, carbonateTransportRhsData );
    carbonateTransport.C.setdata( timeIterationStep, flow.getdata(1));
    hydrogenTransport.gF.setdata( timeIterationStep, hydrogenTransport.gF.getdata(0) + hydrogenTransportRhsData );
    hydrogenTransport.C.setdata( timeIterationStep, flow.getdata(1) );
    magnesiumTransport.gF.setdata( timeIterationStep, magnesiumTransportRhsData );
    magnesiumTransport.C.setdata( timeIterationStep, flow.getdata(1) );
   
 
    SorceScale = zeros(gridHyPHM.numT,1);
    for kT=1:gridHyPHM.numT
        L=0;
        localCount = 0;
        for i=1:3
            if Levels(gridHyPHM.V0E(gridHyPHM.E0T(kT,i),1))>-eps & Levels(gridHyPHM.V0E(gridHyPHM.E0T(kT,i),2))>-eps
                L=L+gridHyPHM.areaE(gridHyPHM.E0T(kT,i)) / gridHyPHM.areaT(kT);
                localCount =  localCount +1;
            end
        end
	if localCount>0
           SorceScale(kT)=L/sqrt( localCount );
        end
    end
    
EpsME = 10*eps; 
    nonlinearFunc = cell( 4, 1 );
    nonlinearFunc{1} = @(x,y,z,w) 1 ...
        .* ( sqrt(z) * rateCoefficientTST ) ...
        .* ( 1 - x .* y .* y .* w./ equilibriumRateConstant ).* SorceScale;
    nonlinearFunc{2} = @(x,y,z,w) 2 * nonlinearFunc{1}( x, y, z, w);
    nonlinearFunc{3} = @(x,y,z,w) 0; 
    nonlinearFunc{4} = @(x,y,z,w) nonlinearFunc{1}( x, y, z, w);
    
    nonlinearJacFunc = cell( 4, 4 );
    nonlinearJacFunc{1,1} = @(x,y,z,w) 1 ...
        .* (sqrt(z) * rateCoefficientTST ) ...
        .* ( -y .* y .* w ./ equilibriumRateConstant ).* SorceScale;;
    nonlinearJacFunc{1,2} = @(x,y,z,w) 1 ...
        .* ( sqrt(z)*rateCoefficientTST ) ...
        .* ( - 2 * x .*y .* w ./ equilibriumRateConstant ).* SorceScale;;
    nonlinearJacFunc{1,3} = @(x,y,z,w) 1 ...
        .* ( 0.5* 1./ sqrt(max(z,EpsME))*rateCoefficientTST ) ...
        .* ( 1 - x .* y .* y .* w ./ equilibriumRateConstant ).* SorceScale;;
    nonlinearJacFunc{1,4} = @(x,y,z,w) 1 ...
        .* ( sqrt(z)*rateCoefficientTST ) ...
        .* ( - x .*y .* y ./ equilibriumRateConstant ).* SorceScale;;
    
    nonlinearJacFunc{2,1} = @(x,y,z,w) 2 * nonlinearJacFunc{1,1}( x, y, z, w );
    nonlinearJacFunc{2,2} = @(x,y,z,w) 2 * nonlinearJacFunc{1,2}( x, y, z, w );
    nonlinearJacFunc{2,3} = @(x,y,z,w) 2 * nonlinearJacFunc{1,3}( x, y, z, w );
    nonlinearJacFunc{2,4} = @(x,y,z,w) 2 * nonlinearJacFunc{1,4}( x, y, z, w );
    
    nonlinearJacFunc{3,1} = @(x,y,z,w) 0; 
    nonlinearJacFunc{3,2} = @(x,y,z,w) 0; 
    nonlinearJacFunc{3,3} = @(x,y,z,w) 0; 
    nonlinearJacFunc{3,4} = @(x,y,z,w) 0; 
    
    nonlinearJacFunc{4,1} = @(x,y,z,w) nonlinearJacFunc{1,1}( x, y, z, w);
    nonlinearJacFunc{4,2} = @(x,y,z,w) nonlinearJacFunc{1,2}( x, y, z, w);
    nonlinearJacFunc{4,3} = @(x,y,z,w) nonlinearJacFunc{1,3}( x, y, z, w);
    nonlinearJacFunc{4,4} = @(x,y,z,w) nonlinearJacFunc{1,4}( x, y, z, w );

    speciesCells = { calciumTransport; carbonateTransport; hydrogenTransport; magnesiumTransport };
    newtonIteration( speciesCells, nonlinearFunc, nonlinearJacFunc, 20 );
    
  %  calciumTransport.U.setdata(timeIterationStep-1, calciumTransport.U.getdata(timeIterationStep-1)-Corrector(:,1));
  %  carbonateTransport.U.setdata(timeIterationStep-1, carbonateTransport.U.getdata(timeIterationStep-1)-Corrector(:,2));
   % hydrogenTransport.U.setdata(timeIterationStep-1, hydrogenTransport.U.getdata(timeIterationStep-1)-Corrector(:,3));
  %  magnesiumTransport.U.setdata(timeIterationStep-1, magnesiumTransport.U.getdata(timeIterationStep-1)-Corrector(:,4));
    
    %  hydrogenTransport.U.setdata(timeIterationStep, max(hydrogenTransport.U.getdata(timeIterationStep),10*eps));
    %  carbonateTransport.U.setdata(timeIterationStep, max(carbonateTransport.U.getdata(timeIterationStep),10*eps));
     % magnesiumTransport.U.setdata(timeIterationStep, max(magnesiumTransport.U.getdata(timeIterationStep),10*eps));
    %  calciumTransport.U.setdata(timeIterationStep, max(calciumTransport.U.getdata(timeIterationStep),10*eps));
     

end % while

save('CylinderImproved2.mat','-v7.3')
save('CalciumOnlyImproved2.mat','calciumTransport','porositiesMICRO','-v7.3')

%calciumConcentration.visualize;
% carbonateConcentration.visualize;
% hydrogenTransport.U.visualize;

% dataFolder = 'data/Molins2017/ConvergenceTests_201805/MacroGridConvergence/Newton/';

% filePrefix = [ dataFolder, 'Physics_Data_', ...
%     '_partsMacro_Peclet_', num2str( pecletNumber ) ];

% dateString = datestr( now, 'yyyymmdd_HHMM' );
% try
%     save( [ dataFolder, 'Physics_Data_', dateString, '.mat' ], ...
%         'dimension', 'spaceScaleFactor', 'diffusionCoefficient', 'molarVolume', ...
%         'rateCoefficientHydrogen', 'rateCoefficientTST', 'equilibriumRateConstant', ...
%         'inletVelocity', 'initialHydrogenConcentration', 'dissolutionReactionRate', ...
%         'lengthXAxis', 'lengthYAxis', 'pecletNumber', 'endTime', '-v7.3' );
% catch
%     disp( 'Physics data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Level_set_', dateString, '.mat' ], ...
%         'initialLevelSetFunc', 'interfaceNormalVelocity', 'levelSet', '-v7.3' );
% catch
%     disp( 'Level set data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Simulation_data_', dateString, '.mat' ], ...
%         'numPartitionsMicroscale', 'macroscaleStepSize', 'numberOfSlices', ...
%         'pecletNumber', 'initialMacroscaleTimeStepSize', 'numMicroscaleTimeSteps', ...
%         'timeStepperType', 'timeSteps', 'timeStepSizeFactor', 'microscaleGrid', ...
%         'gridHyPHM', 'flowStepper', 'transportStepper', '-v7.3' );
% catch
%     disp( 'Simulation data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'General_results_', dateString, '.mat' ], ...
%         'diffusionTensors', 'porosities', 'surfaceArea', 'flow', ...
%         'levelSetEvolutionTime', 'cellProblemTime', 'preprocessingTime', '-v7.3' );
% catch
%     disp( 'General results could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Calcium_', dateString, '.mat' ], ...
%         'calciumConcentration', 'calciumTransport', '-v7.3' );
% catch
%     disp( 'Calcium data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Carbonate_', dateString, '.mat' ], ...
%         'carbonateConcentration', 'carbonateTransport', '-v7.3' );
% catch
%     disp( 'Carbonate data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Hydrogen_', dateString, '.mat' ], ...
%         'hydrogenConcentration', 'hydrogenTransport', '-v7.3' );
% catch
%     disp( 'Hydrogen data could not be saved.' );
% end
% 
% diary off;
% movefile( logFilePath, strcat( dataFolder, 'Log_', dateString ) );

% save( [ filePrefix, '_Calcium_', dateString, '.mat' ], ...
%     'calciumConcentration', 'calciumData', 'calciumDataOld', ...
%     'calciumSliceAverageCells', 'calciumTransport', ...
%     'calciumTransportRhsData' );
% save( [ filePrefix, '_Carbonate_', dateString, '.mat' ], ...
%     'carbonateConcentration', 'carbonateData', 'carbonateDataOld', ...
%     'carbonateSliceAverageCells', 'carbonateTransport', ...
%     'carbonateTransportRhsData' );
% save( [ filePrefix, '_Hydrogen_', dateString, '.mat' ], ...
%     'hydrogenConcentration', 'hydrogenData', 'hydrogenDataOld', ...
%     'hydrogenSliceAverageCells', 'hydrogenTransport', ...
%     'hydrogenTransportRhsData' );
% save( [ dataFolder, 'Data_pH_', num2str( 1 / macroscaleStepSize ), ...
%     '_partsMacro_LevelSet_Peclet_50_', dateString, '.mat' ], ...
%     'levelSet', '-v7.3' );
% save( [ filePrefix, '_Other_', dateString, '.mat' ], ...
%     'cellProblemSystemMatrix', 'cellProblemTime', 'coord', 'coordCell', ...
%     'diffusion', 'diffusionCoefficient', 'diffusionFun', ...
%     'diffusionTensors', 'diffusionTensorsData', 'dimension', ...
%     'dissolutionReactionRate', 'dt_macro', 'dt_micro', 'endTime', ...
%     'equilibriumRateConstant', 'flow', 'flowStepper', 'gridHyPHM', ...
%     'initialData', 'initialHydrogenConcentration', 'inletVelocity', ...
%     'interfaceNormalVelocity', 'isDoF', 'isInArea', 'lengthXAxis', ...
%     'lengthYAxis', 'levelSetData', 'levelSetEvolutionTime', ...
%     'lsf0', 'macroCoordCell', 'macroscaleStepSize', 'microscaleGrid', ...
%     'molarVolume', 'normalSpeed', 'numberOfSlices', ...
%     'numPartitionsMicroscale', 'numTimeSlices', 'phi', 'phiOld', ...
%     'porosities', 'porositiesData', 'porosity', 'porosityFun', ...
%     'preprocessingTime', 'rateCoefficientHydrogen', ...
%     'rateCoefficientTST', 'rateFunction', 'rhs', 'SOL', ...
%     'spaceScaleFactor', 'specificSurfaceArea', 'step', 'surfaceArea', ...
%     'surfaceAreaFun', 't_', 't_micro', 't_old', 'timeCell', 'timeStep', ...
%     'timeSteps', 'transportStepper', 'triangleSurfaces', ...
%     'triangleVolumes', 'velocity', 'X', 'Y' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions

function por = porosityHelperFun( t, x, porosityCells, numberOfSlices, timeStep )

    por = 0;
    for i = 1:numberOfSlices
        por = por + porosityCells{1}(timeStep) ;
    end

%     por = porosityCells{ findSlice( x, numberOfSlices ) }( timeStep );

end

function diff = diffusionHelperFun( t, x, diffusionCells, numberOfSlices, timeStep )

    diff = zeros(2);
    for i = 1:numberOfSlices
        diff = diff + reshape( diffusionCells{1}(:,timeStep), 2, 2 ) ;
    end

%     diff = reshape( diffusionCells{ findSlice( x, numberOfSlices ) }( :, ...
%         timeStep ), 2, 2 ); 

end

function area = areaHelperFun( t, x, surfaceAreaCells, numberOfSlices, timeStep )

    area = surfaceAreaCells{ findSlice( x, numberOfSlices ) }( timeStep );

end

function out = expandNodes(val, grid)
	out = zeros(numel(grid.isFoldedNode),1);
	out(~grid.isFoldedNode ) = val;
	out( grid.isFoldedNode ) = out( grid.foldedIndex( grid.isFoldedNode ) );
end

% function area = areaHelperFun( t, x, grid, levelSetCells, numberOfSlices, timeStep )
% 
% %     area = 0;
% %     for i = 1:numberOfSlices
% %         area = area + surfaceIntegral( grid, levelSetCells{i}( :, timeStep ), ...
% %             @(x) 1 ) * ...
% %             ( x(1) >= (i-1)/numberOfSlices & x(1) < i/numberOfSlices );
% %     end
% 
%     area = surfaceIntegral( grid, ...
%         levelSetCells{ findSlice( x, numberOfSlices ) }( :, timeStep ), @(x) 1 );
% 
% end

function sliceNumber = findSlice( x, numberOfSlices )

    sliceNumber = 1 ;

end

function timeStep = findTimeStep( t, deltaTime )

    timeStep = 1 + fix( t / deltaTime );

end

function out = indizesF(in, coords, scale)
    out = logical(abs(coords(:,1)-in(1))<2*scale& abs(coords(:,2)-in(2))<2*scale);
end
