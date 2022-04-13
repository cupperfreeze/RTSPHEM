% Simple micro-macro simulation with multiple underlying unit-cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General setting variables

global EPS

dim = 2;
n = 32; % Number of partitions in each direction
endTime = 1 - 1 / n;
dt = 0.5 / n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of level set method variables

g = FoldedCartesianGrid(dim, kron(ones(1, dim), [-0.5, 0.5]), n*ones(1, dim));
coord = g.coordinates;
numNodes = g.nodes;

coordCell = mat2cell(coord, ones(1, g.nodes), dim);

lsf0 = @(x) sqrt((1 - 0.76)/pi) - norm(x);
initialData0 = cellfun(lsf0, coordCell);
initialData1 = initialData0;
initialData2 = initialData0;
% lsf1 = @(x) max( ...
%     min( min( x(1) + 0.25, x(2) + 0.25 ), 0.25 - norm( x + 0.25, 1 ) ), ...
%     min( min( -(x(1) - 0.25), -(x(2) - 0.25) ), 0.25 - norm( -(x - 0.25), 1 ) ) );
% initialData1 = cellfun( lsf1, coordCell );
% lsf2 = @(x) max( 0.3 - norm( [2*x(1), x(2)] -[0, -0.5], 5 ), ...
%     0.3 - norm( [2*x(1), x(2)] - [0, 0.5], 5 ) );
% initialData2 = cellfun( lsf2, coordCell );

X = g.reshape(coord(:, 1));
Y = g.reshape(coord(:, 2));
plotInitial2 = g.reshape(initialData2);

stopRadius = 0.15 + g.stepSize(1);

speedFun0 = @(t, x, c) 0.2 * c;
speedFun1 = @(t, x, c) 0.2 * 2 * c;
speedFun2 = @(t, x, c) 0.2 * 4 * c;
% speedFun2 = @(t,x) -0.2 * ( abs( x(1) ) < stopRadius ) + ...
%     0.2 * ( abs( x(1) ) >= stopRadius );

% numFunctions = 22;
% speedFun = cell( numFunctions,1 );
% for i = 1:numFunctions
%     speedFun{i} = @(t,x) -0.2*(i/numFunctions);
% end

velocity = zeros(g.nodes, dim);

t_ = 0;
phi0 = initialData0(:);
phiOld0 = phi0;
levelSet0 = NaN(numel(initialData0), ceil(endTime / dt)+1);
levelSet0(:, 1) = phi0;
phi1 = initialData1(:);
phiOld1 = phi1;
levelSet1 = NaN(numel(initialData1), ceil(endTime / dt)+1);
levelSet1(:, 1) = phi1;
phi2 = initialData2(:);
phiOld2 = phi2;
levelSet2 = NaN(numel(initialData2), ceil(endTime / dt)+1);
levelSet2(:, 1) = phi2;
timeStep = 1;

% phi = repmat( initialData0(:), 1, numFunctions );
% phiOld = phi;
% levelSet = NaN( numel( initialData0 ), ceil( endTime / dt ) + 1, numFunctions );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting of initial interface

fig = figure;
[~, defHandle] = contour(X, Y, plotInitial2, zeros(1, 2), 'LineColor', 'r', 'LineWidth', 3);
% set( gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'fontsize', 20);
axis([-0.5, 0.5, -0.5, 0.5])
axis equal
xlabel('x');
ylabel('y');
hand = copyobj(defHandle, gca);
set(hand, 'Zdata', plotInitial2, 'LineColor', 'b', 'LineWidth', 3);

% print( fig, '-dpng', 'data/MultipleCellProblems/Triangle0.png' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of initial effective parameters

numTimeSlices = size(levelSet0, 2);
diffusionTensors0 = NaN(4, numTimeSlices);
porosities0 = NaN(numTimeSlices, 1);
diffusionTensors1 = NaN(4, numTimeSlices);
porosities1 = NaN(numTimeSlices, 1);
diffusionTensors2 = NaN(4, numTimeSlices);
porosities2 = NaN(numTimeSlices, 1);

[A0, rhs0, isDoF0, triangleVolumes0] = assembleCellProblem(g, levelSet0(:, 1));
SOL0 = solveSystemFE(g, A0, rhs0, isDoF0);
[diffusion0, porosity0] = computeDiffusionTensor(g, SOL0, triangleVolumes0);
diffusionTensors0(:, 1) = diffusion0(:);
porosities0(1) = porosity0;

[A1, rhs1, isDoF1, triangleVolumes1] = assembleCellProblem(g, levelSet1(:, 1));
SOL1 = solveSystemFE(g, A1, rhs1, isDoF1);
[diffusion1, porosity1] = computeDiffusionTensor(g, SOL1, triangleVolumes1);
diffusionTensors1(:, 1) = diffusion1(:);
porosities1(1) = porosity1;

[A2, rhs2, isDoF2, triangleVolumes2] = assembleCellProblem(g, levelSet2(:, 1));
SOL2 = solveSystemFE(g, A2, rhs2, isDoF2);
[diffusion2, porosity2] = computeDiffusionTensor(g, SOL2, triangleVolumes2);
diffusionTensors2(:, 1) = diffusion2(:);
porosities2(1) = porosity2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of HyPHM variables

g_HyPHM = domainRectangle(0, 1, 0, 0.5, 0.02);
if (mod(endTime, dt) < EPS)
    timeSteps = 0:dt:endTime;
else
    timeSteps = [0:dt:endTime, endTime];
end
% timeSteps(end) = [];
st = Stepper(timeSteps);
conc = Variable(g_HyPHM, st, 'Concentration_MultipleCellProblems', 'P0');
% conc.setdata( 0,  @(t, x) 1 * ( norm( x(:) - [0.5; 0.5] ) < 0.05 ) + ...
%     1 * ( norm( x(:) - [1.5; 0.5] ) < 0.05 ) + ...
%     1 * ( norm( x(:) - [2.5; 0.5] ) < 0.05 ) );
conc.setdata(0, @(t, x) 0);
concTransport = Transport(g_HyPHM, st, 'Transport_MultipleCellProblems');
concTransport.id2D = {1, 3};
concTransport.id2N = {2};
concTransport.id2F = {4};
concTransport.U = conc;
concTransport.uD = Variable(g_HyPHM, st, 'uD', 'P0E');
concTransport.uD.setdata(@(t, x) 0);
concTransport.gF = Variable(g_HyPHM, st, 'gF', 'P0E');
concTransport.gF.setdata(@(t, x) -1);
concTransport.A = Variable(g_HyPHM, st, 'Porosity_MultipleCellProblems', 'P0');
concTransport.A.setdata(0, ...
    @(t, x) porosity0*(x(1) >= 0 & x(1) <= 1)+ ...
    porosity1*(x(1) >= 1 & x(1) <= 2)+ ...
    porosity2*(x(1) >= 2 & x(1) <= 3));
concTransport.D = Variable(g_HyPHM, st, 'Diffusion_MultipleCellProblems', 'P0P0P0P0');
concTransport.D.setdata(0, ...
    @(t, x) (diffusion0 * (x(1) >= 0 & x(1) <= 1) + ...
    diffusion1 * (x(1) >= 1 & x(1) <= 2) + ...
    diffusion2 * (x(1) >= 2 & x(1) <= 3))*0.01);
concTransport.C.setdata(@(t, x) [0.1; 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time iteration

%while ( t_ <= endTime + EPS )
while st.next
    disp(['Time step ', num2str(timeStep)]);
    step = st.curstep;

    %     assert( all( phi0 >= phiOld0 - EPS ) );

    if (t_ + dt > endTime - EPS)
        dt = endTime - t_;
        if (abs(dt) < EPS)
            break;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level set evolution step

    disp('Evolution of level set');
    tic;
    coordCell = mat2cell(g.coordinates, ones(1, numNodes), dim);
    timeCell = num2cell(t_*ones(numNodes, 1));

    concData = conc.getdata(timeStep-1);
    isInArea0 = g_HyPHM.baryT(:, 1) < 0.2;
    concentration = sum(concData(isInArea0).*g_HyPHM.areaT(isInArea0)) / ...
        sum(g_HyPHM.areaT(isInArea0));

    concCell = num2cell(concentration*ones(numNodes, 1));
    normalSpeed0 = cellfun(speedFun0, timeCell, coordCell, concCell);
    normalSpeed1 = cellfun(speedFun1, timeCell, coordCell, concCell);
    normalSpeed2 = cellfun(speedFun2, timeCell, coordCell, concCell);

    phi0 = levelSetEquationTimeStep(t_+dt, t_, phiOld0, g, ...
        normalSpeed0);
    phi1 = levelSetEquationTimeStep(t_+dt, t_, phiOld1, g, ...
        normalSpeed1);
    phi2 = levelSetEquationTimeStep(t_+dt, t_, phiOld2, g, ...
        normalSpeed2);


    t_old = t_;
    t_ = t_ + dt;
    phiOld0 = phi0;
    phiOld1 = phi1;
    phiOld2 = phi2;
    timeStep = timeStep + 1;

    levelSet0(:, timeStep) = phi0;
    levelSet1(:, timeStep) = phi1;
    levelSet2(:, timeStep) = phi2;

    %     for i = 1:numFunctions
    %         normalSpeed = cellfun( speedFun{i}, timeCell, coordCell );
    %         phi(:,i) = levelSetEquationTimeStep( t_, phi(:,i), t_old, phiOld(:,i), ...
    %             g, normalSpeed );
    %     end

    %     phiOld = phi;
    %     levelSet(:,timeStep,:) = phi;
    toc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting of interface

    PHI = g.reshape(phi2);
    set(hand, 'Zdata', PHI, 'LineColor', 'b');
    drawnow;
    % %     saveas( fig, [ 'data/MultipleCellProblems/Triangle', ...
    % %              num2str( timeStep ), '.png' ] );
    %     if ( 60 < timeStep && timeStep <= 70 )
    %         print( fig, '-dpng', [ 'data/MultipleCellProblems/Triangle', ...
    %              num2str( timeStep ), '.png' ] );
    %     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation of effective parameters

    disp('Calculation of effective parameters');
    tic;
    %     disp( 'Cell problem' );
    %     tic;
    [A0, rhs0, isDoF0, triangleVolumes0] = assembleCellProblem(g, levelSet0(:, timeStep));
    %     toc
    %     tic;
    SOL0 = solveSystemFE(g, A0, rhs0, isDoF0);
    %     toc
    [diffusion0, porosity0] = computeDiffusionTensor(g, SOL0, triangleVolumes0);
    %     tic;
    diffusionTensors0(:, timeStep) = diffusion0(:);
    porosities0(timeStep) = porosity0;
    %     toc

    %     disp( 'Cell problem' );
    %     tic;
    [A1, rhs1, isDoF1, triangleVolumes1] = assembleCellProblem(g, levelSet1(:, timeStep));
    %     toc
    %     tic;
    SOL1 = solveSystemFE(g, A1, rhs1, isDoF1);
    %     toc
    %     tic;
    [diffusion1, porosity1] = computeDiffusionTensor(g, SOL1, triangleVolumes1);
    diffusionTensors1(:, timeStep) = diffusion1(:);
    porosities1(timeStep) = porosity1;
    %     toc

    %     disp( 'Cell problem' );
    %     tic;
    [A2, rhs2, isDoF2, triangleVolumes2] = assembleCellProblem(g, levelSet2(:, timeStep));
    %     toc
    %     tic;
    size(isDoF2)
    SOL2 = solveSystemFE(g, A2, rhs2, isDoF2);
    %     toc
    %     tic;
    [diffusion2, porosity2] = computeDiffusionTensor(g, SOL2, triangleVolumes2);
    %     toc
    %     tic;
    diffusionTensors2(:, timeStep) = diffusion2(:);
    porosities2(timeStep) = porosity2;
    %     toc

    %     for i = 1:numFunctions
    % %         disp( 'Cell problem' );
    % %         tic;
    %         [ A, rhs, isDoF, triangleVolumes ] = assembleSystemFE( g, levelSet(:,timeStep,i) );
    % %         toc
    % %         tic;
    %         SOL = solveSystemFE( g, A, rhs, isDoF );
    % %         toc
    % %         tic;
    %         [ diffusion, porosity ] = computeDiffusionTensor( g, SOL, triangleVolumes );
    % %         toc
    %     end
    toc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Macroscopic transport step
    concTransport.A.setdata(step, ...
        @(t, x) porosity0*(x(1) >= 0 & x(1) <= 1)+ ...
        porosity1*(x(1) >= 1 & x(1) <= 2)+ ...
        porosity2*(x(1) >= 2 & x(1) <= 3));
    concTransport.D.setdata(step, ...
        @(t, x) (diffusion0 * (x(1) >= 0 & x(1) <= 1) + ...
        diffusion1 * (x(1) >= 1 & x(1) <= 2) + ...
        diffusion2 * (x(1) >= 2 & x(1) <= 3))*0.01);
    concTransport.computeLevel;

end % while

%conc.visualize;

% levelSet0 = solveLevelSetEquation( g, initialData0, speedFun0, velocity, ...
%     endTime, dt, 'ContourHandle', hand, 'SaveTime', 'all', 'ContourFigure', fig );
% levelSet1 = solveLevelSetEquation( g, initialData1, speedFun1, velocity, ...
%     endTime, dt, 'SaveTime', 'all' );
% levelSet2 = solveLevelSetEquation( g, initialData2, speedFun2, velocity, ...
%     endTime, dt, 'SaveTime', 'all' );

% numTimeSlices = size( levelSet0, 2 );
% diffusionTensors0 = NaN( 4, numTimeSlices );
% porosities0 = NaN( numTimeSlices, 1 );
% diffusionTensors1 = NaN( 4, numTimeSlices );
% porosities1 = NaN( numTimeSlices, 1 );
% diffusionTensors2 = NaN( 4, numTimeSlices );
% porosities2 = NaN( numTimeSlices, 1 );

% for ts = 1:numTimeSlices
%     disp( [ 'Time step ', num2str(ts) ] );
%
%     [ A0, rhs0, isDoF0, triangleVolumes0 ] = assembleCellProblem( g, levelSet0(:,ts) );
%     SOL0 = solveSystemFE( g, A0, rhs0, isDoF0 );
%
%     [ diffusion, porosity ] = computeDiffusionTensor( g, SOL0, triangleVolumes0 );
%     diffusionTensors0(:,ts) = diffusion(:);
%     porosities0(ts) = porosity;
%
%     [ A1, rhs1, isDoF1, triangleVolumes1 ] = assembleCellProblem( g, levelSet1(:,ts) );
%     SOL1 = solveSystemFE( g, A1, rhs1, isDoF1 );
%
%     [ diffusion, porosity ] = computeDiffusionTensor( g, SOL1, triangleVolumes1 );
%     diffusionTensors1(:,ts) = diffusion(:);
%     porosities1(ts) = porosity;
%
%     [ A2, rhs2, isDoF2, triangleVolumes2 ] = assembleCellProblem( g, levelSet2(:,ts) );
%     SOL2 = solveSystemFE( g, A2, rhs2, isDoF2 );
%
%     [ diffusion, porosity ] = computeDiffusionTensor( g, SOL2, triangleVolumes2 );
%     diffusionTensors2(:,ts) = diffusion(:);
%     porosities2(ts) = porosity;
% end
