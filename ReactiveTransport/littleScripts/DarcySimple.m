% Test Darcy equation on unit square with obstacle


dimension = 2;                                                              %dimension level-set grid
numPartitions = 150;                                                        %fineness level-set grid
cellGrid = CartesianGrid(dimension, ...                                     %initialize level-set grid
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));


coord = cellGrid.coordinates;
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);

%Define and sample level-set function
initialLevelSetFunc = @(x) 0.2 - norm([x(1), x(2)], 2);
levelSet = cellfun(initialLevelSetFunc, coordCell);

%Initialize Grid for Darcy
gridHyPHM = Grid(cellGrid.coordinates, cellGrid.triangles);

%Align Stokes mesh edges with geometry
[gridHyPHM, levelSet] = localMeshRefinementUnFolded(gridHyPHM, ...
    levelSet);
    
%Identify Edges (left:4 - right:2 - bottom:1- top:3)
gridHyPHM.idE(gridHyPHM.baryE(:, 1) < -0.5+eps) = 4;
gridHyPHM.idE(gridHyPHM.baryE(:, 1) > 0.5-eps) = 2;
gridHyPHM.idE(gridHyPHM.baryE(:, 2) < -0.5 + eps) = 1;
gridHyPHM.idE(gridHyPHM.baryE(:, 2) > 0.5-eps) = 3;

st = Stepper(0:1);                                                          %time stepper
Darcy = TransportLEVEL(gridHyPHM, st, 'Darcy');                             %Initialize Darcy problem
Darcy.L.setdata(levelSet);                                                  %Set geometry variable
Darcy.id2F = {4, 3, 1};                                                     %Set boundary conditions
Darcy.id2D = {2};
Darcy.gF.setdata(@(t, x) -1.*(1-4*(x(2)).^2)*(x(1) < -0.5+eps));
Darcy.uD.setdata(@(t, x) 0);

st.next;                                                                    %advance stepper
Darcy.computeLevel('s');                                                    %solve Darcy


%Plot pressure and velocity data
Darcy.U.visualize();
Darcy.Q.visualize();
