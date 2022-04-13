% Compute Darcy velocity field on pemreability data of SPE10, cf. [2]

load('PermDataSPE.mat');

dimension = 2;
cellGrid = CartesianGrid(dimension, ...
    [1, 220, 1, 60], ...
    [220 - 1, 60 - 1]);
coord = cellGrid.coordinates;

gridHyPHM = Grid(coord, cellGrid.triangles);
gridHyPHM.idE(gridHyPHM.baryE(:, 1) < 1+1000*eps) = 4;
gridHyPHM.idE(gridHyPHM.baryE(:, 1) > 220-1000*eps) = 2;
gridHyPHM.idE(gridHyPHM.baryE(:, 2) < 1+1000*eps) = 1;
gridHyPHM.idE(gridHyPHM.baryE(:, 2) > 60-1000*eps) = 3;
%gridHyPHM = domainRectangle( 0, 2,  0, 1, 0.01 );
%gridHyPHM.visualize('nuE')

Xdir = reshape(Permeability(1:2200, :)', 60, 220)';
Ydir = reshape(Permeability(4401:6600, :)', 60, 220)';
flowStepper = Stepper(0:1);
XauxVar = Variable(gridHyPHM, flowStepper, 'Xaux', 'P1');
XauxVar.setdata(0, Xdir(:));
YauxVar = Variable(gridHyPHM, flowStepper, 'Yaux', 'P1');
YauxVar.setdata(0, Ydir(:));
Xaux = XauxVar.getTSI(0);
Yaux = YauxVar.getTSI(0);

phead = Variable(gridHyPHM, flowStepper, 'pressure head', 'P0');
darcyVelocity = Variable(gridHyPHM, flowStepper, 'Darcy velocity', 'RT0');
phead.setdata(0, @(t, x) 0.0);
darcyVelocity.setdata(0, @(t, x) 0.0);

darcy = Transport(gridHyPHM, flowStepper, 'Darcy problem');
darcy.id2D = {2};
darcy.id2F = {1, 3, 4};
darcy.uD.setdata(@(t, x)0);
darcy.gF.setdata(@(t, x) -1*(x(1) < 1 + 1000 * eps));

darcy.D = Variable(gridHyPHM, flowStepper, 'Permeability', 'P0P0P0P0');
darcy.D.setdata(@(t, x) [Xaux([x(1), x(2)]), 0; 0, Yaux([x(1), x(2)])]);
%darcy.D.setdata( @(t,x) [1 0; 0,1]  );

darcy.F.setdata(@(t, x) 0);
darcy.Q = darcyVelocity;
darcy.U = phead;

flowStepper.next;
darcy.computeLevel;

darcy.Q.visualize();
darcy.U.visualize();
