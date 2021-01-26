% Data related to Figure 12
% Script to visualize the difference between pore and 2-scale approach

try
    load('MacroSimulationData')
catch
    error('Probably out of memory')
end
MGrid = calciumTransport.grid;
Mcalcium = calciumTransport.U;

try
    load('MicroSimulationData')
catch
    error('Probably out of memory')
end
mGrid = calciumTransport.grid;
mcalcium = calciumTransport.U;

%%% Visualize normailzed difference between both solutions

% coarse = Variable( MGrid, Stepper(0:1), ...
%     'Ca^(2+)grob', 'P0' );
%  DIFF = Variable( mGrid, Stepper(0:1), ...
%      'Ca^(2+)Diff', 'P0' );
% fine = Variable( mGrid, Stepper(0:1), ...
%     'Ca^(2+)fine', 'P0' );
% coarse.setdata(Mcalcium.getdata(53));
% fine.setdata(mcalcium.getdata(53));
% fineData = fine.getdata(1);
% wfun = grob.getTSI(0);
% wfun2 = fine.getTSI(0);
% interpol = zeros(MGrid.numT,1);
% for i = 1: mGrid.numT
%       interpol(i) = wfun(mGrid.baryT(i,1),mGrid.baryT(i,2))
% end
% DIFF.setdata(abs(fine.getdata(0)-interpol)./max(abs(fine.getdata(0)-interpol)));
% DIFF.visualize();
% out = distanceL1(mGrid,mcalcium.getdata(5) , @(x)  wfun(x(1),x(2)),0, 0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
fineness = 27;
MacroL1 = zeros(1, fineness);
MicroL1 = zeros(1, fineness);
timeStep = 53;

for i = 1:fineness

    Left = 0.1 / fineness * (i - 1);
    Right = 0.1 / fineness * i;
    MicroL1(i) = ((mGrid.baryT(:, 1) < Right) .* (mGrid.baryT(:, 1) >= Left) .* mGrid.areaT)' * mcalcium.getdata(timeStep);
    MacroL1(i) = ((MGrid.baryT(:, 1) < Right) .* (MGrid.baryT(:, 1) >= Left) .* MGrid.areaT)' * Mcalcium.getdata(timeStep) * porosities{i}(timeStep, 1);

end
MicroL1(end) = MicroL1(end-1); %Correction at outflow boundary
figure
plot(0.0:0.1/fineness:0.1, [0, MacroL1], 'linewidth', 2)
hold on
plot(0.0:0.1/fineness:0.1, [0, MicroL1], 'linewidth', 2)
%title('L1 norm slice wise, ts=53');
xlabel('x [dm]')
ylabel('calcium [mol]')
legend('Multiscale', 'Micro');
legend('Location', 'southeast')
set(gca, 'FontSize', 15)
