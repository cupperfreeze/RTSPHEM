% load data file
try
    load('ChenRelease')
catch
    error('unable to load file')
end

%Calculate total fluxes up to time t as \int_0^t \int_\Gamma j(x,t) \cdot
%\nu d\sigma ds for all three species
AtotalFlux{1} = 0;
endStep = 94;

for i = 1:endStep
    indizes = AconcentrationSave.grids{i+1}.baryE(:, 1) < eps;
    AtotalFlux{i+1} = AtotalFlux{i} + (transportStepperSave.timepts(i + 1) - transportStepperSave.timepts(i)) * sum(AFlux{i}(indizes).*AconcentrationSave.grids{i + 1}.areaE(indizes));
end
AtotalFlux;

BtotalFlux{1} = 0;
for i = 1:endStep
    indizes = BconcentrationSave.grids{i+1}.baryE(:, 1) < eps;
    BtotalFlux{i+1} = BtotalFlux{i} + (transportStepperSave.timepts(i + 1) - transportStepperSave.timepts(i)) * sum(BFlux{i}(indizes).*AconcentrationSave.grids{i + 1}.areaE(indizes));
end
BtotalFlux;

CtotalFlux{1} = 0;

for i = 1:endStep
    indizes = CconcentrationSave.grids{i+1}.baryE(:, 1) < eps;
    CtotalFlux{i+1} = CtotalFlux{i} + (transportStepperSave.timepts(i + 1) - transportStepperSave.timepts(i)) * sum(CFlux{i}(indizes).*AconcentrationSave.grids{i + 1}.areaE(indizes));
end
CtotalFlux;


% calculate mass conservation over time
for i = 1:endStep
    totalAFinal(i) = AconcentrationSave.getdata(i)' * AconcentrationSave.grids{i+1}.areaT - sum(saveXi(:, i + 1) == 2) / sum(saveXi(:, i + 1) > 0) * sum(AconcentrationSave.grids{2}.areaT) / molarVolume + AtotalFlux{i+1};
    totalBFinal(i) = BconcentrationSave.getdata(i)' * BconcentrationSave.grids{i+1}.areaT + sum(saveXi(:, i + 1) == 2) / sum(saveXi(:, i + 1) > 0) * sum(BconcentrationSave.grids{2}.areaT) / molarVolume + BtotalFlux{i+1} + ...
        sum(saveXi(:, i + 1) == 3) / sum(saveXi(:, i + 1) > 0) * sum(BconcentrationSave.grids{2}.areaT) / molarVolume * 4;
    totalCFinal(i) = CconcentrationSave.getdata(i)' * CconcentrationSave.grids{i+1}.areaT + sum(saveXi(:, i + 1) == 3) / sum(saveXi(:, i + 1) > 0) * sum(CconcentrationSave.grids{2}.areaT) / molarVolume * 4 + CtotalFlux{i+1};
end


%plot data over time
plot(transportStepperSave.timepts(2:endStep + 1), (totalAFinal')/totalAFinal(1));
hold on

plot(transportStepperSave.timepts(2:endStep + 1), (totalBFinal')/totalBFinal(1));
hold on

plot(transportStepperSave.timepts(2:endStep + 1), (totalCFinal')/totalCFinal(1));
legend('A', 'B', 'C')