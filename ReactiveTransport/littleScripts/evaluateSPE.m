% Script to evaluate transport simulation based on SPE 10 data

try
    load('SPE10.mat')
catch
    error('unable to load file')
end

%% Level-Set mass preservation


TotalFlux(1, 1:4) = 0;
endStep = 440;
indizes = magnesiumTransport.grid.baryE(:, 1) > lengthXAxis - 10 * eps;
for i = 1:endStep
    for j = 1:numel(speciesCells)
        temp = speciesCells{j}.Q.getdata(i-1);
        TotalFlux(i+1, j) = TotalFlux(i, j) + (transportStepper.timepts(i + 1) - transportStepper.timepts(i)) * ...
            sum(temp(indizes).*speciesCells{j}.grid.areaE(indizes));
    end
end


totalInitial(1:4, 1) = 0;
for i = 1:numberOfSlices
    totalInitial(1) = totalInitial(1) + (outX(i) * outY(i) * 2 / molarVolume{2} + outX(i) * outY(i) * 2 / molarVolume{1}) / numberOfSlices * lengthXAxis * lengthYAxis;
    totalInitial(2) = totalInitial(2) + (2 * outX(i) * outY(i) * 2 / molarVolume{2} + outX(i) * outY(i) * 2 / molarVolume{1}) / numberOfSlices * lengthXAxis * lengthYAxis;
    totalInitial(3) = totalInitial(3) + -outX(i) * outY(i) * 2 / molarVolume{1} / numberOfSlices * lengthXAxis * lengthYAxis;
    totalInitial(4) = totalInitial(4) + outX(i) * outY(i) * 2 / molarVolume{2} / numberOfSlices * lengthXAxis * lengthYAxis;
end
totalInitial(3) = totalInitial(3) + speciesCells{3}.U.getdata(0)' * speciesCells{3}.grid.areaT;

totalTime = zeros(endStep+1, numel(speciesCells));
totalTime(1, 1:4) = totalInitial;

for i = 1:endStep
    totalTime(i+1, 1) = sum(MassOfD(:, i)) + sum(MassOfP(:, i)) + TotalFlux(i, 1) + ...
        speciesCells{1}.U.getdata(i)' * speciesCells{1}.grid.areaT;

    totalTime(i+1, 2) = 2 * sum(MassOfD(:, i)) + sum(MassOfP(:, i)) + TotalFlux(i, 2) + ...
        speciesCells{2}.U.getdata(i)' * speciesCells{2}.grid.areaT;

    totalTime(i+1, 3) = -sum(MassOfP(:, i)) + TotalFlux(i, 3) ...
        -lengthYAxis * inletVelocity * initialHydrogenConcentration * transportStepper.timepts(i+1) + ...
        speciesCells{3}.U.getdata(i)' * speciesCells{3}.grid.areaT;

    totalTime(i+1, 4) = sum(MassOfD(:, i)) + TotalFlux(i, 4) + ...
        speciesCells{4}.U.getdata(i)' * speciesCells{4}.grid.areaT;
end

figure
for i = [1, 2, 4]
    plot(transportStepper.timepts(2:endStep + 1), totalTime(2:end, i)/totalTime(2, i), 'LineWidth', 2)
    hold on
end
legend('Ca^{2+}', 'CO_{3}^{2-}', 'Mg^{2+}', 'Location', 'southwest')
axis([0, endTime, 0.98, 1.01])
set(gca, 'Fontsize', 13)
xlabel('t in s');
ylabel('M_{rel}')

%%%%%

%% Darcy velocity distribution
endStep = timeIterationStep;
minVelocity = zeros(endStep, 2);
maxVelocity = zeros(endStep, 2);
avVelocity = zeros(endStep, 2);
for i = 1:endStep
    velocities = RT0.RT0toP0P0slice(magnesiumTransport.C.grid, magnesiumTransport.C.getdata(i));
    minVelocity(i, :) = min(velocities, [], 1);
    maxVelocity(i, :) = max(velocities, [], 1);
    avVelocity(i, :) = abs(velocities') * magnesiumTransport.C.grid.areaT / (sum(magnesiumTransport.C.grid.areaT));
end

figure
plot(transportStepper.timepts(2:endStep + 1), minVelocity(:, 1)/inletVelocity, 'LineWidth', 2);
hold on
plot(transportStepper.timepts(2:endStep + 1), minVelocity(:, 2)/inletVelocity, 'LineWidth', 2);
hold on
plot(transportStepper.timepts(2:endStep + 1), maxVelocity(:, 1)/inletVelocity, 'LineWidth', 2);
hold on
plot(transportStepper.timepts(2:endStep + 1), maxVelocity(:, 2)/inletVelocity, 'LineWidth', 2);

axis([0, endTime, -6, 11])
set(gca, 'Fontsize', 13)
xlabel('t in s');
ylabel('v_{rel}')
legend('minX', 'minY', 'maxX', 'maxY', 'Location', 'northwest')
