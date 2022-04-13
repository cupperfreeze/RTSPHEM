% Data related to [4] Figure 7
% Script to evaluate permeability evolution in 'Composite Mineral' example

try
    load('CompositeMineral.mat');
catch
    error('unable to load mat file')
end

microscaleGrid = FoldedCartesianGrid(2, [0, 1, 0, 1], 64*[1, 1]);

permeabilityTensors(:, 1) = permeabilityTensors(:, 2);
for i = 1:numel(timeSteps)
    temp(i, :) = permeabilityTensors([1, 4], i);
    [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
        = assembleCellProblem(microscaleGrid, levelSet(:, i));
    porositiesNew(i) = sum(triangleVolumes) / 0.01;
    surfacesNew(i) = sum(triangleSurfaces) / (0.01 - sum(triangleVolumes));
end


figure
imshow(reshape(structure, numPartitions + 1, numPartitions + 1)', 'InitialMagnification', 'fit');

figure
for i = [1, 2]
    plot(timeSteps, temp(:, i), 'lineWidth', 3);
    hold on
end
hold off
%title('Permeability tensor','fontsize',20);
%legend('entry (1,1)', 'entry (1,2)', 'entry (2,2)','Location', 'northwest','fontsize',16)
%axis([-100, 8000, 0, 0.7*10^(-7)])
xlabel('t [s]', 'fontsize', 16);
ylabel('Permeability', 'fontsize', 16);
set(gca, 'fontsize', 15);
legend('K_x', 'K_y');

figure
contour(X, Y, reshape(levelSet(:, 1), numPartitions + 1, numPartitions + 1), [0, 1], 'lineWidth', 3, 'LineColor', 'r');
hold on
contour(X, Y, reshape(levelSet(:, 78), numPartitions + 1, numPartitions + 1), [0, 1], 'lineWidth', 3, 'LineColor', 'b');
axis equal

xlabel('x')
ylabel('y')
%title('shape','fontsize',16);
%legend('initial shape','shape end-time','Location', 'northwest','fontsize',16);

set(gca, 'fontsize', 16);

figure
C = permeabilityTensors(4, 1) * (1 - porositiesNew(1))^2 / porositiesNew(1)^3 * surfacesNew(1)^2;
Kozeny = C * ((1 - porositiesNew).^(-2) .* porositiesNew.^3) ./ surfacesNew.^2;
for i = [1, 2]
    plot(porositiesNew, temp(:, i), 'lineWidth', 3);
    hold on
end
plot(porositiesNew, Kozeny, 'lineWidth', 3);
xlabel('Porosity');
ylabel('Permeability');
legend('K_x', 'K_y', 'Kozeny-Carman', 'location', 'northwest');
set(gca, 'fontsize', 15);