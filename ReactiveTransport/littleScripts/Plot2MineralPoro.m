% Data related to Figure 7
% Script to evaluate permeability evolution in 'Composite Mineral' example

try
    load('CompositeMineral.mat');
catch
    error('unable to load mat file')
end

permeabilityTensors(:, 1) = permeabilityTensors(:, 2);
for i = 1:numel(timeSteps)
    temp(i, :) = permeabilityTensors([1, 4], i);
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