% Evaluate results from simulation 'PermIrregularSoilPeriodic'
%--> obtain porosity and permeability over time as well as permeability
%over porosity


close all
num = size(permeabilityTensors, 2);
for i = 1:num
    Perm(i, 1) = permeabilityTensors(1, i);
    Perm(i, 2) = permeabilityTensors(4, i);
end
figure
plot(linspace(0, 1.2, num), Perm(:, 1), 'lineWidth', 3)
hold on
plot(linspace(0, 1.2, num), Perm(:, 2), 'lineWidth', 3)
%title('Permeability tensor','fontsize',20);
%axis([-100, 8000, 0, 0.7*10^(-7)])
xlabel('t [s]', 'fontsize', 15);
ylabel('Permeability', 'fontsize', 15);
legend('K_x', 'K_y', 'location', 'northwest')
set(gca, 'fontsize', 15);

figure
plot(linspace(0, 1.2, num), porosityValues, 'lineWidth', 3)
%title('Permeability tensor','fontsize',20);
%axis([-100, 8000, 0, 0.7*10^(-7)])
xlabel('t [s]', 'fontsize', 15);
ylabel('Porosity', 'fontsize', 15);
set(gca, 'fontsize', 15);

figure
plot(porosityValues, Perm(:, 1), 'lineWidth', 3);
hold on
plot(porosityValues, Perm(:, 2), 'lineWidth', 3);
hold on

surfaceValues2 = surfaceValues ./ (1 - porosityValues);
C = Perm(1, 2) * (1 - porosityValues(1))^2 / porosityValues(1)^3 * surfaceValues2(1)^2;
Kozeny = C * ((1 - porosityValues').^(-2) .* porosityValues'.^3) ./ surfaceValues2'.^2;

plot(porosityValues, Kozeny, 'lineWidth', 3);

%title('Permeability tensor','fontsize',20);
%axis([-100, 8000, 0, 0.7*10^(-7)])
xlabel('Porosity', 'fontsize', 15);
ylabel('Permeability', 'fontsize', 15);
legend('K_x', 'K_y', 'Kozeny-Carman', 'location', 'northwest')
set(gca, 'fontsize', 15);