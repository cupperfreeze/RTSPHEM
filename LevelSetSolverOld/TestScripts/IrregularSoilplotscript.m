% script to plot eigenvalues of diffusion tensors
% cf. [5] Figure 7
%data=load('data.mat');

steps = size(diffusionTensors, 2);
eigenValues = NaN(steps, 2);

for i = 1:steps
    eigenValues(i, :) = eig(reshape(diffusionTensors(:, i), 2, 2));
end

figure
time = (1:steps) * endTime / steps;
upper = ones(1, steps) * eigenValues(end, 2);
lower = zeros(1, steps);
plot(time, eigenValues(:, 1), 'r', time, eigenValues(:, 2), 'b');
xlabel('Time');
ylabel('Eigenvalues');
title('plot');
