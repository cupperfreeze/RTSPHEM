% Data related to Figure 11
% Script to evaluate the adaptivity scheme in the 'random porosity' example

try
    load('RandomPorosityData.mat');
catch
    error('unable to load mat file')
end

a = cell2mat(porosities);
k = numel(porosities{1}(:));
close all
for i = 25:25
    figure;
    duda = heatmap(flipud(reshape(a(i + (0:numberOfSlices - 1) * k), numberOfSlicesX, numberOfSlicesY)'));
    duda.Title = 'Porosity map, step 25';
    set(gca, 'FontSize', 15);
end

for i = 1:1
    figure;
    duda = heatmap(flipud(reshape(a(i + (0:numberOfSlices - 1) * k), numberOfSlicesX, numberOfSlicesY)'));
    duda.Title = 'Porosity map, step 0';
    set(gca, 'FontSize', 15);
end

a = recalculated{17};
figure;
duda = heatmap(flipud(reshape(a, numberOfSlicesX, numberOfSlicesY)'));
duda.Title = 'Recalculations Cell Problems, step 17';
set(gca, 'FontSize', 15);

a = recalculated{20};
figure;
duda = heatmap(flipud(reshape(a, numberOfSlicesX, numberOfSlicesY)'));
duda.Title = 'Recalculations Cell Problems, step 20';
set(gca, 'FontSize', 15);

a = recalculated{30};
figure;
duda = heatmap(flipud(reshape(a, numberOfSlicesX, numberOfSlicesY)'));
duda.Title = 'Recalculations Cell Problems, step 30';
set(gca, 'FontSize', 15);

a = recalculated{40};
figure;
duda = heatmap(flipud(reshape(a, numberOfSlicesX, numberOfSlicesY)'));
duda.Title = 'Recalculations Cell Problems, step 40';
set(gca, 'FontSize', 15);