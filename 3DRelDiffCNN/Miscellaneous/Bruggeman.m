% Fit Bruggeman equation to data

load('SampleSat2Data.mat');
load('DiffSample2.mat');

DiffData = reshape(DiffData, 5000, 6);
idx = ((DiffData(:, 6) > 0.01) & (DiffData(:, 6) < 0.3));
DiffData = DiffData(idx, :);

x = DiffData(:, 1);
y = DiffData(:, 6);
g = fittype('0.5.^a .*x');
f0 = fit(x, y, g);
fit = f0.a;

z = 0.5.^fit * x;
plot(x, y, '*')
hold on
plot(x, z, '.')

R_squred = 1 - sum((z-(y)).^2) / (sum(((y)-mean((y))).^2))