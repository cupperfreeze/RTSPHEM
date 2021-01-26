function [outX, outY] = GetSpeRectangles(coarseningFactor)
%obtain side - lengths of rectangles that approximately match the
%permeability distribution provided by SPE10 permeability map
%coarseningFactor^2 data points will be averaged to a single one

global EPS
EPS = eps;
load('PermDataSPE.mat')
Xlength = linspace(0.01, 0.5, 50);
fine = linspace(0.01, 0.5, 200);
[InterpolX, InterpolY] = createLookUp(Xlength, fine);

temp = Permeability(1:2200, :);
readX = reshape(temp', 60, 220)'; %result 220x60 matrix (first component~x-Axis)
temp = Permeability(2201:4400, :);
readY = reshape(temp', 60, 220)';

coarseX = zeros(size(readX)/coarseningFactor);
coarseY = zeros(size(readY)/coarseningFactor);

%compute mean value over each batch
for i = 1:(220 / coarseningFactor)
    for j = 1:(60 / coarseningFactor)
        temp = readX(((i - 1)*coarseningFactor+1):i*coarseningFactor, ...
            ((j - 1) * coarseningFactor + 1):j*coarseningFactor);
        coarseX(i, j) = mean(temp(:));

        temp = readY(((i - 1)*coarseningFactor+1):i*coarseningFactor, ...
            ((j - 1) * coarseningFactor + 1):j*coarseningFactor);
        coarseY(i, j) = mean(temp(:));
    end
end

%normalize result
coarseX = coarseX / 20000;
coarseY = coarseY / 20000;

outX = zeros(size(coarseX));
outY = zeros(size(coarseY));

% map rectangle side-lengths to averaged permeability value by nearest
% neighbor
for i = 1:size(coarseX, 1)
    for j = 1:size(coarseX, 2)
        [~, index] = min((coarseX(i, j)-InterpolX(:)).^2+(coarseY(i, j) - InterpolY(:)).^2);
        temp = meshgrid(fine);
        outX(i, j) = temp(index);
        outY(i, j) = temp(index);
    end
end
%linear projection to regularize
outX = 0.1 + 0.6 * outX;
outY = 0.1 + 0.6 * outY;

porosity = mean(1-4*outX(:).*outY(:))
end
