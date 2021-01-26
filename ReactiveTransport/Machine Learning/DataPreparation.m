% Script to obtain Training and Validation data sets from CT-Images

global EPS
EPS = eps;
numberOut = 4000; % number of samples to extract
TotalSize = 575; % side-length of CT-scan
AimSize = 64; % output size samples
innerIterations = 5; % Level-Set adaptations per original sample

% in order to use commented code below to obatin training and validation
% data sets on natural geometries download the following data:
% DOI: https://doi.org/10.6084/m9.figshare.1153794.v2

% read in data file
%  file = fopen('Berea_4.52.raw','r');
%  image = reshape(fread(file,TotalSize^3,'uint8'),TotalSize,TotalSize,TotalSize);
%  fclose(file);
%
% [TrainingIm,TrainingData, TrainingPhy] = getTrainingData(numberOut,AimSize, innerIterations, image);
% [ValidationIm,ValidationData, ValidationPhy] = getTrainingData(floor(numberOut/10),AimSize, innerIterations, image);
%
% save('CNNPermGeneral.mat', 'TrainingIm', 'TrainingData', 'ValidationIm', 'ValidationData', 'TrainingPhy', 'ValidationPhy', '-v7.3')

% obtain traning and validation data
tic
[TrainingIm, TrainingData, TrainingPhy] = rectangleTrain(5*numberOut, AimSize);
[ValidationIm, ValidationData, ValidationPhy] = rectangleTrain(5*floor(numberOut / 10), AimSize);
toc
save('DataRectangle.mat', 'TrainingIm', 'TrainingData', 'ValidationIm', 'ValidationData', 'TrainingPhy', 'ValidationPhy', '-v7.3')