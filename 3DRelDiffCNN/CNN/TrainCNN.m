%Train the neural Network

p = genpath('./VoxelData/SampleBENTHEIMER');
addpath(p);
% set up file paths, functions, load data
load('DiffSample2.mat');
imSize = [100,100,100];
PercentageOfTrainImages = 0.90;												%Share training data wrt. to $complete data set 
InidzesOfData = 1:30000;													%Indizes of data files
TrafoPermtoLabel = @(x) log(10.*x);
TrafoLabeltoPerm = @(x) exp(x)./10;
trainingIm = zeros([imSize(1), imSize(2), imSize(3),1, numel(InidzesOfData )/6], 'uint8');
miniBatchSize = 8;                                                         %Mini-batch Size for network training


% read in Data

disp( 'read in data' );
for i = InidzesOfData(1:5000) 
    name = strcat('out_connected', num2str(i-1),['_',num2str(imSize(1)),'x',num2str(imSize(2)),'x',num2str(imSize(3))]);
    im = rwd2mat(name);
    trainingIm(:, :, :, 1, i) = im;
    if mod(i,200)==0
	disp([num2str(i/numel(InidzesOfData ) *100), ' Procent done'])
    end
end
disp( 'read in data done' );


trainingData = reshape(DiffData(InidzesOfData),numel(InidzesOfData)/6,6);
trainingData = trainingData(:,6);

% split Data

numberOfRows = PercentageOfTrainImages*numel(InidzesOfData)/6;

numberOfData = size(trainingIm,5);
validationIm = trainingIm(:, :, :, :, (numberOfRows +1):end);
validationData = trainingData( (numberOfRows +1):end ,:);

trainingIm = trainingIm(:, :, :, :, 1:numberOfRows);
trainingData = trainingData(1:numberOfRows,:);

idx = ((trainingData >0.01) & (trainingData <0.3));
trainingIm = trainingIm(:,:,:,1,idx);
trainingData = trainingData (idx,:);

idx = ((validationData >0.01) & (validationData <0.3));
validationIm = validationIm(:,:,:,1,idx);
validationData = validationData (idx,:);



% Define network layers
layers = [ ...
    image3dInputLayer([imSize(1), imSize(2), imSize(3), 1],'Name', 'input1'); ...
    convolution3dLayer(5 * [1, 1, 1], 32, 'Padding', 2*[1,1,1], 'Name', 'conv1'); ...
batchNormalizationLayer('Name', 'B1'); ...
   leakyReluLayer(0.1,'Name', 'relu1'); ...
    maxPooling3dLayer(5, 'Stride', 5,'Name', 'MP1'); ...
    convolution3dLayer(5 * [1, 1, 1], 64, 'Padding', 2*[1,1,1],'Name', 'conv2'); ...
batchNormalizationLayer('Name', 'B2'); ...
    leakyReluLayer(0.1,'Name', 'relu2'); ...
    maxPooling3dLayer(4, 'Stride', 4,'Name', 'MP2'); ...
    convolution3dLayer(3 * [1, 1, 1], 100, 'Padding', 1*[1,1,1], 'Name', 'conv3'); ...
batchNormalizationLayer('Name', 'B3'); ...
   leakyReluLayer(0.1,'Name', 'relu3'); ...

%maxPooling3dLayer(2, 'Stride', 2,'Name', 'MP3'); ...
      fullyConnectedLayer(64, 'Name', 'fc2'); ...
    leakyReluLayer(0.1,'Name', 'relu5'); ...
     fullyConnectedLayer(32, 'Name', 'fc4'); ...
  leakyReluLayer(0.1,'Name', 'relu6'); ...
      fullyConnectedLayer(1,'Name', 'fc3'); ...
    regressionLayer('Name', 'out')];

layers = [ ...
    image3dInputLayer([imSize(1), imSize(2), imSize(3), 1],'Name', 'input1'); ...
    convolution3dLayer(5 * [1, 1, 1], 16, 'Padding', 2*[1,1,1], 'Name', 'conv1'); ...
batchNormalizationLayer('Name', 'B1'); ...
   leakyReluLayer(0.1,'Name', 'relu1'); ...
    maxPooling3dLayer(5, 'Stride', 5,'Name', 'MP1'); ...
    convolution3dLayer(5 * [1, 1, 1], 32, 'Padding', 2*[1,1,1],'Name', 'conv2'); ...
batchNormalizationLayer('Name', 'B2'); ...
    leakyReluLayer(0.1,'Name', 'relu2'); ...
    maxPooling3dLayer(4, 'Stride', 4,'Name', 'MP2'); ...
    convolution3dLayer(3 * [1, 1, 1], 50, 'Padding', 1*[1,1,1], 'Name', 'conv3'); ...
batchNormalizationLayer('Name', 'B3'); ...
   leakyReluLayer(0.1,'Name', 'relu3'); ...
 convolution3dLayer(3 * [1, 1, 1], 150, 'Padding', 1*[1,1,1], 'Name', 'conv4'); ...
batchNormalizationLayer('Name', 'B4'); ...
   leakyReluLayer(0.1,'Name', 'relu4'); ...

%maxPooling3dLayer(2, 'Stride', 2,'Name', 'MP3'); ...
      fullyConnectedLayer(64, 'Name', 'fc2'); ...
    leakyReluLayer(0.1,'Name', 'relu5'); ...
     fullyConnectedLayer(32, 'Name', 'fc4'); ...
  leakyReluLayer(0.1,'Name', 'relu6'); ...
      fullyConnectedLayer(1,'Name', 'fc3'); ...
    regressionLayer('Name', 'out')];


%layers = [ ...
%    image3dInputLayer([imSize(1), imSize(2), imSize(3), 1],'Name', 'input1'); ...
%    convolution3dLayer(5 * [1, 1, 1], 20, 'Padding', 2*[1,1,1], 'Name', 'conv1'); ...
% %   batchNormalizationLayer('Name', 'B1'); ...
%    leakyReluLayer(0.1,'Name', 'relu1'); ...
%    maxPooling3dLayer(2, 'Stride', 2,'Name', 'MP1'); ...
%    
%convolution3dLayer(5 * [1, 1, 1], 40, 'Padding', 2*[1,1,1],'Name', 'conv2'); ...
%batchNormalizationLayer('Name', 'B2'); ...
%    leakyReluLayer(0.1,'Name', 'relu2'); ...
%    maxPooling3dLayer(2, 'Stride', 4,'Name', 'MP2'); ...
%
%    convolution3dLayer(5 * [1, 1, 1], 60, 'Padding', 1*[1,1,1], 'Name', 'conv3'); ...
%    dropoutLayer(0.5);...
%batchNormalizationLayer('Name', 'B3'); ...
%   leakyReluLayer(0.1,'Name', 'relu3'); ...
%
% convolution3dLayer(5 * [1, 1, 1], 80, 'Padding', 1*[1,1,1], 'Name', 'conv4'); ...
%    dropoutLayer(0.5);...
%    fullyConnectedLayer(32, 'Name', 'fc4'); ...
%  leakyReluLayer(0.1,'Name', 'relu6'); ...
%      fullyConnectedLayer(1,'Name', 'fc3'); ...
%   regressionLayer('Name', 'out')];



%Group pore-space images, physics input and permeability labels in Datastores

a1=arrayDatastore(trainingIm,"IterationDimension",5);
a3=arrayDatastore(TrafoPermtoLabel(trainingData),"IterationDimension",1);
training = combine(a1,a3);

a1=arrayDatastore(validationIm,"IterationDimension",5);
a3=arrayDatastore(TrafoPermtoLabel(validationData),"IterationDimension",1);
validation = combine(a1,a3);


%training options

opts = trainingOptions('sgdm', ...
    'MaxEpochs', 13, ...
'L2Regularization',0.02,...
    'InitialLearnRate', 0.001, ...
    'MiniBatchSize', miniBatchSize , ...	
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 4, ...
    'LearnRateDropFactor', 0.4, ...
    'ValidationFrequency', 500, ...
    'ExecutionEnvironment', 'gpu', ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress', 'ValidationData', validation, ...
    'Verbose', true);


%Train the network
net = trainNetwork(training, layers, opts);

% evaluate training success

% Compute network prediction on validation data
predicted = zeros(size(validationIm,5),1);
for i=1:size(validationIm,5)
    a1=arrayDatastore(validationIm(:,:,:,1,i),"IterationDimension",5);

    predicted(i,:) = TrafoLabeltoPerm(net.predict(a1));
end

% Compute network prediction on training data
predictedTrain = zeros(size(trainingIm,5),1);
for i=1:size(trainingIm,5)
    a1=arrayDatastore(trainingIm(:,:,:,1,i),"IterationDimension",5);

   predictedTrain(i,:) = TrafoLabeltoPerm(net.predict(a1));
end

valuesTrain = (trainingData);
save('CNNDiffPlain05Small.mat','net','predicted','validationData', 'predictedTrain', 'trainingData','-v7.3');


% Regression plots and R^2 value computation

figure
plot(validationData, predicted, '.');
hold on
plot(validationData, validationData, 'LineWidth', 2);
xlabel('Perm data');
ylabel('Perm prediction');
set(gca, 'FontSize', 14);
title('Validation');

figure
plot(trainingData,predictedTrain,'*')
hold on
plot(trainingData,trainingData, 'LineWidth', 2)
xlabel('Perm data');
ylabel('Perm prediction');
set(gca, 'FontSize', 14);
title('Training');

R_squredLOG = 1 - sum((log10(predicted)-log10(validationData)).^2) / (sum((log10(validationData)-mean(log10(validationData))).^2))
R_squred = 1 - sum((predicted-(validationData)).^2) / (sum(((validationData)-mean((validationData))).^2))


R_squredLOG = 1 - sum((log10(predictedTrain)-log10(trainingData)).^2) / (sum((log10(trainingData)-mean(log10(trainingData))).^2))
R_squred = 1 - sum((predictedTrain-(trainingData)).^2) / (sum(((trainingData)-mean((trainingData))).^2))

