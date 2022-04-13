%Train the neural Network


% set up file paths, functions, load data
p=genpath('./ComputedData');
addpath(p);
p=genpath('./3DPermCNN.git');
addpath(p);
p=genpath('./VoxelData/Sample1');
addpath(p);

TrafoPermtoLabel = @(x) 2*(log10(x)+6) ;
TrafoLabeltoPerm = @(x) 10.^(x/2-6);
TrafoMinCuts = @(x) TrafoPermtoLabel(x.^1.407*10^(-8.183));
imSize = [100,100,100];
load('PermDataBENTHEIMEROrd1.mat');
load('SampleDataBENTHEIMER.mat');

PercentageOfTrainImages = 0.90;												%Share training data wrt. to complete data set 
InidzesOfData = 1:10000;													%Indizes of data files

trainingIm = zeros([imSize(1), imSize(2), imSize(3),1, numel(InidzesOfData )], 'uint8');
miniBatchSize = 16;                                                         %Mini-batch Size for network training


% read in Data

disp( 'read in data' );
for i = InidzesOfData 
    name = strcat('out_connected', num2str(i-1),['_',num2str(imSize(1)),'x',num2str(imSize(2)),'x',num2str(imSize(3))]);
    im = rwd2mat(name);
    trainingIm(:, :, :, 1, i) = im;
    if mod(i,200)==0
	disp([num2str(i/numel(InidzesOfData ) *100), ' Procent done'])
    end
end
disp( 'read in data done' );


trainingData = PermData(InidzesOfData );
trainingMincut = MinCut(InidzesOfData );

% restrict to certain value range
upperBound = 10^(-3);
lowerBound = 10^(-6);
indexIm = find(trainingData(1:end, 1) < upperBound & trainingData(1:end, 1) > lowerBound);

trainingIm = trainingIm(:, :, :, :, indexIm);
trainingData = trainingData(indexIm);
trainingMincut = trainingMincut(indexIm);


% split Data
numberOfData = size(trainingIm,5);
validationIm = trainingIm(:, :, :, :, ceil(PercentageOfTrainImages*numberOfData): end);
validationData = trainingData(ceil(PercentageOfTrainImages*numberOfData): end,1);
validationMincut = trainingMincut(ceil(PercentageOfTrainImages*numberOfData): end,1);

trainingIm = trainingIm(:, :, :, :, 1:ceil(PercentageOfTrainImages*numberOfData-1));
trainingData = trainingData(1:ceil(PercentageOfTrainImages*numberOfData-1),1);
trainingMincut = trainingMincut(1:ceil(PercentageOfTrainImages*numberOfData-1),1);


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
    maxPooling3dLayer(5, 'Stride', 5,'Name', 'MP2'); ...
    convolution3dLayer(3 * [1, 1, 1], 128, 'Padding', 1*[1,1,1], 'Name', 'conv3'); ...
    batchNormalizationLayer('Name', 'B3'); ...
    leakyReluLayer(0.1,'Name', 'relu3'); ...
    maxPooling3dLayer(2, 'Stride', 2,'Name', 'MP3'); ...
    fullyConnectedLayer(64,'Name', 'fc1'); ...
    depthConcatenationLayer(2,'Name','add');...
    leakyReluLayer(0.1,'Name', 'relu4'); ...
    fullyConnectedLayer(32, 'Name', 'fc2'); ...
    leakyReluLayer(0.1,'Name', 'relu5'); ...
    fullyConnectedLayer(1,'Name', 'fc3'); ...
    regressionLayer('Name', 'out')];


% Connect second input (physics input) to the first dense network layer

additionalLayers = [image3dInputLayer([1,1,1,1],'Name','physics'); fullyConnectedLayer(64,'Name','physics2','WeightsInitializer','ones','Biasinitializer','zero', 'WeightLearnRateFactor',0,'BiasLearnRateFactor',0)];
lgraph = layerGraph(layers);
lgraph = addLayers(lgraph,additionalLayers);
lgraph = connectLayers(lgraph,'physics2','add/in2');

%Group pore-space images, physics input and permeability labels in Datastores

a1=arrayDatastore(trainingIm,"IterationDimension",5);
a2=arrayDatastore(TrafoMinCuts(trainingMincut),"IterationDimension",1);
a3=arrayDatastore(TrafoPermtoLabel(trainingData),"IterationDimension",1);
training = combine(a1,a2,a3);

a1=arrayDatastore(validationIm,"IterationDimension",5);
a2=arrayDatastore(TrafoMinCuts(validationMincut),"IterationDimension",1);
a3=arrayDatastore(TrafoPermtoLabel(validationData),"IterationDimension",1);
validation = combine(a1,a2,a3);


%training options

opts = trainingOptions('sgdm', ...
    'MaxEpochs', 15, ...
    'L2Regularization',0.0001,...
    'InitialLearnRate', 0.002, ...
    'MiniBatchSize', miniBatchSize , ...	
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 4, ...
    'LearnRateDropFactor', 0.4, ...
    'ValidationFrequency', 500, ...
    'ExecutionEnvironment', 'multi-gpu', ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress', 'ValidationData', validation, ...
    'Verbose', true);


%Train the network
net = trainNetwork(training, lgraph, opts);

% evaluate training success

% Compute network prediction on validation data
predicted = zeros(size(validationIm,5),1);
for i=1:size(validationIm,5)
    a1=arrayDatastore(validationIm(:,:,:,1,i),"IterationDimension",5);
    a2=arrayDatastore(TrafoMinCuts(validationMincut(i)),"IterationDimension",1);
    a3=arrayDatastore(TrafoPermtoLabel(validationData(i)),"IterationDimension",1);
    temp = combine(a1,a2,a3);

    predicted(i) = TrafoLabeltoPerm(net.predict(temp));
end

% Compute network prediction on training data
predictedTrain = zeros(size(trainingIm,5),1);
for i=1:size(trainingIm,5)
    a1=arrayDatastore(trainingIm(:,:,:,1,i),"IterationDimension",5);
    a2=arrayDatastore(TrafoMinCuts(trainingMincut(i)),"IterationDimension",1);
    a3=arrayDatastore(TrafoPermtoLabel(trainingData(i)),"IterationDimension",1);
    temp = combine(a1,a2,a3);

   predictedTrain(i) = TrafoLabeltoPerm(net.predict(temp));
end

%save('CNNPhys2.mat','net','predicted','validationData', 'predictedTrain', 'trainingData','validationMincut','trainingMincut','-v7.3');


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
