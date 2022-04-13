% Evaluate prediction time by CNN

tic

% Setup timers          %1000 samples		10000 samples
SetupTime = 0;          %15.8197            9.4544
ReadInTime = 0;         %7.5627             51.3034
PhysicsTime = 0;        %87.8383            865.7195
PredictionTime = 0;     %45.5129            367.8259
EvalTime = 0;           %0.0319             0.0774


% Preparation (add file pathsa, load data, prepare variables)
p=genpath('./VoxelData/SampleBENTHEIMER');
addpath(p);
p=genpath('RawConversion');
addpath(p);
p=genpath('ComputedData');
addpath(p);
p=genpath('CNN');
addpath(p);


load('CNNPhys2.mat','net');

numberOfDataGlobal = 10000;                                                 % total number of samples to predict
numberOfData = 1000;                                                        % number of samples to load per iteration (in case of RAM restrictions)
minibatchSize = 10;                                                         % batch size for GPU inferencing
TrafoPermtoLabel = @(x) 2*(log10(x)+6) ;									% permeability label transformations
TrafoLabeltoPerm = @(x) 10.^(x/2-6);
TrafoMinCuts = @(x) TrafoPermtoLabel(x.^1.407*10^(-8.183));
imSize = [100,100,100];                                                     % dimensions of subsample
predictedTrainGlobal = zeros(numberOfDataGlobal,1);

[X,Y,Z] = meshgrid(1:imSize(1),1:imSize(2),1:imSize(3));
grid = CartesianGrid3D(3,[0,1,0,1,0,1], imSize-1);


SetupTime = SetupTime + toc;
for globCount = 1:ceil(numberOfDataGlobal / numberOfData )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in Data
    tic
    disp( 'read in data' );
    if globCount == ceil(numberOfDataGlobal / numberOfData )
        numberOfData = numberOfDataGlobal - (globCount -1)*numberOfData; 
    end
    
    inputIm = zeros([imSize(1), imSize(2), imSize(3), 1, numberOfData], 'uint8');
    MinCut = ones(numberOfData,1);
    
    parfor i = 1:numberOfData
        name = strcat('out_connected', num2str((globCount - 1)*numberOfData +i-1),['_',num2str(imSize(1)),'x',num2str(imSize(2)),'x',num2str(imSize(3))]);
        %name = strcat('VaryPoro', num2str((globCount - 1)*numberOfData +i-1),['_',num2str(imSize(1)),'x',num2str(imSize(2)),'x',num2str(imSize(3))]);

        im = rwd2mat(name);
        inputIm(:, :, :, 1, i) = im;
        if mod(i,200)==0
            disp([num2str(i/numberOfData *100), ' Procent done'])
        end
    end
    disp( 'read in data done' );
    ReadInTime = ReadInTime + toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform physics precomputations (graph flow problems)

    parfor i=1:numberOfData

        A=inputIm(:,:,:,1,i);

        nodesFluidLogical = (A(:)>0);
        nodesFluid = find(nodesFluidLogical);
        numFluid = numel(nodesFluid);

        %Indizes of nodes on left and right boundary
        leftIdx = find(A(:)>0 & grid.coordinates(:,1)<eps);
        rightIdx = find(A(:)>0 & grid.coordinates(:,1)>1-eps);


        %in the end, the pair (nodes,edges) contains the vertices of all edges of
        %the graph
        neighborIndices = uint32([grid.pull{1, 1, 1}, grid.pull{2, 1, 1}, grid.pull{1, 1, 2}, grid.pull{2, 1, 2},grid.pull{1, 1, 3}, grid.pull{2, 1, 3}]);
        neighborIndices(neighborIndices<(1:size(neighborIndices,1))')=0;
        neighborIndices = neighborIndices(nodesFluid,:)';
        edges = neighborIndices(:);
        nodes = kron(nodesFluid,ones(6,1));

        inverseOrder = zeros(prod(imSize),1);
        inverseOrder(nodesFluid) = 1:numel(nodesFluid);
        idx = (edges>0);
        nodes = nodes(idx);
        edges = edges(idx);

        idx = nodesFluidLogical(edges);
        nodes = inverseOrder(nodes(idx));
        edges = inverseOrder(edges(idx));
        weights = ones(numel(nodes),1);

        % add one node connecting to all nodes on the left side, and other
        % connecting the the right hand side
        numNodes = max([nodes;edges]);
        nodes = uint32([nodes; (numNodes+1)*ones(numel(leftIdx),1); (numNodes+2)*ones(numel(rightIdx),1)  ]);
        edges = uint32([edges; inverseOrder(leftIdx); inverseOrder(rightIdx)]);

        % Setup graph and and calculate maxflow (=minimal cut)
         G=graph(nodes, edges, [weights;ones(numel(leftIdx)+numel(rightIdx),1)]);
         mf = maxflow(G,numNodes+1, numNodes+2) ;
         MinCut(i)= mf;
    end

    PhysicsTime = PhysicsTime + toc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform CNN inferencing

    tic
    predicted = zeros(size(inputIm,5),1);

    for i=1:ceil(size(inputIm,5)/minibatchSize )
        idx=(((i-1)*minibatchSize +1):min(i*minibatchSize, numberOfData)); 
        a1=arrayDatastore(inputIm(:,:,:,1,idx),"IterationDimension",5);
        a2=arrayDatastore(TrafoMinCuts(MinCut(idx)),"IterationDimension",1);
        temp = combine(a1,a2);

       predicted(idx) = TrafoLabeltoPerm(net.predict(temp,'ExecutionEnvironment','gpu'));
    end
    PredictionTime = PredictionTime + toc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    predictedTrainGlobal(((globCount -1)*numberOfData+1):min(globCount *numberOfData,numberOfDataGlobal)) = predicted;
end
% Compare CNN prediction against ground-truth permeabiity data

tic
load('PermDataBENTHEIMEROrd1.mat');
validationData = PermData(1:numberOfDataGlobal);
R_squred = 1 - sum((predictedTrainGlobal-(validationData)).^2) / (sum(((validationData)-mean((validationData))).^2))
R_squredLog = 1 - sum((log10(predictedTrainGlobal)-log10(validationData)).^2) / (sum((log10(validationData)-mean(log10(validationData))).^2))
EvalTime = toc;

SetupTime 
ReadInTime 
PhysicsTime 
PredictionTime 
EvalTime 
