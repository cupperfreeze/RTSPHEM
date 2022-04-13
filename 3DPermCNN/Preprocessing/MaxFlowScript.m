% Interpret connected fluid domain as undirected graph (Voxel == nodes, Edges == neighboring relations). 
% This script computes the minimal number of edges to be deleted to make
% the graph disconnected and other interesting quantities

p=genpath('VoxelData/VaryPoro');
addpath(p);
p=genpath('ComputedData');
addpath(p);
p=genpath('RawConversion');
addpath(p);

numberOfSamples = 35;
MiniCubeSize = [100,100,100];
[X,Y,Z] = meshgrid(1:MiniCubeSize(1),1:MiniCubeSize(2),1:MiniCubeSize(3));
grid = CartesianGrid3D(3,[0,1,0,1,0,1],MiniCubeSize-1);
MinCut = nan(numberOfSamples,1);
shortest = nan(numberOfSamples,1);
MinSpan  = nan(numberOfSamples,1);
Central = nan(numberOfSamples,1);
area = nan(numberOfSamples,1);
porosities = nan(numberOfSamples,1);

CentralRed = nan(numberOfSamples,1);
shortestRed = nan(numberOfSamples,1);
SpaceUsed = nan(numberOfSamples,1);
areaRed = nan(numberOfSamples,1);

for i=1:numberOfSamples

    A=rwd2mat(['VaryPoro',num2str(i-1),'_100x100x100']);
    %A=rwd2mat(['out_connected',num2str(i-1),'_100x100x100']);
    %A=rwd2mat('Analytics6x4_100x100x100');
    
    porosities(i) = sum(reshape(A, prod(MiniCubeSize),1))/ prod(MiniCubeSize);

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

    inverseOrder = zeros(prod(MiniCubeSize),1);
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
    [mf,GF,compnent1,component2] = maxflow(G,numNodes+1, numNodes+2) ;

    temp=unique(reshape(table2array(GF.Edges(:,1)),2*size(table2array(GF.Edges(:,1)),1),1));
    temp = temp(1:(end-2));
    C=zeros(100,100,100,'uint8');
    C(nodesFluid(temp))=1;
    SpaceUsed(i) = numel(nodesFluid(temp));

    MinCut(i)= mf;
    [~,d]=shortestpath(G,numNodes+1, numNodes+2);
    shortest(i) = d;

    [~,d]=shortestpath(GF,numNodes+1, numNodes+2);
    shortestRed(i) = d;

    temp = minspantree(G);
    MinSpan(i)  = numel(temp.Edges);

    temp=centrality(G,'pagerank');
    Central(i) = min(mean(temp(inverseOrder(leftIdx))),mean(temp(inverseOrder(rightIdx))));

    temp=centrality(GF,'pagerank');
    leftIdxRed = find(C(:)>0 & grid.coordinates(:,1)<eps);
    rightIdxRed = find(C(:)>0 & grid.coordinates(:,1)>1-eps);
    CentralRed(i) = min(mean(temp(inverseOrder(leftIdxRed))),mean(temp(inverseOrder(rightIdxRed))));




    out = isosurface(X,Y,Z,A);
    B = zeros(size(out.faces,1), 3,3);
    for j=1:size(out.faces,1)
        B(j,:,:) = out.vertices(out.faces(j,:),:);
    end

    B(:,2,:) = B(:,2,:) - B(:,1,:);
    B(:,3,:) = B(:,3,:) - B(:,1,:);
    area(i) = 0.5*sum(vecnorm(cross(squeeze(B(:,2,:)), squeeze(B(:,3,:)),2),2,2));

    out = isosurface(X,Y,Z,C);
    B = zeros(size(out.faces,1), 3,3);
    for j=1:size(out.faces,1)
        B(j,:,:) = out.vertices(out.faces(j,:),:);
    end

    B(:,2,:) = B(:,2,:) - B(:,1,:);
    B(:,3,:) = B(:,3,:) - B(:,1,:);
    areaRed(i) = 0.5*sum(vecnorm(cross(squeeze(B(:,2,:)), squeeze(B(:,3,:)),2),2,2));
    

end
 
save('SampleDataVaryPoro.mat','MinCut','shortest','MinSpan','Central','CentralRed','shortestRed','SpaceUsed','area','areaRed','porosity');
