% Perform CNN inferencing on a given data set and compare to simulated Perm
% values

% R^2 values obtained (standard, log):
% BEREA: 	Ord0:   0.8878	0.9075    Ord1:    0.8943	0.9413
% CASLTE:	Ord0:   0.9402	0.9222	  Ord1:    0.9455   0.9451	  


% file paths 	
p=genpath('./VoxelData/SampleCASTLE');
addpath(p);
p=genpath('RawConversion');
addpath(p);
p=genpath('CNN');
addpath(p);
p=genpath('ComputedData');
addpath(p);

numberOfData = 200;                                                         %number of elements in data set
TrafoPermtoLabel = @(x) 2*(log10(x)+6) ;
TrafoLabeltoPerm = @(x) 10.^(x/2-6);
TrafoMinCuts = @(x) TrafoPermtoLabel(x.^1.407*10^(-8.183));
imSize = [100,100,100];

%load ground-truth data, physics data and CNN

load('PermDataCASTLEOrd1.mat');
load('SampleDataCASTLE.mat');
load('CNNPhys2.mat','net');


inputIm = zeros([imSize(1), imSize(2), imSize(3), 1, numberOfData], 'uint8');

% read in Data
disp( 'read in data' );
for i = 1:numberOfData
    name = strcat('out_connected', num2str(i-1),['_',num2str(imSize(1)),'x',num2str(imSize(2)),'x',num2str(imSize(3))]);
    %name = strcat('VaryPoro', num2str(i-1),['_',num2str(imSize(1)),'x',num2str(imSize(2)),'x',num2str(imSize(3))]);
	
    im = rwd2mat(name);
    inputIm(:, :, :, 1, i) = im;
    if mod(i,200)==0
	disp([num2str(i/numberOfData *100), ' Procent done'])
    end
end
disp( 'read in data done' );

validationData = PermData(1:numberOfData);
inputMincut = MinCut(1:numberOfData);


predicted = zeros(size(inputIm,5),1);
tic

% Perform CNN inferencing

for i=1:size(inputIm,5)
    a1=arrayDatastore(inputIm(:,:,:,1,i),"IterationDimension",5);
    a2=arrayDatastore(TrafoMinCuts(inputMincut(i)),"IterationDimension",1);
    temp = combine(a1,a2);

   predicted(i) = TrafoLabeltoPerm(net.predict(temp));
end
toc

% Compute R^2 value in natural and logarithmic scale

R_squred = 1 - sum((predicted-(validationData)).^2) / (sum(((validationData)-mean((validationData))).^2))
R_squredLog = 1 - sum((log10(predicted)-log10(validationData)).^2) / (sum((log10(validationData)-mean(log10(validationData))).^2))

save('CompareCASTLE.mat','predicted','validationData')
