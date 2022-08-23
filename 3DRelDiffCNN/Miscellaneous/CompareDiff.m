% Apply CNN1 (trained on fully saturated samples only) to partially
% saturated samples 

p=genpath('./VoxelData/SampleSat2');
%p=genpath('./VoxelData/SampleBENTHEIMER');

addpath(p);

for offset = 25000%[5000,10000,15000,20000,25000];

load('DiffSample2.mat');
load('CNNDiffPlain2.mat','net');
DiffData=DiffData(offset+(4501:5000));
imSize = [100,100,100];												%Indizes of data files
TrafoPermtoLabel = @(x) log(10.*x);
TrafoLabeltoPerm = @(x) exp(x)./10;
storage = zeros([imSize(1), imSize(2), imSize(3),1, 5], 'uint8');
predicted = nan(500,1);                                                        %Mini-batch Size for network training

for big = 1:50
	if mod(big,5)==1
		disp('check');
	end
	for i = 1:10
    		name = strcat('out_connected', num2str(offset+4500+(big-1)*10+i-1),['_',num2str(imSize(1)),'x',num2str(imSize(2)),'x',num2str(imSize(3))]);
    		im = rwd2mat(name);
    		storage(:, :, :, 1, i) = im;
   	end
	predicted((big-1)*10+(1:10)) = TrafoLabeltoPerm(net.predict(storage));
end


R_squredLOG = 1 - sum((log10(predicted)-log10(DiffData)).^2) / (sum((log10(DiffData)-mean(log10(DiffData))).^2))
R_squred = 1 - sum((predicted-(DiffData)).^2) / (sum(((DiffData)-mean((DiffData))).^2))
end


%CNNDiffPlain2:
%5000:   82.14%		91.17%
%10000:	 83.16%         91.31%
%15000:	 79.39%         86.67%
%20000:  59.84%		69.23%
%25000:  17.13%         34.88%			
