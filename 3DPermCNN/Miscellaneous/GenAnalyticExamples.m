% Produce different analytical examples with rectangular flow channels

p=genpath('./RawConversion');
addpath(p);

imSize=[100,100,100];                                                       %Size of subsamples

precission = 10;                                                            %Parameter for precission of analytical solution
A=zeros(imSize,'uint8');
x = [2, 5, 10, 10, 6, 6,6];                                                 %pairs of channel width and height in voxel
y = [2, 5, 1, 2, 3, 4,6];

correct = zeros(numel(x),1);


% Compute analytical permeability value on reference cube and compute error bound
for j = 1:numel(x);
    b=x(j)/100;
    h=y(j)/100;
    K1=0;
    for i=precission:-1:1
        K1=K1+(2*i-1)^(-5)*192/pi^5 *min(b,h)/max(b,h)*tanh((2*i-1)*pi/2 *max(b,h)/min(b,h));
    end
    K(j) = 1-K1;
end
error = 1/8 *(2*precission)^(-4)


%K=[0.4218,0.4218, 0.9370, 0.874, 0.6861, 0.5873];

%Write files and scale analytical results
for i=1:numel(x)
    B=A;
    B(:, 50:(50+x(i)-1), 50:(50+y(i)-1)) = 1;
    mat2rawrwd(B, ['Analytics',num2str(x(i)),'x',num2str(y(i))], 1, [0,0,0], [1,1,1], './Analytics/');
    correct(i) = K(i) * min(x(i)/imSize(2), y(i)/imSize(3))^3 *  max(x(i)/imSize(2), y(i)/imSize(3)) /12;
end
    
correct

