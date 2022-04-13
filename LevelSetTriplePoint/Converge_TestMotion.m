% Evaluate order of convergence for Voronoi Implicit Interface Method
% for normal advection
% results of this script are saved in Data file 'LevelConv2eps.mat'
% cf. [2] Section 4.1

% k = 7; %number of refinement levels
% containerXi = cell(k, 1);
% containerPhi = cell(k, 1);
% containerDistFunctions = cell(k, 1);
% containerInterfaceLength = cell(k, 1);
% containerCoordTripel = cell(k, 1);
% containerTime = cell(k, 1);
% 
% for l = 1:k
%     %Full test case
%     close all
%     %parameter
%     StepSizeTime = 1 / 40 / 2^(l - 1);
%     numTimeSteps = 10 * 2^(l - 1);
% 
%     %setup grid
%     dimension = 2;
%     numPartitions = 20 * 2^(l - 1);
%     cellGrid = CartesianGrid(dimension, ...
%         kron(ones(1, dimension), [-0.5, 0.5]), ...
%         numPartitions*ones(1, dimension));
%     coord = cellGrid.coordinates;
%     temp = numPartitions + 1; %temporary variable for plotting
% 
%     % Initial indicator function for different subdomains
%     Xi = ones(size(coord, 1), 1);
%     Xi((coord(:, 2) < 0) & (coord(:, 1) <= 0)) = 2;
%     Xi((coord(:, 2) < 0) & (coord(:, 1) > 0)) = 3;
% 
%     % Velocities at interfaces
%     % interfaceVelocities(i,j) = velocity between Xi=i and Xi=j
%     interfaceVelocities = [0, -0, -1; 0, 0, 0; 0, 0, 0];
%     interfaceVelocities = interfaceVelocities - interfaceVelocities'; %add lower half
% 
%     %Initial configuration Phi
%     epsilon = 2 / numPartitions;
%     initialPhiFunc = @(x) max((epsilon-norm([x(1) + eps], 2))*(x(2) < 0)-10*(x(2) >= 0), epsilon-norm([x(2) + eps], 2));
%     coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
%     Phi = cellfun(initialPhiFunc, coordCell);
% 
%     %Calculate distance function d^i only up to this value
%     restrictDist = 2 * epsilon;
%     tic
%     for i = 1:numTimeSteps
%         %Calculate distance function from each subdomain
%         [distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);
% 
% 
%         %Reconstruct Vonoroi Interface and velocity extension
%         [dV, S] = step234(distFunctions, cellGrid, Xi, interfaceVelocities, 15/numPartitions);
% 
%         %Update Phi, but only every 10 steps
%         if mod(i, 10) == 0
%             Phi = epsilon - dV;
%         end
% 
%         %Evolve Phi
%         Phi = levelSetEquationTimeStep(StepSizeTime, 0, Phi, ...
%             cellGrid, S', 2);
% 
%     end
%     containerTime{l} = toc;
%     containerTime{l}
%     %plot evolved functions
%     [distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);
%     [interfaceLength, coordTripel] = evaluateInterface(cellGrid, Xi, distFunctions, false)
% 
% 
%     containerXi{l} = Xi;
%     containerPhi{l} = Phi;
%     containerDistFunctions{l} = distFunctions;
%     containerInterfaceLength{l} = interfaceLength;
%     containerCoordTripel{l} = coordTripel;
%     hold on
%     Circle = 0.25 - sqrt(coord(:, 1).^2+coord(:, 2).^2);
%     contour(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Circle, temp, temp), [-0.1, 0, 1]);
% end

for l=1:k
dimension = 2;
numPartitions= 20*2^(l-1); 
cellGrid = CartesianGrid( dimension, ...
kron( ones(1,dimension), [-0.5, 0.5] ), ...
numPartitions*ones(1,dimension) );
[interfaceLength, coordTripel,varargout]=evaluateInterfaceReference(cellGrid, containerXi{l}, containerDistFunctions{l}, true,1);
containerError{l,1}=varargout;
end

figure
subplot(2,2,1)


dimension = 2;
numPartitions= 20*2^(k-1)-1; 
cellGrid = CartesianGrid( dimension, ...
kron( ones(1,dimension), [-0.5, 0.5] ), ...
numPartitions*ones(1,dimension) );
coord = cellGrid.coordinates;
temp = numPartitions+1; epsilon = 2/numPartitions;
initialPhiFunc = @(x) max((epsilon - norm([x(1)+eps],2))*(x(2)<0) -10*(x(2)>=0),epsilon - norm([x(2)+eps],2));    
coordCell = mat2cell( coord, ones(1,cellGrid.nodes), dimension );
Phi = cellfun( initialPhiFunc, coordCell );
restrictDist = 3*epsilon;  
Xi = ones(size(coord,1),1);
Xi((coord(:,2)<0) & (coord(:,1)<=0)) = 2;
Xi((coord(:,2)<0) & (coord(:,1)>0)) = 3;
[dist,Xi]=step1(Xi, Phi, cellGrid, restrictDist);
[interfaceLength, coordTripel] = evaluateInterfaceReference(cellGrid, Xi, dist, false)

w=2;
x = linspace(0,0.5);
y=0*x;
line(x,y,'Color',[0.5,0,0],'LineWidth',w)

x = linspace(-0.5,0,10);
y=0.*x;
line(x,y,'Color',[0,0,0.5],'LineWidth',w)

y = linspace(-0.5,0,100);
x=0.*y;
line(x,y,'Color',[0.5,0.5,0],'LineWidth',w)

plot(0,0,'bo','MarkerSize',5,'MarkerFaceColor',0.5*[0.466,0.674, 0.188],'MarkerEdgeColor',0.5*[0.466,0.674, 0.188])

xlim([-0.5,0.5])
ylim([-0.5,0.5])

xlabel('x')
ylabel('y')
set(gca,'FontSize',15);


subplot(2,2,3)
temp = cell2mat(containerError);
lul = [containerCoordTripel{:}];
plot(1:k,log2(abs(-0.25-lul(1:2:2*k))),'--o','LineWidth',3);
hold on
plot(1:k,log2(sum(temp,2)),'--o','LineWidth',3);
hold on
plot(1:k,-(1:k)-1,'LineWidth',3,'Color',[0,0,0])
legend('triple point','interfaces','linear rate')
xlabel('refinement level');
ylabel('log2 distance');
set(gca,'FontSize',15);
set(gca,'FontSize',15);

subplot(2,2,2)

l=7;
dimension = 2;
numPartitions= 20*2^(l-1); 
cellGrid = CartesianGrid( dimension, ...
kron( ones(1,dimension), [-0.5, 0.5] ), ...
numPartitions*ones(1,dimension) );
[interfaceLength, coordTripel]=evaluateInterfaceReference(cellGrid, containerXi{l}, containerDistFunctions{l}, false);

hold on

w=2;
x = [0,0.5];
y=0*x+0.25;
line(x,y,'Color',[0.5,0,0],'LineWidth',w)

x = linspace(-0.5,-0.25,10);
y=0*x;
line(x,y,'Color',[0,0,0.5],'LineWidth',w)

x = linspace(-0.25,0,100);
y=sqrt(0.25^2-x.^2);
line(x,y,'Color',[0.5,0,0],'LineWidth',w)

x = [0,0];
y=[-0.5,0];
line(x,y,'Color',[0.5,0.5,0],'LineWidth',w)

x = [-0.25,0];
y=[0,0];
line(x,y,'Color',[0.5,0.5,0],'LineWidth',w)


plot(-0.25,0,'bo','MarkerSize',5,'MarkerFaceColor',0.5*[0.466,0.674, 0.188],'MarkerEdgeColor',0.5*[0.466,0.674, 0.188])

xlim([-0.5,0.5])
ylim([-0.5,0.5])

xlabel('x')
ylabel('y')
set(gca,'FontSize',15);

subplot(2,2,4)
plot(log2(abs(cell2mat(containerInterfaceLength)-[0.25,0.5+0.25*pi/2,0.75])),'--o','LineWidth',3)
hold on
plot(log2(1./2.^(1:7))-1,'LineWidth',3,'Color',[0,0,0])
legend('interface','interface','interface','linear rate')
xlabel('refinement level');
ylabel('log2 difference length');
set(gca,'FontSize',15);
set(gca,'FontSize',15);