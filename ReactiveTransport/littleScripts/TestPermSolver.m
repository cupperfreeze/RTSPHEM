% Script to fool around with permeability solver

global EPS
EPS = eps;

for i = 1:1
    dimension = 2;
    numPartitions = ceil(50*2^((i - 1) / 2));
    cellGrid = FoldedCartesianGrid(dimension, ...
        kron(ones(1, dimension), [-0.5, 0.5]), ...
        numPartitions*ones(1, dimension));
    helpGridHyPHM = Grid(cellGrid.coordinates, cellGrid.triangles);
    coord = cellGrid.coordinates;
    alpha = pi / 6;
    rot = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];

    % circle     0.019903637086084  -0.000000371117165  -0.000000371117165   0.019903637086084
    %initialLevelSetFunc = @(x) 0.25-norm(x-[0.0,0.0],2) ;%radius-norm(x,2)-norm(x-[0.2,0.1],3);

    % square       0.013036823857188  -0.000003698234083  -0.000003698234083   0.013036823857187
    % initialLevelSetFunc = @(x) 0.25-norm(x-[0.0,0.0],inf) ;

    % plus 0.017293803828765  -0.000003935984802   -0.000003935984802   0.017293803828764
    % initialLevelSetFunc = @(x) max(0.25-norm([2*x(1),x(2)]-[0.0,0.0],inf),0.25-norm([x(1),2*x(2)]-[0.0,0.0],inf)) ;

    % flower     0.024138621159248  -0.000000529198813  -0.000000529198813   0.024138621159248
    % initialLevelSetFunc = @(x) max(0.25-norm([2*x(1),x(2)]-[0.0,0.0],2),0.25-norm([x(1),2*x(2)]-[0.0,0.0],2)) ;

    % tilted ellipse   0.031019829978579  -0.006224203434641  -0.006224203453558   0.044629816625863
    % initialLevelSetFunc = @(x) 0.25-norm([2;1].*(rot*[x(1);x(2)]),2) ;

    % tilted ellipse shifted   0.031015912959508  -0.006222807775769  -0.006222807810015   0.044625250832434

    initialLevelSetFunc = @(x) 0.01 - norm([x(1), x(2) - 0.00], inf);
    %
    coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
    levelSet = cellfun(initialLevelSetFunc, coordCell);
    tic
    for j = 1:1
        computePermeabilityTensor(helpGridHyPHM, levelSet)
    end
    toc
end
radius = 0.15;
Testset = zeros(64*64, 800);
TestIm = zeros(64, 64, 1, 800);

%  for i=1:800
%    initialLevelSetFunc = @(x) radius - norm([x(1),x(2)-0.00],2);
%
%     coordCell = mat2cell( coord, ones(1,cellGrid.nodes), dimension );
%     levelSet = cellfun( initialLevelSetFunc, coordCell );
%     Testset(:,i) = reshape(TrainingIm(:,:,1,i),4096,1)-0.5;
%     TestIm(:,:,1,i) = reshape(levelSet>0,64,64);
%     radius = radius +0.00125*i;
%  end
%
%  tic
% parfor i=1:800
%      computePermeabilityTensor(helpGridHyPHM,Testset(:,i),'Bubble');
%  end
%  toc
%
%  TestIm =cat(4,TestIm,TestIm);
%  tic
%  net.predict(TestIm);
%  toc
%     [ gridNew, levelSetNew ] = localMeshRefinement(helpGridHyPHM, ...
%     levelSet );


% for i=1:1
%
% % out=cell(1,40);
% %
% % %spmd, mpiprofile('on'); end
% tic
% for i=1:1
% a=computePermeabilityTensor(helpGridHyPHM,levelSet,'noRefine','Bubble')
% end
% toc
% tic
% for i=1:1
% b=computePermeabilityTensor(helpGridHyPHM,levelSet,'noRefine')
% end
%
% toc
%          norm(a...
%                 - b )/...
%                 norm(a)
% %spmd, mpiprofile('viewer'); end
% end

% [a,b]=meshgrid(-0.5:1/numPartitions:0.5,-0.5:1/numPartitions:0.5);
%     scatter(a(levelSet>0),b(levelSet>0),'filled', 'MarkerFaceColor','b')
%         axis([-0.5 0.5 -0.5 0.5])
%         axis equal
%         hold on
%         contour(a,b,reshape(levelSet, numPartitions+1, numPartitions+1),[0,0], 'linewidth',5, 'color', 'b')
%  set(gca,'xtick',[]);
%  set(gca,'ytick',[]);
%tocBytes(gcp);
