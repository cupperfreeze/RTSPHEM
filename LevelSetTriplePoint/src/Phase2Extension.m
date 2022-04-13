% Extend normal velocity from interface to whole domain in 2 phase scenario

function [dV,S]=Phase2Extension(grid, gridHyPHM, levelSet, normalSpeed)

Xi = levelSet>0;

 numNodes = grid.nodes;
 isInitial = false( numNodes, 1 );
 dVInitial = inf(numNodes,1);
 val = NaN( 3, numNodes );
 val(1,:) = Xi';
 
    
%Label all Interface points by checking wheater forward/backward steps in
%any dimension change Xi    
    for d = 1:grid.dimension
        
        forwardIndex = grid.getForward( 1:numNodes, d );
        val( 2, ~isnan(forwardIndex) ) = ...
            Xi( forwardIndex(~isnan(forwardIndex)) );
        
        backwardIndex = grid.getBackward( 1:numNodes, d );
        val( 3, ~isnan(backwardIndex) ) = ...
            Xi( backwardIndex(~isnan(backwardIndex)) );
        
        %when Xi changes in forward step
        changeIndex = val(1,:)~=val(2,:) & ~isnan(val(1,:)) & ~isnan(val(2,:));
        isInitial(changeIndex) = true;
    
        %when Xi changes in backward step
        changeIndex = val(1,:)~=val(3,:) & ~isnan(val(1,:)) & ~isnan(val(3,:));
        isInitial(changeIndex) = true;

    end
    
    interfaceVelocitiesNodes = zeros(numNodes,1);
    BoundaryTriangles = abs(sum(levelSet(gridHyPHM.V0T(:,:))>0,2) - 2)<eps; 
    indexBoundaryTriangles = find(BoundaryTriangles);
    baryBoundaryTriangles = gridHyPHM.baryT(BoundaryTriangles,:);
     
    
     interfaceVelocitiesNodes = nan(numNodes,1);
    isInitial(gridHyPHM.numV+1:end)=[]; 

    for i = find(isInitial)'
       [~,idx] =min( 2*abs(baryBoundaryTriangles(:,1) - gridHyPHM.coordV(i,1)) + ...
                         abs (baryBoundaryTriangles(:,2) - gridHyPHM.coordV(i,2))); 
       interfaceVelocitiesNodes(i) = normalSpeed(indexBoundaryTriangles(idx)); 
    
    end
    
dVInitial(isInitial) = abs(levelSet(isInitial));
[dV,S] = reinitializeLevelSetWithS( grid, dVInitial,inf,  interfaceVelocitiesNodes(:,1)');
S(isinf(S))=0;

end

function out = expandNodes(val, grid)
	out = zeros(numel(grid.isFoldedNode),1);
	out(~grid.isFoldedNode ) = val;
	out( grid.isFoldedNode ) = out( grid.foldedIndex( grid.isFoldedNode ) );
end
