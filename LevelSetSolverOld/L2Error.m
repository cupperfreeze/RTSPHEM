function [ l2Error, errorValues ] = L2Error( coarseGrid, coarseValues, fineGrid, fineValues )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % assert: 
    %   . ( nodesPerDimension of fineGrid - 1 ) / ( nodesPerDimension of
    %   coarseGrid - 1) is an array of natural numbers
    %   . coarseValues is a function defined on coarseGrid
    %   . fineValues is a function defined on fineGrid

    nodeRatios = ( fineGrid.nodesPerDimension - 1 ) ./ ...
        ( coarseGrid.nodesPerDimension - 1 );
    
    l2Error = 0;
    errorValues = zeros( size( coarseValues ) );
    
    for coarseIndex = 1:coarseGrid.nodes
        
        % Find corresponding node on fine grid.
        coarseSubscripts = coarseGrid.ind2sub( coarseIndex, : );
        fineSubscripts = nodeRatios .* ( coarseSubscripts - 1 ) + 1;
        fineSubscripts = num2cell( fineSubscripts );
        fineIndex = fineGrid.sub2ind( fineSubscripts{:} );
        
        error = abs( coarseValues( coarseIndex ) - fineValues( fineIndex ) );
        errorValues( coarseIndex ) = error;
        
        l2Error = l2Error + error^2;
        
    end
    
    l2Error = sqrt( l2Error ) * coarseGrid.stepSize(1);

end

