function [ ] = plotGrid( grid, fig )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    assert( isa( grid, 'CartesianGrid' ) );
    figure(fig);
    holdState = ishold;
    hold on;

    coord = grid.coordinates;
    grey = 0.8 * ones(1,3);
    
    xCoordIndex = find( abs( coord(:,2) - 0.5 ) < eps );
    yCoordIndex = find( abs( coord(:,1) - 0.5 ) < eps );
    
    for i = 1:numel( xCoordIndex )
        plot( coord( xCoordIndex(i), 1 ) * ones(1,2), [-0.5 0.5], 'Color', grey );
    end
    
    for i = 1:numel( yCoordIndex )
        plot( [-0.5 0.5], coord( yCoordIndex(i), 2 ) * ones(1,2), ...
            'Color', grey );
    end
    
    if ( ~holdState )
        hold off;
    end

end

