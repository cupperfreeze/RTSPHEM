function [ output_args ] = plotLogical( grid, logicalField, fig )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    figure(fig);
    holdState = ishold;
    hold on;
    
    coord = grid.coordinates;
    xCoord = coord(:,1);
    xCoord = reshape( xCoord, grid.size )';
    xCoord(end,:) = [];
    xCoord(:,end) = [];
    yCoord = coord(:,2);
    yCoord = reshape( yCoord, grid.size )';
    yCoord(end,:) = [];
    yCoord(:,end) = [];
    
    for i = 1:numel( logicalField )
        
        if ( logicalField(i) )
            rectangle( 'Position', [ xCoord(i) yCoord(i) grid.stepSize ], ...
                'FaceColor', 'r' );
        end
        
    end

    if ( ~holdState )
        hold off;
    end
    
end

