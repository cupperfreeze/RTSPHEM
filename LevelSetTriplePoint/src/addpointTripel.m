function [deactNode, add] = addpoint(vertices, coord, Levels)

%add nodes near interface if splitting a cell is stable, otherwise
%deactivate close neighbour

posVertices = round(sum(Levels > eps));
add = [100, 100]; % dummy values, sorted out later
deactNode = 0;
delta = 0.15;

switch posVertices
    case 2
        return

    case 0
        return

    case 1

        if abs(Levels(1)-Levels(2)) > 100 * eps
            alpha = Levels(2) / (Levels(2) - Levels(1));
            if alpha < delta
                alpha = delta;
            end
            if alpha > 1 - delta
                alpha = 1 - delta;
            end

            %    if abs(alpha-0.5)<0.5 - delta                                                %regularisation
            add = [alpha * coord(1, :) + (1 - alpha) * coord(2, :)];
            %    end

            if ((alpha > 1-delta) & (Levels(1) < 0)) %shut off points too close to boundary
                deactNode = vertices(1);
            end

            if ((alpha < delta) & (Levels(2) < 0))
                deactNode = vertices(2);
            end

        end


end


end