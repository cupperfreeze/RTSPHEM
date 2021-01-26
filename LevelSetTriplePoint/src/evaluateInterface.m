% Calculate interface and triple point position at a single time instance
function [interfaceLength, coordTripel] = evaluateInterface(grid, Xi, distFunctions, slt)


edges = grid.edges;
pointsOnInterface = inf(size(edges));
identifyInterface = inf(size(edges, 1), 1);
temp = grid.nodesPerDimension(1);

%the following functions are almost linear around the respective interfaces
%and change their sign at the interface
d1 = distFunctions(:, 1) - distFunctions(:, 2);
d2 = distFunctions(:, 1) - distFunctions(:, 3);
d3 = distFunctions(:, 2) - distFunctions(:, 3);

CheckEdges = find(Xi(edges(:, 1)) ~= Xi(edges(:, 2)));
for i = CheckEdges'
    %Check for each edge wheather it crosses interface
    %If yes, calculate alpha and calculate approximate coordinates
    %of interface point

    if ((Xi(edges(i, 1)) == 1 & Xi(edges(i, 2)) == 2) | (Xi(edges(i, 1)) == 2 & Xi(edges(i, 2)) == 1))
        coords = grid.coordinates(edges(i, :), :);
        alpha = d1(edges(i, 2)) / (d1(edges(i, 2)) - d1(edges(i, 1)));
        pointsOnInterface(i, :) = coords(1, :) * alpha + coords(2, :) * (1 - alpha);
        identifyInterface(i) = 1;
    end

    if ((Xi(edges(i, 1)) == 3 & Xi(edges(i, 2)) == 2) | (Xi(edges(i, 1)) == 2 & Xi(edges(i, 2)) == 3))
        coords = grid.coordinates(edges(i, :), :);
        alpha = d3(edges(i, 2)) / (d3(edges(i, 2)) - d3(edges(i, 1)));
        pointsOnInterface(i, :) = coords(1, :) * alpha + coords(2, :) * (1 - alpha);
        identifyInterface(i) = 3;
    end

    if ((Xi(edges(i, 1)) == 3 & Xi(edges(i, 2)) == 1) | (Xi(edges(i, 1)) == 1 & Xi(edges(i, 2)) == 3))
        coords = grid.coordinates(edges(i, :), :);
        alpha = d2(edges(i, 2)) / (d2(edges(i, 2)) - d2(edges(i, 1)));
        pointsOnInterface(i, :) = coords(1, :) * alpha + coords(2, :) * (1 - alpha);
        identifyInterface(i) = 2;
    end
end

%plot
if ~slt
    for j = 1:3
        scatter(pointsOnInterface(identifyInterface == j, 1), pointsOnInterface(identifyInterface == j, 2));
        hold on
    end
    axis equal
    %  title('Interface');
    axis([-0.5, 0.5, -0.5, 0.5]);

    set(gca, 'FontSize', 15)
    hold on
end
%position 1 corresponds to transition (Xi=1, X=2)
%position 2 corresponds to transition (Xi=1, X=3)
%position 3 corresponds to transition (Xi=2, X=3)
interfaceLength = [0, 0, 0];
coordTripel = nan(0, 2);

CheckTriangles = find(sum(~isinf(identifyInterface(grid.triangleEdges)), 2));
for i = CheckTriangles'
    indizes = grid.triangleEdges(i, :);

    cutEdges = ~isinf(identifyInterface(indizes));

    if sum(cutEdges) == 2
        localIndex = indizes(find(cutEdges));
        lengthPart = norm(pointsOnInterface(localIndex(1), :)- ...
            pointsOnInterface(localIndex(2), :));

        interfaceLength(identifyInterface(localIndex(1))) = ...
            interfaceLength(identifyInterface(localIndex(1))) + lengthPart;

    end

    if sum(cutEdges) == 3
        %First approximation of tripel point position
        coordTripelApprox = mean(pointsOnInterface(indizes, :), 1);

        %within the area around this approximation, calculate distance
        %functions on finer grid via linear interpolation

        [Xl, Yl] = meshgrid(linspace(-0.5, 0.5, temp), linspace(-0.5, 0.5, temp));
        [Xs, Ys] = meshgrid(linspace(coordTripelApprox(1) - 1 / temp, coordTripelApprox(1) + 1 / temp, 100), linspace(coordTripelApprox(2) - 1 / temp, coordTripelApprox(2) + 1 / temp, 100));
        d1New = interp2(Xl, Yl, reshape(distFunctions(:, 1), temp, temp)', Xs, Ys);
        d2New = interp2(Xl, Yl, reshape(distFunctions(:, 2), temp, temp)', Xs, Ys);
        d3New = interp2(Xl, Yl, reshape(distFunctions(:, 3), temp, temp)', Xs, Ys);

        %Choose tripel point as Minimizer of abs(d^1-d^2) + abs(d^3-d^2)
        [~, k] = min(abs(d1New(:) - d2New(:))+abs(d3New(:) - d1New(:)));
        coordTripelApprox = [Xs(k), Ys(k)];
        coordTripel = [coordTripel; coordTripelApprox];

        %Add corresponding lengths
        for j = 1:3
            lengthPart = norm(pointsOnInterface(indizes(j), :)-coordTripel(end, :));
            interfaceLength(identifyInterface(indizes(j))) = ...
                interfaceLength(identifyInterface(indizes(j))) + lengthPart;
        end
    end
end

if size(coordTripel, 1) < 0.5 & ~slt
    display('Tripel point not found');
end


%plot tripel points
if ~slt
    scatter(coordTripel(:, 1), coordTripel(:, 2), 100, 'MarkerFaceColor', [0, .7, .7]);
    hold on
    temp = grid.nodesPerDimension(1);
    coord = grid.coordinates;
    out = double(Xi);
    out(Xi == 1) = 0;
    out(Xi == 2) = -0.5;
    out(Xi == 3) = -1;
    %surf(reshape(coord(:,1),temp,temp),reshape(coord(:,2),temp,temp),reshape(out,temp,temp),'EdgeColor','none');
    colormap('gray');
    set(gca, 'Layer', 'top')
    scatter(coordTripel(:, 1), coordTripel(:, 2), 100, 'MarkerFaceColor', [0, .7, .7]);
end
end
