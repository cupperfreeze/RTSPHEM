% isolate connected fluid domain with respect to axis 'direction'
function [Out, deadNodes, connected] = fillholes3D(numV, coords, neighborIndices, image, direction)

connected = false;
Out = image;

length = (round(double(numV)^(1 / 3))); %side length of cube
ActiveGlobal = zeros(numV, 1); %1 if related index contributues to connected domain

%Iterate over number of connected domains
while true
    nodesVisited = zeros(numV, 1, 'logical'); %1 if node has already been considered
    activeNodes = zeros(numV, 1); %lists nodes sorted by number of propagation front they belong to
    idx = find(image < 0.5, 1); %find index belonging to fluid domain
    if numel(idx) < 1  %no fluid domain remaining --> done
        deadNodes = ~logical(ActiveGlobal);
        Out(deadNodes) = 1; %output = input with all disconnected porespace turned into solid
        return
    end

    Start = 1;
    End = 1;
    activeNodes(1) = idx;
    nodesVisited(idx) = 1;

    %keep track of coordinates of all nodes within thihs connected domain
    coordsTraveled = zeros(uint32(length), 3);

    %Iterate over front propagation steps
    while true

        localcount = 1;
        for m = 1:3
            coordsTraveled(round(coords(idx, m) * (length - 1) + 1), m) = 1;
        end

        %Iterate over all nodes within current front
        for i = Start:End
            neighbors = neighborIndices(activeNodes(i), :);
            %Iterate over all neighbors of current front node
            for j = 1:6
                if neighbors(j) == 0
                    continue
                end

                %if neighbor belongs to fluid domain and hasn't been
                %regarded yet, make it an active node
                if nodesVisited(neighbors(j)) == 0 && image(neighbors(j), 1) < 0.5
                    activeNodes(End+localcount) = neighbors(j);
                    nodesVisited(neighbors(j), 1) = 1;

                    localcount = localcount + 1;
                    for m = 1:3
                        coordsTraveled(round(coords(neighbors(j), m) * (length - 1) + 1), m) = 1;
                    end
                end
            end
        end

        %Calculate loop range such that all nodes in the current front are included
        Start = End + 1;
        End = End + localcount - 1;

        % if loop is empty, all nodes within connected domain visited -->
        % done
        if End < Start
            break;
        end

    end


    temp = sum(coordsTraveled);
    if temp(direction) == length %connected component connects both surfaces on the group along 'direction' axis
        connected = true;
        image(nodesVisited(:, 1)) = 1; %turn connected component to solid
        ActiveGlobal(nodesVisited(:, 1)) = 1; %mark nodes of this connected domain as active

    else
        image(nodesVisited(:, 1)) = 1; %turn connected component to solid
    end
end

end