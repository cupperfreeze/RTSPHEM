function [out, correction] = Continuation(speciesCells, timestep, Levels, depth)

% Compute a continuous extension of the concentration fields into the solid
% phase for a cell array of transport instances. The depth parameter controls
% the thickness of the extension stripe. Out consists of transport
% instances with the extended concentrations, correction contains the added
% values.

grid = speciesCells{1}.grid;
g = grid.T0E;

% old concentration values
old = zeros(grid.numT, numel(speciesCells));
for i = 1:numel(speciesCells)
    old(:, i) = speciesCells{i}.U.getdata(timestep-1);
end


correction = zeros(grid.numT, numel(speciesCells));
for sp = 1:numel(speciesCells)
    for step = 1:depth
        c = speciesCells{sp}.U.getdata(timestep-1);
        normalize = max(c);
        indizes = find((g(:, 1) > 0).*(g(:, 2) > 0)); %           indizes of edges connecting 2 triangles
        indizes = indizes(logical(abs(c(g(indizes, 1)) .* c(g(indizes, 2))) <= 2 * eps * normalize)); % indizes of edges such that one adjacet triangle has concentration 0

        %         correction(g(indizes,1),sp) = max([correction(g(indizes,1),sp), c(g(indizes,1)),c(g(indizes,2)) ],[],2);
        %         correction(g(indizes,2),sp) = max([correction(g(indizes,2),sp), c(g(indizes,1)),c(g(indizes,2)) ],[],2);

        % copy concentration values over the edges into the solid
        for i = indizes'
            if Levels(grid.V0E(i, 1)) > -eps && Levels(grid.V0E(i, 2)) > -eps
                if abs(c(g(i, 1))) >= 2 * eps;
                    correction(g(i, 2), sp) = max(correction(g(i, 2), sp), c(g(i, 1)));
                else
                    correction(g(i, 1), sp) = max(correction(g(i, 1), sp), c(g(i, 2)));
                end
            end
        end

        % correction(logical(correction(:,sp)>0),sp) = mean(correction(logical(correction(:,sp)>0),sp));
        speciesCells{sp}.U.setdata(timestep-1, max(c, correction(:, sp)));
        correction(:, sp) = max(c, correction(:, sp)) - old(:, sp);
    end
    %norm(speciesCells{sp}.U.getdata(timestep-1) - old(:,sp));

end
out = speciesCells;
end