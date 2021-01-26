%> @file stencilwalk.m Seek the triangle @f$T_\mathrm{aim}@f$ where the point @f$\vec{x}_\mathrm{seed}@f$ is located.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Seek the triangle @f$T_\mathrm{aim}@f$ where the point @f$\vec{x}_\mathrm{seed}@f$ is located.  The first search
%> starts at @f$T_\mathrm{seed}@f$.  If the triangle is not found, it is repeated
%> @c tries times with random seed.  After that, if @f$T_\mathrm{aim}@f$ is still not
%> located (maybe the requested point is outside the domain), @f$-1@f$ is
%> returned.

function Taim = stencilwalk(grid, xseek, Tseed, isVis, tries)

assert(isa(grid, 'AbstractGrid'))
assert(isequal(size(xseek), [2, 1]), 'HyPHM: xseek has to be a 2x1 vector on the considered domain.')
assert(isscalar(Tseed))

g = grid;


if isVis
    hold on
    pause(2)
    g.visualize('numT')
    tt = squeeze(g.coordV0T(Tseed, :, :));
    line([tt(1, 1), tt(2, 1), tt(3, 1), tt(1, 1)], [tt(1, 2), tt(2, 2), tt(3, 2), tt(1, 2)], 'Color', 'k', 'LineWidth', 3, 'LineStyle', '--')
    text(xseek(1), xseek(2), 'X', 'Color', 'r', 'FontSize', 18)
end

for Tseed = [Tseed, ceil(rand(1, tries - 1) * g.numT)]
    Taim = oneStencilWalk(g, xseek, Tseed, isVis);
    if Taim ~= -1 % if the correct triangle is found, return its number
        break
    end
end


end % function

function nextT = getnextT(g, curT, locE)
psblTs = g.T0E(g.E0T(curT, locE), :);
if psblTs(1) == curT, nextT = psblTs(2);
else nextT = psblTs(1);
end
end


function Taim = oneStencilWalk(g, xseek, Tseed, isVis)

assert(Tseed > 0 && Tseed <= g.numT)
barycoord = zeros(1, 3); % initialization
Tlist = [Tseed, Tseed]; % growing list of all travelled triangles (needs an extra entry since Tlist(end-1) will be requested)
while true

    % compute barycentric coordinates
    coordV = squeeze(g.coordV0T(Tlist(end), :, :));
    D = [coordV, ones(3, 1)];
    detD = det(D);
    for kV = 1:3
        C = D;
        C(kV, 1:2) = xseek;
        barycoord(kV) = det(C);
    end
    barycoord = barycoord / detD;

    if all(barycoord >= 0)
        Taim = Tlist(end);
        if isVis
            printline(1, 'Found T %d!', Taim)
            tt = squeeze(g.coordV0T(Taim, :, :));
            line([tt(1, 1), tt(2, 1), tt(3, 1), tt(1, 1)], [tt(1, 2), tt(2, 2), tt(3, 2), tt(1, 2)], 'Color', 'k', 'LineWidth', 3)
            %       title(sprintf('Found Triangle %d!', Taim))
        end
        break
    end

    if isVis, printline(1, 'Point [%.3f; %.3f] has on T %d the barycentric coordinates [%.3f, %.3f, %.3f].', ...
            xseek(1), xseek(2), Tlist(end), barycoord(1), barycoord(2), barycoord(3)); end

    % decide which edge is the preferred
    [~, idxs] = sort(barycoord);

    Told = Tlist(end);
    for kE = idxs
        % get new triangle in preferred direction
        tmpT = getnextT(g, Tlist(end), kE); %#ok<*AGROW>
        if isVis, printline(1, 'propose T %d for new way...', tmpT); end
        if tmpT == 0 || tmpT == Tlist(end-1) % if the proposed way is back or over boundary choose other way
            continue
        elseif tmpT == 0 || any(tmpT == Tlist(1:end - 3)) % if we run in a loop, xseek is in an obstacle
            Taim = -1;
            if isVis, printline(1, 'Loop end at T %d!', Tlist(end)); end
            %       if isVis, title('Loop detected!'); end
            return
        else % allowed travelling
            Tlist(end+1) = tmpT;
            break
        end
    end
    if Told == Tlist(end)
        error('HyPHM: Stencilwalk algorithm run into a dead end at T %d!', Tlist(end))
    end

    if isVis
        line(g.baryT([Tlist(end -1), Tlist(end)], 1), g.baryT([Tlist(end -1), Tlist(end)], 2), 'Color', 'b', 'LineWidth', 2)
        printline(1, 'Travel over local edge %d to T %d.', Tlist(end -1), Tlist(end))
    end
end

end
