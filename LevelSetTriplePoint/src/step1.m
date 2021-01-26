% First step of VIIM: calculate distance functions to respective connected
% components of \Gamma_\epsilon

function [distFunctions, Xi] = step1(Xi, Phi, grid, restrictDist)

%Initialize 3 subcells with Phi
subs = cell(1, 3);
[subs{:}] = deal(Phi);
distFunctions = zeros(size(Phi, 1), 3);

%Eliminate all interface parts not belonging to Xi=i by setting to 1 in
%these areas

for i = 1:3
    for j = 1:3
        if i ~= j
            subs{i}(Xi == j) = 1;
        end
    end
    distFunctions(:, i) = reinitializeLevelSet(grid, subs{i}, false, restrictDist);

end

XiChangeIndex = logical(~isinf(distFunctions(:, 1))+ ~isinf(distFunctions(:, 2))+ ~isinf(distFunctions(:, 3)));
[~, XiNew] = min(distFunctions, [], 2); %Update characteristic function Xi

Xi(XiChangeIndex) = XiNew(XiChangeIndex);
end