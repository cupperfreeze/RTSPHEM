% Same as 'step1' for arbitrary number of phases
function [distFunctions, Xi] = step1General(Xi, Phi, grid, restrictDist)

%Initialize n subcells with Phi for n phases
numberOfPhases = max(Xi);
subs = cell(1, numberOfPhases);
[subs{:}] = deal(Phi);
distFunctions = zeros(size(Phi, 1), numberOfPhases);

%Eliminate all interface parts not belonging to Xi=i by setting to 1 in
%these areas

for i = 1:numberOfPhases
    subs{i}(Xi ~= i) = 1;
    distFunctions(:, i) = reinitializeLevelSet(grid, subs{i}, false, restrictDist);
end

XiChangeIndex = logical(sum(~isinf(distFunctions(:, :)), 2));
[~, XiNew] = min(distFunctions, [], 2); %Update characteristic function Xi

Xi(XiChangeIndex) = XiNew(XiChangeIndex);
end