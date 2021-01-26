%> @file tp_indexsets.m Definition of index sets.
%>
%> Assumption: Edge types are the same for each species.
%>
%> @retval markDirE   [#E x 1 logical]
%> @retval markNeumE  [#E x 1 logical]
%> @retval markFluxE  [#E x 1 logical]
%> @retval markFreeE  Dirichlet + inner edges [#E x 1 logical]
%> @retval idxDirE
%> @retval idxNeumE
%> @retval idxFluxE
%> @retval idxFreeE


function [markDirE, markNeumE, markFluxE, markFreeE, idxDirE, idxNeumE, idxFluxE, idxFreeE, markDirV] = tp_indexsets(d)

g = d.grid;


% Dirichlet Edges
markDirE = false(g.numE, 1); % [#E x 1]
for k = 1:length(d.id2D) % fetch edge ids acc. to Dirichlet boundary
    markDirE = logical(markDirE+(d.id2D{k} == g.idE));
end


% Neumann Edges
markNeumE = false(g.numE, 1); % [#E x 1]
for k = 1:length(d.id2N)
    markNeumE = logical(markNeumE+(d.id2N{k} == g.idE));
end
% Flux Edges
markFluxE = false(g.numE, 1); % [#E x 1]
for k = 1:length(d.id2F)
    markFluxE = logical(markFluxE+(d.id2F{k} == g.idE));
end

markFreeE = ~(markFluxE + markNeumE);

% Boundary has to be well-defined
assert(sum(g.idE ~= 0) == sum(markDirE)+sum(markNeumE)+sum(markFluxE), ...
    'HyPHM: Some boundary edges were not assigned to a type.  Please set id2D, id2N and/or id2F.');

idxDirE = find(markDirE)'; % [1 x #E]  % ordered  Dirichlet edge numbers
idxFluxE = find(markFluxE)'; %    "       Flux          "
idxNeumE = find(markNeumE)'; %    "     Neumann         "
idxFreeE = find(markFreeE)'; %    "   inner and Dir     "


end
