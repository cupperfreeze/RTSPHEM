%> @file Grid/ascii.m This method pipes the grid properties into an ascii file.

%%% GRID INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ascii(obj)

% delete old file
delete 'gridProperties.txt'

p = properties(obj);
for iProp = 1:length(p)
    if ~isequal(p{iProp}, 'rectGrid') && ~isequal(p{iProp}, 'mapV') ...
            && ~isequal(p{iProp}, 'mapE') && ~isequal(p{iProp}, 'mapT')
        % print name of the property
        dlmwrite('gridProperties.txt', p{iProp}, '-append', 'delimiter', ' ');
        % print property data
        dlmwrite('gridProperties.txt', full(obj.(p{iProp})), '-append', 'delimiter', '\t');
    end
end

printline(3, 'Grid information piped to gridProperties.txt')

end % print
