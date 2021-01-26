%> @file Grid/export.m Export a HyPHM Grid to an M++ grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%>
%> The data is piped into the M++ file ./[name].geo.@n
%> @b Example: @code grid.export('mymesh') @endcode
%>
%> @param   this    Grid        [Grid]
%> @param   name    description [string]

function export(this, name)

assert(ischar(name))
printline(1, 'Exporting grid to M++ format.')

% open and pipe to .geo file
file = fopen(['./', name, '.geo'], 'wt');

fprintf(file, 'POINTS: %d\n', this.numV);
for iV = 1:this.numV
    fprintf(file, '%f %f\n', this.coordV(iV, 1), this.coordV(iV, 2));
end

%% CELLS (point numbers are from 0 to size(this.nct, 1) - 1)
fprintf(file, 'CELLS: %d\n', this.numT);
for iT = 1:size(this.V0T, 1)
    fprintf(file, '3 %d %d %d %d\n', 1000, this.V0T(iT, 1)-1, this.V0T(iT, 2)-1, this.V0T(iT, 3)-1);
end

%% FACES/EDGES (point numbers are from 0 to size(this.numV, 1) - 1)
fprintf(file, 'FACES: %d\n', this.numE);
for kE = 1:this.numE
    fprintf(file, '2 %d %d %d\n', this.idE(kE), this.V0E(kE, 1)-1, this.V0E(kE, 2)-1);
end
fclose(file);

printline(3, ['Grid exported to ', name, '.geo.'])

end % export
