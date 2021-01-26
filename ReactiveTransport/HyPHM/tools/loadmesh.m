%> @file loadmesh.m Load and convert a mesh generated via gmsh in the mesh format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param filename File in the .mesh format.
%> @retval coordV  cf Grid.coordV.
%> @retval V0T     cf Grid.V0T.
%> @retval V0Ebdry Vertices of (unoriented) edges only on @f$\partial\Omega@f$, cf Grid.V0E.
%> @retval idEbdry Referring edge IDs to V0Ebdry, cf Grid.idE.
function [coordV, V0T, V0Ebdry, idEbdry] = loadmesh(filename)
msg = 'HyPHM: loadmesh(filename) where filename refers to a .mesh file.';
assert(ischar(filename), msg), assert(nargin == 1, msg)
[~, ~, ext] = fileparts(filename);
assert(strcmp(ext, '.mesh'), [msg, '  Extension string mismatch, .mesh asserted.'])

printline(1, 'Loading mesh file %s', filename)

try
    fid = fopen(filename, 'r');
    tline = fgets(fid);
    while ischar(tline)
        % VERTICES
        if strfind(tline, 'Vertices')
            numV = fscanf(fid, '%d', [1, 1]);
            coordV = reshape(fscanf(fid, '%f'), 4, numV)';
            coordV(:, 3:4) = [];
        end
        if strfind(tline, 'Edges')
            numE_bdry = fscanf(fid, '%d', [1, 1]);
            tmp = reshape(fscanf(fid, '%f'), 3, numE_bdry)';
            V0Ebdry = tmp(:, 1:2);
            idEbdry = tmp(:, 3);
            clear tmp;
        end
        if strfind(tline, 'Triangles')
            numT = fscanf(fid, '%d', [1, 1]);
            V0T = reshape(fscanf(fid, '%f'), 4, numT)';
            V0T(:, 4) = [];
        end
        tline = fgets(fid);
    end

    fclose(fid);
catch
    error('HyPHM: Something went wrong, the .mesh file is damaged.  Please check the Gmsh log message.')
end
printline(3, '...done.')

printline(2, 'Checking for unconnected points')
for kV = 1:size(coordV, 1)
    if ~any(any(V0T == kV))
        msg = sprintf('HyPHM: Point %d (and maybe others) is not a vertex from a triangle.  This DOES make (Navier-)Stokes useless!', kV);
        printline(3, msg)
        if strcmp(input('continue nevertheless? [y,N] ', 's'), 'y')
            break
        else
            error('see above')
        end
    end
end
printline(3, '...done.')

end
