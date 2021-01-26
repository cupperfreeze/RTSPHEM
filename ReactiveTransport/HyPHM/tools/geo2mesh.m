%> @file geo2mesh.m Calls the mesh generator gmsh to generate a .mesh file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> The output is a <tt>.mesh</tt> file which will be written to the same directory as the <tt>.geo</tt> file.
%> @sa loadmesh.m.
%>
%> @param filename The name of the <tt>.geo</tt> file including file extension.
%> @param scaling Scaling factor which is multiplied to the mesh size @f$h@f$.


function geo2mesh(filename, scaling)

printline(1, 'Calling mesh generator gmsh for %s', filename)


% fetching pathname, filename and extension
filename = which(filename);
[pathstr, name, ext] = fileparts(filename);
assert(strcmp(ext, '.geo'), 'HyPHM: geo2mesh expects a .geo file, however, extension was %s.', ext)

% removing old mesh file, if any
printline(3, 'removing old mesh file %s, if exists', [pathstr, '/', name, '.mesh'])
delete([pathstr, '/', name, '.mesh'])

% building system call command
try
    systemcall = sprintf('gmsh -2 -format mesh -clscale %f -o "%s" "%s"', scaling, [pathstr, '/', name, '.mesh'], filename);
end

% information to workspace
printline(2, 'The following command is executed as system call')
printline(3, systemcall)

% this call requires gmsh to be in the path variable (in particular under windows)
unix(systemcall);

printline(3, 'If the above stuff is working the mesh was written to %s', [pathstr, '/', name, '.mesh'])


end
