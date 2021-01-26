% VTKTRISURF % VTKQUIVER creates a vtk-file containing a plot similar to
%    Matlab's trisurf which can be called by a vtk-visualization tool.
%
%    VTKTRISURF(TRI,X,Y,Z,VALS,VARNAME,FILENAME)
%    See Matlab documentation acc. to TRISURF(TRI,X,Y,Z,C).
%
%    VTKTRISURF(TRI,X,Y,Z,VARNAME,FILENAME)
%    See Matlab documentation acc. to TRISURF(TRI,X,Y,Z).
%
%    The (optional) argument VALS should contain the scalar data, which can
%    be given on the vertices - OR - on the cells/triangles of the triangu-
%    lation - OR - on the vertices of each triangle.
%    The type is determined automatically.  If VALS is not
%    specified, Z is used as point data.  The argument VARNAME should
%    contain the description of the visualized variable which is required
%    by the vtk file format. FILENAME is the name of the .vtu-file to store
%    the data.  If the string does not end with '.vtu', this extension is
%    attached.
%
% Example 1
% 1   [X,Y] = meshgrid(-2:0.25:2,-1:0.2:1);
% 2   Z = X.* exp(-X.^2 - Y.^2);
% 3   TRI = delaunay(X,Y);
% 4   vtktrisurf(TRI,X,Y,Z,'pressure','vtktrisurf.vtu')
% Now open `vtktrisurf.vtu' via Paraview or Mayavi2.
%
% Example 2
% If you do not want to use the Z-data as Z-coordinate you can write
% 1   [Lines 1 to 3 from Example 1]
% 2   vtktrisurf(TRI,X,Y,zeros(size(Z)),Z,'pressure','vtktrisurf.vtu')
% If you use 'Warp by Scalar' in Paraview the data is lifted in Z-direction again.
%
% For a minimal example you can generate a grid with only 2 triangles using
% '[X,Y] = meshgrid([0,1]);' instead.
%
%
% See also trisurf, vtkquiver
%
% Copyright see license.txt
%
% Author:        Florian Frank
% eMail:         frank@math.fau.de
% Version/Date:  20140409

function vtktrisurf(tri, x, y, z, argin5, argin6, argin7)

% INITIALIZATION
if nargin == 6
    vals = z;
    varname = argin5;
    filename = argin6;
elseif nargin == 7
    vals = argin5;
    varname = argin6;
    filename = argin7;
else
    error('Wrong number of input arguments.')
end

% ASSERTIONS
assert(ischar(varname) && ischar(filename))
numC = size(tri, 1); % number of cells
numP = length(x(:)); % number of points
assert(numP == length(y(:)) && numP == length(z(:)))

% DETERMINE DATA TYPE
if length(vals(:)) == numP
    datatype = 'pointdata';
elseif length(vals(:)) == numC
    datatype = 'celldata';
elseif isequal(size(vals), [numC, 3])
    datatype = 'discpointdata'; % discontinuous point data
else
    error('Input argument VALS has wrong dimensions.')
end

% PREPROCESSING
nodelist = reshape(tri', numC*3, 1);

% AVOID NAN'S BY PRINTING SOME RANDOM FIELD (to reveal the error)
if any(isnan(vals(:)))
    vals = -ones(size(vals));
end

% OPEN FILE
if ~strcmp(filename(end -3:end), '.vtu') % append file extension if not specified yet
    filename = [filename, '.vtu'];
end

file = fopen(filename, 'wt'); % if this file exists, it is overwritten

% HEADER
fprintf(file, '<?xml version="1.0"?>\n');
fprintf(file, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n');
fprintf(file, '  <UnstructuredGrid>\n');
fprintf(file, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n', numC*3, numC);

% POINTS
fprintf(file, '      <Points>\n');
fprintf(file, '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
fprintf(file, '          %.3e %.3e %.3e\n', [x(nodelist), y(nodelist), z(nodelist)]');
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </Points>\n');

% CELLS
fprintf(file, '      <Cells>\n');
fprintf(file, '        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
fprintf(file, '           ');
fprintf(file, '%d ', 0:numC*3-1);
fprintf(file, '\n        </DataArray>\n');
fprintf(file, '        <DataArray type="Int32" Name="offsets" format="ascii">\n');
fprintf(file, '           %d\n', 3:3:3*numC);
fprintf(file, '        </DataArray>\n');
fprintf(file, '        <DataArray type="UInt8" Name="types" format="ascii">\n');
fprintf(file, '           %d\n', 5*ones(numC, 1));
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </Cells>\n');

% CELL DATA
if strcmp(datatype, 'celldata')
    fprintf(file, '      <CellData Scalars="%s">\n', varname); % def of std value
    fprintf(file, '        <DataArray type="Float32" Name="%s" format="ascii">\n', varname);
    fprintf(file, '           %.3e\n', vals);
    fprintf(file, '        </DataArray>\n');
    fprintf(file, '      </CellData>\n');
end

% POINT DATA
if strcmp(datatype, 'pointdata')
    fprintf(file, '      <PointData Scalars="%s">\n', varname); % def of std value
    fprintf(file, '        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n', varname);
    fprintf(file, '          %.3e\n', vals(nodelist));
    fprintf(file, '        </DataArray>\n');
    fprintf(file, '      </PointData>\n');
end
if strcmp(datatype, 'discpointdata')
    fprintf(file, '      <PointData Scalars="%s">\n', varname); % def of std value
    fprintf(file, '        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n', varname);
    fprintf(file, '          %.3e\n', vals');
    fprintf(file, '        </DataArray>\n');
    fprintf(file, '      </PointData>\n');
end

% FOOTER
fprintf(file, '    </Piece>\n');
fprintf(file, '  </UnstructuredGrid>\n');
fprintf(file, '</VTKFile>\n');

fclose(file);

return

end
