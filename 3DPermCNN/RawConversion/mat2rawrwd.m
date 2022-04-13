% Convert 3D Matlab matrix to raw file
% Parameters: A: 		matrix
%	      filename: 	filename (string)
%	      bits:		number of bits in raw format
%	      Mincoord: 	Triple of real numbers, minimal coordinate
%             Maxcoord: 	Triple of real numbers, maximal coordinate
%	      relativePath:	Add relative path to file location	

function mat2rawrwd(A, filename, bits, Mincoord, Maxcoord, relativePath)

[filepath, name, ext] = fileparts(filename); 
A = round(A); % Convert to binary.

%%%%%%%%%%%%%%%%%%%%%%%
%% Modify file name. %%
%%%%%%%%%%%%%%%%%%%%%%%
spaceDim = length(size(A));
switch spaceDim
  case 3
    [Nx, Ny, Nz] = size(A);
    name = sprintf('%s_%dx%dx%d', name, Nx, Ny, Nz);
  otherwise
    error('Unexpected data format.')
end

%%%%%%%%%%%%%%%%%%%%%%%
%% Create .raw file. %%
%%%%%%%%%%%%%%%%%%%%%%%
rawfilename = [name, '.raw'];
fid = fopen(rawfilename, 'w+');
switch bits
  case 1
    fwrite(fid, A, 'ubit1');
  otherwise
    warning('Input parameter ''bits'' not specified -- writing 1 bit.')
    fwrite(fid, A, 'ubit1');
end
fclose(fid);
disp(['Data written to ''', rawfilename, '''.'])

%%%%%%%%%%%%%%%%%%%%%%%
%% Create .rwd file. %%
%%%%%%%%%%%%%%%%%%%%%%%
assert(spaceDim==numel(Mincoord) & spaceDim==numel(Maxcoord),'Min/Max coordinate do not fit dimension');

assert(isequal(min((Maxcoord-Mincoord)./size(A)), max((Maxcoord-Mincoord)./size(A))), 'Voxels are not regular');


% Open .rwd file.
rwdfilename = [name, '.rwd'];
fid = fopen(rwdfilename, 'wt');
if fid == -1
  error('Cannot open file: %s.', rwdfilename);
end

% Writing the data.
fprintf(fid, 'RAW-FILE DESCRIPTOR VERSION 1\n'); % To know what content to expect.
fprintf(fid, 'FILE %s\n', [relativePath,rawfilename]);      % The target raw file.
fprintf(fid, 'BITS %i\n', bits);             % Unsigned integer, uint8, uint16, etc.
fprintf(fid, 'DIMENSION %i\n', spaceDim);    % Space dimension (3).

switch spaceDim
    case 3
       fprintf(fid, 'SIZE %i %i %i\n', Nx, Ny, Nz); % Voxels in each space direction.
       fprintf(fid, 'MINCOORDINATE %f %f %f\n',Mincoord(1), Mincoord(2), Mincoord(3)); % Coordinates of point with smallest indices.
       fprintf(fid, 'MAXCOORDINATE %f %f %f\n',Maxcoord(1), Maxcoord(2), Maxcoord(3)); % Coordinates of point with largest indices.
end

fclose(fid);

% Output information.
disp(['Data written to ''', rwdfilename, '''.'])


end
