% Convert raw file 'filename' to Matlab 3D-matrix
% Data required to be of dimension 'siz' (triple of integers) and in format of 'bits' bit 

function A = raw2mat(filename, siz, bits)

assert(numel(siz) == 3)

% Add '.raw' if filename does not contain the extension.
[~,~,ext] = fileparts(filename);
if ~strcmp(ext, '.raw')
  filename = [filename, '.raw'];
end

%% Read raw file and store voxel data in A.
% A has a value unequal to zero at each index (i,j,k) that is
% associated with fluid.
fid = fopen(filename, 'r');
if fid == -1
  error('Cannot open file: %s.', filename);
end
switch bits
  case 1
    A = fread(fid, siz(1)*siz(2)*siz(3), 'ubit1=>uint8');
  case 8
    A = fread(fid, siz(1)*siz(2)*siz(3), 'uint8=>uint8');
  otherwise
    warning('Input parameter ''bits'' not specified -- reading 8 bits.')
    A = fread(fid, siz(1)*siz(2)*siz(3), 'uint8=>uint8');
end

fclose(fid);
A = reshape(A, [siz(1), siz(2), siz(3)]);

end % function
