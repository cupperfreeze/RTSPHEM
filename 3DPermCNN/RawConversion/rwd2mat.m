% Read rwd file 'filename' to convert associate raw file to Matlab 3D-matrix

function [A,varargout] = rwd2mat(filename)

% Add '.rwd' if filename does not contain the extension.
[~,~,ext] = fileparts(filename);
if ~strcmp(ext, '.rwd')
  filename = [filename, '.rwd'];
end

% Initialization.
BITS = [];
SIZE = [];
DIMENSION = [];
MINCOORDINATE = [];
MAXCOORDINATE = [];

% Read content of raw descriptor file (two loops, since the second one depends on the first one).
fid = fopen(filename);                                                     % Open text file.
tline = fgetl(fid);                                                        % Get first line of text file.
while ischar(tline)                                                        % While text line is not empty...
  if isVariableInLine(tline, 'DIMENSION')
    DIMENSION = sscanf(tline, 'DIMENSION %i');
  end
  tline = fgetl(fid);                                                      % Read next line.
end
fclose(fid);                                                               % Close text file.

fid = fopen(filename);                                                     % Open text file.
tline = fgetl(fid);                                                        % Get first line of text file.
while ischar(tline)                                                        % While text line is not empty...
  if isVariableInLine(tline, 'BITS')
    BITS = sscanf(tline, 'BITS %i');
  end
  if isVariableInLine(tline, 'MINCOORDINATE')
    MINCOORDINATE = sscanf(tline, 'MINCOORDINATE %f %f %f');
  end
  if isVariableInLine(tline, 'MAXCOORDINATE')
    MAXCOORDINATE = sscanf(tline, 'MAXCOORDINATE %f %f %f');
  end
  if isVariableInLine(tline, 'SIZE')
    switch DIMENSION
      case 3
        SIZE = sscanf(tline, 'SIZE %i %i %i');
      otherwise
        error('DIMENSION must be 3.');
    end % switch
  end
  tline = fgetl(fid);                                                      % Read next line.
end
fclose(fid);                                                               % Close text file.

% Load data from the actual raw file.
filename(end-3:end) = '.raw';
A = raw2mat(filename, SIZE, BITS);
varargout{1} = SIZE';

end

% Returns true if the variable name exists in the line (at the beginning).
function ret = isVariableInLine(textline, varname)
assert(ischar(textline))
assert(ischar(varname))
n = length(varname);
ret = length(textline) > n && isequal(textline(1:n), varname);
end % function


