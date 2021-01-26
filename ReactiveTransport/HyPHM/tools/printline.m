%> @file printline.m This is the HyPHM print function to standard out.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> Identification Numbers (example):
%>
%>   -  1:  |Visualization... (brown)
%>   -  2:  |  Using the following paramters: (green)
%>   -  3:  |    approximation  (black)
%>   -  0:  |[Silent mode, do nothing]
%>   - -1:  |Warning or Request: (orange)

function printline(id, s, varargin)

if id == 0
    return
end

persistent lastused % last used type of line, required for newlines

msg = 'HyPHM: Usage is printline(formatid, formatstring, formatvars)';
assert(nargin >= 2, msg)
assert(isnumeric(id), msg)
assert(ischar(s), msg)

% color definition
prettycolors = {[0.5, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.2, 0.0]};


switch id

    %% HEADING
    case 1
        savecprintf(prettycolors{1}, ['\n', s, '\n'], varargin{:});
        lastused = 1;

        %% SUBHEADING
    case 2
        if ~isempty(lastused) && (lastused == 2 || lastused == 3)
            fprintf('\n');
        end
        savecprintf(prettycolors{2}, ['  ', s, '\n'], varargin{:});
        lastused = 2;

        %% DATA/INFORMATION
    case 3
        if lastused == 1
            fprintf('\n');
        end
        savecprintf(prettycolors{3}, ['    ', s, '\n'], varargin{:});
        lastused = 3;

        %% WARNING / REQUEST INFORMATION etc
    case -1
        savecprintf([1.0, 0.4, 0.2], ['\n', s, '\n'], varargin{:});
        lastused = 1;
    otherwise
        error('HyPHM: id for printline not implemented.')
end

end


function savecprintf(varargin) % fallback to fprintf
try
    cprintf(varargin{:})
catch
    fprintf(varargin{2:end})
end
end
