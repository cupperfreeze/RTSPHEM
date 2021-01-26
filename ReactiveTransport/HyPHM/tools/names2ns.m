%> @file names2ns.m Converts a list of strings to a so-called `name space matrix' required by the pdetool function @c decsg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> <h2>Example</h2>
%> Definition of a square @c SQ1 with a round perforation @c C1 (circle).
%> @code
%>     gd = [ 3.0000    1.0000             % geometry description
%>            4.0000         0
%>           -1.0000         0
%>            1.0000    0.4000
%>            1.0000         0
%>           -1.0000         0
%>           -1.0000         0
%>           -1.0000         0
%>            1.0000         0
%>            1.0000         0 ];
%>     sf = 'SQ1-C1';                      % set formula
%>     ns = names2ns('SQ1', 'C1');         % name space
%>     [p,e,t] = initmesh(decsg(gd,sf,ns), 'Hmax', 0.2);
%> @endcode
%>
%> @sa decsg, pdetool

function ns = names2ns(varargin)
if iscell(varargin{1}) % if 1st arg is already a cell
    varargin = varargin{1};
end
num = length(varargin); % number of words
ns = zeros(1, num); % first dim unknown yet
for k = 1:num
    assert(ischar(varargin{k}), 'list{%d} does not contain a string.', k)
    for ell = 1:length(varargin{k})
        ns(ell, k) = double(varargin{k}(ell));
    end
end

end
