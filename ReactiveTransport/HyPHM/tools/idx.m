%> @file idx.m Extracts the considered indices of the argument, even it is a vector-valued answer of a function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @b Example:
%> @code
%> fun = @(x) [sin(x(2)), sin(x(1)); 3, cos(x(2)); cos(x(1)), pi];
%> display(fun)
%> display(fun([3;2]))
%> % fun([3;2])(1,:) won't work, so simply do
%> idx(fun([3;2]), 1, ':')
%> @endcode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param funout Matrix-/ or vector-valued data or function answer.
%> @param varargin One or two indices, see Example.
%> @retval ret funout(k,l), funout(:,l) or funout(k,:).

function ret = idx(funout, varargin)

% NO INDEXING
if nargin == 1
    ret = funout; % funout
    % SINGLE INDEX
elseif nargin == 2
    arg1 = varargin{1};
    fun = @(x, k) x(k);
    if isequal(arg1, ':')
        ret = fun(funout, 1:length(funout(:))); % funout(:)
    else % eg arg1 = 3
        ret = fun(funout, arg1); % funout(k)
    end
    % DOUBLE INDEX
elseif nargin == 3
    arg1 = varargin{1};
    arg2 = varargin{2};
    fun = @(x, k, ell) x(k, ell);
    if isequal(arg1, ':') && isequal(arg2, ':')
        ret = funout; % funout(:, :) == funout
    elseif isequal(arg1, ':') && ~isequal(arg2, ':')
        ret = fun(funout, 1:size(funout, 1), arg2); % funout(:, ell)
    elseif ~isequal(arg1, ':') && isequal(arg2, ':')
        ret = fun(funout, arg1, 1:size(funout, 2)); % funout(k, :)
    else
        ret = fun(funout, arg1, arg2); % funout(k, ell)
    end
else
    error('HyPHM: RTFM.')
end


end
