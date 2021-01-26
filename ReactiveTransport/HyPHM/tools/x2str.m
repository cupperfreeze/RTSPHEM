%> @file x2str.m
%> @brief Converts the input argument into a string.
%>
%> @b Example: @n
%> The call
%> @code
%>   fun = @(x) x.^2 + 4;
%>   x2str(fun)
%> @endcode
%> returns <tt>'@(x)x.^2+4'</tt> as string.

%> param in Arbitrary input argument to be converted to a string.
function ret = x2str(in)
if isempty(in)
    ret = '[]';
elseif isa(in, 'function_handle')
    ret = func2str(in);
elseif isnumeric(in);
    ret = num2str(in);
elseif ischar(in)
    ret = in;
elseif isa(in, 'Unknown')
    ret = '[Unknown]';
elseif isa(in, 'Variable')
    ret = sprintf('Variable of type %s', in.type);
else
    error('HyPHM: Unknown type.');
end
end
