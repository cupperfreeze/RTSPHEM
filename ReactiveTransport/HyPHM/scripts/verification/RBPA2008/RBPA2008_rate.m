%> @file RBPA2008_rate.m Rate acc to @ref RBPA2008 to be mounted as <tt>function_handle</tt>.
%> @f$2A+D \rightarrow \emptyset @f$, with a rate constant equal to one.

function dc = RBPA2008_rate(c)

numSpec = length(c);
assert(numSpec == 2)
dc = zeros(numSpec, 1); % a column vector

dc(1) = -2 * c(1)^2 * c(2);
dc(2) = -1 * c(1)^2 * c(2);


end
