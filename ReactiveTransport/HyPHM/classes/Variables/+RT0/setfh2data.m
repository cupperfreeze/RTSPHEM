%> @file +RT0/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$, @f$f:\mathbb{R}\times\mathbb{R}^2\rightarrow\mathbb{R}^2@f$.

%> @brief Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$, @f$f:\mathbb{R}\times\mathbb{R}^2\rightarrow\mathbb{R}^2@f$.
%>
%> The degrees of freedom are on the barycenters of edges.
function ret = setfh2data(time, g, fun2)

ret = P2P2.P2P2toRT0slice(g, P2P2.setfh2data(time, g, fun2));

% ret = zeros(g.numE, 1);
% if isa(g, 'FoldedGrid') % this works with Grid also, but is slower
%   for kT = 1 : g.numT
%     for k = 1 : 3
%       ret(g.E0T(kT, k)) = dot(fun2(time, g.baryE0T(kT,k, :)'), g.nuE(g.E0T(kT, k),:)) ;
%     end
%   end
%
% elseif isa(g, 'Grid')
%   for kE = 1 : g.numE
%     ret(kE) = dot(fun2(time, g.baryE(kE,:)'), g.nuE(kE,:)) ;
%   end
% end

end
