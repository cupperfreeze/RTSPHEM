function [SOL] = solveSystemFE(grid, A, rhs, isDoF)
%Solve system of equations given by system matrix A and rhs for the
%unknowns given by isDoF

RHS = rhs(isDoF, :);
%     if( isa( grid, 'FoldedCartesianGrid' ) )
%         RHS( end + 1, : ) = 0;
%     end
if (isa(grid, 'FoldedCartesianGrid'))
    ANEW = A(1:end-1, :);
else ANEW = A;
end
%Remark: last line contained normalization condition --> solution now
%defined uniquely up to constant. PCG will chosen one that is normalized
%later

% sol = A\[RHS; 0 0];
r = symrcm(ANEW);
r = 1:size(ANEW, 1);
ANEW = ANEW(r, r);
RHS = RHS(r, :);
opts.type = 'ict';
opts.droptol = 0.0005;
L = ichol(ANEW, opts);
temp(r) = 1:size(ANEW, 1);

[sol2, flag2] = pcg(ANEW, RHS(:, 1), 10^(-10), 1000, L, L');
sol2 = sol2 - sum(sol2) / size(sol2, 1);
sol2 = sol2(temp);
[sol3, flag3] = pcg(ANEW, RHS(:, 2), 10^(-10), 1000, L, L');
sol3 = sol3 - sum(sol3) / size(sol3, 1);
sol3 = sol3(temp);


assert(flag2+flag3 == 0, 'Diffusion solver failed');
%    norm(sol2-sol(:,1))/ norm(sol(:,1)) + norm(sol3-sol(:,2))/ norm(sol(:,2))
sol = [sol2, sol3];
SOL = NaN(grid.nodes, 2);
SOL(isDoF, :) = sol;

if (isa(grid, 'FoldedCartesianGrid'))
    SOL(:, 1) = grid.synchronizeValues(SOL(:, 1));
    SOL(:, 2) = grid.synchronizeValues(SOL(:, 2));
end

end
