%> @file +CCCC/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$.
%> This file does only thow a commented error, to prevent the user to make the mistake
%> using a space-dependant function handle here.

function ret = setfh2data(~, ~, ~)

error('HyPHM: Please specify the constant tensor directly by a constant 2x2-matrix.')

end
