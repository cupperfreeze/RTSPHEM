% exact integration
% returns int_T lam1^a*lam2^b*lam3^c / 2 / |T|
function ret = iex(a, b, c)
str = sprintf('%d!*%d!*%d!/(%d+%d+%d+2)!', a, b, c, a, b, c);
ret = sym(str);
end