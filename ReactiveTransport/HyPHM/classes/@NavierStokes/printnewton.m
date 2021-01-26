function printnewton(obj)

printline(2, 'Data for Newton iteration (convection term)');
printline(3, 'Residual bound: %e', obj.maxRes);
printline(3, 'Maximum number of iterations: %d', obj.maxIter);

end
