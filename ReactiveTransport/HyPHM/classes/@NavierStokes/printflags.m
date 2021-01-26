function printflags(obj)

printline(2, 'Problem specification');

if obj.isEvolution
    msgevo = 'transient';
else
    msgevo = 'stationary';
end

if obj.isConvection
    msgconv = 'Navier-Stokes';
else
    msgconv = 'Stokes';
end

if obj.isSymmetric
    msgsym = 'symmetric';
else
    msgsym = 'non-symmetric';
end

printline(3, 'Problem is the %s %s equation with %s diffusion tensor.', msgevo, msgconv, msgsym)

printline(3, 'Mean pressure contraint: %s', x2str(obj.balanceP));


end
