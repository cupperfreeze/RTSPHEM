%> @file isVariable.m Check for class Variable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @argin arg anything.
%> @retval ret {true, false}

function ret = isVariable(arg)

ret = isa(arg, 'Variable');

end
