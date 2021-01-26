%> @file AbstractProblem.m All problem classes should inherit this class.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> All problem classes should inherit this class.  Listed <code>properties</code>
%> and <code>methods</code> of AbstractProblem with the attribute <code>Abstract</code>
%> <i>must</i> be implemented by inheriting classes with <i>the very same attributes</i>
%> but <i>without the attribute <code>Abstract</code></i>.  Inherited <code>methods</code>
%> (i.e. functon definitions within the block '<code>methods ... end</code>'
%> may have different arguments.
%>
%> @b Example
%> Consider the following snippet from the @c AbstractProblem class:
%> @code
%> classdef AbstractProblem < hgsetget & handle
%>   properties (Abstract, SetAccess = private)
%>     grid
%>   end
%>   methods (Abstract, Hidden)
%>     solve(dummy)
%>   end
%> end
%> @endcode
%> An inheriting class MyProblem may then be implemented as
%> @code
%> classdef MyProblem < AbstractProblem
%>   properties (SetAccess = private)
%>     grid
%>   end
%>   methods (Hidden)
%>     function out1 = solve(arg1, arg2)
%>     [...]
%>     end
%>   end
%> end
%> @endcode
classdef AbstractProblem < hgsetget & handle & matlab.mixin.Copyable

    properties (Abstract, SetAccess = private)
        %> <i>(Abstract, SetAccess = private)</i> Instance of Grid or FoldedGrid (should be set in construction).
        grid
        %> <i>(Abstract, SetAccess = private)</i> Instance of Stepper (should be set in construction).
        stepper
        %> <i>(Abstract, SetAccess = private)</i> Name of the problem [string] (should be set in construction).
        name
    end

    methods (Abstract)
        %> <i>(Abstract)</i> Iterate this.stepper and compute the unknown(s) of the new time level.
        computeLevel(dummy)
        %> <i>(Abstract)</i> Print properties in a nice way to the command window.
        print(dummy)
    end

    methods (Abstract, Hidden)
        %> <i>(Abstract, Hidden)</i> Assemble the system matrix and the right-hand side.
        assembleSystem(dummy)
        %> <i>(Abstract, Hidden)</i> Solve the linear system.
        solveSystem(dummy)
    end

end
