%> @file Stepper.m Time-management modul that returns current time discretization data like step size or step number.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Time-management modul that returns current time discretization data like step size or step number.
%>
%> @code st = Stepper(timepts) @endcode
%> creates an instance st, where the vector timepts contains the decomposition
%> of the time interval.  Stepper has the following properties:
%>
%>   - @c curstep  current step number (step @f$k@f$ is from @f$t_{k-1}@f$ to @f$t_k@f$)@n
%>   - @c curtau   current step size (@f$\tau_k = t_k - t_k-1@f$)@n
%>   - @c curtime  current time (this is the time where the solution has to be computed)@n
%>   - @c initime  first time point@n
%>   - @c endtime  last time point@n
%>   - @c numsteps total number of steps (equal to time points minus one)
%>
%> The stepper can be iterated with the method @c next, which returns @c false
%> if the current step is the final one/if there is no next step,
%> otherwise true.  You can also go back one time step via method @c prev,
%> which analogously returns true if the start iterate has not been
%> reached yet.
%>

%% Example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @section steppersection1 Example
%> @code
%>   st = Stepper([1.0, 1.2, 1.5, 1.6]);
%>   st.print
%>   fprintf('curstep: %d, curtime: %.1f, curtau: %.1f\n', st.curstep, st.curtime, st.curtau);
%>   while st.next
%>     fprintf('curstep: %d, curtime: %.1f, curtau: %.1f\n', st.curstep, st.curtime, st.curtau);
%>   end
%> @endcode
%>

classdef Stepper < handle

    %% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        %> <i>(read only)</i> current step number (step @f$k@f$ is from @f$t_{k-1}@f$ to @f$t_k@f$)@n
        curstep
        %> <i>(read only)</i> current step size (@f$\tau_k = t_k - t_k-1@f$)@n
        curtau
        %> <i>(read only)</i> current time (this is the time where the solution has to be computed)
        curtime
        %> <i>(read only)</i> time of first time point
        initime
        %> <i>(read only)</i> time of last time point
        endtime
        %> <i>(read only)</i> total number of steps (equal to time points minus one)
        numsteps
    end

    properties (Access = public, Hidden)
        %> time points (index shift: first entry referres to initial, (numsteps+1)th referres to final)
        timepts
    end

    properties (Access = private, Hidden)
        %> <i>(private)</i> vector of steps (length = length(timepts) minus one)
        stepvec
        %> the waitbar handle (optional, void if no waitbar requested)
        wbarh
    end

    methods

        %% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> Constructor.
        %> @param timepts Vector of time points, e.g., <code>0:.1:10</code>.
        %> @retval Instance of Stepper.
        function this = Stepper(timepts)
            this.timepts = unique(sort(timepts(:))); % time points/time grid
            this.numsteps = length(this.timepts) - 1; % total number of time steps

            if this.numsteps == 0
                error('HyPHM: Stepper has to be defined with a least one time step.')
            end

            this.stepvec = zeros(this.numsteps, 1);
            for k = 1:this.numsteps
                this.stepvec(k) = this.timepts(k+1) - this.timepts(k);
            end

            this.initime = this.timeofstep(0);
            this.endtime = this.timeofstep(this.numsteps);

            this.curtau = 0.0;
            this.curtime = this.initime;
            this.curstep = 0;

            this.wbarh = []; % no waitbar handle until it is requested
        end

    end

    methods

        %% ITERATE TIMER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function bool = next(this)
            %> Iterate stepper, i.e., go to the next step.
            if this.curstep < this.numsteps % if end time is has not reached yet
                this.curstep = this.curstep + 1;
                this.curtime = this.timepts(this.curstep+1);
                this.curtau = this.stepvec(this.curstep);
                bool = true;
            else
                bool = false; % return 0 if we reach end time
            end
            this.updatewbar;
        end

        %% Update TimeStepSize
        % Overwrite length of next timestep
        function [] = setTimeStepSize(this, num, val)
            this.timepts(num+1) = this.timepts(num) + val;
            this.curtau = val;
            this.curtime = this.timepts(num+1);
        end

        %% DE-ITERATE TIMER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> De-iterate stepper, i.e., go to the previous step.
        function bool = prev(this)
            if this.curstep > 0 % if end start is has not reached yet
                this.curstep = this.curstep - 1;
                this.curtime = this.timepts(this.curstep+1);
                if this.curstep == 0
                    this.curtau = 0.0;
                else
                    this.curtau = this.stepvec(this.curstep);
                end
                bool = true;
            else
                bool = false; % return 0 if we reach start time
            end
            this.updatewbar;
        end

        %% REVERT TIMER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> Revert stepper, i.e., go to the initial step.
        function revert(this)
            this.curtau = 0.0;
            this.curtime = this.timeofstep(0);
            this.curstep = 0;
        end

        %% TIME POINT TO A STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> Returns the time value of the requested step
        %> @param this Instance of Stepper.
        %> @param stepno Number of requested step.
        function ret = timeofstep(this, stepno)
            assert(stepno >= 0 && stepno <= this.numsteps, ...
                'HyPHM: Step number out of range [0, %d]', this.numsteps)
            ret = this.timepts(stepno+1);
        end

        %% PRINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> Print the Stepper information to terminal.
        function print(this)
            printline(2, 'The stepper contains the following data')
            printline(3, 'Number of steps: %d', this.numsteps)
            s = sprintf('Time grid:       ');
            for k = 1:min(5, this.numsteps) % first 5 entries
                s = [s, sprintf('%.1e  ', this.timepts(k))]; %#ok<AGROW>
            end
            if this.numsteps > 5
                s = [s, ' . . .   ']; % dots if entries omitted
            end
            s = [s, sprintf('%.1e', this.timepts(end))]; % last entry
            printline(3, s);
        end

        %% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> Plot the time points vs time in a diagram.
        function visualize(this)
            plot(this.timepts)
        end

        %% STATUS BAR -- GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> Open a status bar displaying the current position of the Stepper.
        function status(this)
            if isempty(this.wbarh)
                this.wbarh = 'dummy';
                progressbar('Progress...') % initialize bar (globally defined, does not need to be stored/pointed to)
            end
            this.updatewbar;
        end
    end


    methods (Hidden = true, Access = private)

        %% STATUR BAR -- UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> Update the status bar (does not happen automatically).
        function updatewbar(this)
            if ~isempty(this.wbarh)
                % msg = sprintf('Step %d of %d', this.curstep, this.numsteps);
                progressbar(this.curstep/this.numsteps);
            end
        end
    end


end
