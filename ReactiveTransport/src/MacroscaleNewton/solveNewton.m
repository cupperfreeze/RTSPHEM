function sol = solveNewton(M, b, iteration, numTransportVariables, varargin)

%Solve linear systems arising in Newton iterations
global Solver
persistent PREC PREC1 PREC2
persistent options pl pr Dl Dr
persistent levs

if numel(Solver) == 0 | ~ismember(Solver, {'StandardDirect', 'ilupack', 'BlockPrec', 'StandardIterative', 'ilukIterative'})
    Solver = 'StandardDirect'; %default choice
    disp('Use standard direct solver');
end

if  ismember(Solver, { 'ilupack', 'BlockPrec'})
    try
       AMGinit(sparse(1));
    catch
        warning('Solver option requested requires ilupack, switching to Standard solver');
        Solver = 'StandardDirect';
    end
end

if  ismember(Solver, {'BlockPrec','ilukIterative'})
    try
       [~,~]=iluk(sparse(1),1);
    catch
        warning('Solver option requested requires iluk, switching to Standard solver');
        Solver = 'StandardDirect';
    end
end
% Choose different kind of solvers

switch Solver
    % Matlab direct solver
    case 'StandardDirect'
        sol = M \ b;
        return

        % ilupack package
    case 'ilupack'
        if length(varargin) > 0
            PREC = AMGdelete(PREC);
            options = 0;
            sol = 0;
            return
        end

        % biuld up preconditioner
        if iteration == 0
            if numel(PREC) > 1
                PREC = AMGdelete(PREC);
            end
            options = AMGinit(M);
            options.mixedprecision = 1;
            options.isdefinite = 0;
            options.droptol = 10^(-1.7);
            %options.amg='amli';
            %options.ncoarse=2;
            %options.solver='bicgstab';
            m = size(M, 1) / numTransportVariables;
            cells = cell(numTransportVariables, 1);
            for i = 1:numTransportVariables
                cells{i} = M((i-1)*m+1:i*m, (i - 1)*m+1:i*m);
            end
            temp = blkdiag(cells{:});
            [PREC, options] = AMGfactor(temp, options);


        end

        % solve system
        options.nrestart = 100;
        options.maxit = 3000;
        options.restol = 1e-3;
        [sol, options] = AMGsolver(M, PREC, options, b);
        options.niter;

        clear M
        clear b


        %Use block structure
    case 'BlockPrec'
        m = size(M, 1) / numTransportVariables;

        if length(varargin) > 0
            sol = 0;
            return
        end


        if iteration == 0
            PREC1 = cell(numTransportVariables, 1);
            PREC2 = cell(numTransportVariables, 1);
            pl = cell(1, numTransportVariables);
            pr = cell(1, numTransportVariables);
            Dl = cell(numTransportVariables, 1);
            Dr = cell(numTransportVariables, 1);
            tic

            for i = 1:numTransportVariables
                %resort and rescale each submatrix
                temp = M((i-1)*m+1:i*m, (i - 1)*m+1:i*m);
                [pl{i}, pr{i}, Dl{i}, Dr{i}] = mwmamd(temp);
                temp2 = Dl{i} * temp(pl{i}, pr{i}) * Dr{i};
                pl{i} = pl{i} + (i - 1) * m;
                pr{i} = pr{i} + (i - 1) * m;

                % biuld up preconditioner for each submatrix
                [PREC1{i}, PREC2{i}] = iluk(temp2, 2);

            end
            toc
            pl = cell2mat(pl);
            pr = cell2mat(pr);
            Dl = blkdiag(Dl{:});
            Dr = blkdiag(Dr{:});
        end
        solvePREChelp = @(rhs) solvePREC(PREC1, PREC2, rhs, m);
        %[sol] = gmres(Dl*M(pl,pr)*Dr,Dl*b(pl),100,10^(-5),3,solvePREChelp);

        [sol] = bicgstab(Dl*M(pl, pr)*Dr, Dl*b(pl), 10^(-5), 300, solvePREChelp);

        back(pr) = 1:(m * numTransportVariables);
        sol = Dr * sol;
        sol = sol(back);


    case 'StandardIterative'

        [PREC1, PREC2] = ilu(M);
        [sol] = bicgstab(M, b, 2*10^(-2), 300, PREC1, PREC2);
        
    case 'ilukIterative'
        if iteration == 0
            [l, u, levs] = iluk(M, 2);
        else
            [l, u] = iluk(M, [], levs); 
        end
        sol = bicgstab(M, b, 4*10^(-2), 400, l, u);
end

    function help = solvePREC(PREC1, PREC2, rhs, m)
        num = size(rhs, 1) / m;
        for i = 1:num
            [x{i, 1}] = PREC2{i} \ (PREC1{i} \ rhs((i-1) * m + 1:i * m));
        end
        help = cell2mat(x);
    end
end
