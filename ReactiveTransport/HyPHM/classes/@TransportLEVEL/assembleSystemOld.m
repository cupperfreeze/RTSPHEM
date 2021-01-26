%> @file TransportLEVEL/assembleSystem.m Assemble the linear system for the Transport problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Assemble the linear system for the Transport problem.
%>
%> [Qnew, Unew] = assembleSystem(d, Upasts, BDFq, curTau, curTime, dataA, dataAold, dataAvold, dataB, dataBold, dataC, dataD, dataE, dataF, datauD, datagF)
%>
%> @sa HyPHM, TransportLEVEL, Grid


%> @param   this       [ TransportLEVEL ]  Data required by solver.
%> @param   Upasts     [ cell{BDFq} ]  Scalar unknowns of old time steps [#T].
%> @param   BDFq       [ scalar    ]  Order of BDF time discretization.
%> @param   curTau     [ scalar    ]  Time step length.
%> @param   curTime    [ scalar    ]  Evaluation time (of new solution).
%> @param   dataA      [ #T        ]  Coordinates of A  in @f$\mathbb{P}_0(\mathcal{T})@f$ at current time.
%> @param   dataAold   [ #T        ]  Coordinates of A  in @f$\mathbb{P}_0(\mathcal{T})@f$ at previous time.
%> @param   dataC      [ #E        ]  Coordinates of C  in @f$\mathbb{RT}_0(\mathcal{T})@f$ at current time.
%> @param   dataD      [ #T x 4    ]  Coordinates of D  in @f$\mathbb{P}_0(\mathcal{T})^{4}@f$ at current time.
%> @param   dataE      [ #E        ]  Coordinates of E  in @f$\mathbb{RT}_0(\mathcal{T})@f$ at current time.
%> @param   datauD     [ #E        ]  Coordinates of uD in @f$\mathbb{P}_0(\mathcal{E})@f$ on edges at current time (def on ALL edges, but maybe NaN on non-Dirichlet ones).
%> @param   datagF     [ #E        ]  Coordinates of gF in @f$\mathbb{P}_0(\mathcal{E})@f$ on edges at current time (def on ALL edges, but maybe NaN on non-flux ones).
%> @param   isSlt      [ logical   ]  Suppress output to stdout (is silent)?
%>
%> @retval   Qnew       [ #E ]  flux solution in @f$\mathbb{RT}_0(\mathcal{T})@f$ basis
%> @retval   Unew       [ #T ]  scalar solution in @f$\mathbb{P}_0(\mathcal{T})@f$
%>


function [B, C, D, E, bQ, bU, Y] = assembleSystem(this, Upasts, BDFq, curTau, ...
    dataA, dataAold, dataAvold, dataB, dataBold, dataC, dataD, isDstationary, dataE, dataF, ...
    datauD, datagF, markDirE, markNeumE, markFluxE, idxNeumE, isSlt) %#ok<*SPRIX>

g = this.grid;
numT = g.numT;
numE = g.numE;
assert(BDFq == 1 || BDFq == 2, 'HyPHM: Requested BDF order of %d not implemented.', BDFq)
assert(iscell(Upasts)), assert(length(Upasts) == BDFq)

%%%%%%%%%%%%%%

%% Assembly %%
%%%%%%%%%%%%%%

tAssembly = tic;
printline(~isSlt*2, 'Assembly of linear terms for the problem ''%s''', this.name)

%% Assembly of H's
if isDstationary && ~isempty(this.globH)
    % load
    printline(~isSlt*3, 'Coefficient D did not change in step %d.  Loading matrix of step %d instead of a new assembly.', this.stepper.curstep, this.stepper.curstep-1)
    HlocT = this.globH;
else % D is instationary or matrix hasn't been assembled, yet
    % assemble
    HlocT = zeros(numT, 3, 3);
    dataDinv2 = [dataD(:, 4), -dataD(:, 2), -dataD(:, 3), dataD(:, 1)] ./ (dataD(:, 1) .* dataD(:, 4) - dataD(:, 2) .* dataD(:, 3));
    for kT = 1:numT %improved SG
        dataDinv = reshape(dataDinv2(kT, :), 2, 2);
        HlocT(kT, :, :) = locMatAlg(this, kT, dataDinv);
    end
    % storage
    if isDstationary
        this.globH = HlocT;
    end
end

%% Assembly of B
if isDstationary && ~isempty(this.globB) % load
    B = this.globB;
else % assemble (optimized!)
    B = sparse(g.E0T(:, [1, 2, 3, 1, 2, 3, 1, 2, 3]), ...
        g.E0T(:, [1, 1, 1, 2, 2, 2, 3, 3, 3]), ...
        -HlocT(:, :), numE, numE); % note that the second index (:, 1:9) stands for (:, 1:3, 1:3)
    % storage
    if isDstationary
        this.globB = B;
    end
end
% Assembly of B (SLOW VERSION)
%   B = sparse(numE, numE);
%   for kT = 1 : numT
%     idxsE = g.E0T(kT, :); % edge numbers of kTth element
%     B(idxsE, idxsE) = B(idxsE, idxsE) - HlocT{kT};
%   end

%% Assembly of stationary parts of D
if ~isempty(this.globD) % load
    D = this.globD;
else % assemble
    D = sparse(repmat((1:g.numT)', 1, 3), ...
        g.E0T, ...
        g.sigE0T.*g.areaE(g.E0T), g.numT, g.numE);
    this.globD = D; % storage
end
% Assembly of stationary parts of D (SLOW VERSION)
% D = sparse(g.numT, g.numE); % [#T x #E]
% for kT = 1 : g.numT  % sigma_EkTm |E_k|
%   D(kT, g.E0T(kT, :)) = g.sigE0T(kT, :)' .* g.areaE(g.E0T(kT, :));
% end

%% Calculation of local Peclet Numbers.
Pe = tp_pecletnumbers(this, dataC, dataD); % [#T x 1]
printline(~isSlt*3, 'Maximum Peclet number: Pe_max = %.3e', max(Pe))

%% Assembly of C (optionally with upwinding)
C = sparse(numE, numT); %[#E x #T]
switch this.isUpwind
    case 'none' % no upwinding
        alpha = @(z) ones(size(z));
    case 'full' % full upwinding
        alpha = @(z) 0.5 * (sign(z) + 1);
    case 'exp' % exponential upwinding
        alpha = @(z) alpha_exp(z);
    case 'alt' % exponential like upwinding, but only if PeT > 2
        alpha = @(z) alpha_alt(z);
    otherwise
        error('HyPHM: Unknown upwind flag (Transport.isUpwind).  Admissable upwind flags are ''none'', ''full'', ''exp'', ''alt''.')
end

idxE = g.E0T(:, :);
dataCET = dataC(idxE) .* g.sigE0T(:, :);
alphas = alpha(sign(dataCET(:, :)).*Pe(:));

marknonneumE = ~markNeumE(g.E0T(:, :));

kEs = zeros(3*numT, 1);
kTs = zeros(3*numT, 1);
ads = zeros(3*numT, 1);
value1 = zeros(3*numT, 1);
value2 = zeros(3*numT, 1);

for kT = 1:numT
    locmarknonneumE = marknonneumE(kT, :);
    %idxE    = g.E0T(kT, :);               % edge numbers of T_k
    idxEwoN = g.E0T(kT, locmarknonneumE); %          ''          exclusive Neum edges
    % locC = squeeze(HlocT(kT, locmarknonneumE, :))' * dataC(idxEwoN);
    locC = reshape(HlocT(kT, locmarknonneumE, :), sum(locmarknonneumE), 3)' * dataC(idxEwoN); %improved SG

    %dataCET = dataC(idxE) .* g.sigE0T(kT, :)'; % [3 x 1] flow over boundary (negative value: inflow, pos. value: outflow)


    for iET = 1:3 % loop over local edges
        kE = g.E0T(kT, iET); % current edge (global index)
        % define the index of the adjacent triangle adjT
        if all(g.T0E(g.E0T(kT, iET), :)) % upwinding if E is inflow edge and there exists an adjacent inflow triangle
            adjTidx = ones(3, 1);
            adjTidx = adjTidx + (g.sigE0T(kT, :) > 0)'; % column-index in g.T0E (ordering depds. on edge signs)
            adjT = g.T0E(g.E0T(kT, iET), adjTidx(iET)); % printline(3, 'Upwinding: T: %d:  E = %d --> take u_T from T == %d', kT, kE, upT)
        else
            adjT = kT; % if current edge is boundary edge -> no upwinding possible
        end
        %C(kE, kT)   = C(kE, kT)   +   alphas(kT,iET)   * locC(iET); % if upT==kT or if alpha = 0 we have the non-stabilized version
        %C(kE, adjT) = C(kE, adjT) + (1-alphas(kT,iET)) * locC(iET);

        index = 3 * (kT - 1) + iET;
        kEs(index) = kE;
        kTs(index) = kT;
        ads(index) = adjT;
        value1(index) = alphas(kT, iET) * locC(iET);
        value2(index) = (1 - alphas(kT, iET)) * locC(iET);
    end
end
C = sparse([kEs; kEs], [kTs; ads], [value1, value2]);

%% Assembly of E
if BDFq == 1
    E = sparse(1:numT, 1:numT, g.areaT.*dataA);
elseif BDFq == 2
    E = sparse(1:numT, 1:numT, 3/2*g.areaT.*dataA);
end
for kE = idxNeumE % update of E due to Neumann boundary conditions     % loop over all Neumann edges is cheaper than over triangles
    jT = g.T0E(kE, 1); % find triangle containing E_N
    idxNE = g.E0T(jT, markNeumE(g.E0T(jT, :))); % indices of Neumann edges on jT (may be 1 or 2)
    E(jT, jT) = E(jT, jT) + curTau * dot(dataC(idxNE), D(jT, idxNE));
end

%% Assembly of bQ [#E]
bQ = zeros(numE, 1);
bQ(markDirE) = g.areaE(markDirE) .* datauD(markDirE);
for kT = 1:numT
    idxE = g.E0T(kT, :);
    bQ(idxE) = bQ(idxE) - reshape(HlocT(kT, :, :), 3, 3)' * dataE(idxE); %improved SG
end

%% Assembly of bU [#T]
if BDFq == 1
    bU = g.areaT .* (curTau * dataF - dataB + dataBold + Upasts{1} .* dataAold);
elseif BDFq == 2
    bU = g.areaT .* (curTau * dataF - dataB + dataBold + 2 * Upasts{1} .* dataAold - 1 / 2 * Upasts{2} .* dataAvold);
end

% Flux and Neumann conditions to bQ and bU and computing flux data on
% non-free boundary edges in y (requires B and D) (cf remark 2)
Y = zeros(numE, 1); % [#E x 1]
Y(markFluxE) = datagF(markFluxE);
Y(markNeumE) = dataE(markNeumE);

bQ = bQ - B * Y;
bU = bU - curTau * D * Y;

clear HlocT

printline(~isSlt*3, '                      ...done [%.3f sec]', toc(tAssembly)) % stop timer assembly linear part


% coefficient for exponential upwinding (requires if-else construct which
% is forbidden in a function_handle)
%   function ret = alpha_exp(z)
%     if z ~= 0
%       ret = 1 - (1/z)*(1 - z/(exp(z) - 1));
%     else
%       ret = 0.5;
%     end
%   end


    function ret = alpha_exp(z)
        ret = 0.5 * ones(size(z));
        idx = find(z);
        h = z(idx);
        ret(idx) = 1 - (1 ./ h) .* (1 - h ./ (exp(h) - 1));
    end

% coefficient for the alternative upwinding (requires if-else construct which
% is forbidden in a function_handle)
%   function ret = alpha_alt(z)
%     if z ~= 0
%       ret =  0.5*(1 + sign(z) * max([0, 1 - 2/abs(z)]));
%     else
%       ret = 0.5;
%     end
%   end

    function ret = alpha_alt(z)
        ret = 0.5 * ones(size(z));
        idx = find(z);
        h = z(idx);
        ret(idx) = 0.5 * (1 + sign(h) * max([0, 1 - 2 ./ abs(h)]));
    end

end


%%%%%%%%%%%%%%%

%% Footnotes %%
%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remark 2: Variable Y
% Incorporation of flux data on non-free edges (won't be touched by solver)
% reminder: the degrees of freedom wrt to the flux lie on the inner and
% Dirichlet edges whereas the flux data on Neumann and flux edges are given
% by respective boundary conditions.  This data is written to the (yet
% uncomputed) solution vector.  The solver works on a subindex set and does
% not change these values.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remark 5: Calculating additional Flux Terms for Neumann Edges
% See also remark 2.  The flux-solution on Neumann edge depends directly of the
% solution U itself.  The dependency is neglected within the solution and
% postprocessed by delay.
