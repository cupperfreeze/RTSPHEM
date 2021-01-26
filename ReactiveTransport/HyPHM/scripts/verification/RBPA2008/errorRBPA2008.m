%> @file errorRBPA2008.m Error estimation for @ref RBPA2008, each evaluated at end time.
%> Requires some stuff.


printline(1, 'Calculation of L2 Error')

sol_u = @RBPA2008_u; % exact scalar solution
sol_q = @RBPA2008_q; % exact flux solution
%sol_u   = @BC2005_u; % exact scalar solution
%sol_q   = @BC2005_q; % exact flux solution

endTime = st.endtime;
lastStep = st.numsteps;

scalUnk = {conc1, conc2};
fluxUnk = {mflux1, mflux2};

numSpecs = length(scalUnk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% L2 and P0 Errors of Scalars %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err_u_L2 = cell(numSpecs, 1);
err_u_P0 = cell(numSpecs, 1);

u_apprx = cell(numSpecs, 1);
for kSpec = 1:numSpecs
    u_apprx{kSpec} = scalUnk{kSpec}.getdata(lastStep);
end

parfor kSpec = 1:numSpecs


    err_u_L2{kSpec} = 0.0;
    for kT = 1:g.numT
        fun = @(x) abs(sol_u(endTime, x, kSpec)-u_apprx{kSpec}(kT))^2;
        err_u_L2{kSpec} = err_u_L2{kSpec} + intT(g, kT, fun, '612');
    end
    err_u_L2{kSpec} = err_u_L2{kSpec}^(1 / 2);

    err_u_P0{kSpec} = 0.0;
    for kT = 1:g.numT
        baryT = g.baryT(kT, :)';
        areaT = g.areaT(kT);
        err_u_P0{kSpec} = err_u_P0{kSpec} + areaT * (abs(sol_u(endTime, baryT, kSpec) - u_apprx{kSpec}(kT)))^2;
    end
    err_u_P0{kSpec} = err_u_P0{kSpec}^(1 / 2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% L2 and RT0 Errors of Fluxes %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err_q_L2 = cell(numSpecs, 1);
err_q_RT0 = cell(numSpecs, 1);

q_apprx = cell(numSpecs, 1);
for kSpec = 1:numSpecs
    q_apprx{kSpec} = fluxUnk{kSpec}.getdata(lastStep);
end

parfor kSpec = 1:numSpecs
    err_q_L2{kSpec} = 0.0;
    for kT = 1:g.numT
        edgeIdxs = g.E0T(kT, :);
        fun = @(x) norm(sol_q(endTime, x, kSpec)-RT0.getCartesianCoords(g, q_apprx{kSpec}(edgeIdxs)', kT, x))^2;
        err_q_L2{kSpec} = err_q_L2{kSpec} + intT(g, kT, fun, '612');
    end
    err_q_L2{kSpec} = err_q_L2{kSpec}^(1 / 2);

    err_q_RT0{kSpec} = 0.0;
    for kT = 1:g.numT
        areaT = g.areaT(kT);
        for iEdge = g.E0T(kT, :)
            baryE = g.baryE(iEdge, :)';
            nuE = g.nuE(iEdge, :)';
            err_q_RT0{kSpec} = err_q_RT0{kSpec} + areaT / 3 * (abs(sol_q(endTime, baryE, kSpec)' * nuE - q_apprx{kSpec}(iEdge)))^2;
        end
    end
    err_q_RT0{kSpec} = err_q_RT0{kSpec}^(1 / 2);

end

for kSpec = 1:numSpecs
    printline(3, 'L2 error of %s at t = %.1f: %.3e', scalUnk{kSpec}.name, endTime, err_u_L2{kSpec});
    printline(3, 'P0 error of %s at t = %.1f: %.3e', scalUnk{kSpec}.name, endTime, err_u_P0{kSpec});
    printline(3, 'L2 error of %s at t = %.1f: %.3e', fluxUnk{kSpec}.name, endTime, err_q_L2{kSpec});
    printline(3, 'RT0 error of %s at t = %.1f: %.3e', fluxUnk{kSpec}.name, endTime, err_q_RT0{kSpec});
end
