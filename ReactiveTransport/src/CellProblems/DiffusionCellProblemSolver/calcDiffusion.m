%Compute Diffusion tensor using a mixed discretization

function DiffusionTensor = calcDiffusion(grid, LevelSet, diffCoeff)

DiffusionTensor = zeros(0);
st = Stepper(0:1);
problem = TransportLEVEL(grid, st, 'Diffusion');
problem.L.setdata(LevelSet);
problem.C.setdata(@(t, x) [0, 0]);
problem.E.setdata(@(t, x) [0, 0]);
problem.D.setdata(eye(2));
problem.gF.setdata(-grid.nuE(:, 1));

fluidNodes = (LevelSet < -eps);
fluidTriangles = false(grid.numT, 1);
for tri = 1:grid.numT
    fluidTriangles(tri) = (sum(LevelSet(grid.V0T(tri, :))) < 2.5);
end


st.next;
problem.computeLevel('s');
firstout = problem.Q.getdata(1);
st.prev;
problem.gF.setdata(-grid.nuE(:, 2));
st.next;
problem.computeLevel('s');
secondout = problem.Q.getdata(1);


fluidArea = grid.areaT' * fluidTriangles;

q = zeros(grid.numT, 2);
for kT = 1:grid.numT
    if fluidTriangles(kT)
        q(kT, :) = getCartesianCoords(grid, firstout(grid.E0T(kT, :))', kT, grid.baryT(kT, :)');
    end
end

DiffusionTensor(1, 1) = fluidArea + dot(q(:, 1), grid.areaT);
DiffusionTensor(2, 1) = dot(q(:, 2), grid.areaT);

q = zeros(grid.numT, 2);
for kT = 1:grid.numT
    if fluidTriangles(kT)
        q(kT, :) = getCartesianCoords(grid, secondout(grid.E0T(kT, :))', kT, grid.baryT(kT, :)');
    end
end

DiffusionTensor(1, 2) = dot(q(:, 1), grid.areaT);
DiffusionTensor(2, 2) = fluidArea + dot(q(:, 2), grid.areaT);
end
