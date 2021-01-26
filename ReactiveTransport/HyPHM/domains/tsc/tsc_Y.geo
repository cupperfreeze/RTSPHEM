//! @file tsc_Y.geo Definition of the domain @f$Y@f$.

Include "tsc_geofuns.geo";   // geometry
Include "tsc_epsilon.geo";   // scaling parameter eps
Include "tsc_meshwidth.geo"; // the mesh width h 

// Line Loop for outer surface: Line Loop(theloops[0]) (has to be first Line Loop in Plane Surface generation)
xmin = 0.0;
xmax = 1.0;
ymin = 0.0;
ymax = 1.0;
Call Rectangle; // will have ID 1 to 4

// Line Loop for inner surfaces: Line Loop(theloops[1])
t = 1;
xtrans = 0.0;
eps = 1; 
Call Obstacles;

ps0 = newreg; Plane Surface(ps0) = {theloops[]};
