//! @file tsc_OMe.geo Definition of the domain @f$\Omega_\varepsilon@f$.
//!

Include "tsc_geofuns.geo";   // geometry
Include "tsc_epsilon.geo";   // scaling parameter eps
Include "tsc_meshwidth.geo"; // the mesh width h 

// Line Loop for outer surface: Line Loop(theloops[0]) (has to be first Line Loop in Plane Surface generation)
xmin = 0.0;
xmax = 1.0;
ymin = 0.0;
ymax = eps;
Call Rectangle; // will have ID 1 to 4

// Line Loops for inner surfaces: Line Loop(theloops[t]), t = 1, 2, ...
For t In {1/eps:1:-1} // the expression `xtrans = 0' has to go last (that's a Bug (?) workaround)
  xtrans = (t-1)*eps;
  Printf("------------------------------");
  Printf("Calling Obstacles for t: %g.", t);
  Call Obstacles;
  Printf("Called Obstacles for t: %g.", t);
  Printf("------------------------------");
EndFor

ps0 = newreg; Plane Surface(ps0) = {theloops[]};

