//! @file tsc_OM.geo Definition of the domain @f$\Omega@f$.

Include "tsc_geofuns.geo";   // geometry
Include "tsc_epsilon.geo";   // scaling parameter eps
Include "tsc_meshwidth.geo"; // the mesh width h 

// Line Loop for outer surface: Line Loop(theloops[0]) (has to be first Line Loop in Plane Surface generation)
xmin = 0.0;
xmax = 1.0;
ymin = 0.0;
ymax = eps;
Call Rectangle;

ps0 = newreg; Plane Surface(ps0) = {theloops[]};
