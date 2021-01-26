//! @file stdRectangle.geo Definition of a rectangular bounded domain @f$\Omega = [x_\mathrm{min},x_\mathrm{max}]\times[y_\mathrm{min},y_\mathrm{max}]@f$.
//!

h = 1;
xmin = -1/2;
xmax = 1/2;
ymin = -1/2;
ymax = 1/2;

Point(0) = {xmin, ymin, 0, h};
Point(1) = {xmax, ymin, 0, h};
Point(2) = {xmax, ymax, 0, h};
Point(3) = {xmin, ymax, 0, h};

Line(1)  = {0, 1};
Line(2)  = {1, 2};
Line(3)  = {2, 3};
Line(4)  = {3, 0};

Line Loop(9) = {1,2,3,4}; // edge IDs will be 1 to 4

Plane Surface(10) = {9};
