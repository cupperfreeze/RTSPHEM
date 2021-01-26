//! @file RBPA2008.geo Definition of a rectangular bounded domain @f$\Omega = [0,2]\times[0,3]@f$ for the test scenario @ref RBPA2008.
//! The mesh size is @f$h=1@f$ and can be scaled via scaling option in Grid generation.

h = 1;

Point(0) = {0.0, 0.0, 0, h};
Point(1) = {2.0, 0.0, 0, h};
Point(2) = {2.0, 3.0, 0, h};
Point(3) = {0.0, 3.0, 0, h};

Line(1)  = {0, 1};
Line(2)  = {1, 2};
Line(3)  = {2, 3};
Line(4)  = {3, 0};

Line Loop(9) = {1,2,3,4}; // edge IDs will be 1 to 4

Plane Surface(10) = {9};
