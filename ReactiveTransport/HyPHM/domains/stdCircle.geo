//! @file stdCircle.geo Definition of a circular bounded domain @f$\Omega@f$.
//!

h = 1;
cx = 0.0;
cy = 0.0;
r = 0.5;

Point(0) = {cx, cy,   0, h};
Point(1) = {cx-r, cy, 0, h};
Point(2) = {cx, cy+r, 0, h};
Point(3) = {cx+r, cy, 0, h};
Point(4) = {cx, cy-r, 0, h};

Circle(5) = {1,0,2};
Circle(6) = {2,0,3};
Circle(7) = {3,0,4};
Circle(8) = {4,0,1};

Line Loop(9) = {5,6,7,8};

Plane Surface(10) = {9};
