//! @file stdEllipse.geo Definition of a elliptic bounded domain @f$\Omega@f$.
//!

h = 0.1;
cx = 0.5;
cy = 1.5;
ra = 1.0;
rb = 0.6;

Point(0) = {cx, cy,   0, h};
Point(1) = {cx-ra, cy, 0, h};
Point(2) = {cx, cy+rb, 0, h};
Point(3) = {cx+ra, cy, 0, h};
Point(4) = {cx, cy-rb, 0, h};

Ellipse(5) = {1,0,1,2};
Ellipse(6) = {2,0,2,3};
Ellipse(7) = {3,0,3,4};
Ellipse(8) = {4,0,4,1};

Line Loop(9) = {5,6,7,8};

Plane Surface(10) = {9};
