//! @file stdCircle.geo Definition of a circular bounded domain @f$\Omega@f$.
//!

h = 0.1;
cx = 0.5;
cy = 1.5;
r = 1.0;

Point(0) = {cx, cy,   0, h};
Point(1) = {cx-r, cy, 0, h/3};
Point(2) = {cx, cy+r, 0, h/3};


Circle(6) = {2,0,1};

Point(11) = {cx-2*r, cy, 0, h};
Point(12) = {cx, cy+2*r, 0, h};
Point(13) = {cx+2*r, cy, 0, h};
Point(14) = {cx, cy-2*r, 0, h};

Circle(15) = {11,0,12};
Circle(16) = {12,0,13};
Circle(17) = {13,0,14};
Circle(18) = {14,0,11};


Line Loop(10) = {15,16,17,18};
Line Loop(20) = {6};

Plane Surface(19) = {10,20};

