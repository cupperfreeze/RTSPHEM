//! @file BC2005.geo Definition of a L-shaped domain @f$\Omega = \,]-1,1[^2\setminus [0,1]\times [-1,0]@f$.
//!

h = 1; // will be scaled by compiling option

Point(0) = {-1, -1, 0, h};
Point(1) = {0,  -1, 0, h};
Point(2) = {0,   0, 0, h}; // Point(2) = {0,   0, 0, h/10};
Point(3) = {1,   0, 0, h};
Point(4) = {1,   1, 0, h};
Point(5) = {-1,  1, 0, h};

Line(1)  = {0, 1};
Line(2)  = {1, 2};
Line(3)  = {2, 3};
Line(4)  = {3, 4};
Line(5)  = {4, 5};
Line(6)  = {5, 0};

Line Loop(9) = {1,2,3,4,5,6}; // edge IDs will be 1 to 6

Plane Surface(10) = {9};
