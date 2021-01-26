
a = 1.5; // width
b = 3.0; // heigth
h = 0.05;
h2 =  0.04;
Point(1) = {0, 0, 0, h};
Point(2) = {a, 0, 0, h};
Point(3) = {a, b, 0, h};
Point(4) = {0, b, 0, h};
Line(101) = {1,2} ;
Line(102) = {2,3} ;
Line(103) = {3,4} ;
Line(104) = {4,1} ;

// M1
cx = 0.60;
cy = 0.35;
k = 5;
Point(k)      = {cx-0.5, cy-0.05, 0, h2};
Point(k+1)    = {cx+0.5, cy-0.05, 0, h2};
Point(k+2)    = {cx+0.5, cy+0.05, 0, h2};
Point(k+3)    = {cx-0.5, cy+0.05, 0, h2};
Rotate { {0, 0, 1}, { cx, cy, 0},  0.0 } {Point{k,k+1,k+2,k+3}; }
Line(100+k)   = {k+3,k+2};
Line(100+k+1) = {k+2,k+1};
Line(100+k+2) = {k+1,k};
Line(100+k+3) = {k,  k+3};


// M2
cx = 0.98;
cy = 0.79;
k = 9;
Point(k)   = {cx-0.5, cy-0.05, 0, h2};
Point(k+1) = {cx+0.5, cy-0.05, 0, h2};
Point(k+2) = {cx+0.5, cy+0.05, 0, h2};
Point(k+3) = {cx-0.5, cy+0.05, 0, h2};
Rotate { {0, 0, 1}, { cx, cy, 0},  -0.698131700797732 } {Point{k,k+1,k+2,k+3}; }
Line(100+k)   = {k+3,k+2};
Line(100+k+1) = {k+2,k+1};
Line(100+k+2) = {k+1,k};
Line(100+k+3) = {k,  k+3};

// M3
cx = 0.55;
cy = 1.21;
k = 13;
Point(k)   = {cx-0.5, cy-0.05, 0, h2};
Point(k+1) = {cx+0.5, cy-0.05, 0, h2};
Point(k+2) = {cx+0.5, cy+0.05, 0, h2};
Point(k+3) = {cx-0.5, cy+0.05, 0, h2};
Rotate { {0, 0, 1}, { cx, cy, 0},  0.872664625997165 } {Point{k,k+1,k+2,k+3}; }
Line(100+k)   = {k+3,k+2};
Line(100+k+1) = {k+2,k+1};
Line(100+k+2) = {k+1,k};
Line(100+k+3) = {k,  k+3};

// M4
cx = 0.60;
cy = 1.86;
k = 17;
Point(k)   = {cx-0.5, cy-0.05, 0, h2};
Point(k+1) = {cx+0.5, cy-0.05, 0, h2};
Point(k+2) = {cx+0.5, cy+0.05, 0, h2};
Point(k+3) = {cx-0.5, cy+0.05, 0, h2};
Rotate { {0, 0, 1}, { cx, cy, 0},  -0.174532925199433  } {Point{k,k+1,k+2,k+3}; }
Line(100+k)   = {k+3,k+2};
Line(100+k+1) = {k+2,k+1};
Line(100+k+2) = {k+1,k};
Line(100+k+3) = {k,  k+3};

// M5
cx = 1.12;
cy = 2.09;
k = 21;
Point(k)   = {cx-0.5, cy-0.05, 0, h2};
Point(k+1) = {cx+0.5, cy-0.05, 0, h2};
Point(k+2) = {cx+0.5, cy+0.05, 0, h2};
Point(k+3) = {cx-0.5, cy+0.05, 0, h2};
Rotate { {0, 0, 1}, { cx, cy, 0}, -1.047197551196598 } {Point{k,k+1,k+2,k+3}; }
Line(100+k)   = {k+3,k+2};
Line(100+k+1) = {k+2,k+1};
Line(100+k+2) = {k+1,k};
Line(100+k+3) = {k,  k+3};

// M6
cx = 0.59;
cy = 2.60;
k = 25;
Point(k)   = {cx-0.5, cy-0.05, 0, h2};
Point(k+1) = {cx+0.5, cy-0.05, 0, h2};
Point(k+2) = {cx+0.5, cy+0.05, 0, h2};
Point(k+3) = {cx-0.5, cy+0.05, 0, h2};
Rotate { {0, 0, 1}, { cx, cy, 0},  0.349065850398866 } {Point{k,k+1,k+2,k+3}; }
Line(100+k)   = {k+3,k+2};
Line(100+k+1) = {k+2,k+1};
Line(100+k+2) = {k+1,k};
Line(100+k+3) = {k,  k+3};

Line Loop(200) = {101,102,103,104, 105,106,107,108,  109,110,111,112, 113,114,115,116,  117,118,119,120,  121,122,123,124, 125,126,127,128} ;

Plane Surface(300) = {200} ;



