
ta = 2/8;

//+
Point(1) = {-1, -1, 0, ta};
//+
Point(2) = {1, -1, 0, ta};
//+
Point(3) = {1, 1, 0, ta};
//+
Point(4) = {-1, 1, 0, ta};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve(1) = {1};
//+
Physical Curve(2) = {2};
//+
Physical Curve(3) = {3};
//+
Physical Curve(4) = {4};
//+
Physical Surface(1000) = {1};
//+
Physical Point(1) = {1};
//+
Physical Point(2) = {2};
//+
Physical Point(3) = {3};
//+
Physical Point(4) = {4};
