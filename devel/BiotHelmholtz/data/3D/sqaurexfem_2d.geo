nbele = 5;
ta = 0.25/nbele;

//+
Point(1) = {-0.75, 0, 0, ta};
//+
Point(2) = {0.25-0.002, 0, 0, ta};
//+
Point(3) = {0.25-0.002, 0.2, 0, ta};
//+
Point(4) = {-0.75, 0.2, 0, ta};

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
Physical Curve(101) = {4};
//+
Physical Curve(102) = {1};
//+
Physical Curve(103) = {2};
//+
Physical Curve(104) = {3};
//+
Physical Surface(5) = {1};

Physical Point(8) = {1};
