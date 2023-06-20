
/*number of elements*/
nbele = 11;
ta = 0.2/nbele;
tf = 0.0006/1;
t2 = 0.0004/1;
Point(1) = {-0.100, 0, 0, ta};
//+
Point(2) = {0.1000, 0, 0, ta};
//+
Point(3) = {0.1000, 0.05, 0, ta};
//+
Point(4) = {-0.100, 0.05, 0, ta};
//+
Point(9) = {-0.001, 0, 0, tf};
Point(10) = {-0.0, 0, 0, tf};
Point(11) = {-0.001, 0.05, 0, tf};
Point(12) = {-0.0, 0.05, 0, tf};

Point(13) = {-0.0007, 0, 0, t2};
Point(14) = {-0.0, 0, 0, t2};
Point(15) = {-0.0007, 0.05, 0, t2};
Point(16) = {-0.0, 0.05, 0, t2};

Point(17) = {-0.0003, 0, 0, t2};
Point(18) = {-0.0, 0, 0, t2};
Point(19) = {-0.0003, 0.05, 0, t2};
Point(20) = {-0.0, 0.05, 0, t2};

//+
Line(1) = {1, 9};
//+
Line(2) = {9, 13};
//+
Line(3) = {13, 17};
//+
Line(4) = {17, 10};
//+
Line(5) = {10, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 12};
//+
Line(8) = {12, 19};
//+
Line(9) = {19, 15};
//+
Line(10) = {15, 11};
//+
Line(11) = {11, 4};
//+
Line(12) = {4, 1};
//+
Curve Loop(1) = {12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
//+
Line(13) = {9, 11};
//+
Line(14) = {15, 13};
//+
Line(15) = {17, 19};
//+
Line(16) = {10, 12};
//+
Curve Loop(2) = {11, 12, 1, 13};
//+
Plane Surface(1) = {2};
//+
Curve Loop(3) = {13, -10, 14, -2};
//+
Plane Surface(2) = {3};
//+
Curve Loop(4) = {14, 3, 15, 9};
//+
Plane Surface(3) = {4};
//+
Curve Loop(5) = {15, -8, -16, -4};
//+
Plane Surface(4) = {5};
//+
Curve Loop(6) = {16, -7, -6, -5};
//+
Plane Surface(5) = {6};
//+
Physical Curve(101) = {12};
//+
Physical Curve(103) = {6};
//+
Physical Surface(110) = {1};
//+
Physical Surface(111) = {2};
//+
Physical Surface(112) = {3};
//+
Physical Surface(113) = {4};
//+
Physical Surface(114) = {5};
