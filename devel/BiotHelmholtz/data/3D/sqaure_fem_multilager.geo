
/*number of elements*/
nbele = 5;
ta = 0.25/nbele;
tf = 0.0005/1;
t2 = 0.001/1;
Point(1) = {-0.75, 0, 0, ta};
//+
Point(2) = {0.250, 0, 0, ta};
//+
Point(3) = {0.250, 0.20, 0, ta};
//+
Point(4) = {-0.75, 0.20, 0, ta};
//+
Point(9) = {0.0006, 0, 0, tf};
Point(10) = {-0.0, 0, 0, tf};
Point(11) = {0.2006, 0.20, 0, tf};
Point(12) = {0.2, 0.20, 0, tf};

Point(13) = {0.0014, 0, 0, tf};
Point(14) = {-0.0, 0, 0, tf};
Point(15) = {0.2014, 0.20, 0, tf};
Point(16) = {0.2, 0.20, 0, tf};

Point(17) = {0.002, 0, 0, tf};
Point(18) = {-0.0, 0, 0, tf};
Point(19) = {0.202, 0.20, 0, tf};
Point(20) = {0.2, 0.20, 0, tf};


//+
Line(1) = {1, 10};
//+
Line(2) = {10, 9};
//+
Line(3) = {9, 13};
//+
Line(4) = {13, 17};
//+
Line(5) = {17, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 19};
//+
Line(8) = {19, 15};
//+
Line(9) = {15, 11};
//+
Line(10) = {11, 12};
//+
Line(11) = {12, 4};
//+
Line(12) = {4, 1};
//+
Line(13) = {12, 10};
//+
Line(14) = {11, 9};
//+
Line(15) = {15, 13};
//+
Line(16) = {19, 17};
//+
Curve Loop(1) = {11, 12, 1, -13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, 2, -14, 10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {14, 3, -15, 9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, 15, 4, -16};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {16, 5, 6, 7};
//+
Plane Surface(5) = {5};
//+
Physical Curve(101) = {12};
//+
Physical Curve(102) = {6};
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
