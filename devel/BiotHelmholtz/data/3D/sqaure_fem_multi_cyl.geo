
/*number of elements*/
nbele = 5;
ta = 0.25/nbele;
tf = 0.0008/1;
t2 = 0.0004/1;
d = 0.1;

Point(1) = {-0.75, 0, 0, ta};
//+
Point(2) = {0.250, 0, 0, ta};
//+
Point(3) = {0.250, 0.20, 0, ta};
//+
Point(4) = {-0.75, 0.20, 0, ta};


//+
Point(5) = {0.1+d, 0.1, 0, tf};
//+
Point(6) = {-0.+d, 0.2, 0, tf};
//+
Point(7) = {-0.+d, 0, 0, tf};
//+
Circle(1) = {6, 5, 7};
//+
Point(8) = {0.0006+d, 0.2, 0, tf};
//+
Point(9) = {0.0006+d, 0, 0, tf};
//+
Circle(2) = {8, 5, 9};

Point(10) = {0.0014+d, 0.2, 0, tf};
//+
Point(11) = {0.0014+d, 0, 0, tf};
//+
Circle(3) = {10, 5, 11};

Point(12) = {0.002+d, 0.2, 0, tf};
//+
Point(13) = {0.002+d, 0, 0, tf};
//+
Circle(4) = {12, 5, 13};

//+
Line(5) = {1, 7};
//+
Line(6) = {7, 9};
//+
Line(7) = {9, 11};
//+
Line(8) = {11, 13};
//+
Line(9) = {13, 2};
//+
Line(10) = {2, 3};
//+
Line(11) = {3, 12};
//+
Line(12) = {12, 10};
//+
Line(13) = {10, 8};
//+
Line(14) = {8, 6};
//+
Line(15) = {6, 4};
//+
Line(16) = {4, 1};
//+
Curve Loop(1) = {16, 5, -1, 15};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 6, -2, 14};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, 7, -3, 13};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 8, -4, 12};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 9, 10, 11};
//+
Plane Surface(5) = {5};
//+
Physical Curve(101) = {16};
//+
Physical Curve(102) = {10};
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
