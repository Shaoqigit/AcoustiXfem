
/*number of elements*/
nbele = 5;
ta = 0.25/nbele;
tf = 0.0008/1;
t2 = 0.0004/1;
Point(1) = {-0.75, 0, 0, ta};
//+
Point(2) = {0.250, 0, 0, ta};
//+
Point(3) = {0.250, 0.20, 0, ta};
//+
Point(4) = {-0.75, 0.20, 0, ta};



//+
Point(5) = {0.25, 0.1, 0, tf};
//+
Point(6) = {-0., 0.2, 0, tf};
//+
Point(7) = {-0., 0, 0, tf};
//+
Circle(1) = {6, 5, 7};

Point(8) = {0.25, 0.1, 0, tf};
//+
Point(9) = {0.002, 0.2, 0, tf};
//+
Point(10) = {0.002, 0, 0, tf};
//+
Circle(2) = {9, 8, 10};

/*Point(11) = {0.25-0.2693, 0.1, 0, tf};*/
//+
Line(3) = {1, 7};
//+
Line(4) = {7, 10};
//+
Line(5) = {10, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 9};
//+
Line(8) = {9, 6};
//+
Line(9) = {6, 4};
//+
Line(10) = {4, 1};
//+
Curve Loop(1) = {10, 3, -1, 9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 4, -2, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, 2, 5, 6};
//+
Plane Surface(3) = {3};
//+
Physical Curve(101) = {10};
//+
Physical Curve(102) = {6};
//+
Physical Surface(103) = {1};
//+
Physical Surface(104) = {2};
//+
Physical Surface(105) = {3};

