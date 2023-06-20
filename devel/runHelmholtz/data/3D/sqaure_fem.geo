
/*number of elements*/
nbele = 11;
ta = 0.2/nbele;
tf = 0.001/1;

Point(1) = {-0.1, 0, 0, ta};
//+
Point(2) = {0.1, 0, 0, ta};
//+
Point(3) = {0.1, 0.05, 0, ta};
//+
Point(4) = {-0.1, 0.05, 0, ta};
//+
Point(5) = {-0.1, 0., 0.05, ta};
//+
Point(6) = {0.1, 0, 0.05, ta};
//+
Point(7) = {0.1, 0.05, 0.05, ta};
//+
Point(8) = {-0.1, 0.05, 0.05, ta};

Point(9) = {-0.001, 0, 0, tf};
Point(10) = {-0.0, 0, 0, tf};
Point(11) = {-0.001, 0.05, 0, tf};
Point(12) = {-0.0, 0.05, 0, tf};

Point(13) = {-0.001, 0, 0.05, tf};
Point(14) = {-0.0, 0, 0.05, tf};
Point(15) = {-0.001, 0.05, 0.05, tf};
Point(16) = {-0.0, 0.05, 0.05, tf};

//+
Line(1) = {5, 13};
//+
Line(2) = {13, 14};
//+
Line(3) = {14, 6};
//+
Line(4) = {6, 7};
//+
Line(5) = {7, 16};
//+
Line(6) = {16, 15};
//+
Line(7) = {15, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {1, 9};
//+
Line(10) = {9, 10};
//+
Line(11) = {10, 2};
//+
Line(12) = {2, 3};
//+
Line(13) = {3, 12};
//+
Line(14) = {12, 11};
//+
Line(15) = {11, 4};
//+
Line(16) = {4, 1};
//+
Line(17) = {1, 5};
//+
Line(18) = {4, 8};
//+
Line(19) = {13, 15};
//+
Line(20) = {15, 11};
//+
Line(21) = {11, 9};
//+
Line(22) = {9, 13};
//+
Line(23) = {14, 16};
//+
Line(24) = {16, 12};
//+
Line(25) = {12, 10};
//+
Line(26) = {10, 14};
//+
Line(27) = {7, 3};
//+
Line(28) = {6, 2};
//+
Curve Loop(1) = {18, 8, -17, -16};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 19, 7, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, 23, 6, -19};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {23, -5, -4, -3};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 27, -12, -28};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {12, 13, 25, 11};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {25, -10, -21, -14};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {21, -9, -16, -15};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {7, -18, -15, -20};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {20, -14, -24, 6};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {24, -13, -27, 5};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {17, 1, -22, -9};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {22, 2, -26, -10};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {26, 3, 28, -11};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {20, 21, 22, 19};

Plane Surface(15) = {15};
//+
Curve Loop(16) = {23, 24, 25, 26};
//+
Plane Surface(16) = {16};
//+
Surface Loop(1) = {1, 9, 2, 12, 8, 15};
//+
Physical Surface(1) = {1};
//+
Physical Surface(2) = {5};

Volume(1) = {1};
//+
Surface Loop(2) = {15, 16, 3, 13, 7, 10};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {5, 4, 11, 6, 14, 16};
//+
Volume(3) = {3};
//+
Physical Volume(1) = {1};
//+
Physical Volume(2) = {2};
//+
Physical Volume(3) = {3};
