nbele = 5;
ta = 0.2/nbele;

//+
Point(1) = {-0.1+0.001, 0, 0, ta};
//+
Point(2) = {0.1, 0, 0, ta};
//+
Point(3) = {0.1, 0.05, 0, ta};
//+
Point(4) = {-0.1+0.001, 0.05, 0, ta};
//+
Point(5) = {-0.1+0.001, 0., 0.05, ta};
//+
Point(6) = {0.1, 0, 0.05, ta};
//+
Point(7) = {0.1, 0.05, 0.05, ta};
//+
Point(8) = {-0.1+0.001, 0.05, 0.05, ta};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {2, 6};
//+
Line(11) = {3, 7};
//+
Line(12) = {4, 8};
//+
Curve Loop(1) = {12, 8, -9, -4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {11, -6, -10, 2};
//+
Plane Surface(2) = {2};
//+
//+
Curve Loop(3) = {8, 5, 6, 7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 4, 1, 2};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {3, 12, -7, -11};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {1, 10, -5, -9};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 5, 4, 6, 2, 3};
//+
Volume(1) = {1};

Physical Surface(101) = {1};
//+
Physical Surface(102) = {2};

Physical Volume(121) = {1};


//+
Transfinite Curve {7, 3, 5, 1} = 10 Using Progression 1;
//+
Transfinite Curve {12, 9, 8, 4} = 3 Using Progression 1;
//+
Transfinite Curve {11, 2, 10, 6} = 3 Using Progression 1;
