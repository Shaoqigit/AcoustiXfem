
/*number of elements*/
nbele = 2;
ta = 0.25/nbele;
tf = 0.003/1;

Point(1) = {-0.2525, 0, 0, ta};
//+
Point(2) = {0.25, 0, 0, ta};
//+
Point(3) = {0.25, 0.5, 0, ta};
//+
Point(4) = {-0.2525, 0.5, 0, ta};
//+
Point(9) = {-0.0025, 0, 0, tf};
Point(10) = {-0.0, 0, 0, tf};
Point(11) = {-0.0025, 0.5, 0, tf};
Point(12) = {-0.0, 0.5, 0, tf};

//+
Line(1) = {1, 9};
//+
Line(2) = {9, 10};
//+
Line(3) = {10, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 12};
//+
Line(6) = {12, 11};
//+
Line(7) = {11, 4};
//+
Line(8) = {4, 1};
//+
Line(9) = {9, 11};
//+
Line(10) = {10, 12};
//+
Curve Loop(1) = {8, 1, 9, 7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, -6, -10, -2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 4, 5, -10};
//+
Plane Surface(3) = {3};
//+
Physical Curve(101) = {8};
//+
Physical Curve(102) = {4};
//+
Physical Surface(103) = {1};
//+
Physical Surface(104) = {2};
//+
Physical Surface(105) = {3};
//+
/*Transfinite Curve {7, 1, 5, 3} = 5 Using Progression 1;
//+
Transfinite Curve {8, 9, 10, 4} = 3 Using Progression 1;*/
