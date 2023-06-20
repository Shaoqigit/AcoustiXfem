/*40 : F*/
nbele = 4;
nbptx = nbele+1;
nbpty = 51;
ta = 2/nbele;

//+
Point(1) = {-0.1, 0, -0, ta};
//+
Point(2) = {0.1, 0, -0, ta};
//+
Point(3) = {0.1, 0.2, -0, ta};
//+
Point(4) = {-0.1, 0.2, -0, ta};
//+
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};


//+
Point(5) = {0, 0.1, 0, ta};
//+
Point(6) = {0, 0.06, 0, ta};
//+
Point(7) = {0, 0.14, 0, ta};
//+
Point(8) = {0, 0.059, 0, ta};
//+
Point(9) = {0, 0.141, 0, ta};
//+

//+
Circle(5) = {9, 5, 8};
//+
Circle(6) = {8, 5, 9};
//+
Circle(7) = {7, 5, 6};
//+
Circle(8) = {6, 5, 7};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5, 6};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {7, 8};
//+
Plane Surface(2) = {2, 3};
//+
Plane Surface(3) = {3};
//+
Physical Curve(1) = {4, 1, 3, 2};
//+

Physical Surface(2) = {1};
//+
Physical Surface(3) = {3};
Physical Surface(4) = {2};
//+
MeshSize {4, 1, 2, 3} = 0.02;
//+
Transfinite Curve {5, 7, 8, 6} = 40 Using Progression 1;
//+
Transfinite Curve {5, 6} = 60 Using Progression 1;
//+
/*Transfinite Curve {7, 8} = 60 Using Progression 1;*/
//+
MeshSize {5} = 0.02;
//+
MeshSize {7, 6} = 0.005;
