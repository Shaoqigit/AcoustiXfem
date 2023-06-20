nbele = 30;

SetFactory("OpenCASCADE");

Circle(1) = {0, -0, 0, 3.0, 2*Pi};

Line Loop(2) = {1};
//+
Plane Surface(1) = {2};
Transfinite Line {1} = nbele Using Progression 1;
//+
Physical Line    (1)  = {1};
Physical Surface(121) = {1};

