// Gmsh project created on Thu Nov  5 16:26:34 2020
//+
nbele = 10;
nbptx = nbele+1;
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 3, 0, 2*Pi};
//+
Line Loop(1) = {1};
//+
Plane Surface(21) = {1};
//+
//+
//+
Transfinite Line {1} = nbptx Using Progression 1;

Physical Line    (1)  = {1};
Physical Surface (121)  = {21} ;
