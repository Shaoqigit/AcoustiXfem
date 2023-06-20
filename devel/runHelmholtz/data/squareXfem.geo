/*40 : F*/
nbele = 2;
nbptx = nbele+1;
nbpty = 51;
ta = 2/nbele;
/* Point      1 */
Point(newp) = {-0.1,0.0,0.0,ta};
/* Point      2 */
Point(newp) = {0.1,0.0,0.0,ta};
/* Point      3 */
Point(newp) = {0.1,0.2,0.0,ta};
/* Point      4 */
Point(newp) = {-0.1,0.2,0.0,ta};

//+
Line(9)  = {1,2};
Line(10) = {2,3};
Line(11) = {4,3};
Line(12) = {1,4};
  
Line Loop(18) = {-12,9,10,-11};
/*Ruled Surface(21) = {18};*/
Plane Surface(21) = {18};

Transfinite Line {9,11} = nbptx Using Power 1.0; 
Transfinite Line {10,12} = nbptx Using Power 1.0; 
Transfinite Surface {21} = {1,2,3,4};

/*Recombine Surface {21} = 90;*/

Physical Point   (1)  = {1} ;
Physical Point   (2)  = {2} ;
Physical Point   (3)  = {3} ;
Physical Point   (4)  = {4} ;

Physical Line    (1)  = {9};
Physical Line    (2)  = {10};
Physical Line    (3)  = {11};
Physical Line    (4)  = {12};

Physical Surface (121)  = {21} ;

/*SetFactory("OpenCASCADE");
Circle(13) = {0, 0.1, 0, 0.04, 0, 2*Pi};
//+
Circle(14) = {0, 0.1, 0, 0.041, 0, 2*Pi};*/



