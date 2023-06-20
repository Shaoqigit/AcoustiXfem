/*40 : F*/
nbele = 3;
nbpt = nbele+1;
ta = 1/nbele;
/* Point      1 */
Point(newp) = {-0.5,-0.0,0.0,ta};
/* Point      2 */
Point(newp) = {0.5,-0.0,0.0,ta};
/* Point      3 */
Point(newp) = {0.5,0.5,0.0,ta};
/* Point      4 */
Point(newp) = {-0.5,0.5,0.0,ta};

Line(9)  = {1,2};
Line(10) = {2,3};
Line(11) = {4,3};
Line(12) = {1,4};
  
Line Loop(18) = {-12,9,10,-11};
/*Ruled Surface(21) = {18};*/
Plane Surface(21) = {18};
/*MeshSize {4, 1, 2, 3} = 0.2;
Transfinite Curve {11} = 10 Using Progression 1;
//+
Transfinite Curve {9} = 5 Using Progression 1;*/
//+
Transfinite Line {9,11} = nbpt Using Power 1.0; 
Transfinite Line {10,12} = nbpt/2+1 Using Power 1.0;
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

/*Physical Point   (1)  = {1} ;
Physical Point   (2)  = {2} ;
Physical Point   (3)  = {3} ;
Physical Point   (4)  = {4} ;
Physical Point   (5)  = {5} ;
Physical Point   (6)  = {6} ;
Physical Point   (7)  = {7} ;
Physical Point   (8)  = {8} ;
Physical Point   (9)  = {9} ;
Physical Point   (10)  = {10} ;*/


