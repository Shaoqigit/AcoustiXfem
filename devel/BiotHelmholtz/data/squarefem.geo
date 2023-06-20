/*40 : F*/
nbele = 259;
nbptx = nbele+1;
nbpty = nbptx;
ta = 2/nbele;
/* Point      1 */
Point(newp) = {-0.2,-0.2,0.0,ta};
/* Point      2 */
Point(newp) = {0,-0.2,0.0,ta};
/* Point      3 */
Point(newp) = {0.2,-0.2,0.0,ta};
/* Point      4 */
Point(newp) = {0.2,0.2,0.0,ta};
/* Point      5 */
Point(newp) = {0,0.2,0.0,ta};
/* Point      6 */
Point(newp) = {-0.2,0.2,0.0,ta};

Line(9)  = {1,2};
Line(10) = {2,3};
Line(11) = {3,4};
Line(12) = {4,5};
Line(13) = {5,6};
Line(14) = {6,1};
Line(15) = {2,5};
  
Line Loop(18) = {9,15,13,14};
Line Loop(19) = {10,11,12,-15};
/*Ruled Surface(21) = {18};*/
Plane Surface(21) = {18};
Plane Surface(22) = {19};

Transfinite Line {9,10, 13, 12} = nbptx/2+1 Using Power 1.0; 
Transfinite Line {14,15,11} = nbpty Using Power 1.0; 
Transfinite Surface {21} = {1,2,5,6};
Transfinite Surface {22} = {2,3,4,5};
/*Recombine Surface {21} = 90;*/

Physical Point   (1)  = {1} ;
Physical Point   (2)  = {2} ;
Physical Point   (3)  = {3} ;
Physical Point   (4)  = {4} ;

Physical Line    (1)  = {14};
Physical Line    (2)  = {15};
Physical Line    (3)  = {11};
Physical Line    (4)  = {9};
Physical Line    (5)  = {10};
Physical Line    (6)  = {13};
Physical Line    (7)  = {12};

Physical Surface (121)  = {21} ;
Physical Surface (122)  = {22} ;










