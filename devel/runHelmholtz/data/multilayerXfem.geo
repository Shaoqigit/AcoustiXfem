
/*40 : F*/
nbele = 200;
nbptx = nbele+1;
nbpty = 51;
ta = 2/nbele;
/* Point      1 */
Point(newp) = {-0.038-0.000225, 0, 0,ta};
/* Point      2 */
Point(newp) = {0.000225+0.038, 0, 0,ta};
/* Point      3 */
Point(newp) = {0.000225+0.038,0.01,0.0,ta};
/* Point      4 */
Point(newp) = {-0.038-0.000225,0.01,0.0,ta};



Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};



Line Loop(1) = {4,1,2,3};
Plane Surface(21) = {1};


//+
Physical Line(1) = {4};
//+
Physical Line(2) = {1};
//+
Physical Line(3) = {2};
//+
Physical Line(4) = {3};


//+
Physical Surface(121) = {21};
//+


Transfinite Line {1,3} = (nbele+1) Using Progression 1;
Transfinite Line {2,4} = (nbele/50+1) Using Progression 1;
/*Transfinite Line {6,7,3} = (nbele/50) Using Progression 1;*/

Transfinite Surface {21} = {1,2,3,4};




