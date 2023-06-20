
/*40 : F*/
nbele = 200;
nbptx = nbele+1;
nbpty = 51;
ta = 2/nbele;
/* Point      1 */
Point(newp) = {-0.038, 0, 0,ta};
/* Point      2 */
Point(newp) = {0,0,0.0,ta};
/* Point      3 */
Point(newp) = {0.00045, 0, 0,ta};
/* Point      4 */
Point(newp) = {0.00045+0.038, 0, 0,ta};
/* Point      5 */
Point(newp) = {0.00045+0.038,0.01,0.0,ta};
/* Point      6 */
Point(newp) = {0.00045,0.01,0.0,ta};
/* Point      7 */
Point(newp) = {0,0.01,0.0,ta};
/* Point      8 */
Point(newp) = {-0.038,0.01,0.0,ta};



Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {2,7};
Line(10) = {3,6};


Line Loop(1) = {8,1,9,7};
Plane Surface(21) = {1};

Line Loop(2) = {-9,2,10,6};
Plane Surface(22) = {2};

Line Loop(3) = {-10,3,4,5};
Plane Surface(23) = {3};


//+
Physical Line(1) = {8};
//+
Physical Line(2) = {1,2,3};
//+
Physical Line(3) = {4};
//+
Physical Line(4) = {5,6,7};


//+
Physical Surface(121) = {21};
//+
Physical Surface(122) = {22};
//+
Physical Surface(123) = {23};

Transfinite Line {1,7} = (nbele+1) Using Progression 1;
Transfinite Line {2,6} = (nbele/10) Using Progression 1;
Transfinite Line {3,5} = (nbele+1) Using Progression 1;
Transfinite Line {8,4,9,10} = (nbele/50) Using Progression 1;

Transfinite Surface {21} = {1,2,7,8};
Transfinite Surface {22} = {2,3,6,7};
Transfinite Surface {23} = {3,4,5,6};


