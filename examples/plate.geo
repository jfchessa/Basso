cl__1 = 1e+22;
cl__2 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
Physical Line(8) = {3};
Physical Line(9) = {1};
Physical Surface(7) = {6};
