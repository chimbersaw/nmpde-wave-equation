SetFactory("OpenCASCADE");

h = 0.0625;
Point(1) = {0, 0, 0, h};
Point(2) = {5, 0, 0, h};
Point(3) = {5, 5, 0, h};
Point(4) = {0, 5, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};
Physical Curve(1) = {1, 2, 3, 4};
