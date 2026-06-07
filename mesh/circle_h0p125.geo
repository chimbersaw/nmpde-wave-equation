SetFactory("OpenCASCADE");

h = 0.125;
Point(1) = {2.5, 2.5, 0, h};
Point(2) = {5, 2.5, 0, h};
Point(3) = {2.5, 5, 0, h};
Point(4) = {0, 2.5, 0, h};
Point(5) = {2.5, 0, 0, h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};
Physical Curve(1) = {1, 2, 3, 4};
