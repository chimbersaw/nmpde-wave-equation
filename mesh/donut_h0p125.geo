SetFactory("OpenCASCADE");

h = 0.125;
Point(1) = {2.5, 2.5, 0, h};
Point(2) = {5, 2.5, 0, h};
Point(3) = {2.5, 5, 0, h};
Point(4) = {0, 2.5, 0, h};
Point(5) = {2.5, 0, 0, h};
Point(6) = {3.4, 2.5, 0, h};
Point(7) = {2.5, 3.4, 0, h};
Point(8) = {1.6, 2.5, 0, h};
Point(9) = {2.5, 1.6, 0, h};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};
Physical Surface(1) = {1};
Physical Curve(1) = {1, 2, 3, 4, 5, 6, 7, 8};
