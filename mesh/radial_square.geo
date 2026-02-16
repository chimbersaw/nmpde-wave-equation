lc = 0.2;
n = 81;

Point(1) = {-5.0, -5.0, 0.0, lc};
Point(2) = {5.0, -5.0, 0.0, lc};
Point(3) = {5.0, 5.0, 0.0, lc};
Point(4) = {-5.0, 5.0, 0.0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve {1, 2, 3, 4} = n;
Transfinite Surface {1};
Recombine Surface {1};

Physical Curve("boundary") = {1, 2, 3, 4};
Physical Surface("domain") = {1};

Mesh.MshFileVersion = 2.2;
