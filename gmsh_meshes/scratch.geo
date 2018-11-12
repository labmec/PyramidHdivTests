Mesh.Algorithm3D=4; //1-Delaunay, 4-Frontal, 5-Frontal Delaunay 6-Frontal hex, 7-MMG3D, 9-RTree

lc = 2;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Merge "InitialMeshFalcon.stl";
//Merge "SurfaceMeshUniform.stl";
//Merge "SurfaceMeshAnisoCurvature.stl";

CreateTopology;

z0 = 0.272342;
zmax = 11;
xmin = -20;
xmax = 3;
ymin = -3.0;
ymax = 6.0;

Point(20001) = {xmin, ymin, z0};
Point(20002) = {xmax, ymin, z0};
Point(20003) = {xmax, ymax, z0};
Point(20004) = {xmin, ymax, z0};
Point(20005) = {xmin, ymin, zmax};
Point(20006) = {xmax, ymin, zmax};
Point(20007) = {xmax, ymax, zmax};
Point(20008) = {xmin, ymax, zmax};

Line(10) = {20001, 20002};
Line(11) = {20002, 20003};
Line(12) = {20003, 20004};
Line(13) = {20004, 20001};
Line(14) = {20005, 20006};
Line(15) = {20006, 20007};
Line(16) = {20007, 20008};
Line(17) = {20008, 20005};
Line(18) = {20001, 20005};
Line(19) = {20002, 20006};
Line(20) = {20003, 20007};
Line(21) = {20004, 20008};

Line Loop(22) = {10,19,-14,-18};
Line Loop(23) = {11,20,-15,-19};
Line Loop(24) = {12,21,-16,-20};
Line Loop(25) = {13,18,-17,-21};
Line Loop(26) = {14,15,16,17};
Line Loop(27) = {1,10,11,12,13};

Plane Surface(2) = {22};
Plane Surface(3) = {23};
Plane Surface(4) = {24};
Plane Surface(5) = {25};
Plane Surface(6) = {26};
Plane Surface(7) = {27};

Surface Loop(8) = {1,2,3,4,5,6,7};
Volume(1) = {8};

Physical Surface(1) = {1};
Physical Surface(2) = {5, 4, 3, 2, 6};
Physical Surface(3) = {7};
Physical Volume(4) = {1};
