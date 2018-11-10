
// ---- Gmsh Script ----
// ---- Channel with a sphere with Hexahedra mesh  ----
// Created 01/02/2018 by Omar Durán
// Labmec, University of Campinas
// --------------------------------------------


Geometry.Surfaces=0;
Mesh.SurfaceFaces=0;


////////////////////////////////////////////////////////////////
// Geometry Parameters
////////////////////////////////////////////////////////////////

outer_r = 5.0; // reservoir radius
inner_r = 0.1; // wellbore radius
h = 1; // vertical lenght
s = 5.0; // amplification factor for the wellbore box

////////////////////////////////////////////////////////////////
// Mesh Parameters
///////////////////////////////////////////////////////////////

n_sphere = 4;
n_structured = 3;
n_radial = 2;
radial_progression = 1.5;

x_length = 1.0;
y_length = 1.0;
z_length = 1.0;

o_x_length = 5.0;
o_y_length = 5.0;
o_z_length = 5.0;

r = 0.25;
x = 0;
y = 0;
z = 0;


////////////////////////////////////////////////////////////////
// Mesh Type
///////////////////////////////////////////////////////////////
// // mesh_type = 1; // Tetrahedra dominated
// // mesh_type = 2; // Hexahedra dominated
// // mesh_type = 3; // Hybrid {Pyramids,Hexahdra,Tetrahedra}

mesh_type = 2; 



// Spherical hole
h=0.5;
s=1.0/Sqrt(2.0);


spc = newp; Point(spc) = {x,  y,  z} ;
sp1 = newp; Point(sp1) = {x + r*h,  y + r*h,  z + r*s} ;
sp2 = newp; Point(sp2) = {x + r*h,  y - r*h,  z + r*s} ;
sp3 = newp; Point(sp3) = {x - r*h,  y - r*h,  z + r*s} ;
sp4 = newp; Point(sp4) = {x - r*h,  y + r*h,  z + r*s} ;
sp5 = newp; Point(sp5) = {x + r*h,  y + r*h,  z - r*s} ;
sp6 = newp; Point(sp6) = {x + r*h,  y - r*h,  z - r*s} ;
sp7 = newp; Point(sp7) = {x - r*h,  y - r*h,  z - r*s} ;
sp8 = newp; Point(sp8) = {x - r*h,  y + r*h,  z - r*s} ;


hl1 = newl; Circle(hl1) = {sp1,spc,sp2};
hl2 = newl; Circle(hl2) = {sp2,spc,sp3};
hl3 = newl; Circle(hl3) = {sp3,spc,sp4};
hl4 = newl; Circle(hl4) = {sp4,spc,sp1};

hl5 = newl; Circle(hl5) = {sp5,spc,sp6};
hl6 = newl; Circle(hl6) = {sp6,spc,sp7};
hl7 = newl; Circle(hl7) = {sp7,spc,sp8};
hl8 = newl; Circle(hl8) = {sp8,spc,sp5};

hl9 = newl; Circle(hl9) = {sp5,spc,sp1};
hl10 = newl; Circle(hl10) = {sp6,spc,sp2};
hl11 = newl; Circle(hl11) = {sp7,spc,sp3};
hl12 = newl; Circle(hl12) = {sp8,spc,sp4};

hll1  = newll; Line Loop(hll1) = {hl1,  hl2,   hl3, hl4}; // Bottom
hll2  = newll; Line Loop(hll2) = {hl5,  hl6,   hl7, hl8}; // Top
hll3  = newll; Line Loop(hll3) = {hl1, -hl10, -hl5, hl9}; // South
hll4  = newll; Line Loop(hll4) = {hl2, -hl11, -hl6, hl10}; // East
hll5  = newll; Line Loop(hll5) = {hl3, -hl12, -hl7, hl11}; // North
hll6  = newll; Line Loop(hll6) = {hl4, -hl9,  -hl8, hl12}; // West


hs1  = news; Ruled Surface(hs1) = {hll1}; // Bottom unstructured region
hs2  = news; Ruled Surface(hs2) = {hll2}; // Top unstructured region
hs3  = news; Ruled Surface(hs3) = {hll3}; // South unstructured region
hs4  = news; Ruled Surface(hs4) = {hll4}; // East unstructured region
hs5  = news; Ruled Surface(hs5) = {hll5}; // North unstructured region
hs6  = news; Ruled Surface(hs6) = {hll6}; // West unstructured region

h_B[] = {hs1};
h_T[] = {hs2};
h_S[] = {hs3};
h_E[] = {hs4};
h_N[] = {hs5};
h_W[] = {hs6};

spherical_ribs[] = {hl1,hl2,hl3,hl4,hl5,hl6,hl7,hl8,hl9,hl10,hl11,hl12};
spherical_hole[] = {h_B[],h_T[],h_S[],h_E[],h_N[],h_W[]};


// Common meshing controls
Transfinite Line {spherical_ribs[]} = n_sphere;
//Transfinite Line {ibox_ribs[],obox_ribs[]} = n_structured;
//Transfinite Line {radial_ribs[]} =  n_radial Using Progression radial_progression;
Transfinite Surface {spherical_hole[]};


// Meshing directives for surfaces
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
Recombine Volume "*";

//Mesh.RecombinationAlgorithm = 0;
//Mesh.Algorithm=6; // mesh algorithm
// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
//Mesh.Algorithm3D = 4;


isl1 = newsl; Surface Loop(isl1) = {spherical_hole[]};
iv1  = newv; Volume(iv1) = {isl1};
structured[] = {iv1};

Physical Volume("domain") = {structured[]};
Physical Surface("surface_bc") = {spherical_hole[]};



