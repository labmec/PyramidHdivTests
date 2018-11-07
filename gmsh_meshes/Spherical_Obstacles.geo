
// ---- Gmsh Script ----
// ---- Channel with a spherical obstacle wih hybrid 3D mesh  ----
// Created 01/02/2018 by Omar Durán
// Labmec, University of Campinas
// --------------------------------------------




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

n_radial = 8;
n_azimuthal = 8; 
n_vertical = 1;
n_wellbore = 5;
radial_progression = 1.25;

////////////////////////////////////////////////////////////////
// Mesh Type
///////////////////////////////////////////////////////////////
// // mesh_type = 0; // Tetrahedra dominated
// // mesh_type = 1; // Hexahedra dominated
// // mesh_type = 2; // Prism dominated
// // mesh_type = 3; // Hybrid {Pyramids,Hexahdra,Tetrahedra}

mesh_type = 3; 

n_sphere = 10;
n_unstructured = 5;

x_length = 1.0;
y_length = 1.0;
z_length = 1.0;

o_x_length = 2.0;
o_y_length = 2.0;
o_z_length = 2.0;

r = 0.25;
x = 0;
y = 0;
z = 0;

// Spherical hole

p1 = newp; Point(p1) = {x,  y,  z} ;
p2 = newp; Point(p2) = {x+r,y,  z} ;
p3 = newp; Point(p3) = {x,  y+r,z} ;
p4 = newp; Point(p4) = {x,  y,  z+r};
p5 = newp; Point(p5) = {x-r,y,  z} ;
p6 = newp; Point(p6) = {x,  y-r,z} ;
p7 = newp; Point(p7) = {x,  y,  z-r};

lc1 = newl; Circle(lc1) = {p2,p1,p7}; lc2 = newl; Circle(lc2) = {p7,p1,p5};
lc3 = newl; Circle(lc3) = {p5,p1,p4}; lc4 = newl; Circle(lc4) = {p4,p1,p2};
lc5 = newl; Circle(lc5) = {p2,p1,p3}; lc6 = newl; Circle(lc6) = {p3,p1,p5};
lc7 = newl; Circle(lc7) = {p5,p1,p6}; lc8 = newl; Circle(lc8) = {p6,p1,p2};
lc9 = newl; Circle(lc9) = {p7,p1,p3}; lc10 = newl; Circle(lc10) = {p3,p1,p4};
lc11 = newl; Circle(lc11) = {p4,p1,p6}; lc12 = newl; Circle(lc12) = {p6,p1,p7};

// We need non-plane surfaces to define the spherical holes. Here we use ruled
// surfaces, which can have 3 or 4 sides:

ll1 = newll; Line Loop(ll1) = {lc5,lc10,lc4};    s1 = news; Surface(s1) = {ll1};
ll2 = newll; Line Loop(ll2) = {lc9,-lc5,lc1};    s2 = news; Surface(s2) = {ll2};
ll3 = newll; Line Loop(ll3) = {lc12,-lc8,-lc1};  s3 = news; Surface(s3) = {ll3};
ll4 = newll; Line Loop(ll4) = {lc8,-lc4,lc11};   s4 = news; Surface(s4) = {ll4};
ll5 = newll; Line Loop(ll5) = {-lc10,lc6,lc3};   s5 = news; Surface(s5) = {ll5};
ll6 = newll; Line Loop(ll6) = {-lc11,-lc3,lc7};  s6 = news; Surface(s6) = {ll6};
ll7 = newll; Line Loop(ll7) = {-lc2,-lc7,-lc12}; s7 = news; Surface(s7) = {ll7};
ll8 = newll; Line Loop(ll8) = {-lc6,-lc9,lc2};   s8 = news; Surface(s8) = {ll8};

spherical_ribs[] = {lc1,lc2,lc3,lc4,lc5,lc6,lc7,lc8,lc9,lc10,lc11,lc12};
spherical_hole[] = {s1,s2,s3,s4,s5,s6,s7,s8};

// Innner volume

ip1 = newp; Point(ip1) = {-x_length/2.0, -y_length/2.0, -z_length/2.0};
ip2 = newp; Point(ip2) = { x_length/2.0, -y_length/2.0, -z_length/2.0};
ip3 = newp; Point(ip3) = { x_length/2.0,  y_length/2.0, -z_length/2.0};
ip4 = newp; Point(ip4) = {-x_length/2.0,  y_length/2.0, -z_length/2.0};

ip5 = newp; Point(ip5) = {-x_length/2.0, -y_length/2.0, z_length/2.0};
ip6 = newp; Point(ip6) = { x_length/2.0, -y_length/2.0, z_length/2.0};
ip7 = newp; Point(ip7) = { x_length/2.0,  y_length/2.0, z_length/2.0};
ip8 = newp; Point(ip8) = {-x_length/2.0,  y_length/2.0, z_length/2.0};

il1 = newl; Line(il1) = {ip1,ip2};
il2 = newl; Line(il2) = {ip2,ip3};
il3 = newl; Line(il3) = {ip3,ip4};
il4 = newl; Line(il4) = {ip4,ip1};

il5 = newl; Line(il5) = {ip5,ip6};
il6 = newl; Line(il6) = {ip6,ip7};
il7 = newl; Line(il7) = {ip7,ip8};
il8 = newl; Line(il8) = {ip8,ip5};

il9  = newl; Line(il9)  = {ip5,ip1};
il10 = newl; Line(il10) = {ip6,ip2};
il11 = newl; Line(il11) = {ip7,ip3};
il12 = newl; Line(il12) = {ip8,ip4};

ill1  = newll; Line Loop(ill1) = {il1,  il2,   il3, il4}; // Bottom
ill2  = newll; Line Loop(ill2) = {il5,  il6,   il7, il8}; // Top
ill3  = newll; Line Loop(ill3) = {il1, -il10, -il5, il9}; // South
ill4  = newll; Line Loop(ill4) = {il2, -il11, -il6, il10}; // East
ill5  = newll; Line Loop(ill5) = {il3, -il12, -il7, il11}; // North
ill6  = newll; Line Loop(ill6) = {il4, -il9,  -il8, il12}; // West

is1  = news; Plane Surface(is1) = {ill1}; // Bottom unstructured region
is2  = news; Plane Surface(is2) = {ill2}; // Top unstructured region
is3  = news; Plane Surface(is3) = {ill3}; // South unstructured region
is4  = news; Plane Surface(is4) = {ill4}; // East unstructured region
is5  = news; Plane Surface(is5) = {ill5}; // North unstructured region
is6  = news; Plane Surface(is6) = {ill6}; // West unstructured region

ibox_B[] = {is1};
ibox_T[] = {is2};
ibox_S[] = {is3};
ibox_E[] = {is4};
ibox_N[] = {is5};
ibox_W[] = {is6};

ibox_ribs[] = {il1,il2,il3,il4,il5,il6,il7,il8,il9,il10,il11,il12};
ibox_boundaries[] = {ibox_B[],ibox_T[],ibox_S[],ibox_E[],ibox_N[],ibox_W[]};


p1 = newp; Point(p1) = {-o_x_length/2.0, -o_y_length/2.0, -o_z_length/2.0};
p2 = newp; Point(p2) = { o_x_length/2.0, -o_y_length/2.0, -o_z_length/2.0};
p3 = newp; Point(p3) = { o_x_length/2.0,  o_y_length/2.0, -o_z_length/2.0};
p4 = newp; Point(p4) = {-o_x_length/2.0,  o_y_length/2.0, -o_z_length/2.0};

p5 = newp; Point(p5) = {-o_x_length/2.0, -o_y_length/2.0, o_z_length/2.0};
p6 = newp; Point(p6) = { o_x_length/2.0, -o_y_length/2.0, o_z_length/2.0};
p7 = newp; Point(p7) = { o_x_length/2.0,  o_y_length/2.0, o_z_length/2.0};
p8 = newp; Point(p8) = {-o_x_length/2.0,  o_y_length/2.0, o_z_length/2.0};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

l5 = newl; Line(l5) = {p5,p6};
l6 = newl; Line(l6) = {p6,p7};
l7 = newl; Line(l7) = {p7,p8};
l8 = newl; Line(l8) = {p8,p5};

l9  = newl; Line(l9)  = {p5,p1};
l10 = newl; Line(l10) = {p6,p2};
l11 = newl; Line(l11) = {p7,p3};
l12 = newl; Line(l12) = {p8,p4};

ll1  = newll; Line Loop(ll1) = {l1,  l2,   l3, l4}; // Bottom
ll2  = newll; Line Loop(ll2) = {l5,  l6,   l7, l8}; // Top
ll3  = newll; Line Loop(ll3) = {l1, -l10, -l5, l9}; // South
ll4  = newll; Line Loop(ll4) = {l2, -l11, -l6, l10}; // East
ll5  = newll; Line Loop(ll5) = {l3, -l12, -l7, l11}; // North
ll6  = newll; Line Loop(ll6) = {l4, -l9,  -l8, l12}; // West

s1  = news; Plane Surface(s1) = {ll1}; // Bottom unstructured region
s2  = news; Plane Surface(s2) = {ll2}; // Top unstructured region
s3  = news; Plane Surface(s3) = {ll3}; // South unstructured region
s4  = news; Plane Surface(s4) = {ll4}; // East unstructured region
s5  = news; Plane Surface(s5) = {ll5}; // North unstructured region
s6  = news; Plane Surface(s6) = {ll6}; // West unstructured region

obox_B[] = {s1};
obox_T[] = {s2};
obox_S[] = {s3};
obox_E[] = {s4};
obox_N[] = {s5};
obox_W[] = {s6};
obox_boundaries[] = {obox_B[],obox_T[],obox_S[],obox_E[],obox_N[],obox_W[]};


// Creating another entities

rl1  = newl; Line(rl1)  = {ip1,p1};
rl2  = newl; Line(rl2)  = {ip2,p2};
rl3  = newl; Line(rl3)  = {ip3,p3};
rl4  = newl; Line(rl4)  = {ip4,p4};

rl5  = newl; Line(rl5)  = {ip5,p5};
rl6  = newl; Line(rl6)  = {ip6,p6};
rl7  = newl; Line(rl7)  = {ip7,p7};
rl8  = newl; Line(rl8)  = {ip8,p8};

//radial_lines[] = {rl1,rl2,rl3,rl4,rl5,rl6,rl7,rl8};

rll1  = newll; Line Loop(rll1) = {il1, rl2, -l1,  -rl1};
rll2  = newll; Line Loop(rll2) = {il2, rl3, -l2,  -rl2};
rll3  = newll; Line Loop(rll3) = {il3, rl4, -l3,  -rl3};
rll4  = newll; Line Loop(rll4) = {il4, rl1, -l4,  -rl4};

rll5  = newll; Line Loop(rll5) = {il5, rl6, -l5,  -rl5};
rll6  = newll; Line Loop(rll6) = {il6, rl7, -l6,  -rl6};
rll7  = newll; Line Loop(rll7) = {il7, rl8, -l7,  -rl7};
rll8  = newll; Line Loop(rll8) = {il8, rl5, -l8,  -rl8};

rll9  = newll; Line Loop(rll9) = {-rl1, -il9, rl5, l9};
rll10  = newll; Line Loop(rll10) = {-rl2, -il10, rl6, l10};
rll11  = newll; Line Loop(rll11) = {-rl3, -il11, rl7, l11};
rll12  = newll; Line Loop(rll12) = {-rl4, -il12, rl8, l12};

rs1  = news; Plane Surface(rll1) = {rll1};
rs2  = news; Plane Surface(rll2) = {rll2};
rs3  = news; Plane Surface(rll3) = {rll3};
rs4  = news; Plane Surface(rll4) = {rll4};

rs5  = news; Plane Surface(rll5) = {rll5};
rs6  = news; Plane Surface(rll6) = {rll6};
rs7  = news; Plane Surface(rll7) = {rll7};
rs8  = news; Plane Surface(rll8) = {rll8};

rs9  = news; Plane Surface(rll9) = {rll9};
rs10  = news; Plane Surface(rs10) = {rll10};
rs11  = news; Plane Surface(rs11) = {rll11};
rs12  = news; Plane Surface(rs12) = {rll12};

obox1_boundaries[] = {rs1,rs9,rs5,rs10};
obox2_boundaries[] = {rs2,rs10,rs6,rs11};
obox3_boundaries[] = {rs3,rs11,rs7,rs12};
obox4_boundaries[] = {rs4,rs12,rs8,rs9};

Transfinite Line {spherical_ribs[]} = n_sphere;
Transfinite Line {ibox_ribs[]} = n_unstructured;
Transfinite Surface {ibox_boundaries[],obox_boundaries[]};

isl1 = newsl; Surface Loop(isl1) = {ibox_boundaries[],spherical_hole[]};
iv1  = newv; Volume(iv1) = {isl1} ;
unstructured[] = {iv1};


sl1 = newsl; Surface Loop(sl1) = {obox1_boundaries[]};
sl2 = newsl; Surface Loop(sl2) = {obox2_boundaries[]};
sl3 = newsl; Surface Loop(sl3) = {obox3_boundaries[]};
sl4 = newsl; Surface Loop(sl4) = {obox4_boundaries[]};

v1  = newv; Volume(v1) = {sl1};
v2  = newv; Volume(v2) = {sl2};
v3  = newv; Volume(v3) = {sl3};
v4  = newv; Volume(v4) = {sl4};


structured[] = {v1,v2,v3,v4};


If(mesh_type == 3)

// Meshing directives for surfaces
Recombine Surface{obox_boundaries[]};


// Meshing directives for volumes
Transfinite Volume{structured[]}; // Regular partition for the reservoir region
TransfQuadTri {structured[],unstructured[]}; // Directive to force the pyramids between volumes : 1 (quads) and 5 (tri)

// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
Mesh.Algorithm3D = 4;

EndIf


Physical Volume("domain") = {unstructured[],structured[]};
Physical Surface("outer_bc") = {obox_boundaries[]};
Physical Surface("inner_bc") = {spherical_hole[]};
Physical Surface("obox1_boundaries") = {rs1,rs5};

//Physical Surface("non_flux_bc") = {top_bottom_wellbore_region_bc[], top_bottom_reservoir_bc[]};



