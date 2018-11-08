
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

n_sphere = 4;
n_structured = 4;
n_radial = 8;
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

mesh_type = 3; 



// Spherical hole
h=0.5;
s=1.0/Sqrt(2.0);


pc = newp; Point(pc) = {x,  y,  z} ;
p1 = newp; Point(p1) = {x + r*h,  y + r*h,  z + r*s} ;
p2 = newp; Point(p2) = {x + r*h,  y - r*h,  z + r*s} ;
p3 = newp; Point(p3) = {x - r*h,  y - r*h,  z + r*s} ;
p4 = newp; Point(p4) = {x - r*h,  y + r*h,  z + r*s} ;
p5 = newp; Point(p5) = {x + r*h,  y + r*h,  z - r*s} ;
p6 = newp; Point(p6) = {x + r*h,  y - r*h,  z - r*s} ;
p7 = newp; Point(p7) = {x - r*h,  y - r*h,  z - r*s} ;
p8 = newp; Point(p8) = {x - r*h,  y + r*h,  z - r*s} ;


hl1 = newl; Circle(hl1) = {p1,pc,p2};
hl2 = newl; Circle(hl2) = {p2,pc,p3};
hl3 = newl; Circle(hl3) = {p3,pc,p4};
hl4 = newl; Circle(hl4) = {p4,pc,p1};

hl5 = newl; Circle(hl5) = {p5,pc,p6};
hl6 = newl; Circle(hl6) = {p6,pc,p7};
hl7 = newl; Circle(hl7) = {p7,pc,p8};
hl8 = newl; Circle(hl8) = {p8,pc,p5};

hl9 = newl; Circle(hl9) = {p5,pc,p1};
hl10 = newl; Circle(hl10) = {p6,pc,p2};
hl11 = newl; Circle(hl11) = {p7,pc,p3};
hl12 = newl; Circle(hl12) = {p8,pc,p4};

hll1  = newll; Line Loop(hll1) = {hl1,  hl2,   hl3, hl4}; // Bottom
hll2  = newll; Line Loop(hll2) = {hl5,  hl6,   hl7, hl8}; // Top
hll3  = newll; Line Loop(hll3) = {hl1, -hl10, -hl5, hl9}; // South
hll4  = newll; Line Loop(hll4) = {hl2, -hl11, -hl6, hl10}; // East
hll5  = newll; Line Loop(hll5) = {hl3, -hl12, -hl7, hl11}; // North
hll6  = newll; Line Loop(hll6) = {hl4, -hl9,  -hl8, hl12}; // West


hs1  = news; Surface(hs1) = {hll1}; // Bottom unstructured region
hs2  = news; Surface(hs2) = {hll2}; // Top unstructured region
hs3  = news; Surface(hs3) = {hll3}; // South unstructured region
hs4  = news; Surface(hs4) = {hll4}; // East unstructured region
hs5  = news; Surface(hs5) = {hll5}; // North unstructured region
hs6  = news; Surface(hs6) = {hll6}; // West unstructured region

h_B[] = {hs1};
h_T[] = {hs2};
h_S[] = {hs3};
h_E[] = {hs4};
h_N[] = {hs5};
h_W[] = {hs6};

spherical_ribs[] = {hl1,hl2,hl3,hl4,hl5,hl6,hl7,hl8,hl9,hl10,hl11,hl12};
spherical_hole[] = {h_B[],h_T[],h_S[],h_E[],h_N[],h_W[]};

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

obox_ribs[] = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12};
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

rs1  = news; Plane Surface(rs1) = {rll1};
rs2  = news; Plane Surface(rs2) = {rll2};
rs3  = news; Plane Surface(rs3) = {rll3};
rs4  = news; Plane Surface(rs4) = {rll4};

rs5  = news; Plane Surface(rs5) = {rll5};
rs6  = news; Plane Surface(rs6) = {rll6};
rs7  = news; Plane Surface(rs7) = {rll7};
rs8  = news; Plane Surface(rs8) = {rll8};

rs9  = news; Plane Surface(rs9) = {rll9};
rs10  = news; Plane Surface(rs10) = {rll10};
rs11  = news; Plane Surface(rs11) = {rll11};
rs12  = news; Plane Surface(rs12) = {rll12};

obox1_boundaries[] = {rs1,rs9,rs5,rs10};
obox2_boundaries[] = {rs2,rs10,rs6,rs11};
obox3_boundaries[] = {rs3,rs11,rs7,rs12};
obox4_boundaries[] = {rs4,rs12,rs8,rs9};

obox5_boundaries[] = {rs1,rs2,rs3,rs4};
obox6_boundaries[] = {rs5,rs6,rs7,rs8};

radial_ribs[] = {rl1,rl2,rl3,rl4,rl5,rl6,rl7,rl8};
radial_planes[] = {rs1,rs2,rs3,rs4,rs5,rs6,rs7,rs8,rs9,rs10,rs11,rs12};

isl1 = newsl; Surface Loop(isl1) = {ibox_boundaries[],spherical_hole[]};
iv1  = newv; Volume(iv1) = {isl1};
unstructured[] = {iv1};

sl1 = newsl; Surface Loop(sl1) = {obox1_boundaries[],ibox_S[],obox_S[]};
sl2 = newsl; Surface Loop(sl2) = {obox2_boundaries[],ibox_E[],obox_E[]};
sl3 = newsl; Surface Loop(sl3) = {obox3_boundaries[],ibox_N[],obox_N[]};
sl4 = newsl; Surface Loop(sl4) = {obox4_boundaries[],ibox_W[],obox_W[]};
sl5 = newsl; Surface Loop(sl5) = {obox5_boundaries[],ibox_B[],obox_B[]};
sl6 = newsl; Surface Loop(sl6) = {obox6_boundaries[],ibox_T[],obox_T[]};

v1  = newv; Volume(v1) = {sl1};
v2  = newv; Volume(v2) = {sl2};
v3  = newv; Volume(v3) = {sl3};
v4  = newv; Volume(v4) = {sl4};
v5  = newv; Volume(v5) = {sl5};
v6  = newv; Volume(v6) = {sl6};


structured[] = {v1,v2,v3,v4,v5,v6};


// Common meshing controls

Transfinite Line {spherical_ribs[]} = n_sphere;
Transfinite Line {ibox_ribs[],obox_ribs[]} = n_structured;
Transfinite Line {radial_ribs[]} =  n_radial Using Progression radial_progression;
Transfinite Surface {ibox_boundaries[],obox_boundaries[],radial_planes[]};

// Mesh types

If(mesh_type == 2)

// Meshing directives for surfaces
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
Recombine Volume "*";

// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
Mesh.Algorithm3D = 4;

EndIf

If(mesh_type == 3)

// Meshing directives for surfaces
Recombine Surface{obox_boundaries[],radial_planes[]};
Recombine Volume{structured[]};

// Meshing directives for volumes
Transfinite Volume{structured[]}; // Regular partition for the reservoir region
TransfQuadTri {v1,iv1}; // Directive to force the pyramids between volumes : v1 (quads) and iv1 (tri)
TransfQuadTri {v2,iv1}; // Directive to force the pyramids between volumes : v2 (quads) and iv1 (tri)
TransfQuadTri {v3,iv1}; // Directive to force the pyramids between volumes : v3 (quads) and iv1 (tri)
TransfQuadTri {v4,iv1}; // Directive to force the pyramids between volumes : v4 (quads) and iv1 (tri)
TransfQuadTri {v5,iv1}; // Directive to force the pyramids between volumes : v5 (quads) and iv1 (tri)
TransfQuadTri {v6,iv1}; // Directive to force the pyramids between volumes : v6 (quads) and iv1 (tri)

// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
Mesh.Algorithm3D = 4;

EndIf


Physical Volume("domain") = {unstructured[],structured[]};
Physical Surface("outer_bc") = {obox_boundaries[]};
Physical Surface("inner_bc") = {spherical_hole[]};
//Physical Surface("non_flux_bc") = {top_bottom_wellbore_region_bc[], top_bottom_reservoir_bc[]};



