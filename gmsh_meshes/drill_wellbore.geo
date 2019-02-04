
// ---- Gmsh Macro ----
// ---- Vertical wellbore wih hybrid 3D mesh  ----
// Created 01/02/2018 by Omar Durán
// Labmec, University of Campinas
// --------------------------------------------

Macro MakeVerticalWellbore

l = 1.0;
r1 = inner_r/Sqrt(2.0);
r2 = outer_r/Sqrt(2.0);

p1 = newp; Point(p1) = {0,0,-h/2,l};

// interior circle
p2 = newp; Point(p2) = {r1,r1,-h/2,l};
p3 = newp; Point(p3) = {-r1,r1,-h/2,l};
p4 = newp; Point(p4) = {-r1,-r1,-h/2,l};
p5 = newp; Point(p5) = {r1,-r1,-h/2,l};

l1 = newl; Circle(l1) = {p3,p1,p2};
l2 = newl; Circle(l2) = {p2,p1,p5};
l3 = newl; Circle(l3) = {p5,p1,p4};
l4 = newl; Circle(l4) = {p4,p1,p3};

// exterior circle
pe2 = newp; Point(pe2) = {r2,r2,-h/2,l};
pe3 = newp; Point(pe3) = {-r2,r2,-h/2,l};
pe4 = newp; Point(pe4) = {-r2,-r2,-h/2,l};
pe5 = newp; Point(pe5) = {r2,-r2,-h/2,l};
le1 = newl; Circle(le1) = {pe3,p1,pe2};
le2 = newl; Circle(le2) = {pe2,p1,pe5};
le3 = newl; Circle(le3) = {pe5,p1,pe4};
le4 = newl; Circle(le4) = {pe4,p1,pe3};

pm2 = newp; Point(pm2) = {s*r1,s*r1,-h/2,l};
pm3 = newp; Point(pm3) = {-s*r1,s*r1,-h/2,l};
pm4 = newp; Point(pm4) = {-s*r1,-s*r1,-h/2,l};
pm5 = newp; Point(pm5) = {s*r1,-s*r1,-h/2,l};

ilr1  = newl; Line(ilr1)  = {p2,pm2};
ilr2  = newl; Line(ilr2)  = {p3,pm3};
ilr3  = newl; Line(ilr3)  = {p4,pm4};
ilr4  = newl; Line(ilr4)  = {p5,pm5};

lr1  = newl; Line(lr1)  = {pm2,pe2};
lr2  = newl; Line(lr2)  = {pm3,pe3};
lr3  = newl; Line(lr3)  = {pm4,pe4};
lr4  = newl; Line(lr4)  = {pm5,pe5};

lbox1  = newl; Line(lbox1)  = {pm2,pm3};
lbox2  = newl; Line(lbox2)  = {pm3,pm4};
lbox3  = newl; Line(lbox3)  = {pm4,pm5};
lbox4  = newl; Line(lbox4)  = {pm5,pm2};

// Hard coded
Line Loop(1) = {17, -10, 1, 9};
Plane Surface(1) = {1};
Line Loop(2) = {10, 18, -11, 4};
Plane Surface(2) = {2};
Line Loop(3) = {3, 11, 19, -12};
Plane Surface(3) = {3};
Line Loop(4) = {2, 12, 20, -9};
Plane Surface(4) = {4};
Line Loop(5) = {13, -5, -14, -17};
Plane Surface(5) = {5};
Line Loop(6) = {18, 15, 8, -14};
Plane Surface(6) = {6};
Line Loop(7) = {19, 16, 7, -15};
Plane Surface(7) = {7};
Line Loop(8) = {20, 13, 6, -16};
Plane Surface(8) = {8};

base[] = {1,2,3,4,5,6,7,8};

If(mesh_type == 0)
Extrude {0, 0, h} {
  Surface{base[]};
  Layers{1}; 
}
Else
Extrude {0, 0, h} {
  Surface{base[]};
}
EndIf

// Line Groups
wellbore_v_edges[] = {71,36,32,58};
wellbore_h_edges[] = {1,2,3,4,66,24,88,47};
box_h_edges[] = {17,18,19,20,45,68,90,22};
box_v_edges[] = {54,80,27,28};
outer_h_edges[] = {5,6,7,8,134,156,178,111};
outer_v_edges[] = {142,120,116,164};
box_radial_edges[]= {9,10,11,12,-46,-69,25,-23};
radial_edges[]= {13,14,15,16,133,155,110,-112};



// Surface Groups
mid_surf_box[]= {38, 104, 82, 60};
mid_surf_rad[]= {143,165,117,125};
mid_surf[] = {mid_surf_rad[],mid_surf_box[]};
top_bottom_reservoir_bc[]={5,6,7,8,174,130,152,196};
top_bottom_wellbore_region_bc[]={1,2,3,4,42,64,86,108};
wellbore_bc[]= {73,37,63,95};
reservoir_bc[]= {169,147,191,121};

// Volume Groups
wellbore_region[]={1,2,3,4};
reservoir[]={5,6,7,8};


// Mesh types

If(mesh_type == 0)

// Meshing directives for lines
Transfinite Line {box_radial_edges[]} = n_radial_wb Using Progression radial_progression_wb; // Radial control for wellbore
Transfinite Line {radial_edges[]} = n_radial Using Progression radial_progression; // Radial control
Transfinite Line {outer_h_edges[],box_h_edges[]} = n_azimuthal; // Azimuthal control
Transfinite Line {box_v_edges[],wellbore_v_edges[]} = n_vertical; // Vertical control
//Transfinite Surface {top_bottom_reservoir_bc[],mid_surf_box[],mid_surf_rad[]};
Transfinite Surface "*";
//Transfinite Volume "*";
//Transfinite Volume {reservoir[]};

// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
Mesh.Algorithm3D = 4;

EndIf


If(mesh_type == 1)

// Meshing directives for lines
Transfinite Line {box_radial_edges[]} = n_radial_wb Using Progression radial_progression_wb; // Radial control for wellbore
Transfinite Line {radial_edges[]} = n_radial Using Progression radial_progression; // Radial control
Transfinite Line {outer_h_edges[],box_h_edges[],wellbore_h_edges[]} = n_azimuthal; // Azimuthal control
Transfinite Line {outer_v_edges[],box_v_edges[],wellbore_v_edges[]} = n_vertical; // Vertical control

// Meshing directives for surfaces
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
Recombine Volume "*";

// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
Mesh.Algorithm3D = 4;

EndIf

If(mesh_type == 2)

// Meshing directives for lines
Transfinite Line {box_radial_edges[]} = n_radial_wb Using Progression radial_progression_wb; // Radial control for wellbore
Transfinite Line {radial_edges[]} = n_radial Using Progression radial_progression; // Radial control
Transfinite Line {outer_h_edges[],box_h_edges[],wellbore_h_edges[]} = n_azimuthal; // Azimuthal control
Transfinite Line {outer_v_edges[],box_v_edges[],wellbore_v_edges[]} = n_vertical; // Vertical control

// Meshing directives for surfaces
//Transfinite Surface{mid_surf[],reservoir_bc[],wellbore_bc[]};
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface{mid_surf[],reservoir_bc[],wellbore_bc[]};


EndIf


If(mesh_type == 3)

// Meshing directives for lines
Transfinite Line {box_radial_edges[]} = n_radial_wb Using Progression radial_progression_wb; // Radial control for wellbore
Transfinite Line {radial_edges[]} = n_radial Using Progression radial_progression; // Radial control
Transfinite Line {outer_h_edges[],box_h_edges[]} = n_azimuthal; // Azimuthal control
Transfinite Line {outer_v_edges[],box_v_edges[]} = n_vertical; // Vertical control

// Meshing directives for surfaces
//Transfinite Surface{mid_surf[],reservoir_bc[],top_bottom_reservoir_bc[]};
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface{mid_surf_rad[],top_bottom_reservoir_bc[],reservoir_bc[],reservoir_bc[]}; 


// Meshing directives for volumes
Transfinite Volume{reservoir[]}; // Regular partition for the reservoir region
TransfQuadTri {5,1}; // Directive to force the pyramids between volumes : 5 (quads) and 1 (tri)
TransfQuadTri {6,2}; // Directive to force the pyramids between volumes : 6 (quads) and 2 (tri)
TransfQuadTri {7,3}; // Directive to force the pyramids between volumes : 7 (quads) and 3 (tri)
TransfQuadTri {8,4}; // Directive to force the pyramids between volumes : 8 (quads) and 4 (tri)


// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
Mesh.Algorithm3D = 4;

EndIf


Transfinite Line {wellbore_v_edges[]} = n_v_wellbore; // Wellbore control
Transfinite Line {wellbore_h_edges[]} = n_wellbore; // Wellbore control
Transfinite Surface {wellbore_bc[]};

// Tagging
Physical Volume("reservoir") = {reservoir[],wellbore_region[]};
//Physical Volume("wellbore") = {wellbore_region[]};
Physical Surface("outer_bc") = {reservoir_bc[]};
Physical Surface("inner_bc") = {wellbore_bc[]};
Physical Surface("non_flux_bc") = {top_bottom_wellbore_region_bc[], top_bottom_reservoir_bc[]};


Another_entities = 1;
If(Another_entities)
Physical Surface("non_flux_res") = {top_bottom_reservoir_bc[]};
Physical Surface("non_flux_well") = {top_bottom_wellbore_region_bc[]};
Physical Surface("mid_box") = {mid_surf_box[]};
Physical Surface("mid_rad") = {mid_surf_rad[]};
Physical Line("box_h_edges") = {box_h_edges[]};
Physical Line("box_v_edges") = {box_v_edges[]};
Physical Line("wellbore_v_edges") = {wellbore_v_edges[]};
Physical Line("wellbore_h_edges") = {wellbore_h_edges[]};
Physical Line("outer_h_edges") = {outer_h_edges[]};
Physical Line("outer_v_edges") = {outer_v_edges[]};
Physical Line("radial_edges") = {radial_edges[]};
Physical Line("box_radial_edges") = {box_radial_edges[]};
EndIf



Return
