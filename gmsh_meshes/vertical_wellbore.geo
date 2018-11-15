
// ---- Gmsh Script ----
// ---- Vertical wellbore wih hybrid 3D mesh  ----
// Created 01/02/2018 by Omar Durán
// Labmec, University of Campinas
// --------------------------------------------

Include "drill_wellbore.geo";


////////////////////////////////////////////////////////////////
// Geometry Parameters
////////////////////////////////////////////////////////////////

outer_r = 5.0; // reservoir radius
inner_r = 0.1; // wellbore radius
h = 0.1; // vertical lenght
s = 5.0; // amplification factor for the wellbore box

////////////////////////////////////////////////////////////////
// Mesh Parameters
///////////////////////////////////////////////////////////////

n = 0;

n_radial = 6 + n;
n_azimuthal = 4 + n; 
n_vertical = 1;
n_wellbore = 4 + n;
n_v_wellbore = 1;
radial_progression = 1.75;

////////////////////////////////////////////////////////////////
// Mesh Type
///////////////////////////////////////////////////////////////
// // mesh_type = 0; // Tetrahedra dominated
// // mesh_type = 1; // Hexahedra dominated
// // mesh_type = 2; // Prism dominated
// // mesh_type = 3; // Hybrid {Pyramids,Hexahdra,Tetrahedra}

mesh_type = 1; 

Call MakeVerticalWellbore;

Coherence;
Geometry.Tolerance=1e-05;
Coherence Mesh;

// optimize the mesh
//Mesh  3;
//Mesh.Optimize=1;
//Mesh.OptimizeNetgen=1;


//If(mesh_type == 0)
//Save "./vertical_wellbore_Te.msh";
//EndIf

//If(mesh_type == 1)
//Save "./vertical_wellbore_He.msh";
//EndIf

//If(mesh_type == 2)
//Save "./vertical_wellbore_Pe.msh";
//EndIf

//If(mesh_type == 3)
//Save "./vertical_wellbore_hybrid.msh";
//EndIf


