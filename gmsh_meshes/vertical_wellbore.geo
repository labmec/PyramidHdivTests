
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
s = 4.0; // amplification factor for the wellbore box

////////////////////////////////////////////////////////////////
// Mesh Parameters
///////////////////////////////////////////////////////////////

n_radial = 6;
n_azimuthal = 6; 
n_vertical = 1;
n_wellbore = 4;
n_v_wellbore = 1;
radial_progression = 1.75;

////////////////////////////////////////////////////////////////
// Mesh Type
///////////////////////////////////////////////////////////////
// // mesh_type = 0; // Tetrahedra dominated
// // mesh_type = 1; // Hexahedra dominated
// // mesh_type = 2; // Prism dominated
// // mesh_type = 3; // Hybrid {Pyramids,Hexahdra,Tetrahedra}

mesh_type = 0; 

Call MakeVerticalWellbore;

