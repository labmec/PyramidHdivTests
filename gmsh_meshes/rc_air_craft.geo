
// ---- Gmsh Script ----
// ---- Vertical wellbore wih hybrid 3D mesh  ----
// Created 01/02/2018 by Omar Durán
// Labmec, University of Campinas
// --------------------------------------------

SetFactory("OpenCASCADE");

Mesh.Algorithm=1; // 2D mesh algorithm  (1) MeshAdapt (default)(5) Delaunay (6) Frontal
Mesh.ElementOrder=1; // 1=linear elements, N (<6) = elements of higher order
Geometry.OCCSewFaces = 1;
Merge "rc_aircraft_half_V.IGS";
CreateTopology;
Coherence;

//Line{78,81} In Surface{6};