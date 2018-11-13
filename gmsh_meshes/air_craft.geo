
// ---- Gmsh Script ----
// ---- Vertical wellbore wih hybrid 3D mesh  ----
// Created 01/02/2018 by Omar Durán
// Labmec, University of Campinas
// --------------------------------------------

SetFactory("OpenCASCADE");

Mesh.Algorithm=1; // 2D mesh algorithm  (1) MeshAdapt (default)(5) Delaunay (6) Frontal
Mesh.ElementOrder=1; // 1=linear elements, N (<6) = elements of higher order
Geometry.OCCSewFaces = 1;
Merge "base.igs";
Merge "benchmarks_3d_Falcon_SurfaceMeshAnisoCurvature.stl";
//Merge "Falcon_InitialMeshFalcon.msh";
CreateTopology;
Coherence;

//Falconl[] = Line "*";
//Compound Line {Falconl[]};

//Falcons[] = Surface "*";
//Compound Surface {Falcons[]};


Extrude {0, 0, 10} {
 Line{1,2,3,4};
 Layers{1};
}

lid_ribs[] = {206,208,210,211};
lidll = newsl; Line Loop(lidll) = {lid_ribs[]};
lid  = news; Plane Surface(lid) = {lidll}; // To close the extrusion

falcon_region[] = {1,2,3,4,5,213};
isl1 = newsl; Surface Loop(isl1) = {falcon_region[]};
iv1  = newv; Volume(iv1) = {isl1};
unstructured[] = {iv1};


