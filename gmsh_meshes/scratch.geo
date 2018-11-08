
SetFactory("OpenCASCADE");

Mesh.Algorithm=7; //(1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)

Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthMin = 0.005;
Mesh.CharacteristicLengthMax = 2.0;
Mesh.LcIntegrationPrecision=1.e-5;
Mesh.MinimumCirclePoints=50;
Mesh.CharacteristicLengthExtendFromBoundary=0;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.CharacteristicLengthFromPoints = 0;

Merge "InitialMeshFalcon.msh";


CreateTopology;
ll[] = Line "*";
Compound Line {ll[]};

ss[] = Surface "*";
Compound Surface {ss[]};




//+
BooleanFragments{ Surface{1}; Surface{13}; Surface{14}; Surface{15}; Surface{16}; Delete; }{ }
