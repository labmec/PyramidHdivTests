
SetFactory("OpenCASCADE");

//Mesh.RemeshParametrization=7; //0=harmonic_circle, 1=conformal, 2=rbf, 3=harmonic_plane, 4=convex_circle, //5=convex_plane, 6=harmonic square
//Mesh.RemeshAlgorithm=1; //(0) nosplit (1) automatic (2) split only with metis
Mesh.Algorithm=4; //(1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)

csm = 0.01;
cs = 0.1;

Mesh.CharacteristicLengthMin = csm;
Mesh.CharacteristicLengthMax = cs;
Mesh.LcIntegrationPrecision = 1.e-5;
Mesh.MinimumCirclePoints = 50;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthFromPoints = 0;


Merge "Boeing_737_800.IGS";
CreateTopology;

B737P[] = Point "*";

B737l[] = Line "*";
//Compound Line {ll[]};

B737s[] = Surface "*";
//Compound Surface {ss[]};

righthalfB737s[] = {};
rc=0;
For i In {5 : 47}
  righthalfB737s[rc] = B737s[i];
  rc=rc+1;
EndFor


lefthalfB737s[] = {};
lc=0;
For i In {0 : 4}
  lefthalfB737s[lc] = B737s[i];
  lc=lc+1;
EndFor
For i In {49 : 83}
  lefthalfB737s[lc] = B737s[i];
  lc=lc+1;
EndFor

For i In {0 : #righthalfB737s[]-1}
  Recursive Delete{ Surface{righthalfB737s[i]};}
EndFor


l1 = newl; Line(l1) = {3414,3415};
ll1  = newll; Line Loop(ll1) = {l1, 3424};
lid1  = news; Surface(lid1) = {ll1};

l2 = newl; Line(l2) = {3411,3412};
ll2  = newll; Line Loop(ll2) = {l2, 3421};
lid2  = news; Surface(lid2) = {ll2};


Physical Surface("BoeingL") = {lefthalfB737s[]};
Physical Surface("lid1") = {lid1};
Physical Surface("lid2") = {lid2};
Physical Surface("BoeingR") = {righthalfB737s[]};

