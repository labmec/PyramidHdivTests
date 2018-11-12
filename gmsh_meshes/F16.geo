
SetFactory("OpenCASCADE");

//Mesh.RemeshParametrization=7; //0=harmonic_circle, 1=conformal, 2=rbf, 3=harmonic_plane, 4=convex_circle, //5=convex_plane, 6=harmonic square
//Mesh.RemeshAlgorithm=1; //(0) nosplit (1) automatic (2) split only with metis
Mesh.Algorithm=4; //(1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)


csm = 0.1;
cs = 5.0;

Mesh.CharacteristicLengthMin = csm;
Mesh.CharacteristicLengthMax = cs;
Mesh.LcIntegrationPrecision = 1.e-5;
Mesh.MinimumCirclePoints = 50;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthFromPoints = 0;


Merge "f16.STEP";
CreateTopology;

F16p[] = Point "*";

F16l[] = Line "*";
//Compound Line {F16l[]};

F16s[] = Surface "*";
//Compound Surface {F16s[]};

F16v[] = Volume "*";



// Mid plane

y_dim[] = {-40.0,190};
z_dim[] = {-710,10.0};

mp1 = newp; Point(mp1) = {0.0, y_dim[0], z_dim[0]};
mp2 = newp; Point(mp2) = {0.0, y_dim[1], z_dim[0]};
mp3 = newp; Point(mp3) = {0.0, y_dim[1], z_dim[1]};
mp4 = newp; Point(mp4) = {0.0, y_dim[0], z_dim[1]};

ml1 = newl; Line(ml1) = {mp1,mp2};
ml2 = newl; Line(ml2) = {mp2,mp3};
ml3 = newl; Line(ml3) = {mp3,mp4};
ml4 = newl; Line(ml4) = {mp4,mp1};

mll  = newll; Line Loop(mll) = {ml1,  ml2, ml3, ml4};
mplane  = news; Plane Surface(mplane) = {mll}; // Bottom unstructured region

mplane_ribs[] = {ml1,  ml2, ml3, ml4};
mid_plane[] = {mplane};

VolumeFragments[] = BooleanFragments { Volume{F16v[]};}{ Surface{mid_plane[]}; };
Recursive Delete{ Volume{2};}
Delete{ Volume{1};}

surfaces_to_deletion[] = {2,15,16,17,25,26,27,31,32,37,40,45,46,47,48,49,50,55,63,64,65,66,67,68,69,70,71,77,78,80};
For i In {0 : #surfaces_to_deletion[]-1}
  Recursive Delete{ Surface{surfaces_to_deletion[i]};}
EndFor

surfaces_to_divide[] = {1,3,6,7,8,9,10,11,12,13,20,22,33,34,35,36,39,43,44,56,57,81,82,83};
SurfaceFragments[] = BooleanFragments { Surface{surfaces_to_divide[]}; }{  Surface{mid_plane[]};};


For i In {0 : #surfaces_to_divide[]-1}
  Recursive Delete{ Surface{surfaces_to_divide[i]};}
EndFor

Recursive Delete {
  Surface{348}; Surface{352}; Surface{364}; Surface{347}; Surface{350}; Surface{355}; Surface{367}; Surface{333}; Surface{359}; Surface{356}; Surface{334}; Surface{336}; Surface{340}; Surface{329}; Surface{327}; Surface{368}; Surface{361}; Surface{362}; Surface{338}; Surface{345}; Surface{330}; Surface{343}; 
}

//Cleaning mid planes
Delete{ Surface{mid_plane[]};}

F16lines[] = Line "*";
Compound Line {F16lines[]};

F16surf[] = Surface "*";
Compound Surface {F16surf[]};

HalfofF16s[] =  Surface "*";
Physical Surface("F16") = {HalfofF16s[]};


