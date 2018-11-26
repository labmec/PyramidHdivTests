//
//  TSimulationControl.h
//  PyramidHdivTests
//
//  Created by Omar Durán on 10/19/18.
//

#ifndef TSimulationControl_h
#define TSimulationControl_h

#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>


/// Enumerate that defines the type of approximation space
enum EApproxSpace {EPyramid, EDividedPyramid, EDividedPyramidIncreasedOrder, EDividedPyramid4, EDividedPyramidIncreasedOrder4, ETetrahedra, EHexaHedra};

/// Enumerate that defines the type of geometry type
enum EGeometryType {EAcademic, EVerticalWellHePyTe, EVerticalWellTe, EVerticalWellHe, ESphericalBarrierHePyTe, ESphericalBarrierTe, ESphericalBarrierHe};

extern EGeometryType gCurrentRun;

// Class that defines all the simulation controls
class TSimulationControl {
    
public:
    
    /// Type of approximation space
    EApproxSpace m_run_type = EDividedPyramidIncreasedOrder4;
    
    /// Type of geometry description
    EGeometryType m_geometry_type = EAcademic;
    
    /// Number of h refinements for H convergence
    int m_h_levels = 3;
    
    /// Directive to create a red black mixed mesh (red pyramids, black hexahedra)
    bool m_red_black_stride_Q = true;
    
    /// Maximum valure of polynomial approximation order
    int m_p_levels = 1;
    
    /// Directive for the use of augmented pressure accuracy
    bool m_Hdiv_plusplus_Q = false;
    
    /// Directive for the generation of vtk file
    bool m_draw_vtk_Q = false;
    
    /// Directive to compute only the size of the system of equations
    bool m_dry_run = false;
    
    /// Directive to run the hybrid formulation or not
    bool m_hybrid = false;
    
public:
    
    /// Default constructor
    TSimulationControl();
    
    /// Destructor
    ~TSimulationControl();
    
    /// Constructor based on char *argv[]
    TSimulationControl(char *argv[]);
    
    /// Copy Constructor
    TSimulationControl(const TSimulationControl &other);
    
    /// Copy Constructor
    TSimulationControl & operator=(const TSimulationControl &other);
    
    /// Print object attributes
    void Print(std::ostream &out = std::cout);
    
    /// Convert the run type into a string
    std::string RunType()
    {
        switch(m_run_type)
        {
            case EPyramid:
                return "Pyramid";
            case EDividedPyramid:
                return "DividedPyramid2Tetrahedra";
            case EDividedPyramid4:
                return "DividedPyramid4Tetrahedra";
            case EDividedPyramidIncreasedOrder:
                return "DividedPyramid2TetrahedraIncreasedOrder";
            case EDividedPyramidIncreasedOrder4:
                return "DividedPyramid4TetrahedraIncreasedOrder";
            case ETetrahedra:
                return "Tetrahedra";
            case EHexaHedra:
                return "HexaHedra";
            default:
                return "UNKNOWN RUN TYPE EXPECT TROUBLE";
        }
    }
    
    /// Convert the geometry of the mesh into a string
    std::string MeshGeometry()
    {
        switch(m_geometry_type)
        {
            case EAcademic:
                return "AcademicMesh";
            case EVerticalWellHe:
                return "VerticalWellHexahedraMesh";
            case EVerticalWellHePyTe:
                return "VerticalWellMixedMesh";
            case EVerticalWellTe:
                return "VerticalWellTetrahedraMesh";
            case ESphericalBarrierHe:
                return "SphereHexahedraMesh";
            case ESphericalBarrierTe:
                return "SphereTetrahedraMesh";
            case ESphericalBarrierHePyTe:
                return "SphereMixedMesh";
            default:
                return "UNKNOWN MESH EXPECT TROUBLE";
        }
    }
};

#endif /* TSimulationControl_h */
