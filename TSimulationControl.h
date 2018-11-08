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
enum EGeometryType {EAcademic, EVerticalWellHePyTe, EVerticalWellTe, EVerticalWellHe, ESphericalBarrierHePyTe, ESphericalBarrierTe};


// Class that defines all the simulation controls
class TSimulationControl {
    
public:
    
    /// Type of approximation space
    EApproxSpace m_run_type;
    
    /// Type of geometry description
    EGeometryType m_geometry_type;
    
    /// Number of h refinements for H convergence
    int m_h_levels;
    
    /// Directive to create a red black mixed mesh (red pyramids, black hexahedra)
    bool m_red_black_stride_Q;
    
    /// Maximum valure of polynomial approximation order
    int m_p_levels;
    
    /// Directive for the use of augmented pressure accuracy
    bool m_Hdiv_plusplus_Q;
    
    /// Directive for the generation of vtk file
    bool m_draw_vtk_Q;
    
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
    
};

#endif /* TSimulationControl_h */
