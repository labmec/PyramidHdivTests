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
enum EApproxSpace {ETetrahedra, EPyramid, EDividedPyramid, EDividedPyramidIncreasedOrder, EDividedPyramid4, EDividedPyramidIncreasedOrder4};

/// Enumerate that defines the type of geometry type
enum EGeometryType {EAcademic, EVerticalWellbore};


// Class that defines all the simulation controls
class TSimulationControl {
    
public:
    
    /// Type of approximation space
    EApproxSpace m_run_type;
    
    /// Type of geometry description
    EGeometryType m_geometry_type;
    
    /// Number of h refinements for H convergence
    int m_h_levels;
    
    /// Number of elements in x, y and z
    int m_n_elements;
    
    /// Directive to create a red black mixed mesh (red pyramids, black hexahedra)
    bool m_red_black_stride_Q;
    
    /// Polynomial approximation order
    int m_approx_order;
    
    /// Directive for the use of augmented pressure accuracy
    bool m_Hdiv_plusplus_Q;
    
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
    
    void Print();
    
};

#endif /* TSimulationControl_h */
