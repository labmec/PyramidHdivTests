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
enum EApproxSpace {ETetrahedra, EPyramid,EDividedPyramid, EDividedPyramidIncreasedOrder, EDividedPyramid4, EDividedPyramidIncreasedOrder4};

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
    
    /// Stride size to define where the refine patterns are applied
    int m_cartesian_stride;
    
    /// Polynomial approximation order
    int m_approx_order;
    
    /// Directive for the use of augmented pressure accuracy
    bool m_Hdiv_plusplus_Q;
    
public:
    
    /// Default constructor
    TSimulationControl(){
        
        m_run_type          = ETetrahedra;
        m_geometry_type     = EAcademic;
        m_h_levels          = 0;
        m_n_elements        = 1;
        m_cartesian_stride  = 1;
        m_approx_order      = 1;
        m_Hdiv_plusplus_Q   = 0;
        
    }
    
    /// Destructor
    ~TSimulationControl(){
        
    }
    
    /// Constructor based on char *argv[]
    TSimulationControl(char *argv[]){
        
        m_run_type          = static_cast<EApproxSpace>(atoi(argv[1]));
        m_geometry_type     = static_cast<EGeometryType>(atoi(argv[2]));
        m_h_levels          = atoi(argv[3]);
        m_n_elements        = atoi(argv[4]);
        m_cartesian_stride  = atoi(argv[5]);
        m_approx_order      = atoi(argv[6]);
        m_Hdiv_plusplus_Q   = atoi(argv[7]);
        
    }
    
    /// Copy Constructor
    TSimulationControl(const TSimulationControl &other){
        
        m_run_type          = other.m_run_type;
        m_geometry_type     = other.m_geometry_type;
        m_h_levels          = other.m_h_levels;
        m_n_elements        = other.m_n_elements;
        m_cartesian_stride  = other.m_cartesian_stride;
        m_approx_order      = other.m_approx_order;
        m_Hdiv_plusplus_Q   = other.m_Hdiv_plusplus_Q;
        
    }
    
    /// Copy Constructor
    TSimulationControl & operator=(const TSimulationControl &other){
        
        if (this != &other) {
            m_run_type          = other.m_run_type;
            m_geometry_type     = other.m_geometry_type;
            m_h_levels          = other.m_h_levels;
            m_n_elements        = other.m_n_elements;
            m_cartesian_stride  = other.m_cartesian_stride;
            m_approx_order      = other.m_approx_order;
            m_Hdiv_plusplus_Q   = other.m_Hdiv_plusplus_Q;
        }
        
        return *this;
        
    }
    
    void Print(){
        std::cout << "m_run_type            = " << m_run_type <<std::endl;
        std::cout << "m_geometry_type       = " << m_geometry_type <<std::endl;
        std::cout << "m_h_levels            = " << m_h_levels <<std::endl;
        std::cout << "m_n_elements          = " << m_n_elements <<std::endl;
        std::cout << "m_cartesian_stride    = " << m_cartesian_stride <<std::endl;
        std::cout << "m_approx_order        = " << m_approx_order <<std::endl;
        std::cout << "m_Hdiv_plusplus_Q     = " << m_Hdiv_plusplus_Q <<std::endl;
    }
    
};

#endif /* TSimulationControl_h */
