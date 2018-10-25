//
//  TSimulationControl.cpp
//  PyramidHdivTests
//
//  Created by Omar Durán on 10/19/18.
//

#include "TSimulationControl.h"


TSimulationControl::TSimulationControl(){
    
    m_run_type          = EDividedPyramidIncreasedOrder4;
    m_geometry_type     = EAcademic;
    m_h_levels          = 3;
    m_n_elements        = 1;
    m_red_black_stride_Q  = true;
    m_approx_order      = 1;
    m_Hdiv_plusplus_Q   = false;
    
}

TSimulationControl::~TSimulationControl(){
    
}

TSimulationControl::TSimulationControl(char *argv[]){
    
    m_run_type          = static_cast<EApproxSpace>(atoi(argv[1]));
    m_geometry_type     = static_cast<EGeometryType>(atoi(argv[2]));
    m_h_levels          = atoi(argv[3]);
    m_n_elements        = atoi(argv[4]);
    m_red_black_stride_Q  = atoi(argv[5]);
    m_approx_order      = atoi(argv[6]);
    m_Hdiv_plusplus_Q   = atoi(argv[7]);
    
}

TSimulationControl::TSimulationControl(const TSimulationControl &other){
    
    m_run_type          = other.m_run_type;
    m_geometry_type     = other.m_geometry_type;
    m_h_levels          = other.m_h_levels;
    m_n_elements        = other.m_n_elements;
    m_red_black_stride_Q  = other.m_red_black_stride_Q;
    m_approx_order      = other.m_approx_order;
    m_Hdiv_plusplus_Q   = other.m_Hdiv_plusplus_Q;
    
}

TSimulationControl & TSimulationControl::operator=(const TSimulationControl &other){
    
    if (this != &other) {
        m_run_type          = other.m_run_type;
        m_geometry_type     = other.m_geometry_type;
        m_h_levels          = other.m_h_levels;
        m_n_elements        = other.m_n_elements;
        m_red_black_stride_Q  = other.m_red_black_stride_Q;
        m_approx_order      = other.m_approx_order;
        m_Hdiv_plusplus_Q   = other.m_Hdiv_plusplus_Q;
    }
    
    return *this;
    
}

void TSimulationControl::Print(){
    std::cout << "m_run_type            = " << m_run_type <<std::endl;
    std::cout << "m_geometry_type       = " << m_geometry_type <<std::endl;
    std::cout << "m_h_levels            = " << m_h_levels <<std::endl;
    std::cout << "m_n_elements          = " << m_n_elements <<std::endl;
    std::cout << "m_red_black_stride_Q  = " << m_red_black_stride_Q <<std::endl;
    std::cout << "m_approx_order        = " << m_approx_order <<std::endl;
    std::cout << "m_Hdiv_plusplus_Q     = " << m_Hdiv_plusplus_Q <<std::endl;
}
