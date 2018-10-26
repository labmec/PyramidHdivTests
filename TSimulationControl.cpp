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
    m_h_levels          = 2;
    m_n_elements        = 1;
    m_red_black_stride_Q  = true;
    m_p_levels      = 1;
    m_Hdiv_plusplus_Q   = false;
    m_draw_vtk_Q        = false;
    
}

TSimulationControl::~TSimulationControl(){
    
}

TSimulationControl::TSimulationControl(char *argv[]){
    
    m_run_type          = static_cast<EApproxSpace>(atoi(argv[1]));
    m_geometry_type     = static_cast<EGeometryType>(atoi(argv[2]));
    m_h_levels          = atoi(argv[3]);
    m_n_elements        = atoi(argv[4]);
    m_red_black_stride_Q  = atoi(argv[5]);
    m_p_levels          = atoi(argv[6]);
    m_Hdiv_plusplus_Q   = atoi(argv[7]);
    m_draw_vtk_Q        = atoi(argv[8]);
    
}

TSimulationControl::TSimulationControl(const TSimulationControl &other){
    
    m_run_type          = other.m_run_type;
    m_geometry_type     = other.m_geometry_type;
    m_h_levels          = other.m_h_levels;
    m_n_elements        = other.m_n_elements;
    m_red_black_stride_Q  = other.m_red_black_stride_Q;
    m_p_levels          = other.m_p_levels;
    m_Hdiv_plusplus_Q   = other.m_Hdiv_plusplus_Q;
    m_draw_vtk_Q        = other.m_draw_vtk_Q;
    
}

TSimulationControl & TSimulationControl::operator=(const TSimulationControl &other){
    
    if (this != &other) {
        m_run_type          = other.m_run_type;
        m_geometry_type     = other.m_geometry_type;
        m_h_levels          = other.m_h_levels;
        m_n_elements        = other.m_n_elements;
        m_red_black_stride_Q  = other.m_red_black_stride_Q;
        m_p_levels          = other.m_p_levels;
        m_Hdiv_plusplus_Q   = other.m_Hdiv_plusplus_Q;
        m_draw_vtk_Q        = other.m_draw_vtk_Q;
    }
    
    return *this;
    
}

void TSimulationControl::Print(std::ostream &out){
    
    out << "m_run_type            = " << m_run_type <<std::endl;
    out << "m_geometry_type       = " << m_geometry_type <<std::endl;
    out << "m_h_levels            = " << m_h_levels <<std::endl;
    out << "m_n_elements          = " << m_n_elements <<std::endl;
    out << "m_red_black_stride_Q  = " << m_red_black_stride_Q <<std::endl;
    out << "m_p_levels            = " << m_p_levels <<std::endl;
    out << "m_Hdiv_plusplus_Q     = " << m_Hdiv_plusplus_Q <<std::endl;
    out << "m_draw_vtk_Q          = " << m_draw_vtk_Q <<std::endl;
}
