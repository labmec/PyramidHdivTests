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
    m_h_levels          = 1;
    m_red_black_stride_Q  = true;
    m_p_levels          = 1;
    m_Hdiv_plusplus_Q   = false;
    m_draw_vtk_Q        = true;
    
}

TSimulationControl::~TSimulationControl(){
    
}

TSimulationControl::TSimulationControl(char *argv[]){
    
    EApproxSpace run_type = static_cast<EApproxSpace>(atoi(argv[1]));
    EGeometryType geometry_type     = static_cast<EGeometryType>(atoi(argv[2]));
    int h_levels          = atoi(argv[3]);
    int p_levels          = atoi(argv[4]);
    bool Hdiv_plusplus_Q   = atoi(argv[5]);
    bool draw_vtk_Q        = atoi(argv[6]);

    bool red_black_stride_Q;
    if (run_type == ETetrahedra || geometry_type != EAcademic) {
        red_black_stride_Q = false;
    }else{
        red_black_stride_Q = true;
    }

    
    m_run_type          = run_type;
    m_geometry_type     = geometry_type;
    m_h_levels          = h_levels;
    m_red_black_stride_Q  = red_black_stride_Q;
    m_p_levels          = p_levels;
    m_Hdiv_plusplus_Q   = Hdiv_plusplus_Q;
    m_draw_vtk_Q        = draw_vtk_Q;
    
    
    
}

TSimulationControl::TSimulationControl(const TSimulationControl &other){
    
    m_run_type          = other.m_run_type;
    m_geometry_type     = other.m_geometry_type;
    m_h_levels          = other.m_h_levels;
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
    out << "m_red_black_stride_Q  = " << m_red_black_stride_Q <<std::endl;
    out << "m_p_levels            = " << m_p_levels <<std::endl;
    out << "m_Hdiv_plusplus_Q     = " << m_Hdiv_plusplus_Q <<std::endl;
    out << "m_draw_vtk_Q          = " << m_draw_vtk_Q <<std::endl;
}
