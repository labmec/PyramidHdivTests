//
//  TSimulationControl.cpp
//  PyramidHdivTests
//
//  Created by Omar Durán on 10/19/18.
//

#include "TSimulationControl.h"

EGeometryType gCurrentRun = EAcademic;

TSimulationControl::TSimulationControl(){
   // all variables are initialized in the declaration
    
}

TSimulationControl::~TSimulationControl(){
    
}

TSimulationControl::TSimulationControl(char *argv[]){
    
    EApproxSpace run_type = static_cast<EApproxSpace>(atoi(argv[1]));
    EGeometryType geometry_type     = static_cast<EGeometryType>(atoi(argv[2]));
    int h_level_min          = atoi(argv[3]);
    int h_level_max          = atoi(argv[4]);
    int p_level_min          = atoi(argv[5]);
    int p_level_max          = atoi(argv[6]);
    bool Hdiv_plusplus_Q   = atoi(argv[7]);
    bool draw_vtk_Q        = atoi(argv[8]);
    bool hybrid            = atoi(argv[9]);
    bool dry_run           = atoi(argv[10]);
    bool red_black_stride_Q;
    if (run_type == ETetrahedra || geometry_type != EAcademic) {
        red_black_stride_Q = false;
    }else{
        red_black_stride_Q = true;
    }

    
    m_run_type          = run_type;
    m_geometry_type     = geometry_type;
    m_h_level_min          = h_level_min;
    m_h_level_max          = h_level_max;
    m_red_black_stride_Q  = red_black_stride_Q;
    m_p_level_min          = p_level_min;
    m_p_level_max          = p_level_max;
    m_Hdiv_plusplus_Q   = Hdiv_plusplus_Q;
    m_draw_vtk_Q        = draw_vtk_Q;
    m_hybrid            = hybrid;
    m_dry_run           = dry_run;
    
    
}

TSimulationControl::TSimulationControl(const TSimulationControl &other){
    
    m_run_type          = other.m_run_type;
    m_geometry_type     = other.m_geometry_type;
    m_h_level_min          = other.m_h_level_min;
    m_h_level_max          = other.m_h_level_max;
    m_red_black_stride_Q  = other.m_red_black_stride_Q;
    m_p_level_min          = other.m_p_level_min;
    m_p_level_max          = other.m_p_level_max;
    m_Hdiv_plusplus_Q   = other.m_Hdiv_plusplus_Q;
    m_draw_vtk_Q        = other.m_draw_vtk_Q;
    m_dry_run           = other.m_dry_run;
    m_hybrid            = other.m_hybrid;
    
}

TSimulationControl & TSimulationControl::operator=(const TSimulationControl &other){
    
    if (this != &other) {
        m_run_type          = other.m_run_type;
        m_geometry_type     = other.m_geometry_type;
        m_h_level_min          = other.m_h_level_min;
        m_h_level_max          = other.m_h_level_max;
        m_red_black_stride_Q  = other.m_red_black_stride_Q;
        m_p_level_min          = other.m_p_level_min;
        m_p_level_max          = other.m_p_level_max;
        m_Hdiv_plusplus_Q   = other.m_Hdiv_plusplus_Q;
        m_draw_vtk_Q        = other.m_draw_vtk_Q;
        m_dry_run           = other.m_dry_run;
        m_hybrid            = other.m_hybrid;
    }
    
    return *this;
    
}

void TSimulationControl::Print(std::ostream &out){
    
    out << "m_run_type            = " << m_run_type << " " << RunType() << std::endl;
    out << "m_geometry_type       = " << m_geometry_type << " " << MeshGeometry() << std::endl;
    out << "m_h_levels            = " << m_h_level_min << " " << m_h_level_max <<std::endl;
    out << "m_red_black_stride_Q  = " << m_red_black_stride_Q <<std::endl;
    out << "m_p_levels            = " << m_p_level_min << " " << m_p_level_max <<std::endl;
    out << "m_Hdiv_plusplus_Q     = " << m_Hdiv_plusplus_Q <<std::endl;
    out << "m_draw_vtk_Q          = " << m_draw_vtk_Q <<std::endl;
    out << "m_hybrid              = " << m_hybrid << std::endl;
    out << "m_dry_run             = " << m_dry_run;
    if(m_dry_run == true) out << " results were not computed";
    out << std::endl;
    
}
