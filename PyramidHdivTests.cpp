/**
 * @file
 * @brief Tests for hdiv pyramid
 * @author Nathan Shauer
 * @since 2016
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include <math.h>

#include "tpzgeoelrefpattern.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZTimer.h"
#include "pzshapepiramHdiv.h"
#include "mixedpoisson.h"
#include "pzelctemp.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzgengrid.h"
//#include "pzpoisson3d.h"
#include "TPZDualPoisson.h"
#include "TPZPrimalPoisson.h"
#include "TPZCompMeshTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzelchdiv.h"
#include "pzshapetetra.h"
#include "pzelementgroup.h"
#include "TPZVecL2.h"
#include "pzl2projection.h"
#include "pzmatrix.h"
#include "TPZAcademicGeoMesh.h"
#include "TPZGmshReader.h"
#include "pzcondensedcompel.h"
#include "pzshapetriang.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

// Simulation Control
#include "TSimulationControl.h"

#include "TPZHybridizeHDiv.h"

#include "run_stats_table.h"
#ifdef USING_TBB
#include <tbb/tbb.h>
#endif

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.pyramtests"));
#endif

using namespace pzshape;

using namespace std;


void PrintArrayInMathematica(TPZVec<REAL> &array, std::ofstream &out, std::string arrayName);

void GenerateMathematicaWithConvergenceRates(TPZVec<REAL> &neqVec, TPZVec<REAL> &hSizeVec,
                                             TPZVec<REAL> &h1ErrVec, TPZVec<REAL> &l2ErrVec,
                                             TPZVec<REAL> &semih1ErrVec, TSimulationControl * control);
TPZGeoMesh *MalhaCubo(string &projectpath, const int &nref);
TPZGeoMesh *MalhaQuadrada(int &nelx, int &nely);
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);
void InsertElasticityCubo(TPZCompMesh *mesh);
void InsertBidimensionalPoisson(TPZCompMesh *cmesh, int &dim);
TPZGeoMesh * CreateGeoMesh1Pir();
TPZGeoMesh * CreateGeoMesh1Tet();
TPZGeoMesh * CreateGeoMeshHexaOfPir();
TPZGeoMesh * CreateGeoMeshHexaOfPirTetra();
TPZGeoMesh * CreateGeoMeshPrism();
TPZCompMesh * CreateCmeshPressure(TPZGeoMesh *gmesh, TSimulationControl * control, int p);
TPZCompMesh * CreateCmeshFlux(TPZGeoMesh *gmesh, TSimulationControl * control, int p);
TPZCompMesh * CreateCmeshMulti(TPZVec<TPZCompMesh *> &meshvec, TSimulationControl * control);
void LoadSolution(TPZCompMesh *cpressure);
void ProjectFlux(TPZCompMesh *cfluxmesh);
void GroupElements(TPZCompMesh *cmesh);
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void LaplaceExact(const TPZVec<REAL> &pt, TPZVec<STATE> &f);

void ExactSolution(const TPZVec<REAL> &pt, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);

TPZAutoPointer<TPZRefPattern> PyramidRef();

TPZAutoPointer<TPZRefPattern> PyramidTo4Tetrahedra();

void DividePyramids(TPZGeoMesh &gmesh);

void IncreasePyramidSonOrder(TPZVec<TPZCompMesh *> &meshvec, TSimulationControl * control, int pFlux);

void DivideBoundaryElements(TPZGeoMesh &gmesh, int exceptmatid = 3);

/// verify if the pressure space is compatible with the flux space
void VerifyDRhamCompatibility(TSimulationControl * control);

int gIntegrationOrder = 5;

/// Print Volumetric elements
void PrintGeometryVols(TPZGeoMesh * gmesh, std::stringstream & file_name);

//#define Solution_Sine
//#define Solution_MonoFourthOrder
//#define Solution_MonoCubic
//#define Solution_TriQuadratic
//#define Solution_MonoQuadratic
//#define Solution_MonoLinear
//#define Solution_Dupuit_Thiem
#define Solution_Spherical_Barrier

void Analytic(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &flux_and_f){
    
    flux_and_f.Resize(4, 1);
#ifdef Solution_Sine
    // Sine problem
    u[0] = sin(M_PI*pt[0])*sin(M_PI*pt[1])*sin(M_PI*pt[2]);
    flux_and_f(0,0) = -M_PI*cos(M_PI*pt[0])*sin(M_PI*pt[1])*sin(M_PI*pt[2]);
    flux_and_f(1,0) = -M_PI*cos(M_PI*pt[1])*sin(M_PI*pt[0])*sin(M_PI*pt[2]);
    flux_and_f(2,0) = -M_PI*cos(M_PI*pt[2])*sin(M_PI*pt[0])*sin(M_PI*pt[1]);
    flux_and_f(3,0) = 3*pow(M_PI,2)*sin(M_PI*pt[0])*sin(M_PI*pt[1])*sin(M_PI*pt[2]);
    return;
#endif
    
#ifdef Solution_MonoFourthOrder
    // x^4
    u[0] = pt[0]*pt[0]*pt[0]*pt[0];
    flux_and_f(0,0) = -4*pt[0]*pt[0]*pt[0];
    flux_and_f(1,0) = 0.0;
    flux_and_f(2,0) = 0.0;
    flux_and_f(3,0) = -12.0*pt[0]*pt[0];
    return;
#endif
    
#ifdef Solution_MonoCubic
    // x^3
    u[0] = pt[0]*pt[0]*pt[0];
    flux_and_f(0,0) = -3*pt[0]*pt[0];
    flux_and_f(1,0) = 0.0;
    flux_and_f(2,0) = 0.0;
    flux_and_f(3,0) = -6.0*pt[0];
    return;
#endif

#ifdef Solution_TriQuadratic
    // x^2 + y^2 + z^2
    u[0] = pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2];
    flux_and_f(0,0) = -2*pt[0];
    flux_and_f(1,0) = -2*pt[1];
    flux_and_f(2,0) = -2*pt[2];
    flux_and_f(3,0) = -6.0;
    return;
#endif
    
#ifdef Solution_MonoQuadratic
    // x^2 only
    u[0] = pt[0]*pt[0];
    flux_and_f(0,0) = -2*pt[0];
    flux_and_f(1,0) = 0.;
    flux_and_f(2,0) = 0.;
    flux_and_f(3,0) = -2.0;
    return;
#endif
    
#ifdef Solution_MonoLinear
    // x only
    u[0] = pt[0];
    flux_and_f(0,0) = -1.0;
    flux_and_f(1,0) = 0.;
    flux_and_f(2,0) = 0.;
    flux_and_f(3,0) = 0.0;
    return;
#endif
    
#ifdef Solution_Dupuit_Thiem
    REAL x = pt[0];
    REAL y = pt[1];
    
    flux_and_f.Resize(4,1);
    
    REAL r0 = 100.0;
    REAL r = sqrt(x*x+y*y);
    REAL theta = atan2(y,x);
    
    REAL costheta = cos(theta);
    REAL sintheta = sin(theta);
    
    // Gradient computations
    REAL Radialunitx = costheta;
    REAL Radialunity = sintheta;
    REAL Radialunitz = 0.0;
    
    REAL Thetaunitx = -sintheta;
    REAL Thetaunity = costheta;
    REAL Thetaunitz = 0.0;
    
    u[0] = 20 + log(r/r0);
    
    REAL dfdr = 1.0/r;
    REAL dfdTheta = 0.0;
    
    flux_and_f(0,0) = -1.0*(dfdr * Radialunitx + dfdTheta * Thetaunitx);
    flux_and_f(1,0) = -1.0*(dfdr * Radialunity + dfdTheta * Thetaunity);
    flux_and_f(2,0) = -1.0*(dfdr * Radialunitz + dfdTheta * Thetaunitz);
    
    flux_and_f(3,0) = 0.0;
    return;
#endif
    
#ifdef Solution_Spherical_Barrier
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    flux_and_f.Resize(4,1);
    

    REAL v = 1.0;
    REAL a = 0.1;
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
    REAL phi = atan2(y,x);
    
    REAL costheta = cos(theta);
    REAL sintheta = sin(theta);
    REAL sinphi = sin(phi);
    REAL cosphi = cos(phi);
    
    
    REAL p = -(1.0/(2.0*r*r))*((a*a*a) + 2*(r*r*r))*v*costheta;
    REAL dpdr = ((a*a*a)/(r*r*r) - 1.0)*v*costheta;
    REAL dpdTheta = 0.5*((a*a*a)/(r*r*r) + 2.0)*v*sintheta;
    REAL dpdPhi = 0.0;
    
    // Gradient computations
    REAL Radialunitx = sintheta*cosphi;
    REAL Radialunity = sintheta*sinphi;
    REAL Radialunitz = costheta;

    REAL Thetaunitx = cosphi*costheta;
    REAL Thetaunity = costheta*sinphi;
    REAL Thetaunitz = -sintheta;

    REAL Phiunitx = -sinphi;
    REAL Phiunity = cosphi;
    REAL Phiunitz = 0.0;
    
    u[0] = p;
    
    flux_and_f(0,0) = -1.0*(dpdr * Radialunitx + dpdTheta * Thetaunitx + dpdPhi * Phiunitx);
    flux_and_f(1,0) = -1.0*(dpdr * Radialunity + dpdTheta * Thetaunity + dpdPhi * Phiunity);
    flux_and_f(2,0) = -1.0*(dpdr * Radialunitz + dpdTheta * Thetaunitz + dpdPhi * Phiunitz);
    
    flux_and_f(3,0) = 0.0;

#endif
    
}

void AnalyticPotential(const TPZVec<REAL> &pt, TPZVec<STATE> &u) {
    TPZFMatrix<STATE> flux_and_f;
    Analytic(pt, u, flux_and_f);
    return;
}

void AnalyticFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &flux) {
    TPZManVector<STATE,1> u(1);
    TPZFMatrix<STATE> flux_and_f;
    Analytic(pt, u, flux_and_f);
    flux.Resize(3);
    flux[0] = flux_and_f(0,0);
    flux[1] = flux_and_f(1,0);
    flux[2] = flux_and_f(2,0);
    return;
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &u) {
    TPZFMatrix<STATE> flux_and_f;
    Analytic(pt, u, flux_and_f);
    return;
}

void BodyForcing(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
    TPZVec<STATE> u(1);
    TPZFMatrix<STATE> flux_and_f;
    Analytic(pt, u, flux_and_f);
    f[0]= flux_and_f(3,0);
    return;
}

void FluxFunc(const TPZVec<REAL> &pt, TPZVec<STATE> &flux)
{
    flux[0] = 0.;
}

void ExactSolution(const TPZVec<REAL> &pt, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol)
{
    int dir = 0;
    sol[0] = pt[dir];
    dsol.Zero();
    sol[0] = 1.;
    return;
    dsol(dir,0) = -1.;
    return;
    for (int i=0; i<3; i++) {
        dsol(i,0) = 1.-2.*pt[i];
    }
    for (int i=0; i<3; i++) {
        sol[0] *= pt[i]*(1.-pt[i]);
        for (int j=0; j<3; j++) {
            if (i != j) {
                dsol(j,0) *= pt[i]*(1.-pt[i]);
            }
        }
    }
}

void LaplaceExact(const TPZVec<REAL> &pt, TPZVec<STATE> &f)
{
    f[0] = 0.;
    return;
    for (int i=0; i<3; i++) {
        STATE term = 1.;
        for (int j=0; j<3; j++) {
            if (i!= j) {
                term *= pt[j]*(1.-pt[j]);
            }
        }
        f[0] += 2.*term;
    }
}

void CheckNormalContinuity(TPZCompMesh * cmesh);

void CheckNormalContinuity(TPZCompElSide & small, TPZCompElSide & large);

void ProjectAnalyticalSolution(TPZManVector<TPZCompMesh*,2> & mesh_vector);

void CompareShapeFunctions(TPZCompElSide & cel_hdiv_side, TPZVec<REAL> & xi_hdiv, TPZCompEl * cel_bound, TPZVec<REAL> & xi_bound, TPZVec<REAL> & normal);

int ConvergenceTest(TSimulationControl * control);

void ComputeConvergenceRates(TPZVec<REAL> &error, TPZVec<REAL> &h_size, TPZVec<REAL> &convergence);

int ComputeApproximation(TSimulationControl * sim_control);

TPZGeoMesh * GeometryConstruction(int n_h_ref_level, REAL & h_min, int & n_elements, TSimulationControl * sim_control);

void UnwrapMesh(TPZCompMesh *cmesh);

std::string PyramidApproxSpaceType(TSimulationControl * control);

int integer_power(int base, unsigned int exp){
    
    if (exp == 0)
        return 1;
    int temp = integer_power(base, exp/2);
    if (exp%2 == 0)
        return temp*temp;
    else
        return base*temp*temp;
    
}

void ComputeCharacteristicHElSize(TPZGeoMesh * geometry, REAL & h_min, int & n_elements);

void CreateFlattenGeometry(TPZGeoMesh & gmesh);

void FlipPyramids(TPZGeoMesh * gmesh);

void TestingCondensation();

int main(int argc, char *argv[])
{
    /// Global controls
    HDivPiola = 1;

#ifdef LOG4CXX
//    std::string dirname = PZSOURCEDIR;
//    std::string FileName = dirname;
//    FileName = dirname + "/Projects/PyramidHdivTests/";
//    FileName += "pyramlogfile.cfg";
    InitializePZLOG();
#endif
    
//    TestingCondensation();
//    return 0;
    
    TSimulationControl * sim_control = NULL;
    if(argc != 7){
        sim_control = new TSimulationControl;
        sim_control->m_run_type = EHexaHedra;
        sim_control->m_h_levels = 3;
    }
    else{
        sim_control = new TSimulationControl(argv);
    }
    std::cout << "Simulation control object with parameters " << std::endl;
    sim_control->Print();
    
    ComputeApproximation(sim_control);
    return 0;
}

void ComputeCharacteristicHElSize(TPZGeoMesh * geometry, REAL & h_min, int & n_elements){
    
    h_min   = 1.0;
    n_elements = 0;
    
    REAL h;
    int nel = geometry->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if(!gel){
            DebugStop();
        }
#endif
        if (gel->Dimension() != geometry->Dimension() || gel->HasSubElement() == 1) {
            continue;
        }
        
        n_elements++;
        
        h = gel->CharacteristicSize();
        
        if (h < h_min) {
            h_min = h;
        }
    }
    
}

void TestingCondensation(){
    
    
    TPZMatRed<STATE, TPZFMatrix<STATE> > fCondensed(10,5);
    fCondensed.Zero();
    TPZAutoPointer<TPZMatrix<STATE> > K00_ini = new TPZFMatrix<STATE>(5,5);
    fCondensed.SetK00(K00_ini);
    
    int64_t dim0 = fCondensed.Dim0();
    int64_t dim1 = fCondensed.Dim1();
    int64_t rows = fCondensed.Rows();
    int64_t cols = fCondensed.Cols()+1;
    
    TPZFMatrix<STATE> ef(rows,1,10);
    fCondensed.SetF(ef);
    
    // Prencher the Matrix
    // Devloo
    unsigned int n = fCondensed.Rows();
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            fCondensed(i,j) = (1.0/((i+j+2)*(i+j+2)));
        }
        fCondensed(i,i) += 1.0;
    }
    
    {// Initialization for pivot structure
        TPZFMatrix<STATE> * K00_temp = dynamic_cast<TPZFMatrix<STATE> * >(fCondensed.K00().operator->());
        K00_temp->InitializePivot();
    }

    TPZFMatrix<STATE> KF(rows,cols); //  Local object
    for(int64_t i=0; i<rows;i++) //
        for (int64_t j=0; j<cols; j++)
        {
            if (j<rows)
                KF(i,j) = fCondensed(i,j);
            else
                KF(i,j) = ef(i,j-rows);
        }
    
    for (int64_t i=0; i<dim1; i++) //
        for (int64_t j=0; j<dim0; j++)
            fCondensed.K10().operator()(i,j)=KF(i+dim0,j);
    
    for (int64_t i=0; i<dim1; i++) //
        for (int64_t j=0; j<dim1; j++)
            fCondensed(i+dim0,j+dim0)=KF(i+dim0,j+dim0);
    
    fCondensed.SetF(ef);
    
    // LDLt Decomposition
    for (int64_t i=0; i<rows-dim1; i++)
    {
        for(int64_t j=i+1;j<cols;j++)
        {
            if (j<rows)
            {
                KF(j,i)/=KF(i,i);
                KF(i,j)/=KF(i,i);
            }
            else
                KF(i,j)/=KF(i,i);
        }
        
#ifdef STATEdouble
        cblas_dger (CblasColMajor, rows-i-1, cols-i-1,
                    -KF(i,i), &KF(i+1,i), 1,
                    &KF(i,i+1), rows, &KF(i+1,i+1), rows);
#else
        DebugStop();
#endif
    }
    
    for (int64_t i=dim0; i< rows; i++)
    {
        ef(i,0) = KF.GetVal(i,fCondensed.Rows());
        for (int64_t j=dim0; j< fCondensed.Rows(); j++)
        {
            fCondensed(i,j) = KF.GetVal(i,j);
        }
    }
    
    TPZAutoPointer<TPZMatrix<STATE> > K00 = fCondensed.K00();
    for (int64_t i=0; i<dim0; i++)
        for (int64_t j=0; j<dim0; j++)
            K00->operator()(i, j) = KF(i,j);

    fCondensed.K00()->SetIsDecomposed(ELDLt);
    
//    if(0){
//
//        for (int64_t i=0; i<dim0; i++){ // Substituindo valores obtidos para K01 usando o BLAS
//            for (int64_t j=0; j<dim1; j++){
//                fCondensed.K01().operator()(i,j) = KF(i,j+dim0)/KF(i,i);
//            }
//
//        }
//
//        double scale = 1.0;
//        cblas_dtrsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasUnit,dim0,dim1,scale,&fCondensed.K00().operator->()->operator()(0,0),dim0,&fCondensed.K01()(0,0),dim0);
//    }else{
    fCondensed.K00()->SolveDirect(fCondensed.K01(), ELDLt);
//    }
    fCondensed.SetF(ef);
    fCondensed.SetK01IsComputed(true);
    fCondensed.SetF0IsComputed(true);
    fCondensed.SetReduced();// Directive that instructs the object to behave as reduced matrix.
    fCondensed.Print("Condensed = ",std::cout);
}

int ComputeApproximation(TSimulationControl * sim_control)
{
    
    std::string e_accuracy_type = "Hdiv";
    if(sim_control->m_Hdiv_plusplus_Q){
        e_accuracy_type = "Hdiv++";
    }
    std::ofstream output("convergence_summary.txt",std::ios::app);
    output << std::endl;
    output << "Approximation space type : " << PyramidApproxSpaceType(sim_control) << std::endl;
    output << "Enhanced accuracy type : " << e_accuracy_type << std::endl;
    output << "TSimulation control used : " << std::endl;
    sim_control->Print(output);
    output << std::endl;
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream str;
        str << std::endl;
        str << " Runing Mixed Formulation with Pyramids. " << std::endl;
        LOGPZ_DEBUG(logger,str.str())
    }
#endif

    EApproxSpace run_type = sim_control->m_run_type;
    int n_p_levels = sim_control->m_p_levels;
    int n_h_levels  = sim_control->m_h_levels;
    const int dim = 3;
    
    /// Hard code controls
    bool should_renumber_Q = true;
    bool use_pardiso_Q = true;
    const int n_threads_error = 64;
    const int n_threads_assembly = 64;
    bool keep_lagrangian_multiplier_Q = true;
    bool keep_matrix_Q = false;
    TPZGeoMesh *gmesh = NULL;
    
    
    TPZAutoPointer<TPZRefPattern> pyramid_ref_pattern;
    bool refine_into_4_tetrahedra_Q = sim_control->m_run_type == EDividedPyramid4 || sim_control->m_run_type == EDividedPyramidIncreasedOrder4;
    if(refine_into_4_tetrahedra_Q){
        pyramid_ref_pattern = PyramidTo4Tetrahedra();
    }
    else {
        pyramid_ref_pattern = PyramidRef();
    }

    REAL assemble_time, solving_time, error_time;

    TPZManVector<REAL,20> h_vec(n_h_levels+1,0.);
    TPZManVector<REAL,20> primal_error(n_h_levels+1,0.);
    TPZManVector<REAL,20> dual_error(n_h_levels+1,0.);
    TPZManVector<REAL,20> div_error(n_h_levels+1,0.);
    
    for (int p = 1; p <= n_p_levels; p++) {
        
        output << std::endl;
        output << " Polynomial order  =  " << p << std::endl;
        output << setw(5) <<  " h_level " << setw(10) << " n_elements" << setw(10) << " h" << setw(15) << " ndof" << setw(15) << " ndof_cond" << setw(25) << " assemble_time (msec)" << setw(25) << " solving_time (msec)" << setw(25) << " error_time (msec)" << setw(25) << " Primal l2 error" << setw(25) << " Dual l2 error"  << setw(25) << " Div l2 error" << endl;
    
        for (int i = 0 ; i <= n_h_levels; i++){
            
#ifdef USING_BOOST
            boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
#endif
//            int n_elements = integer_power(2,i);
//            REAL h = 1./REAL(n_elements);
            int n_elements;
            REAL h;
            gmesh = GeometryConstruction(i,h,n_elements,sim_control);
            
            TPZManVector<TPZCompMesh*,2> meshvecOrig(2);
            
            /// Construction for Hdiv (velocity) approximation space
            meshvecOrig[0] = CreateCmeshFlux(gmesh, sim_control, p);
            
            /// Construction for L2 (pressure) approximation space
            meshvecOrig[1] = CreateCmeshPressure(gmesh, sim_control, p);
            
            if (run_type == EDividedPyramidIncreasedOrder || run_type == EDividedPyramidIncreasedOrder4)
            {
                IncreasePyramidSonOrder(meshvecOrig,sim_control,p);
            }
            
//            ProjectAnalyticalSolution(meshvecOrig);
            
#ifdef LOG4CXX
            if (logger->isDebugEnabled())
            {
                std::stringstream sout;
                meshvecOrig[0]->Print(sout);
                meshvecOrig[1]->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            // ------------------ Create CMesh multiphysics -------------------
            TPZCompMesh *cmeshMultOrig = CreateCmeshMulti(meshvecOrig,sim_control);

            TPZCompMesh *cmeshMult = 0;
            TPZManVector<TPZCompMesh *,2> meshvec(2,0);
            bool hybrid = false;
            if(hybrid)
            {
                TPZHybridizeHDiv hybrid;
                std::tuple<TPZCompMesh*, TPZVec<TPZCompMesh*> > chunk;
                chunk = hybrid.Hybridize(cmeshMultOrig, meshvecOrig);
                TPZCompMesh *cmeshMultHybrid = std::get<0>(chunk);
                TPZManVector<TPZCompMesh *,2> meshvecHybrid = std::get<1>(chunk);
//                hybrid.GroupElements(cmeshMultHybrid);
                cmeshMult = cmeshMultHybrid;
                meshvec = meshvecHybrid;
                
                CheckNormalContinuity(meshvec[0]);
                
            }
            else
            {
            
                TPZCompMeshTools::GroupElements(cmeshMultOrig);
                std::cout << "Created grouped elements\n";
                TPZCompMeshTools::CreatedCondensedElements(cmeshMultOrig, keep_lagrangian_multiplier_Q, keep_matrix_Q);
                std::cout << "Created condensed elements\n";
                cmeshMultOrig->CleanUpUnconnectedNodes();
                cmeshMultOrig->ExpandSolution();
                
                cmeshMult = cmeshMultOrig;
                meshvec = meshvecOrig;
            }
#ifdef LOG4CXX
            if (logger->isDebugEnabled())
            {
                std::stringstream sout;
                cmeshMult->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            cmeshMult->CleanUpUnconnectedNodes();

            
            // ------------------ Creating Analysis object -------------------
            TPZAnalysis an(cmeshMult,should_renumber_Q);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            
            if (use_pardiso_Q) {
                TPZSymetricSpStructMatrix sparse(cmeshMult);
                sparse.SetNumThreads(n_threads_assembly);
                an.SetStructuralMatrix(sparse);
                an.SetSolver(step);
            }else{
                TPZSkylineStructMatrix skyl(cmeshMult);
                skyl.SetNumThreads(n_threads_assembly);
                an.SetStructuralMatrix(skyl);
                an.SetSolver(step);
            }
            
            std::cout << "Starting assemble..." << std::endl;
            std::cout << "Nequations = " << cmeshMult->NEquations() << std::endl;
#ifdef USING_BOOST
            boost::posix_time::ptime tass1 = boost::posix_time::microsec_clock::local_time();
#endif
            // ------------------ Assembling -------------------
            an.Assemble();
#ifdef USING_BOOST
            boost::posix_time::ptime tass2 = boost::posix_time::microsec_clock::local_time();
            assemble_time = boost::numeric_cast<REAL>((tass2 - tass1).total_milliseconds());
            std::cout << "Total wall time of Assemble = " << assemble_time << " ms." << std::endl;
#endif
            
            TPZAutoPointer<TPZMatrix<STATE> > mat = an.Solver().Matrix();
            
            std::cout << "Assembled!" << std::endl;
        
            std::cout << "Starting Solve..." << std::endl;
#ifdef USING_BOOST
            boost::posix_time::ptime tsolve1 = boost::posix_time::microsec_clock::local_time();
#endif
            an.Solve();
            
#ifdef USING_BOOST
            boost::posix_time::ptime tsolve2 = boost::posix_time::microsec_clock::local_time();
            solving_time = boost::numeric_cast<REAL>((tsolve2 - tsolve1).total_milliseconds());
            std::cout << "Total wall time of Solve = " << solving_time << " ms." << std::endl;
#endif
            // Getting dof information before unwrap the mesh
            int ndof = meshvecOrig[0]->Solution().Rows()+ meshvecOrig[1]->Solution().Rows();
            int ndof_cond = cmeshMult->NEquations();
            
            UnwrapMesh(cmeshMult);
            an.LoadSolution();
            cmeshMult->Solution() *= -1.0; // Because the material contributions are expressed in residual form
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmeshMult);
            std::cout << "Solved!" << std::endl;
            
            // ------------------ Post Processing VTK -------------------
            if (sim_control->m_draw_vtk_Q) {
                
                std::cout << "Starting Post-processing..." << std::endl;
                TPZStack<std::string> scalnames, vecnames;
                scalnames.Push("p");
                scalnames.Push("div_q");
                vecnames.Push("q");
                std::string plotfile = "Approximated_Solution.vtk";
                an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
                
                int postprocessresolution = 0;
                an.PostProcess(postprocessresolution);
                
            }
            std::cout << "Calculating error..." << std::endl;
            
            an.SetExact(Analytic);
            TPZManVector<REAL,3> errors(3,1);
            an.SetThreadsForError(n_threads_error);
#ifdef USING_BOOST
            boost::posix_time::ptime terr1 = boost::posix_time::microsec_clock::local_time();
#endif
            an.PostProcessError(errors,false);
#ifdef USING_BOOST
            boost::posix_time::ptime terr2 = boost::posix_time::microsec_clock::local_time();
#endif
            
#ifdef USING_BOOST
            error_time = boost::numeric_cast<REAL>((terr2 - terr1).total_milliseconds());
            std::cout << "Computed errors " << errors << std::endl;
            std::cout << "Total wall time of PostProcessError = " << terr2 - terr1 << " ms." << std::endl;
#endif
            
#ifdef USING_BOOST
            boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
            std::cout << "Total wall time of simulation = " << tsim2 - tsim1 << " s" << std::endl;
#endif
            
            int h_level = i;
            REAL p_error = errors[0]; // primal
            REAL d_error = errors[1]; // dual
            REAL h_error = errors[2]; // div
            output << setw(5) <<  h_level << setw(10) << n_elements << setw(15) << h << setw(15) << ndof << setw(15) << ndof_cond << setw(25) << assemble_time << setw(25) << solving_time << setw(25) << error_time << setw(25) << p_error << setw(25) << d_error  << setw(25) << h_error << endl;
            
            // Storage data for h convergence rates
            primal_error[i] = p_error;
            dual_error[i] = d_error;
            div_error[i] = h_error;
            h_vec[i] = h;
            
#ifdef USING_BOOST
            boost::posix_time::ptime deletion_t1 = boost::posix_time::microsec_clock::local_time();
#endif
            
            delete cmeshMult;
            delete meshvec[0];
            int64_t nel = meshvec[1]->NElements();
            gmesh->ResetReference();
            for (int64_t el=0; el<nel; el++) {
                TPZCompEl *cel = meshvec[1]->Element(el);
                delete cel;
            }
            delete meshvec[1];
            delete gmesh;
            
#ifdef USING_BOOST
            boost::posix_time::ptime deletion_t2 = boost::posix_time::microsec_clock::local_time();
            std::cout << "Deletion time = " << deletion_t2 - deletion_t1 << " s" << std::endl;
#endif
            
            
        }// n_h_levels
        
        // Computing approximation rates
        TPZVec<REAL> p_conv(n_h_levels), d_conv(n_h_levels), h_conv(n_h_levels);
        ComputeConvergenceRates(primal_error, h_vec, p_conv);
        ComputeConvergenceRates(dual_error, h_vec, d_conv);
        ComputeConvergenceRates(div_error, h_vec, h_conv);
        
        // print convergence summary
        output << std::endl;
        output << " Convergence rates summary " << std::endl;
        output << " Polynomial order  =  " << p << std::endl;
        output << " Primal convergence rates = " << setw(5) << p_conv << std::endl;
        output << " Dual convergence rates = " << setw(5) << d_conv << std::endl;
        output << " Divergence convergence rates = " << setw(5) << h_conv << std::endl;
        output << std::endl;
        output.flush();
        
    }// n_p_levels
    
    output << std::endl;
    output << " ------------------------------------------------------------------ " << std::endl;
    output.flush();
    
    return 0;

}

std::string PyramidApproxSpaceType(TSimulationControl * control){
    
    std::string type;
    EApproxSpace run_type = control->m_run_type;
    switch (run_type) {
        case EDividedPyramid:
            type = "Non-conformal pyramid divided in two tetrahedra.";
            break;
        case EDividedPyramid4:
            type = "Non-conformal pyramid divided in four tetrahedra.";
            break;
        case EDividedPyramidIncreasedOrder:
            type = "Conformal pyramid divided in two tetrahedra.";
            break;
        case EDividedPyramidIncreasedOrder4:
            type = "Conformal pyramid divided in four tetrahedra.";
            break;
        case ETetrahedra:
            type = "Conformal tetrahedral mesh.";
            break;
        case EHexaHedra:
            type = "Conformal hexahedral mesh.";
            break;
        default:
            DebugStop();
            break;
    }
    return type;
}

void ProjectAnalyticalSolution(TPZManVector<TPZCompMesh*,2> & mesh_vector){
    
    {
        std::ofstream flux_file("flux_cmesh.txt");
        mesh_vector[0]->Print(flux_file);
        
        std::ofstream pressure_file("pressure_cmesh.txt");
        mesh_vector[1]->Print(pressure_file);
    }
    
    int dim = 3;
    int n_threads_assembly = 12;
    {
        int mesh_index = 0; // flux projection
        TPZAnalysis an(mesh_vector[mesh_index],false);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        TPZSymetricSpStructMatrix sparse(mesh_vector[mesh_index]);
        sparse.SetNumThreads(n_threads_assembly);
        an.SetStructuralMatrix(sparse);
        an.SetSolver(step);
        an.Assemble();
        an.Solve();
        an.LoadSolution();
//        std::cout << "Flux Post-processing..." << std::endl;
//        TPZStack<std::string> scalnames, vecnames;
//        vecnames.Push("Solution");
//        std::string plotfile = "Projected_flux.vtk";
//        an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
//
//        int postprocessresolution = 0;
//        an.PostProcess(postprocessresolution);
        
    }

    {
        int mesh_index = 1; // pressure projection
        TPZAnalysis an(mesh_vector[mesh_index],false); // Not working when second arg is true
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        TPZSymetricSpStructMatrix sparse(mesh_vector[mesh_index]);
        sparse.SetNumThreads(n_threads_assembly);
        an.SetStructuralMatrix(sparse);
        an.SetSolver(step);
        an.Assemble();
        an.Solve();
        an.LoadSolution();
//        std::cout << "Pressure Post-processing..." << std::endl;
//        TPZStack<std::string> scalnames, vecnames;
//        scalnames.Push("Solution");
//        std::string plotfile = "Projected_pressure.vtk";
//        an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
//
//        int postprocessresolution = 0;
//        an.PostProcess(postprocessresolution);
    }
    
}

void CheckNormalContinuity(TPZCompMesh * cmesh){
    
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    
    int mesh_dim  = cmesh->Dimension();
    int64_t n_dof = cmesh->Solution().Rows();
    for (int64_t i = 0; i <n_dof; i++) {
        cmesh->Solution()(i,0) = 100.0*rand()/RAND_MAX;
    }
    cmesh->LoadSolution(cmesh->Solution());
    
    int64_t nel = cmesh->NElements();
    for (int64_t icel = 0; icel < nel; icel++) {
        TPZCompEl * cel  = cmesh->Element(icel);
        if (!cel) {
            continue;
        }
        
        
        TPZGeoEl * gel = cel->Reference();
        if (!gel or gel->Dimension()!= mesh_dim) {
            continue;
        }
        
        TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement * >(cel);
        if (!intel) {
            DebugStop();
        }
        
        int n_sides = gel->NSides();
        for (int iside = 0; iside < n_sides; iside++) {
            if (gel->SideDimension(iside)!=mesh_dim-1) {
                continue;
            }
            int64_t con_index = intel->SideConnectIndex(0, iside);
            TPZGeoElSide gel_side(gel,iside);
            TPZGeoElSide neigh = gel_side.Neighbour();
            
            while (neigh not_eq gel_side) {
                TPZCompElSide neigh_cel = neigh.Reference();
                TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement *>(neigh_cel.Element());
                
                if (intel) {
                    int64_t intel_con_index = intel->SideConnectIndex(0, neigh_cel.Side());
                    if (intel_con_index == con_index ) {
                        break;
                    }
                }

                neigh = neigh.Neighbour();
            }
            
            if (neigh == gel_side) {
                DebugStop();
            }
            
            TPZCompElSide cel_neigh = neigh.Reference();
            TPZCompElSide cel_side(cel,iside);
            CheckNormalContinuity(cel_side, cel_neigh);
            
        }
        
    }
    
}

void CheckNormalContinuity(TPZCompElSide & small, TPZCompElSide & large){
    
    int order = 5;
    REAL eps = 1.0e-8;
    
    TPZGeoElSide gel_small(small.Reference());
    TPZGeoElSide gel_large(large.Reference());
    
    TPZTransform<REAL> small_to_large(gel_small.Dimension(),gel_small.Dimension());
    gel_small.SideTransform3(gel_large, small_to_large);
    
    int small_vol_side = gel_small.Element()->NSides()-1;
    TPZTransform<REAL> side_to_vol_small = gel_small.Element()->SideToSideTransform(small.Side(),small_vol_side);
    
    int large_vol_side = gel_large.Element()->NSides()-1;
    TPZTransform<REAL> side_to_vol_large = gel_large.Element()->SideToSideTransform(large.Side(),large_vol_side);
    
    TPZIntPoints * int_rule = gel_small.Element()->CreateSideIntegrationRule(gel_small.Side(), order);
    
    int n_p = int_rule->NPoints();
    TPZManVector<REAL,3> x_hat(gel_small.Dimension());
    TPZManVector<REAL,3> x_hat_large(gel_large.Dimension());
    
    TPZManVector<REAL,3> x_hat_vol(gel_small.Element()->Dimension());
    TPZManVector<REAL,3> x_hat_vol_large(gel_large.Element()->Dimension());
    
    TPZManVector<REAL,3> sol_small, sol_large;
    TPZManVector<REAL,3> x_small(3), x_large(3);
    REAL w;
    
    TPZManVector<REAL,3> normal(3);
    int flux_var = 0;
    for (int ip = 0; ip < n_p; ip++) {
        
        int_rule->Point(ip, x_hat, w);
        
        side_to_vol_small.Apply(x_hat, x_hat_vol);
        small.Element()->Solution(x_hat_vol, flux_var, sol_small);
        gel_small.Element()->X(x_hat_vol, x_small);
        
        small_to_large.Apply(x_hat, x_hat_large);
        side_to_vol_large.Apply(x_hat_large, x_hat_vol_large);
        large.Element()->Solution(x_hat_vol_large, flux_var, sol_large);
        gel_large.Element()->X(x_hat_vol_large, x_large);
        
        gel_small.Normal(x_hat, gel_small.Element(), gel_large.Element(), normal);
        
        REAL normal_flux = 0;
        REAL var_norm = 0;
        REAL x_norm = 0;
        
        for (int k = 0; k < x_large.size(); k++) {
            normal_flux += sol_small[k]*normal[k];
            x_norm += (x_large[k] - x_small[k])*(x_large[k] - x_small[k]);
        }
        
        var_norm = (normal_flux - sol_large[0])*(normal_flux - sol_large[0]);
        
        var_norm = sqrt(var_norm);
        x_norm = sqrt(x_norm);
        
        bool check_Q = (var_norm < eps) and  (x_norm < eps);
        if (not check_Q) {
            std::cout << "ok the flux is not continuous. " << std::endl;
            CompareShapeFunctions(small, x_hat_vol, large.Element(), x_hat_vol_large, normal);
            DebugStop();
        }
        
        
    }
}

void CompareShapeFunctions(TPZCompElSide & cel_hdiv_side, TPZVec<REAL> & xi_hdiv, TPZCompEl * cel_bound, TPZVec<REAL> & xi_bound, TPZVec<REAL> & normal){
    
    REAL eps = 1.0e-10;
    
    TPZCompEl * cel_hdiv = cel_hdiv_side.Element();
    
    TPZMaterialData data_hdiv,data_bound;
    TPZInterpolatedElement * intel_hdiv = dynamic_cast<TPZInterpolatedElement * >(cel_hdiv);
    TPZInterpolatedElement * intel_bound = dynamic_cast<TPZInterpolatedElement * >(cel_bound);
    
    if (!intel_hdiv or !intel_bound) {
        DebugStop();
    }
    
    int int_con_index = intel_hdiv->SideConnectLocId(0, cel_hdiv_side.Side());
    int first_shape = 0;
    for (int ic = 0; ic < int_con_index; ic++) {
        first_shape += intel_hdiv->Connect(ic).NShape();
    }
    
    intel_hdiv->InitMaterialData(data_hdiv);
    intel_bound->InitMaterialData(data_bound);
    
    intel_hdiv->ComputeRequiredData(data_hdiv, xi_hdiv);
    intel_bound->ComputeRequiredData(data_bound, xi_bound);
    
    TPZFMatrix<REAL> fluxes;
    data_hdiv.ComputeFluxValues(fluxes);
//    {
//        TPZFNMatrix<200,REAL> phi(100,1,-1),dphi(3,100);
//        TPZManVector<int,20> order(15,6);
//        TPZManVector<int64_t,4> id(4,4);
//        for (int i = 0; i < 4; i++) {
//            id[i] = cel_hdiv->Reference()->NodePtr(i)->Id();
//        }
//        pzshape::TPZShapeTetra::Shape(xi_hdiv, id, order, phi, dphi);
//        int aka = 0;
//    }
//
//    {
//        TPZFNMatrix<200,REAL> phi(100,1,-1),dphi(3,100);
//        TPZManVector<int,20> order(15,6);
//        TPZManVector<int64_t,4> id(3,3);
//        for (int i = 0; i < 3; i++) {
//            id[i] = cel_bound->Reference()->NodePtr(i)->Id();
//        }
//        pzshape::TPZShapeTriang::Shape(xi_bound, id, order, phi, dphi);
//        int aka = 0;
//    }
    
    // compute the difference
    int n_shape = data_bound.phi.Rows();
    TPZFMatrix<REAL> flux_normal(n_shape,1,0.0);
    TPZFNMatrix<100,REAL> normal_component(n_shape,1,0.0);
    REAL diff_accum = 0;
    TPZFNMatrix<100,REAL> phi_s(n_shape,1);
    for (int i = 0; i < n_shape; i++) {
        int vec_index = data_hdiv.fVecShapeIndex[first_shape + i].first;
        int shape_index = data_hdiv.fVecShapeIndex[first_shape + i].second;
        phi_s(i,0) = data_hdiv.phi(shape_index,0);
        for (int k = 0; k < 3; k++) {
            normal_component(i,0) += data_hdiv.fNormalVec(k,vec_index)* normal[k];
            flux_normal(i,0) += fluxes(k,first_shape + i) * normal[k];
        }

        diff_accum += (flux_normal(i,0) - data_bound.phi(i,0))*(flux_normal(i,0) - data_bound.phi(i,0));
    }
    
    diff_accum/=n_shape;
    diff_accum = sqrt(diff_accum);

    bool check_Q = diff_accum < eps;
    if (!check_Q) {
        std::cout << "Some shape functions are different." <<std::endl;
        flux_normal.Print("Normal flux = ",std::cout);
        data_bound.phi.Print("Normal flux from boundary = ",std::cout);
        (flux_normal - data_bound.phi).Print("Difference = ",std::cout);
    }
    
}

void ComputeConvergenceRates(TPZVec<REAL> &error, TPZVec<REAL> &h_size, TPZVec<REAL> &convergence){
    
    int ndata = error.size();
    for (int i = 1; i < ndata; i++) {
        STATE logerror = log(error[i-1]);
        STATE logerrori = log(error[i]);
        convergence[i-1] = (logerrori - logerror)/(log(h_size[i])-log(h_size[i-1]));
    }
}

TPZGeoMesh * GeometryConstruction(int h_ref_level, REAL & h_min, int & n_elements, TSimulationControl * sim_control){
    
    // ------------------ Creating GeoMesh -------------------
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    EApproxSpace run_type = sim_control->m_run_type;
    std::cout << "Creating geometry description." << std::endl;
    switch (sim_control->m_geometry_type) {
        case EAcademic:{

            const int n_elements = integer_power(2,h_ref_level); // num of hexes in x y and z
            const int matid = 1;
            /// Defining the type o geometry
            TPZAcademicGeoMesh academic;
            if(sim_control->m_run_type == EHexaHedra)
            {
                academic.SetMeshType(TPZAcademicGeoMesh::EHexa);
            }
            else
            {
                academic.SetMeshType(TPZAcademicGeoMesh::EPyramid);
            }
            if (sim_control->m_run_type == ETetrahedra) {
                academic.SetMeshType(TPZAcademicGeoMesh::ETetrahedra);
            }
            TPZManVector<int,6> BCids(6,-1); // ids of the bcs
            academic.SetBCIDVector(BCids);
            academic.SetMaterialId(matid);
            academic.SetNumberElements(n_elements);
            if (run_type != ETetrahedra && run_type != EHexaHedra) {
                if (sim_control->m_red_black_stride_Q) {
                    gmesh = academic.RedBlackPyramidalAndHexagonalMesh();
                }else{
                    gmesh = academic.PyramidalAndTetrahedralMesh();
                }
            }
            else
            {
                gmesh = academic.CreateGeoMesh();
            }
        }
            break;
        case EVerticalWellHe:{
            TPZGmshReader Geometry;
            REAL s = 1.0;
            Geometry.SetfDimensionlessL(s);
            std::string gmsh_file("vertical_wellbore_He.msh");
            gmesh = Geometry.GeometricGmshMesh(gmsh_file);
            
        }
            break;
        case EVerticalWellTe:{
            TPZGmshReader Geometry;
            REAL s = 1.0;
            Geometry.SetfDimensionlessL(s);
            std::string gmsh_file("vertical_wellbore_Te.msh");
            gmesh = Geometry.GeometricGmshMesh(gmsh_file);
            
        }
            break;
        case EVerticalWellHePyTe:{
            TPZGmshReader Geometry;
            REAL s = 1.0;
            Geometry.SetfDimensionlessL(s);
            std::string gmsh_file("vertical_wellbore_hybrid.msh");
            gmesh = Geometry.GeometricGmshMesh(gmsh_file);
        
        }
            break;
        case ESphericalBarrierHePyTe:{
            TPZGmshReader Geometry;
            REAL s = 1.0;
            Geometry.SetfDimensionlessL(s);
            std::string gmsh_file("spherical_obstacle_hybrid.msh");
            gmesh = Geometry.GeometricGmshMesh(gmsh_file);
            
        }
            break;
        case ESphericalBarrierTe:{
            TPZGmshReader Geometry;
            REAL s = 1.0;
            Geometry.SetfDimensionlessL(s);
            std::string gmsh_file("spherical_obstacle_Te.msh");
            gmesh = Geometry.GeometricGmshMesh(gmsh_file);
            
        }
            break;
        case ESphericalBarrierHe:{
            TPZGmshReader Geometry;
            REAL s = 1.0;
            Geometry.SetfDimensionlessL(s);
            std::string gmsh_file("spherical_obstacle_He.msh");
            gmesh = Geometry.GeometricGmshMesh(gmsh_file);
            
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
    if (sim_control->m_geometry_type!=EAcademic) {
        // ------------------ Uniform Refining -------------------
        UniformRefine(gmesh, h_ref_level);
        if (run_type == EDividedPyramid || run_type == EDividedPyramidIncreasedOrder ) {
            CreateFlattenGeometry(*gmesh);
            FlipPyramids(gmesh);
        }
    }
    
#ifdef PZDEBUG
            if (!gmesh)
            {
                std::cout << "The mesh was not generated." << std::endl;
                DebugStop();
            }
#endif

//    {
//        std::stringstream vols_name;
//        vols_name << "geometry_vols.txt";
//        PrintGeometryVols(gmesh, vols_name);
//    }

    // Computing mehs characteristics
    ComputeCharacteristicHElSize(gmesh, h_min, n_elements);
    
    // ------------------ Dividing pyramids in tets -------------------
    if(run_type == EDividedPyramid || run_type == EDividedPyramidIncreasedOrder ||
       run_type == EDividedPyramid4 || run_type == EDividedPyramidIncreasedOrder4)
    {
        DividePyramids(*gmesh);
        DivideBoundaryElements(*gmesh);
    }
    
#ifdef PZDEBUG
    // ------------------ Generating VTK with GMesh -------------------
    {
        std::string geoMeshName = "geometry.vtk";
        std::ofstream outPara(geoMeshName);
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outPara, true);
    }
#endif
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    return gmesh;
}

void CreateFlattenGeometry(TPZGeoMesh & gmesh){
    
    // Create a flatten mesh
    TPZGeoMesh geo_mesh;
    geo_mesh.SetName("Extacted mesh");
    geo_mesh.SetDimension(gmesh.Dimension());
    geo_mesh.NodeVec().Resize(gmesh.NNodes());
    std::map<int64_t,int64_t> node_map;
    ino64_t counter_node = 0;
    for (auto inode : gmesh.NodeVec()) {
        geo_mesh.NodeVec()[counter_node] = inode;
        node_map[inode.Id()] = counter_node;
        counter_node++;
    }
    TPZManVector <int64_t,8> node_indexes;
    int matid;
    MElementType type;
    int64_t index;
    for (auto iel : gmesh.ElementVec()) {
        if(!iel){
            continue;
        }
        if(iel->HasSubElement()){
            continue;
        }
        iel->GetNodeIndices(node_indexes);
        
        int n_nodes = node_indexes.size();
        for (int i = 0; i < n_nodes; i++) {
            node_indexes[i] = node_map[node_indexes[i]];
        }
        type = iel->Type();
        matid = iel->MaterialId();
        geo_mesh.CreateGeoElement(type, node_indexes, matid, index);
    }
    
    geo_mesh.BuildConnectivity();
    gmesh  = geo_mesh;
}

void FlipPyramids(TPZGeoMesh * gmesh){
    
    // Pyramids With Common Quadrilateral Faces
    int geo_dim = gmesh->Dimension();
    int n_el = gmesh->NElements();
    std::vector<int> pyramid_indexes;
    for (unsigned int iel = 0; iel < n_el; iel++) {
        TPZGeoEl * gel = gmesh->Element(iel);
        if (!gel || gel->Dimension() != geo_dim || gel->HasSubElement() || gel->Type() != EPiramide) {
            continue;
        }
        pyramid_indexes.push_back(gel->Index());
    }
    
    std::map<int, int> left_right_pyramid_map;
//    int n_pyramids = pyramid_indexes.size();
    for (auto ipyramid: pyramid_indexes) {
        TPZGeoEl * gel = gmesh->Element(ipyramid);
        TPZGeoElSide gel_side(gel,13);
        TPZStack<TPZGeoElSide> allneigh;
        gel_side.AllNeighbours(allneigh);
        int n_neighs = allneigh.size();
        if(n_neighs != 1){
            continue;
        }
        TPZGeoEl * gel_neigh = allneigh[0].Element();
        if(gel_neigh->Type()!=EPiramide){
            continue;
        }
        TPZManVector<int64_t,5> left_node_indexes,right_node_indexes;
        gel->GetNodeIndices(left_node_indexes);
        gel_neigh->GetNodeIndices(right_node_indexes);
        
        bool is_candidate_Q = left_node_indexes[0] != right_node_indexes[0] && left_node_indexes[2] != right_node_indexes[2];
        if (is_candidate_Q) {
            left_right_pyramid_map[gel->Index()] = gel_neigh->Index();
        }
    
    }
    
//    int n_pairs = left_right_pyramid_map.size();
    // Flip Right Pyramids
    for (auto ipair: left_right_pyramid_map) {
        TPZGeoEl * gel_left     = gmesh->Element(ipair.first);
        TPZGeoEl * gel_right    = gmesh->Element(ipair.second);
        TPZManVector<int64_t,5> left_node_indexes,right_node_indexes;
        gel_left->GetNodeIndices(left_node_indexes);
        gel_right->GetNodeIndices(right_node_indexes);
//        std::vector<int> perm;
//        perm.push_back(0);
//        perm.push_back(3);
//        perm.push_back(2);
//        perm.push_back(1);
        for (int i = 0; i < 4; i++) {
            right_node_indexes[i] = left_node_indexes[i];
            gel_right->SetNodeIndex(i, left_node_indexes[i]);
        }
        
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
    
}

void ApproximationError(int nref, int porder, TPZVec<STATE> &errors, bool hdivmm, TSimulationControl * control);

int ConvergenceTest(TSimulationControl * control)
{
    string projectpath = "/Projects/PyramidHdivTests/";
    
#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + projectpath;
    FileName += "pyramlogfile.cfg";
    InitializePZLOG();
#endif
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream str;
        str << "\nRodando testes de piramede Hdiv" << std::endl;
        LOGPZ_DEBUG(logger,str.str())
    }
#endif
    
    TPZManVector<STATE,3> errors(3,0.);
    int nref = 0;
    int porder = 1;
    bool hdivmm = true;
    ApproximationError(nref, porder, errors,hdivmm, control);
    
    return 0;
}

void PrintGeometryVols(TPZGeoMesh * gmesh, std::stringstream & file_name){
    
    std::ofstream mesh_data(file_name.str().c_str(),std::ios::app);
    
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    TPZStack<int64_t> gel_indexes;
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = gmesh->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->Dimension() != dim || gel->HasSubElement()) {
            continue;
        }
        gel_indexes.Push(gel->Index());
    }
    
    // Writing element nodes
    int64_t n_vols = gel_indexes.size();
    mesh_data << std::setw(5) << n_vols << std::endl;
    TPZFMatrix<REAL> xcoor;
    
    for (int64_t ivol = 0; ivol < n_vols; ivol++) {
        TPZGeoEl * gel = gmesh->Element(gel_indexes[ivol]);
        gel->NodesCoordinates(xcoor);
        int nr = xcoor.Rows();
        int nc = xcoor.Cols();
        for (int c = 0; c < nc; c++) {
            for (int r = 0; r < nr; r++) {
                mesh_data << std::setw(20) << xcoor(r,c) ;
            }
        }
        mesh_data << std::endl;
    }
    mesh_data.flush();
    return;
}

void ApproximationError(int nref, int porder, TPZVec<STATE> &errors, bool hdivmm, TSimulationControl * control)
{
    //   gRefDBase.InitializeAllUniformRefPatterns();
    HDivPiola = 1;
    TPZGeoMesh *gmesh = CreateGeoMesh1Pir();
    //    TPZGeoMesh *gmesh = CreateGeoMeshHexaOfPir();
    //    TPZGeoMesh *gmesh = CreateGeoMeshHexaOfPirTetra();
    //    TPZGeoMesh *gmesh = CreateGeoMesh1Tet();
    //    TPZGeoMesh *gmesh = CreateGeoMeshPrism();
    
    UniformRefine(gmesh, nref);
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled() && nref < 2)
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZManVector<TPZCompMesh*,2> meshvec(2);
    meshvec[1] = CreateCmeshPressure(gmesh, control, porder);
    LoadSolution(meshvec[1]);
    meshvec[0] = CreateCmeshFlux(gmesh, control, porder);
    TPZCompMeshTools::AddHDivPyramidRestraints(meshvec[0]);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled() && nref < 2)
    {
        std::stringstream sout;
        meshvec[0]->Print(sout);
        meshvec[1]->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZCompMesh *cmeshMult = CreateCmeshMulti(meshvec, control);
    TPZMaterial *mat = cmeshMult->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    mat->SetForcingFunction(LaplaceExact, porder);
//    GroupElements(cmeshMult);

    
    TPZAnalysis an(cmeshMult,false);
    an.SetExact(ExactSolution);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    TPZSkylineStructMatrix skyl(cmeshMult);
    skyl.SetNumThreads(0);
    //    TPZFStructMatrix skyl(cmeshMult);
    an.SetStructuralMatrix(skyl);
    an.SetSolver(step);
    
    an.Run();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled() && nref < 2)
    {
        std::stringstream sout;
        cmeshMult->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    std::ofstream out("AccumErrors.txt",ios::app);
	bool storeElements = false;
    an.PostProcessError(errors, storeElements, std::cout);
    
    out << "nref " << nref <<  " h " << 1./(2<<nref) << " porder " << porder << " hdivmm " << hdivmm <<  " neq " << cmeshMult->NEquations() << " errors " << errors << std::endl;
}

TPZGeoMesh * CreateGeoMesh1Pir()
{
    int64_t nodeids[]={0,1,3,2,4};
    REAL coords[5][3]={{-1.,-1.,0},{1.,-1.,0.},{1.,1.,0.},{-1.,1.,0.},{0.,0.,1.}};
    //    REAL coords[5][3]={{-1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.},{-1.,-1.,1.},{1.,-1.,-1.}};
    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    
    // Setando os nohs
    int nnodes = 5;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    int64_t index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 1
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 2
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 3
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 4
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // Criando elemento
    TPZManVector<int64_t,5> topolPyr(5);
    for (int i = 0; i < 5; i++) {
        topolPyr[i] = i;
    }
    
    TPZGeoEl *gel = gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index,0);
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    gel->CreateBCGeoEl(13, bc0); // fundo
    gel->CreateBCGeoEl(14, bc0); // frente (-1,-1 ateh 1,-1)
    gel->CreateBCGeoEl(15, bc0); // direita
    gel->CreateBCGeoEl(16, bc0); // atras
    gel->CreateBCGeoEl(17, bc0); // esquerda
    
    gmesh->BuildConnectivity();
    
    std::ofstream out("1PyrGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
    
}

TPZGeoMesh * CreateGeoMeshPrism()
{
    //    int64_t nodeids[]={0,1,2,3,4};
    int64_t nodeids[]={3,0,1,2,4};
    int64_t nodeidstet[] = {2,3,4,5};
    int64_t boundaryids[7][4] = {{0,4,3},{3,4,5},{0,4,1},{3,5,2},{4,1,2},{4,2,5},{0,1,2,3}};
    //    REAL coords[5][3]={{-1.,-1.,0},{1.,-1.,0.},{1.,1.,0.},{-1.,1.,0.},{0.,0.,1.}};
    REAL coords[6][3]={{-1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.},{-1.,-1.,1.},{1.,-1.,-1.},{1.,-1.,1.}};
    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    
    // Setando os nohs
    int nnodes = 6;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    int64_t index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 1
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 2
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 3
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 4
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 5
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // Criando elemento
    TPZManVector<int64_t,5> topolPyr(5);
    for (int i = 0; i < 5; i++) {
        topolPyr[i] = nodeids[i];
    }
    
    TPZGeoEl *gel = gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index);
    
    TPZManVector<int64_t,4> topolTet(4);
    for (int i=0; i<4; i++) {
        topolTet[i] = nodeidstet[i];
    }
    gel = gmesh->CreateGeoElement(ETetraedro, topolTet, matid, index);
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    for (int el=0; el<6; el++) {
        TPZManVector<int64_t,3> nodes(3);
        for (int i=0; i<3; i++) {
            nodes[i] = boundaryids[el][i];
        }
        int64_t index;
        gmesh->CreateGeoElement(ETriangle, nodes, bc0, index);
    }
    TPZManVector<int64_t,4> nodes(4);
    for (int i=0; i<4; i++) {
        nodes[i] = boundaryids[6][i];
    }
    gmesh->CreateGeoElement(EQuadrilateral, nodes, bc0, index);
    
    gmesh->BuildConnectivity();
    
    std::ofstream out("PrismGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
    
}


TPZGeoMesh * CreateGeoMesh1Tet()
{
    int64_t nodeids[]={0,1,2,3};
    //    REAL coords[4][3]={{0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
    REAL coords[4][3]={{-1.,-1.,1.},{1.,-1.,1.},{1.,1.,1.},{1.,-1.,-1.}};
    
    
    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    
    // Setando os nohs
    int nnodes = 5;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    int64_t index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 1
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 2
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 3
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    
    // Criando elemento
    TPZManVector<int64_t,5> topolTet(4);
    for (int i = 0; i < 4; i++) {
        topolTet[i] = i;
    }
    
    TPZGeoEl *gel = gmesh->CreateGeoElement(ETetraedro, topolTet, matid, index);
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    gel->CreateBCGeoEl(10, bc0); // fundo
    gel->CreateBCGeoEl(11, bc0); // frente (-1,-1 ateh 1,-1)
    gel->CreateBCGeoEl(12, bc0); // direita
    gel->CreateBCGeoEl(13, bc0); // atras
    
    gmesh->BuildConnectivity();
    
    std::ofstream out("../1TetGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
    
}


TPZGeoMesh * CreateGeoMeshHexaOfPir()
{
    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    
    // Setando os nohs
    int nnodes = 9;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    int64_t index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    nodecoord[0] = -1.;
    nodecoord[1] = -1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 1
    nodecoord[0] = 1.;
    nodecoord[1] = -1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 2
    nodecoord[0] = 1.;
    nodecoord[1] = 1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 3
    nodecoord[0] = -1.;
    nodecoord[1] = 1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 4
    nodecoord[0] = -1.;
    nodecoord[1] = -1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 5
    nodecoord[0] = 1.;
    nodecoord[1] = -1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 6
    nodecoord[0] = 1.;
    nodecoord[1] = 1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 7
    nodecoord[0] = -1.;
    nodecoord[1] = 1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 8
    nodecoord[0] = 0.;
    nodecoord[1] = 0.;
    nodecoord[2] = 0.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    
    // Criando elemento
    TPZManVector<int64_t,5> topolPyr(5);
    int myels[6][5] = {{0,1,5,4,8},{1,2,6,5,8},{2,3,7,6,8},{0,3,7,4,8},{0,1,2,3,8},{4,5,6,7,8}};
    //int myels[6][5] = {{0,1,5,4,8},{6,5,1,2,8},{2,3,7,6,8},{7,4,0,3,8},{0,1,2,3,8},{4,5,6,7,8}}; //Sequencia trocada soh para funcionar o AddHDivPyramidRestraints
    for (int iel = 0; iel < 6; iel++) {
        for (int i = 0; i < 5; i++) {
            topolPyr[i] = myels[iel][i];
        }
        gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index);
    }
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    
    const int64_t nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        gmesh->Element(iel)->CreateBCGeoEl(13, bc0);
    }
    
    gmesh->BuildConnectivity();
    
    std::ofstream out("HexaPyrGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
}

TPZGeoMesh * CreateGeoMeshHexaOfPirTetra()
{
    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    
    // Setando os nohs
    int nnodes = 9;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    int64_t index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    nodecoord[0] = -1.;
    nodecoord[1] = -1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 1
    nodecoord[0] = 1.;
    nodecoord[1] = -1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 2
    nodecoord[0] = 1.;
    nodecoord[1] = 1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 3
    nodecoord[0] = -1.;
    nodecoord[1] = 1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 4
    nodecoord[0] = -1.;
    nodecoord[1] = -1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 5
    nodecoord[0] = 1.;
    nodecoord[1] = -1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 6
    nodecoord[0] = 1.;
    nodecoord[1] = 1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 7
    nodecoord[0] = -1.;
    nodecoord[1] = 1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 8
    nodecoord[0] = 0.;
    nodecoord[1] = 0.;
    nodecoord[2] = 0.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    //    ino++;
    gmesh->SetNodeIdUsed(ino);
    
    // Criando elemento
    TPZManVector<int64_t,5> topolPyr(5), topolTet(4), topolTri(3);
    int myelsp[2][5] = {{4,0,2,6,1},{4,0,2,6,7}};
    int myelst[2][4] = {{4,6,5,1},{0,2,3,7}};
    //                          front           right          top             back            left            bottom
    int triangles[12][3] = {{0,1,4},{1,5,4},{1,2,6},{1,6,5},{4,5,6},{4,6,7},{2,6,7},{2,7,3},{0,3,7},{0,7,4},{0,1,2},{0,2,3} };
    //int myels[6][5] = {{0,1,5,4,8},{6,5,1,2,8},{2,3,7,6,8},{7,4,0,3,8},{0,1,2,3,8},{4,5,6,7,8}}; //Sequencia trocada soh para funcionar o AddHDivPyramidRestraints
    for (int iel = 0; iel < 2; iel++) {
        for (int i = 0; i < 5; i++) {
            topolPyr[i] = myelsp[iel][i];
        }
        gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index,0);
    }
    for (int iel = 0; iel < 2; iel++) {
        for (int i = 0; i < 4; i++) {
            topolTet[i] = myelst[iel][i];
        }
        gmesh->CreateGeoElement(ETetraedro, topolTet, matid, index,0);
    }
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    
    for (int64_t iel = 0; iel < 12; iel++) {
        for (int i = 0; i < 3; i++) {
            topolTri[i] = triangles[iel][i];
        }
        gmesh->CreateGeoElement(ETriangle, topolTri, bc0, index,0);
    }
    
    gmesh->BuildConnectivity();
    
    std::ofstream out("../HexaPyrTetGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
}


TPZCompMesh * CreateCmeshPressure(TPZGeoMesh *gmesh, TSimulationControl * control, int p)
{
    const int matid = 1;
    const int dim = 3;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    if (control->m_Hdiv_plusplus_Q) {
        p = p+1;
    }
    cmesh->SetDefaultOrder(p);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    int nstate =1;
    TPZVec<STATE> sol(1);
    TPZL2Projection * mymat = new TPZL2Projection(matid,dim,nstate,sol);
    mymat->SetForcingFunction(AnalyticPotential, gIntegrationOrder);
    cmesh->InsertMaterialObject(mymat);
    cmesh->AutoBuild();
    
    int64_t ncon = cmesh->NConnects();
    for (int64_t ic=0; ic<ncon; ic++) {
        cmesh->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    cmesh->InitializeBlock();
    return cmesh;
}

TPZCompMesh * CreateCmeshFlux(TPZGeoMesh *gmesh, TSimulationControl * control, int p)
{
    const int matid = 1;
    const int dim = 3;
    const int dirichlet = 0;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(p);
    
    TPZVecL2 *mymat = new TPZVecL2(matid);
    mymat->SetDimension(dim);
    mymat->SetForcingFunction(AnalyticFlux, gIntegrationOrder);
    cmesh->InsertMaterialObject(mymat);
    
    TPZFMatrix<> val1(3,3,0.);
    TPZFMatrix<> val2(3,1,0.);
    
    switch (control->m_geometry_type) {
        case EAcademic:
        {
            int bc0 = -1;
            int matbcfake = 3;
            TPZBndCond *bnd = mymat->CreateBC(mymat, bc0, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bnd);
            
            TPZBndCond *bnd3 = mymat->CreateBC(mymat, matbcfake, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bnd3);
            
        }
            break;
        case EVerticalWellHe:
        {
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            int bc_impervious_id = 4;
            TPZBndCond *bc_outer = mymat->CreateBC(mymat, bc_outer_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_outer);
            
            TPZBndCond *bc_inner = mymat->CreateBC(mymat, bc_inner_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_inner);
            
            TPZBndCond *bc_impervious = mymat->CreateBC(mymat, bc_impervious_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_impervious);
            
        }
            break;
        case EVerticalWellTe:
        {
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            int bc_impervious_id = 4;
            TPZBndCond *bc_outer = mymat->CreateBC(mymat, bc_outer_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_outer);
            
            TPZBndCond *bc_inner = mymat->CreateBC(mymat, bc_inner_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_inner);
            
            TPZBndCond *bc_impervious = mymat->CreateBC(mymat, bc_impervious_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_impervious);
            
        }
            break;
        case EVerticalWellHePyTe:
        {
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            int bc_impervious_id = 4;
            TPZBndCond *bc_outer = mymat->CreateBC(mymat, bc_outer_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_outer);
            
            TPZBndCond *bc_inner = mymat->CreateBC(mymat, bc_inner_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_inner);
            
            TPZBndCond *bc_impervious = mymat->CreateBC(mymat, bc_impervious_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_impervious);
            
        }
            break;
            
        case ESphericalBarrierHePyTe:
        {
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            TPZBndCond *bc_outer = mymat->CreateBC(mymat, bc_outer_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_outer);
            
            TPZBndCond *bc_inner = mymat->CreateBC(mymat, bc_inner_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_inner);
            
        }
            break;
        case ESphericalBarrierTe:
        {
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            TPZBndCond *bc_outer = mymat->CreateBC(mymat, bc_outer_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_outer);
            
            TPZBndCond *bc_inner = mymat->CreateBC(mymat, bc_inner_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_inner);
            
        }
            break;
        case ESphericalBarrierHe:
        {
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            TPZBndCond *bc_outer = mymat->CreateBC(mymat, bc_outer_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_outer);
            
            TPZBndCond *bc_inner = mymat->CreateBC(mymat, bc_inner_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc_inner);
            
        }
            break;
        default:{
            DebugStop();
        }
            break;
    }
    
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    
    if (control->m_Hdiv_plusplus_Q)
    {
        cmesh->SetDefaultOrder(p+1);
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if(!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            int nc = intel->NConnects();
            if (nc <=1) {
                continue;
            }
            TPZGeoEl *gel = intel->Reference();
            int ns = gel->NSides();
            int porder = p;
            porder++;
            intel->ForceSideOrder(ns-1, porder);
        }
    }
    cmesh->ExpandSolution();
    return cmesh;
}

TPZCompMesh * CreateCmeshMulti(TPZVec<TPZCompMesh *> &meshvec, TSimulationControl * control)
{
    const int int_p_order = 10;
    //Creating computational mesh for multiphysic elements
    TPZGeoMesh *gmesh = meshvec[0]->Reference();
    int dim = gmesh->Dimension();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int p1 = meshvec[0]->GetDefaultOrder();
    int p2 = meshvec[1]->GetDefaultOrder();
    int p = p1 < p2 ? p2 : p1;
    
    switch (control->m_geometry_type) {
        case EAcademic:
        {
            const int matid = 1, bc0 = -1;
            const int matbcfake = 3;
            const int dirichlet = 0;
            TPZDualPoisson * mat = new TPZDualPoisson(matid);
            {
                TPZDummyFunction<STATE> *force = new TPZDummyFunction<STATE>(BodyForcing,int_p_order);
                force->SetPolynomialOrder(gIntegrationOrder);
                TPZAutoPointer<TPZFunction<STATE> > bodyforce = force;
                
                mat->SetForcingFunction(bodyforce);
            }
            {
                TPZDummyFunction<STATE> *analytic_f = new TPZDummyFunction<STATE>(Analytic,int_p_order);
                analytic_f->SetPolynomialOrder(gIntegrationOrder);
                TPZAutoPointer<TPZFunction<STATE> > analytic = analytic_f;
                mat->SetForcingFunctionExact(analytic);
            }
            //inserindo o material na malha computacional
            mphysics->InsertMaterialObject(mat);
            
            //Criando condicoes de contorno
            TPZFMatrix<> val1(3,3,0.), val2(3,1,0.);
            TPZBndCond * BCond0 = NULL;
            BCond0 = mat->CreateBC(mat, bc0,dirichlet, val1, val2);
            {
                TPZDummyFunction<STATE> *boundforce = new TPZDummyFunction<STATE>(Forcing,int_p_order);
                boundforce->SetPolynomialOrder(gIntegrationOrder);
                TPZAutoPointer<TPZFunction<STATE> > force = boundforce;
                BCond0->SetForcingFunction(0,force);
            }
            mphysics->InsertMaterialObject(BCond0);
            
            // a zero contribution element
            BCond0 = mat->CreateBC(mat, matbcfake,dirichlet, val1, val2);
            mphysics->InsertMaterialObject(BCond0);
            
        }
            break;
        case EVerticalWellHe:
        {
            const int matid = 1;
            TPZDualPoisson * mat = new TPZDualPoisson(matid);
            mphysics->InsertMaterialObject(mat);
            
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            int bc_impervious_id = 4;
            int dirichlet = 0;
            
            TPZDummyFunction<STATE> *boundforce = new TPZDummyFunction<STATE>(Forcing,int_p_order);
            boundforce->SetPolynomialOrder(gIntegrationOrder);
            TPZAutoPointer<TPZFunction<STATE> > force = boundforce;
            
            TPZFMatrix<> val1(3,3,0.), val2(3,1,0.);
            TPZBndCond * bc_outer = mat->CreateBC(mat, bc_outer_id ,dirichlet, val1, val2);
            bc_outer->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_outer);
            TPZBndCond * bc_inner = mat->CreateBC(mat, bc_inner_id ,dirichlet, val1, val2);
            bc_inner->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_inner);
            TPZBndCond * bc_impervious = mat->CreateBC(mat, bc_impervious_id ,dirichlet, val1, val2);
            bc_impervious->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_impervious);
            
        }
            break;
        case EVerticalWellTe:
        {
            const int matid = 1;
            TPZDualPoisson * mat = new TPZDualPoisson(matid);
            mphysics->InsertMaterialObject(mat);
            
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            int bc_impervious_id = 4;
            int dirichlet = 0;
            
            TPZDummyFunction<STATE> *boundforce = new TPZDummyFunction<STATE>(Forcing,int_p_order);
            boundforce->SetPolynomialOrder(gIntegrationOrder);
            TPZAutoPointer<TPZFunction<STATE> > force = boundforce;
            
            TPZFMatrix<> val1(3,3,0.), val2(3,1,0.);
            TPZBndCond * bc_outer = mat->CreateBC(mat, bc_outer_id ,dirichlet, val1, val2);
            bc_outer->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_outer);
            TPZBndCond * bc_inner = mat->CreateBC(mat, bc_inner_id ,dirichlet, val1, val2);
            bc_inner->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_inner);
            TPZBndCond * bc_impervious = mat->CreateBC(mat, bc_impervious_id ,dirichlet, val1, val2);
            bc_impervious->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_impervious);
            
        }
            break;
        case EVerticalWellHePyTe:
        {
            const int matid = 1;
            TPZDualPoisson * mat = new TPZDualPoisson(matid);
            mphysics->InsertMaterialObject(mat);
            
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            int bc_impervious_id = 4;
            int dirichlet = 0;
            
            TPZDummyFunction<STATE> *boundforce = new TPZDummyFunction<STATE>(Forcing,int_p_order);
            boundforce->SetPolynomialOrder(gIntegrationOrder);
            TPZAutoPointer<TPZFunction<STATE> > force = boundforce;
            
            TPZFMatrix<> val1(3,3,0.), val2(3,1,0.);
            TPZBndCond * bc_outer = mat->CreateBC(mat, bc_outer_id ,dirichlet, val1, val2);
            bc_outer->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_outer);
            TPZBndCond * bc_inner = mat->CreateBC(mat, bc_inner_id ,dirichlet, val1, val2);
            bc_inner->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_inner);
            TPZBndCond * bc_impervious = mat->CreateBC(mat, bc_impervious_id ,dirichlet, val1, val2);
            bc_impervious->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_impervious);
            
        }
            break;
        case ESphericalBarrierHePyTe:
        {
            const int matid = 1;
            TPZDualPoisson * mat = new TPZDualPoisson(matid);
            mphysics->InsertMaterialObject(mat);
            
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            int dirichlet = 0;
            
            TPZDummyFunction<STATE> *boundforce = new TPZDummyFunction<STATE>(Forcing,int_p_order);
            boundforce->SetPolynomialOrder(gIntegrationOrder);
            TPZAutoPointer<TPZFunction<STATE> > force = boundforce;
            
            TPZFMatrix<> val1(3,3,0.), val2(3,1,0.);
            TPZBndCond * bc_outer = mat->CreateBC(mat, bc_outer_id ,dirichlet, val1, val2);
            bc_outer->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_outer);
            TPZBndCond * bc_inner = mat->CreateBC(mat, bc_inner_id ,dirichlet, val1, val2);
            bc_inner->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_inner);
            
        }
            break;
        case ESphericalBarrierTe:
        {
            const int matid = 1;
            TPZDualPoisson * mat = new TPZDualPoisson(matid);
            mphysics->InsertMaterialObject(mat);
            
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            int dirichlet = 0;
            
            TPZDummyFunction<STATE> *boundforce = new TPZDummyFunction<STATE>(Forcing,int_p_order);
            boundforce->SetPolynomialOrder(gIntegrationOrder);
            TPZAutoPointer<TPZFunction<STATE> > force = boundforce;
            
            TPZFMatrix<> val1(3,3,0.), val2(3,1,0.);
            TPZBndCond * bc_outer = mat->CreateBC(mat, bc_outer_id ,dirichlet, val1, val2);
            bc_outer->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_outer);
            TPZBndCond * bc_inner = mat->CreateBC(mat, bc_inner_id ,dirichlet, val1, val2);
            bc_inner->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_inner);
            
        }
            break;
        case ESphericalBarrierHe:
        {
            const int matid = 1;
            TPZDualPoisson * mat = new TPZDualPoisson(matid);
            mphysics->InsertMaterialObject(mat);
            
            int bc_outer_id = 2;
            int bc_inner_id = 3;
            int dirichlet = 0;
            
            TPZDummyFunction<STATE> *boundforce = new TPZDummyFunction<STATE>(Forcing,int_p_order);
            boundforce->SetPolynomialOrder(gIntegrationOrder);
            TPZAutoPointer<TPZFunction<STATE> > force = boundforce;
            
            TPZFMatrix<> val1(3,3,0.), val2(3,1,0.);
            TPZBndCond * bc_outer = mat->CreateBC(mat, bc_outer_id ,dirichlet, val1, val2);
            bc_outer->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_outer);
            TPZBndCond * bc_inner = mat->CreateBC(mat, bc_inner_id ,dirichlet, val1, val2);
            bc_inner->SetForcingFunction(0,force);
            mphysics->InsertMaterialObject(bc_inner);
            
        }
            break;
        default:{
            DebugStop();
        }
            break;
    }
    
    mphysics->SetDefaultOrder(p);
    mphysics->SetDimModel(dim);
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->AutoBuild();
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        mphysics->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    return mphysics;
}

// ------------------------ Para testes do assemble -----------------------------

TPZGeoMesh *MalhaQuadrada(int &nelx, int &nely)
{
    const int bcid1 = -1, bcid2 = -2;
    TPZManVector<int,2> nelxy(2,0);
    nelxy[0] = nelx;
    nelxy[1] = nely;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.); // dois pontos nos 2 cantos do quadrado
    x1[0] = 1.;
    x1[1] = 1.;
    TPZGenGrid gengrid(nelxy,x0,x1);
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gengrid.Read(gmesh);
    gengrid.SetBC(gmesh,7,bcid1); // na direita
    gengrid.SetBC(gmesh,5,bcid2); // na esquerda
    
    gmesh->BuildConnectivity();
    return gmesh;
}

void InsertBidimensionalPoisson(TPZCompMesh *cmesh, int &dim)
{
    int matid = 1;
    const int bcidLeft = -1, bcidRight = -2;
    const int diri = 0;
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(matid,dim);
    cmesh->InsertMaterialObject(mat);
    
    TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
    TPZBndCond *bcLeft = mat->CreateBC(mat, bcidLeft, diri, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    val2(0,0) = 2.;
    TPZBndCond *bcRight = mat->CreateBC(mat, bcidRight, diri, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
}


TPZGeoMesh *MalhaCubo(string &projectpath, const int &nref)
{
    int64_t numnodes=-1;
    int64_t numelements=-1;
    
    string FileName, dirname = PZSOURCEDIR;
    FileName = dirname + projectpath;
    FileName += "cube1.txt";
    
    {
        bool countnodes = false;
        bool countelements = false;
        
        ifstream read (FileName.c_str());
        
        while(read)
        {
            char buf[1024];
            read.getline(buf, 1024);
            std::string str(buf);
            if(str == "Coordinates") countnodes = true;
            if(str == "end coordinates") countnodes = false;
            if(countnodes) numnodes++;
            
            if(str == "Elements") countelements = true;
            if(str == "end elements") countelements = false;
            if(countelements) numelements++;
        }
    }
    
    TPZGeoMesh * gMesh = new TPZGeoMesh;
    
    gMesh -> NodeVec().Resize(numnodes);
    
    TPZManVector <int64_t> TopolTetra(4);
    
    const int64_t Qnodes = numnodes;
    TPZVec <TPZGeoNode> Node(Qnodes);
    
    //setting nodes coords
    int64_t nodeId = 0, elementId = 0, matElId = 1;
    
    ifstream read;
    read.open(FileName.c_str());
    
    double nodecoordX , nodecoordY , nodecoordZ ;
    
    char buf[1024];
    read.getline(buf, 1024);
    read.getline(buf, 1024);
    std::string str(buf);
    int64_t in;
    for(in=0; in<numnodes; in++)
    {
        read >> nodeId;
        read >> nodecoordX;
        read >> nodecoordY;
        read >> nodecoordZ;
        Node[nodeId-1].SetNodeId(nodeId);
        Node[nodeId-1].SetCoord(0,nodecoordX);
        Node[nodeId-1].SetCoord(1,nodecoordY);
        Node[nodeId-1].SetCoord(2,nodecoordZ);
        gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
    }
    
    {
        read.close();
        read.open(FileName.c_str());
        
        int64_t l , m = numnodes+5;
        for(l=0; l<m; l++)
        {
            read.getline(buf, 1024);
        }
        
        
        int64_t el;
        int neumann1 = -4, neumann2 = -5;
        //std::set<int> ncoordz; //jeitoCaju
        for(el=0; el<numelements; el++)
        {
            read >> elementId;
            read >> TopolTetra[0]; //node 1
            read >> TopolTetra[1]; //node 2
            read >> TopolTetra[2]; //node 3
            read >> TopolTetra[3]; //node 4
            
            // O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
            TopolTetra[0]--;
            TopolTetra[1]--;
            TopolTetra[2]--;
            TopolTetra[3]--;
            
            int64_t index = el;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
        }
        
        gMesh->BuildConnectivity();
        
        // Colocando as condicoes de contorno
        for(el=0; el<numelements; el++)
        {
            TPZManVector <TPZGeoNode,4> Nodefinder(4);
            TPZManVector <REAL,3> nodecoord(3);
            TPZGeoEl *tetra = gMesh->ElementVec()[el];
            
            // na face x = 1
            TPZVec<int64_t> ncoordzVec(0); int64_t sizeOfVec = 0;
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gMesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (nodecoord[0] == 1.)
                {
                    sizeOfVec++;
                    ncoordzVec.Resize(sizeOfVec);
                    ncoordzVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordzVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordzVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,neumann1);
            }
            
            // Na face x = -1
            ncoordzVec.Resize(0);
            sizeOfVec = 0;
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gMesh->NodeVec()[pos];
                
                Nodefinder[i].GetCoordinates(nodecoord);
                if (nodecoord[0] == -1.)
                {
                    sizeOfVec++;
                    ncoordzVec.Resize(sizeOfVec);
                    ncoordzVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordzVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordzVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,neumann2);
            }
            
        }
        
        TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
        yz[0] = 1.;
        z[2] = -1;
        int bcidxyz = -1, bcidyz = -2, bcidz = -3;
        SetPointBC(gMesh, xyz, bcidxyz);
        SetPointBC(gMesh, yz, bcidyz);
        SetPointBC(gMesh, z, bcidz);
    }
    
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        const int nel = gMesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gMesh->Element(iel);
            if (gel || !gel->HasSubElement() || gel->Dimension() > 0) {
                gel->Divide(sons);
            }
        }
    }
    
    return gMesh;
}

/// Generate a boundary geometric element at the indicated node
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc)
{
    // look for an element/corner node whose distance is close to start
    TPZGeoNode *gn1 = gr->FindNode(x);
    int64_t iel;
    int64_t nelem = gr->ElementVec().NElements();
    TPZGeoEl *gel;
    for (iel = 0; iel<nelem; iel++) {
        gel = gr->ElementVec()[iel];
        if(!gel) continue;
        int nc = gel->NCornerNodes();
        int c;
        for (c=0; c<nc; c++) {
            TPZGeoNode *gn = gel->NodePtr(c);
            if (gn == gn1) {
                break;
            }
        }
        if (c<nc) {
            TPZGeoElBC(gel, c, bc);
            return;
        }
    }
}

void InsertElasticityCubo(TPZCompMesh *mesh)
{
    mesh->SetDimModel(3);
    int nummat = 1, neumann = 1, mixed = 2;
    //	int dirichlet = 0;
    int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5;   //, dirp2 = -6;
    TPZManVector<STATE> force(3,0.);
    //force[1] = 0.;
    
    STATE ElaE = 1000., poissonE = 0.2;   //, poissonV = 0.1, ElaV = 100.;
    
    STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
    lambdaV = 11.3636;
    muV = 45.4545;
    alpha = 1.;
    deltaT = 0.01;
    
    //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
    //viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
    //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, ElaE, poissonE, lambdaV, muV, alphaT, force);
    TPZElasticity3D *viscoelast = new TPZElasticity3D(nummat, ElaE, poissonE, force);
    
    TPZFNMatrix<6> qsi(6,1,0.);
    //viscoelast->SetDefaultMem(qsi); //elast
    //int index = viscoelast->PushMemItem(); //elast
    TPZMaterial * viscoelastauto(viscoelast);
    mesh->InsertMaterialObject(viscoelastauto);
    
    // Neumann em x = 1;
    TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
    val2(0,0) = 1.;
    TPZBndCond *bc4 = viscoelast->CreateBC(viscoelastauto, neumann1, neumann, val1, val2);
    TPZMaterial * bcauto4(bc4);
    mesh->InsertMaterialObject(bcauto4);
    
    // Neumann em x = -1;
    val2(0,0) = -1.;
    TPZBndCond *bc5 = viscoelast->CreateBC(viscoelastauto, neumann2, neumann, val1, val2);
    TPZMaterial * bcauto5(bc5);
    mesh->InsertMaterialObject(bcauto5);
    
    val2.Zero();
    // Dirichlet em -1 -1 -1 xyz;
    val1(0,0) = 1e4;
    val1(1,1) = 1e4;
    val1(2,2) = 1e4;
    TPZBndCond *bc1 = viscoelast->CreateBC(viscoelastauto, dir1, mixed, val1, val2);
    TPZMaterial * bcauto1(bc1);
    mesh->InsertMaterialObject(bcauto1);
    
    // Dirichlet em 1 -1 -1 yz;
    val1(0,0) = 0.;
    val1(1,1) = 1e4;
    val1(2,2) = 1e4;
    TPZBndCond *bc2 = viscoelast->CreateBC(viscoelastauto, dir2, mixed, val1, val2);
    TPZMaterial * bcauto2(bc2);
    mesh->InsertMaterialObject(bcauto2);
    
    // Dirichlet em 1 1 -1 z;
    val1(0,0) = 0.;
    val1(1,1) = 0.;
    val1(2,2) = 1e4;
    TPZBndCond *bc3 = viscoelast->CreateBC(viscoelastauto, dir3, mixed, val1, val2);
    TPZMaterial * bcauto3(bc3);
    mesh->InsertMaterialObject(bcauto3);
}

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
//    int geometry_dim = gmesh->Dimension();
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
//            if (!gel || gel->HasSubElement()) {
//                continue;
//            }
            gel->Divide(filhos);
        }
    }
    
////    gmesh->Ad
//    if(nDiv != 0){
//        int mat_volumetric_id = 1;
//        DivideBoundaryElements(*gmesh,mat_volumetric_id);
//    }
}

void GroupElements(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Type() == EPiramide) {
            std::list<TPZOneShapeRestraint> shape;
            shape = cel->GetShapeRestraints();
            if (shape.size() != 1) {
                DebugStop();
            }
            TPZOneShapeRestraint local = *shape.begin();
            TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!mult) {
                DebugStop();
            }
            TPZCompElHDiv<pzshape::TPZShapePiram> *hdiv = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapePiram> *>(mult->Element(0));
            if (!hdiv) {
                DebugStop();
            }
            int face = hdiv->RestrainedFace();
            TPZGeoElSide gelside(gel,face);
            TPZGeoElSide neighbour = gelside.Neighbour();
            TPZCompEl *celneigh = neighbour.Element()->Reference();
            TPZMultiphysicsElement *multneigh = dynamic_cast<TPZMultiphysicsElement *>(celneigh);
            if (!multneigh) {
                DebugStop();
            }
            TPZCompElHDiv<pzshape::TPZShapeTetra> *hdivneigh = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTetra> *>(multneigh->Element(0));
            if (!hdivneigh) {
                DebugStop();
            }
            int64_t index;
            TPZElementGroup *grp = new TPZElementGroup(*cmesh,index);
            grp->AddElement(cel);
            grp->AddElement(celneigh);
        }
    }
}

void LoadSolution(TPZCompMesh *cpressure)
{
    int64_t nel = cpressure->NElements();
    for (int64_t iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cpressure->Element(iel);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != 3) {
            continue;
        }
        if (gel->Type() == EPiramide) {
            TPZGeoNode *top = gel->NodePtr(4);
            TPZManVector<REAL,3> topco(3),valvec(1);
            top->GetCoordinates(topco);
            Forcing(topco, valvec);
            STATE topval = valvec[0];
            for (int i=0; i<4; i++) {
                TPZConnect &c = cel->Connect(i);
                int64_t seqnum = c.SequenceNumber();
                cpressure->Block()(seqnum,0,1,0) = topval;
                TPZGeoNode *no = gel->NodePtr(i);
                no->GetCoordinates(topco);
                Forcing(topco, valvec);
                STATE nodeval = valvec[0];
                cpressure->Block()(seqnum,0,0,0) = nodeval-topval;
            }
        }
        else if(gel->Type() == ETetraedro)
        {
            for (int i=0; i<4; i++) {
                TPZConnect &c = cel->Connect(i);
                TPZGeoNode *no = gel->NodePtr(i);
                TPZManVector<REAL,3> topco(3);
                no->GetCoordinates(topco);
                TPZManVector<STATE,3> valvec(1);
                Forcing(topco, valvec);
                STATE nodeval = valvec[0];
                int64_t seqnum = c.SequenceNumber();
                cpressure->Block()(seqnum,0,0,0) = nodeval;
            }
            
        }
    }
}

void ProjectFlux(TPZCompMesh *cfluxmesh)
{
    cout << "\n**** Projecting Flux ******" << endl;
    TPZAnalysis an(cfluxmesh,false);
    TPZSkylineStructMatrix str(cfluxmesh);
    an.SetStructuralMatrix(str);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();
}

static int gfluxorder = 3;

/// Generate the L2 matrix of the pressure space and the inner product of the divergence and the pressure shape functions
static void GenerateProjectionMatrix(TPZCompEl *cel, TPZAutoPointer<TPZMatrix<STATE> > L2, TPZFMatrix<STATE> &inner);

/// Given the multiplier coefficients of the pressure space, verify the correspondence of the divergence of the vector function and the L2 projection
static int VerifyProjection(TPZCompEl *cel, TPZFMatrix<STATE> &multiplier);

/// verify if the divergence of each vector function is included in the pressure space
static void CheckDRham(TPZCompEl *cel)
{
    TPZFMatrix<STATE> inner, multiplier;
    TPZAutoPointer<TPZMatrix<STATE> > L2 = new TPZFMatrix<STATE>;
    GenerateProjectionMatrix(cel, L2, inner);
    int porder = cel->GetgOrder();
    std::string filename;
    {
        std::stringstream sout;
        sout << "../matrices" << gfluxorder << ".nb";
        filename = sout.str();
    }
    
    std::ofstream output(filename.c_str());
    output.precision(16);
    {
        std::stringstream sout;
        sout << "L2" << gfluxorder << " = ";
        filename = sout.str();
    }
    L2->Print(filename.c_str(),output, EMathematicaInput);
    {
        std::stringstream sout;
        sout << "PressHDiv" << gfluxorder << " = ";
        filename = sout.str();
    }
    inner.Print(filename.c_str(),output,EMathematicaInput);
    TPZStepSolver<STATE> step(L2);
    step.SetDirect(ELU);
    step.Solve(inner,multiplier);
    {
        std::stringstream sout;
        sout << "multipl" << gfluxorder << " = ";
        filename = sout.str();
    }
    
    multiplier.Print(filename.c_str(),output,EMathematicaInput);
    output.close();
    int nwrong = 0;
    nwrong = VerifyProjection(cel, multiplier);
    if(nwrong)
    {
        std::cout << "Number of points with wrong pressure projection " << nwrong << std::endl;
    }
    //    return nwrong;
    
}

static void VerifyPressureShapeProperties(int order, TPZMaterialData &data, TPZVec<REAL> &pt)
{
    // quadratic bubble function
    REAL bolha = data.phi(0)*data.phi(4)*16.;
    int ibolha = 2*(order+1)*(order+1)+4*(order-1)-2;
    REAL phival = data.phi(ibolha);
    REAL diff = bolha-phival;
    if (fabs(diff) > 1.e-9) {
        std::cout << "Bolha nao identificada\n";
    }
    // lateral function
    int offset = 0;
    bolha = (data.phi(0+offset)*data.phi(2+offset)+data.phi(0+offset)*data.phi(4+offset))*4.;
    ibolha = offset+8;
    phival = data.phi(ibolha);
    diff = bolha-phival;
    if (fabs(diff) > 1.e-9) {
        std::cout << "Bolha nao identificada\n";
    }
    // quadratic singular bubble function
    offset = 1;
    bolha = data.phi(0)*data.phi(4+offset)*16.;
    ibolha = 2*(order+1)*(order+1)+4*(order-1)-2+offset;
    phival = data.phi(ibolha);
    diff = bolha-phival;
    if (fabs(diff) > 1.e-9) {
        std::cout << "Bolha nao identificada\n";
    }
    // lateral singular function
    bolha = (data.phi(0)*data.phi(2+offset)+data.phi(0)*data.phi(4+offset))*4.;
    ibolha = offset+8;
    phival = data.phi(ibolha);
    diff = bolha-phival;
    if (fabs(diff) > 1.e-9) {
        std::cout << "Bolha nao identificada\n";
    }
}

/// Generate the L2 matrix of the pressure space and the inner product of the divergence and the pressure shape functions
static void GenerateProjectionMatrix(TPZCompEl *cel, TPZAutoPointer<TPZMatrix<STATE> > L2, TPZFMatrix<STATE> &inner)
{
    TPZMaterialData dataA,dataB;
    TPZMultiphysicsElement *celMF = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if (!celMF) {
        DebugStop();
    }
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(celMF->Element(0));
    TPZInterpolationSpace *intelP = dynamic_cast<TPZInterpolationSpace *>(celMF->Element(1));
    if (!intel || ! intelP) {
        DebugStop();
    }
    intel->InitMaterialData(dataA);
    intelP->InitMaterialData(dataB);
    int dim = intel->Reference()->Dimension();
    const TPZIntPoints &intrule = intel->GetIntegrationRule();
    int np = intrule.NPoints();
    int npressure = dataB.phi.Rows();
    int nflux = dataA.fVecShapeIndex.NElements();
    L2->Redim(npressure,npressure);
    inner.Redim(npressure,nflux);
    {
        REAL weight;
        TPZManVector<REAL,3> pos(dim);
        int ip = intrule.NPoints()/2;
        intrule.Point(ip, pos, weight);
        //        intel->ComputeShape(pos, dataA.x, dataA.jacobian, dataA.axes, dataA.detjac, dataA.jacinv, dataA.phi, dataA.dphix);
        intel->ComputeRequiredData(dataA, pos);
        intelP->ComputeRequiredData(dataB, pos);
        {
            int order = intel->Mesh()->GetDefaultOrder();
            std::stringstream filename;
            filename << "../phis"<< order << ".nb";
            std::ofstream out(filename.str().c_str());
            std::string phiHDiv, phiPress, VecShape, PosSub, Directions;
            phiHDiv = "phiHDiv" + std::to_string(order) + " = ";
            phiPress = "phiPress" + std::to_string(order) + " = ";
            VecShape = "VecShape" + std::to_string(order) + " = ";
            Directions = "Directions" + std::to_string(order) + " = ";
            PosSub = "PosSub" + std::to_string(order) + " = ";
            
            dataA.phi.Print(phiHDiv.c_str(),out,EMathematicaInput);
            dataB.phi.Print(phiPress.c_str(),out,EMathematicaInput);
            TPZFMatrix<REAL> normalvec;
            dataA.fNormalVec.Transpose(&normalvec);
            normalvec.Print(Directions.c_str(),out,EMathematicaInput);
            out << VecShape << " {\n";
            for (int i=0; i<dataA.fVecShapeIndex.size(); i++) {
                out << "{" << dataA.fVecShapeIndex[i].first+1 << "," << dataA.fVecShapeIndex[i].second+1 << "}";
                if (i != dataA.fVecShapeIndex.size()-1) {
                    out << ",";
                }
                out << std::endl;
            }
            out << "};\n";
            out.precision(16);
            out << PosSub << " { x-> "<< pos[0] << ", y-> " << pos[1] << ", z-> " << pos[2] << "};" << std::endl;
        }
        
        int order6 = intelP->Connect(6).Order();
        if(order6 > 1)
        {
            VerifyPressureShapeProperties(intelP->Connect(6).Order(), dataB,pos);
        }
        if (order6 == 3)
        {
            int vecindex = dataA.fVecShapeIndex[85].first;
            int phiindex = dataA.fVecShapeIndex[85].second;
            std::cout << "vector column 85 phiindex " << phiindex << ' ';
            for (int i=0; i<3; i++) {
                std::cout << dataA.fNormalVec(i,vecindex) << " ";
            }
            std::cout << std::endl;
            vecindex = dataA.fVecShapeIndex[83].first;
            
            phiindex = dataA.fVecShapeIndex[83].second;
            std::cout << "vector column 83 phiindex " << phiindex << ' ';
            for (int i=0; i<3; i++) {
                std::cout << dataA.fNormalVec(i,vecindex) << " ";
            }
            std::cout << std::endl;
        }
    }
    int ip;
    for (ip=0; ip<np; ip++) {
        REAL weight;
        TPZManVector<REAL,3> pos(dim);
        intrule.Point(ip, pos, weight);
        //        intel->ComputeShape(pos, dataA.x, dataA.jacobian, dataA.axes, dataA.detjac, dataA.jacinv, dataA.phi, dataA.dphix);
        intel->ComputeRequiredData(dataA, pos);
        intelP->ComputeShape(pos, dataB);
        int ish,jsh;
        for (ish=0; ish<npressure; ish++) {
            for (jsh=0; jsh<npressure; jsh++) {
                L2->s(ish,jsh) += dataB.phi(ish,0)*dataB.phi(jsh,0)*weight*fabs(dataB.detjac);
            }
            for (jsh=0; jsh<nflux; jsh++) {
                // compute the divergence of the shapefunction
                TPZManVector<REAL,3> vecinner(intel->Dimension(),0.);
                int vecindex = dataA.fVecShapeIndex[jsh].first;
                int phiindex = dataA.fVecShapeIndex[jsh].second;
                int j;
                int d;
                for (d=0; d<dim; d++) {
                    vecinner[d]=0;
                    for (j=0; j<3; j++) {
                        vecinner[d] += dataA.fNormalVec(j,vecindex)*dataA.axes(d,j);
                    }
                }
                REAL divphi = 0.;
                for (d=0; d<dim; d++) {
                    divphi += dataA.dphix(d,phiindex)*vecinner[d];
                }
                if (jsh == 83) {
                    divphi = dataB.phi(27)*pos[0];
                    REAL da27 = -dataB.phi(27);
                    REAL da37 = dataB.phi(37)/2.;
                    REAL verify = divphi-(-dataB.phi(27)+dataB.phi(37)/2.);
                    std::cout << " phi 27 " << dataB.phi(27) << " phi 27 x " << dataB.phi(27)*pos[0] << " verify " << verify << std::endl;
                    //                    divphi *= pos[0];
                    //                    divphi += dataA.phi(phiindex,0);
                }
                inner(ish,jsh) += dataB.phi(ish,0)*divphi*weight*fabs(dataA.detjac);
            }
        }
    }
}

/// Given the multiplier coefficients of the pressure space, verify the correspondence of the divergence of the vector function and the L2 projection
static int VerifyProjection(TPZCompEl *cel, TPZFMatrix<STATE> &multiplier)
{
    TPZMaterialData dataA,dataB;
    TPZMultiphysicsElement *celMF = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if (!celMF) {
        DebugStop();
    }
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(celMF->Element(0));
    TPZInterpolationSpace *intelP = dynamic_cast<TPZInterpolationSpace *>(celMF->Element(1));
    
    if (!intelP || !intel) {
        DebugStop();
    }
    intel->InitMaterialData(dataA);
    intelP->InitMaterialData(dataB);
    int dim = intel->Reference()->Dimension();
    const TPZIntPoints &intrule = intel->GetIntegrationRule();
    int np = intrule.NPoints();
    int npressure = dataB.phi.Rows();
    int nflux = dataA.fVecShapeIndex.NElements();
    TPZFNMatrix<30> pointpos(2,np);
    TPZFNMatrix<30> divergence(np,nflux);
    int ip;
    //std::cout << dataA.fVecShapeIndex << std::endl;
    int nwrong = 0;
    for (ip=0; ip<np; ip++) {
        REAL weight;
        TPZManVector<REAL,3> pos(dim);
        intrule.Point(ip, pos, weight);
        pointpos(0,ip) = pos[0];
        pointpos(1,ip) = pos[1];
        //        intel->ComputeShape(pos, dataA.x, dataA.jacobian, dataA.axes, dataA.detjac, dataA.jacinv, dataA.phi, dataA.dphix);
        intel->ComputeRequiredData(dataA, pos);
        intelP->ComputeShape(pos, dataB);
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Phi's " << dataA.phi<< " dphix's "<< dataA.dphix<<std::endl;
            
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        int ish,jsh;
        for (jsh=0; jsh<nflux; jsh++) {
            // compute the divergence of the shapefunction
            TPZManVector<REAL,3> vecinner(intel->Dimension(),0.);
            int vecindex = dataA.fVecShapeIndex[jsh].first;
            int phiindex = dataA.fVecShapeIndex[jsh].second;
            int j;
            int d;
            for (d=0; d<dim; d++) {
                vecinner[d]=0;
                for (j=0; j<3; j++) {
                    vecinner[d] += dataA.fNormalVec(j,vecindex)*dataA.axes(d,j);
                }
            }
            REAL divphi = 0.;
            for (d=0; d<dim; d++) {
                divphi += dataA.dphix(d,phiindex)*vecinner[d];
            }
            
            if (jsh == 83) {
                //                divphi *= pos[0];
                //                divphi += dataA.phi(phiindex,0);
                divphi = dataB.phi(27)*pos[0];
            }
            //#ifdef LOG4CXX
            //						{
            //								std::stringstream sout;
            //								sout << "Div " << divphi<< std::endl;
            //
            //								LOGPZ_DEBUG(logger,sout.str())
            //						}
            //#endif
            divergence(ip,jsh) = divphi;
            REAL phival = 0;
            
            for (ish=0; ish<npressure; ish++) {
                
                phival += multiplier(ish,jsh)*dataB.phi(ish);
            }
            // the divergence of the vector function should be equal to the value of projected pressure space
            REAL diff = phival-divphi;
#ifdef LOG4CXX
            if (logger->isDebugEnabled())
            {
                std::stringstream sout;
                sout << "phi: " << phival<<" dphi: "<< divphi <<"\n";
                sout << "flux number " << jsh << " diff: "<<diff<< "\n";
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
            if(fabs(diff) > 1.e-6)
            {
                nwrong++;
                std::cout << "flux number " << jsh << " did not project: diff: "<<diff<<"\n";
                //StopError();
                DebugStop();
            }
        }
    }
    
    /*
     int ifl;
     std::ofstream fluxes("fluxes.nb");
     for (ifl=0; ifl<nflux; ifl++) {
     fluxes << "flux" << ifl << " = {\n";
     for (ip=0; ip<np; ip++) {
     fluxes << "{ " << pointpos(0,ip) << " , " << pointpos(1,ip) << " , " << divergence(ip,ifl) << "} ";
     if(ip<np-1) fluxes << "," << std::endl;
     }
     fluxes << " };\n";
     }
     */
    return nwrong;
}



/// verify if the pressure space is compatible with the flux space
void VerifyDRhamCompatibility(TSimulationControl * control)
{
    // generate a mesh
    //   gRefDBase.InitializeAllUniformRefPatterns();
    HDivPiola = 1;
    TPZGeoMesh *gmesh = CreateGeoMesh1Pir();
    int nref = 0;
    UniformRefine(gmesh, nref);
    {
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, "../PyramidGMesh.vtk", true);
    }
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZManVector<TPZCompMesh*,2> meshvec(2);
    int p = 1;
    meshvec[1] = CreateCmeshPressure(gmesh, control, p);
    LoadSolution(meshvec[1]);
    meshvec[0] = CreateCmeshFlux(gmesh, control, p);
    TPZCompMeshTools::AddHDivPyramidRestraints(meshvec[0]);
    //    ProjectFlux(meshvec[0]);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        meshvec[0]->Print(sout);
        meshvec[1]->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZCompMesh *cmeshMult = CreateCmeshMulti(meshvec, control);
    std::ofstream arg1("cmesh.txt");
    cmeshMult->Print(arg1);
    // for each computational element (not boundary) verify if the Div(vecspace) is included in the pressure space
    int nel = cmeshMult->NElements();
    int meshdim = cmeshMult->Dimension();
    int iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmeshMult->ElementVec()[iel];
        TPZMultiphysicsElement *intel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!intel)
        {
            DebugStop();
        }
        if(intel->Reference()->Dimension() != meshdim) continue;
        CheckDRham(intel);
    }
}

TPZAutoPointer<TPZRefPattern> PyramidRef()
{
    TPZGeoMesh gmesh;
    gmesh.NodeVec().Resize(5);
    REAL nodeco[][3] =
    {
        {-1,-1,0},
        {1,-1,0},
        {1,1,0},
        {-1,1,0},
        {0,0,1}
    };
    int64_t nodeindexes[][5] = {
        {0,1,2,3,4},
        {0,1,3,4},
        {1,2,3,4}
    };
    for (int i=0; i<5; i++) {
        TPZManVector<REAL,3> coord(3);
        for (int c=0; c<3; c++) {
            coord[c] = nodeco[i][c];
        }
        gmesh.NodeVec()[i].Initialize(coord, gmesh);
    }
    TPZManVector<int64_t> corners(5);
    for (int64_t i=0; i<5; i++) {
        corners[i] = nodeindexes[0][i];
    }
    int64_t elindex, subels[2];
    gmesh.CreateGeoElement(EPiramide, corners, 1, elindex);
    int64_t fatherindex = elindex;
    for (int64_t i=0; i<4; i++) {
        corners[i] = nodeindexes[1][i];
    }
    gmesh.CreateGeoElement(ETetraedro, corners, 1, subels[0]);
    gmesh.Element(subels[0])->SetFather(fatherindex);
    for (int64_t i=0; i<4; i++) {
        corners[i] = nodeindexes[2][i];
    }
    gmesh.CreateGeoElement(ETetraedro, corners, 1, subels[1]);
    gmesh.Element(subels[1])->SetFather(fatherindex);
    gmesh.BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(gmesh);
    refpat->SetName("PyramidToTetraedra");
    
    gmesh.Element(subels[0])->SetFather(-1);
    gmesh.Element(subels[1])->SetFather(-1);

    if(!gRefDBase.FindRefPattern(refpat))
    {
        gRefDBase.InsertRefPattern(refpat);
        refpat->InsertPermuted();
    }
    return refpat;
}

TPZAutoPointer<TPZRefPattern> PyramidTo4Tetrahedra()
{
    TPZGeoMesh gmesh;
    gmesh.NodeVec().Resize(6);
    REAL nodeco[][3] =
    {
        {-1,-1,0},
        {1,-1,0},
        {1,1,0},
        {-1,1,0},
        {0,0,0},
        {0,0,1}
    };
    int64_t nodeindexes[][5] = {
        {0,1,2,3,5},
        {0,1,4,5},
        {1,2,4,5},
        {2,3,4,5},
        {3,0,4,5},
    };
    for (int i=0; i<6; i++) {
        TPZManVector<REAL,3> coord(3);
        for (int c=0; c<3; c++) {
            coord[c] = nodeco[i][c];
        }
        gmesh.NodeVec()[i].Initialize(coord, gmesh);
    }
    TPZManVector<int64_t,5> corners(5);
    for (int64_t i=0; i<5; i++) {
        corners[i] = nodeindexes[0][i];
    }
    int64_t elindex, subels[4];
    gmesh.CreateGeoElement(EPiramide, corners, 1, elindex);
    int64_t fatherindex = elindex;
    for (int sub=0; sub<4; sub++)
    {
        for (int64_t i=0; i<4; i++) {
            corners[i] = nodeindexes[sub+1][i];
        }
        gmesh.CreateGeoElement(ETetraedro, corners, 1, subels[sub]);
        gmesh.Element(subels[sub])->SetFather(fatherindex);
    }
    gmesh.BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(gmesh);
    refpat->SetName("PyramidToTetraedra");
    
    for(int sub=0; sub<4; sub++)
    {
        gmesh.Element(subels[sub])->SetFather(-1);
    }
    
    if(!gRefDBase.FindRefPattern(refpat))
    {
        gRefDBase.InsertRefPattern(refpat);
        refpat->InsertPermuted();
    }
    return refpat;
}

void DividePyramids(TPZGeoMesh &gmesh)
{
    TPZAutoPointer<TPZRefPattern> refpat = gRefDBase.FindRefPattern("PyramidToTetraedra");
    if (!refpat) {
        DebugStop();
    }
    int64_t nel = gmesh.NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.Element(el);
        if (!gel || gel->Type() != EPiramide || gel->HasSubElement()) {
            continue;
        }
        TPZManVector<TPZGeoEl *,4> subels(4,0);
        gel->SetRefPattern(refpat);
        gel->Divide(subels);
    }
    
}

void DivideBoundaryElements(TPZGeoMesh &gmesh, int exceptmatid)
{
    int64_t nel = gmesh.NElements();
    int meshdim = gmesh.Dimension();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh.Element(el);
        int matid = gel->MaterialId();
        if (gel->Dimension() != meshdim-1 || matid == exceptmatid) {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        bool shoulddivide = false;
        while (neighbour != gelside) {
            int nsidesub = neighbour.Element()->NSideSubElements(neighbour.Side());
            if (nsidesub > 1) {
                shoulddivide = true;
            }
            neighbour = neighbour.Neighbour();
        }
        if (shoulddivide) {
            // find compatible refinement pattern
            TPZAutoPointer<TPZRefPattern> compatible;
            std::list<TPZAutoPointer<TPZRefPattern> > patlist;
            TPZRefPatternTools::GetCompatibleRefPatterns(gel, patlist);
            if (patlist.size() != 1) {
                DebugStop();
            }
            compatible = *patlist.begin();
            gel->SetRefPattern(compatible);
            // Create a geometric element along the quadrilateral side
            TPZManVector<TPZGeoEl *,4> subels(2,0);
            gel->Divide(subels);
        }
    }
}

void IncreasePyramidSonOrder(TPZVec<TPZCompMesh *> &meshvec, TSimulationControl * control, int pFlux)
{
    
    int n_accuracy = 0;
//    if (control->m_Hdiv_plusplus_Q) {
//        n_accuracy++;
//    }
    
    meshvec[0]->Reference()->ResetReference();
    meshvec[0]->LoadReferences();
    int64_t nel = meshvec[0]->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = meshvec[0]->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Type() != ETetraedro || !gel->Father()) {
            continue;
        }
        // we have a tetrahedra with father
        if (gel->Father()->Type() != EPiramide) {
            continue;
        }
        TPZConnect &c = cel->Connect(0);
        bool hasdependency = c.HasDependency();

        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        intel->PRefine(2*pFlux+n_accuracy);
        TPZCompElSide large;
        
        if(hasdependency)
        {
            TPZGeoElSide gelside(gel,10);
            large = gelside.LowerLevelCompElementList2(1);
            if (!large) {
                DebugStop();
            }
            intel->RemoveSideRestraintWithRespectTo(10,large);
            
            if(c.HasDependency())
            {
                DebugStop();
            }
        }
        int nshape = intel->NConnectShapeF(0, 2*pFlux);
        int64_t cindex = cel->ConnectIndex(0);
        c.SetOrder(2*pFlux, cindex);
        c.SetNShape(nshape);
        int64_t seqnum = c.SequenceNumber();
        meshvec[0]->Block().Set(seqnum, nshape);
        if(hasdependency)
        {
            TPZInterpolatedElement *largeintel = dynamic_cast<TPZInterpolatedElement *>(large.Element());
            intel->RestrainSide(10, largeintel, large.Side());
        }

    }
    meshvec[0]->ExpandSolution();
    meshvec[1]->Reference()->ResetReference();
    nel = meshvec[1]->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = meshvec[1]->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Type() != ETetraedro || !gel->Father()) {
            continue;
        }
        // we have a tetrahedra with father
        if (gel->Father()->Type() != EPiramide) {
            continue;
        }
        
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(meshvec[1]->Element(el));
        intel->PRefine(2*pFlux+n_accuracy);
    }
    meshvec[1]->ExpandSolution();
}

void PrintArrayInMathematica(TPZVec<REAL> &array, std::ofstream &out, std::string arrayName)
{
    const int size = array.size();
    out << arrayName << " = {" << array[0];
    for (int i = 1; i < size ; i++) {
        out <<  "," << array[i];
    }
    out << "};" << endl;
}

void GetAllStringsForVariablesInMathematica(int &pFlux, EApproxSpace &runtype, bool &HDivMaisMais, std::string &neqstr,
                                            std::string &hSizestr, std::string &h1Str, std::string &l2str,
                                            std::string &semih1str, std::string &plotnamestr, std::string &convrateh1str,
                                            std::string &convratel2str, std::string &convratesemih1str, std::string &ptsh1str,
                                            std::string &ptsl2str, std::string &ptssemistr){
    
    const int size = 12;
    std::stringstream strVec[size];
    strVec[0] << "neq";
    strVec[1] << "h";
    strVec[2] << "h1";
    strVec[3] << "l2";
    strVec[4] << "semih1";
    strVec[5] << "plot";
    strVec[6] << "convratesh1";
    strVec[7] << "convratesl2";
    strVec[8] << "convratesh1semi";
    strVec[9] << "ptsh1";
    strVec[10] << "ptsl2";
    strVec[11] << "ptssemih1";
    
    if (runtype == ETetrahedra) {
        for(int i = 0 ; i < size ; i++){
            strVec[i] << "Tet";
        }
    }
    else if (runtype == EPyramid) {
        for(int i = 0 ; i < size ; i++){
            strVec[i] << "Pyr";
        }
    }
    else if(runtype == EDividedPyramid)
    {
        for(int i = 0 ; i < size ; i++){
            strVec[i] << "DivPyr";
        }
    }
    else if(runtype == EDividedPyramidIncreasedOrder)
    {
        for(int i = 0 ; i < size ; i++){
            strVec[i] << "DivPyrIncOrd";
        }
    }
    else if(runtype == EDividedPyramid4)
    {
        for(int i = 0 ; i < size ; i++){
            strVec[i] << "DivPyr4";
        }
    }
    else if(runtype == EDividedPyramidIncreasedOrder4)
    {
        for(int i = 0 ; i < size ; i++){
            strVec[i] << "DivPyr4IncOrd";
        }
    }
    else
    {
        DebugStop();
    }
    
    for(int i = 0 ; i < size ; i++){
        strVec[i] << "P" << pFlux;
    }
    
    if(HDivMaisMais){
        for(int i = 0 ; i < size ; i++){
            strVec[i] << "plusplus";
        }
    }
    
    neqstr = strVec[0].str();
    hSizestr = strVec[1].str();
    h1Str = strVec[2].str();
    l2str = strVec[3].str();
    semih1str = strVec[4].str();
    plotnamestr = strVec[5].str();
    convrateh1str = strVec[6].str();
    convratel2str = strVec[7].str();
    convratesemih1str = strVec[8].str();
    ptsh1str = strVec[9].str();
    ptsl2str = strVec[10].str();
    ptssemistr = strVec[11].str();
    
    
}


void GenerateMathematicaWithConvergenceRates(TPZVec<REAL> &neqVec, TPZVec<REAL> &hSizeVec,
                                             TPZVec<REAL> &h1ErrVec, TPZVec<REAL> &l2ErrVec,
                                             TPZVec<REAL> &semih1ErrVec, TSimulationControl * control)
{
    // ---------------- Defining filename ---------------
    int p_order = control->m_p_levels;
    bool hdiv_plus_plus_Q = control->m_Hdiv_plusplus_Q;
    EApproxSpace run_type = control->m_run_type;
    std::string mathematicaFilename = "NoName.nb";
    std::stringstream Mathsout;
    if(run_type == ETetrahedra){Mathsout << "convergenceRatesTetMesh";}
    if(run_type == EPyramid){Mathsout << "convergenceRatesPyrMesh";}
    if(run_type == EDividedPyramid){Mathsout << "convergenceRatesDividedPyrMesh";}
    if(run_type == EDividedPyramidIncreasedOrder){Mathsout << "convergenceRatesDivPyrIncOrdMesh";}
    if(run_type == EDividedPyramid4){Mathsout << "convergenceRatesDividedPyr4Mesh";}
    if(run_type == EDividedPyramidIncreasedOrder4){Mathsout << "convergenceRatesDivPyr4IncOrdMesh";}
    Mathsout << control->m_p_levels;
    if (control->m_Hdiv_plusplus_Q) {
        Mathsout << "plusplus";
    }
    Mathsout << ".nb";
    mathematicaFilename = Mathsout.str();
    std::ofstream out(mathematicaFilename.c_str());
    
    
    // ---------------- Defining plot parameters ---------------
    out << "myAxesSize = 15;" << endl;
    out << "myLabelSize = 12;" << endl;
    out << "myImageSize = 500;" << endl;
    
    // ---------------- Printing vectors ---------------
    std::string neqstr, hSizestr, h1Str, l2str, semih1str, plotnamestr, convrateh1str, convratel2str, convratesemih1str, ptsh1str, ptsl2str, ptssemistr;
    GetAllStringsForVariablesInMathematica(p_order, run_type, hdiv_plus_plus_Q, neqstr, hSizestr, h1Str, l2str, semih1str, plotnamestr, convrateh1str, convratel2str, convratesemih1str, ptsh1str, ptsl2str, ptssemistr);
    PrintArrayInMathematica(neqVec,out,neqstr);
    PrintArrayInMathematica(hSizeVec,out,hSizestr);
    PrintArrayInMathematica(h1ErrVec,out,h1Str);
    PrintArrayInMathematica(l2ErrVec,out,l2str);
    PrintArrayInMathematica(semih1ErrVec,out,semih1str);
    
    // ---------------- Tables with the points to plot ---------------
    out << ptsh1str << " = Table[{" << hSizestr << "[[i]], " << h1Str << "[[i]]}, {i, Length[" << neqstr << "]}];" << endl;
    out << ptsl2str << " = Table[{" << hSizestr << "[[i]], " << l2str << "[[i]]}, {i, Length[" << neqstr << "]}];" << endl;
    out << ptssemistr <<  " = Table[{" << hSizestr << "[[i]], " << semih1str << "[[i]]}, {i, Length[" << neqstr << "]}];" << endl;
    
    // ---------------- ListLogLogPlot ---------------
    out << plotnamestr << " = ListLogLogPlot[{" << ptsh1str << ", " << ptsl2str << ", " << ptssemistr << "}, PlotRange -> All," <<
    "PlotStyle -> Thickness[0.003], ImageSize -> myImageSize 1.5,\n" <<
    "AxesLabel -> {Style[\"h\", myAxesSize, Black],Style[Err, myAxesSize, Black]}," <<
    "LabelStyle -> Directive[myLabelSize, Bold, Black], Joined -> True,\n" <<
    "PlotMarkers -> Automatic," <<
    "PlotLegends -> Placed[{\"Potential\", \"Flux\", \"Divergence\"}, {0.25, 0.85}]]" << endl;
    
    
    // ---------------- Calculating sloples ---------------
    out << convrateh1str << " = Range[Length[" << neqstr << "] - 1];" << endl;
    out << convratel2str << " = Range[Length[" << neqstr << "] - 1];" << endl;
    out << convratesemih1str << " = Range[Length[" << neqstr << "] - 1];" << endl;
    
    out << "For[i = 1, i <= Length[" << neqstr << "] - 1, i++,\n" <<
    "\t" << convrateh1str << "[[i]] = Log[" << h1Str << "[[i]]/" << h1Str << "[[i + 1]]]/Log[" << hSizestr << "[[i]]/" << hSizestr << "[[i + 1]]];\n" <<
    "\t" << convratel2str << "[[i]] = Log[" << l2str << "[[i]]/" << l2str << "[[i + 1]]]/Log[" << hSizestr << "[[i]]/" << hSizestr << "[[i + 1]]];\n" <<
    "\t" << convratesemih1str << "[[i]] = Log[" << semih1str << "[[i]]/" << semih1str << "[[i + 1]]]/Log[" << hSizestr << "[[i]]/" << hSizestr << "[[i + 1]]];\n" <<
    "]" << endl;
    
    // ---------------- Printing slopes ---------------
    out << convrateh1str << endl;
    out << convratel2str << endl;
    out << convratesemih1str << endl;
    
    std::cout << "WRITTEN OUTPUT FILE ********* " << mathematicaFilename << " *****************" << std::endl;

    
}

/// uncondense the elements unwrap the elements
void UnwrapMesh(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    bool change = true;
    while(change)
    {
        change = false;
        for (int64_t el=0; el<nel; el++) {
            
            TPZCompEl *cel = cmesh->Element(el);
            TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
            if (condense) {
                condense->Unwrap();
                change = true;
            }
            cel = cmesh->Element(el);
            TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
            if (elgr) {
                elgr->Unwrap();
                change = true;
            }
        }
    }
}
