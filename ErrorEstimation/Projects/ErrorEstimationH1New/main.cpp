/**
 * @file This file implements an error estimator in space H1.
 */

#include "pzlog.h"
#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"
#include "TPZBndCondT.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "pzgeoelbc.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGmshReader.h"
#include "TPZAnalyticSolution.h"
#include "TPZFrontSym.h"

#include <fstream>
#include <ctime>
#include <cstdio>
#include <cmath>

#include "ProblemConfig.h"
#include "TPZPostProcessError.h"
#include "InputTreatment.h"
#include "Tools.h"
#include "Solver.h"
//#include "pzgengrid.h"
#include "TPZGenGrid2D.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"
#include "Tools.h"

#include <iostream>
// Global variables
const int problemDimension = 2;
const bool readGMeshFromFile = false;

const int matID = 1;

// Functions declarations
TPZGeoMesh *CreateGeoMesh();
TPZGeoMesh *CreateGeoMesh2();
TPZCompMesh *CMeshPressure(struct SimulationCase &sim_case);
void UniformRefinement(int nDiv, TPZGeoMesh *gmesh);
bool SolvePoissonProblem(struct SimulationCase &sim_case);
bool PostProcessProblem(TPZAnalysis &an, TPZGeoMesh * gmesh, TPZCompMesh * pressuremesh);

void SolveH1Problem(TPZCompMesh *cmeshH1, ProblemConfig &config);
TPZCompMesh *CompMeshH1(ProblemConfig &problem);
TPZGeoMesh *GeometricMesh(int nel, TPZVec<int> &bcids);
TPZGeoMesh* CreateLShapeMesh(int nel, TPZVec<int>& bcids);
TPZGeoMesh* CreateSquareShapeMesh(int nel, TPZVec<int>& bcids);
TPZGeoMesh* CreateSquareShapeMesh2(int nel, TPZVec<int>& bcids);

bool PostProcessing(TPZCompMesh * pressuremesh,TPZFMatrix<STATE> true_elerror, TPZFMatrix<STATE> estimate_elerror);

TPZGeoMesh *CreateLCircleGeoMesh();
TPZGeoMesh *CreateGeoCircleMesh();
double Maximum(TPZVec<double> &vect);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.refine"));
#endif

using namespace std;
// Laplace equation on square 1D 2D 3D - Volker John article 2000
double errorseminorm = 0.;

int main(int argc, char *argv[]) {
    
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    for(int ndiv = 0; ndiv < 1; ndiv++){
        
        std::string meshfilename = "../Quad.msh";
        
        TPZGeoMesh *gmesh = nullptr;
        
        if (readGMeshFromFile) {
            TPZGmshReader gmsh;
            
            gmsh.GetDimNamePhysical()[1]["dirichlet"] = -1;
            gmsh.GetDimNamePhysical()[1]["neuman"] = -2;
            gmsh.GetDimNamePhysical()[2]["domain"] = 1;
            
            //        gmsh.SetFormatVersion("4.0");
#ifdef MACOSX
            gmesh = gmsh.GeometricGmshMesh(meshfilename);
#else
            gmesh = gmsh.GeometricGmshMesh(meshfilename);
#endif
            
        }
        else {
            gmesh = CreateGeoMesh2(); //[0,1]x[0,1] quadrilateral
            //gmesh = CreateGeoMesh(); //[-1,1]x[-1,1] quadrilateral
            
            TPZManVector<int, 4> bcids(8, -1);
            //gmesh = CreateLShapeMesh(1, bcids);//CreateGeoCircleMesh();
            //gmesh = CreateSquareShapeMesh2(1, bcids);//[0,1]x[0,1] triangular
            //gmesh = CreateSquareShapeMesh(1, bcids);//[-1,1]x[-1,1] triangular
        }
        
        ProblemConfig config;
        PreConfig pConfig;
        //SimulationCase Case1;
        
        //Case1.nthreads = 0;
        
        pConfig.refLevel = 3;
        //Case1.numinitialrefine = 1;//ndiv;
        
        pConfig.k = 2;
        //Case1.porder = 2;
        
        config.gmesh = gmesh;
        //Case1.gmesh = gmesh;
        
        config.materialids.insert(1);
        //Case1.materialids.insert(1);
        
        config.bcmaterialids.insert(-1);
        //Case1.bcmaterialids.insert(-1);
        
        //    Case1.bcmaterialids.insert(-2);
        //    Case1.bcmaterialids.insert(2);//para sinmark
        
        TLaplaceExample1 example;
        config.exact = new TLaplaceExample1;
        config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
        //Case1.exact.fExact = example.ESinSin;//ESinMark//ESinSin//ESinSinDirNonHom
        
        pConfig.approx = "H1";                 //// {"H1","Hybrid", "Mixed"}
        pConfig.topology = "Quadrilateral";        //// Triangular, Quadrilateral, LQuad, Tetrahedral, Hexahedral, Prism
        
        pConfig.problem = "ESinSin";
        //Case1.problemname = "ESinSin";//ESinMark,EConst,EBubble2D,ESteepWave,EX
        
        config.dir_name = "ESinSin";
        //Case1.dir_name = "ESinSin";
        std::string command = "mkdir -p " +config.dir_name; //+ porder;
        system(command.c_str());
                
        pConfig.numberAdapativitySteps = 3;        //// Maximum number of adapativity refinement steps.
        pConfig.estimateError = true;              //// Wheater Error Estimation procedure is invoked
        pConfig.debugger = true;                   //// Print geometric and computational mesh for the simulation (Error estimate not involved).
        pConfig.vtkResolution = 1;                 //// Vtk resolution. Set 0 to see a paraview mesh equals the  simulation mesh.
        
        if(0){ //If hp-adaptivity
            
            EvaluateEntry(argc,argv,pConfig);
            InitializeOutstream(pConfig,argv);
            
            {
                std::ofstream outgmesh("gmesh.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outgmesh);
                std::ofstream outgmesh2("gmesh.txt");
                gmesh->Print(outgmesh2);
            }
            
            {
                config.division_threshold = 0.85;
                for (config.adaptivityStep = 0; config.adaptivityStep < pConfig.numberAdapativitySteps+1; config.adaptivityStep++) {     //ndiv = 1 corresponds to a 2x2 mesh.
                    pConfig.h = 1./pConfig.exp;
                    
                    Configure(config,pConfig.refLevel,pConfig,argv);
                    
                    Solve(config,pConfig);
                }
            }
        }
        else{
            
            UniformRefinement(pConfig.refLevel,gmesh);
            std::ofstream file("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
            std::ofstream outgmesh("gmesh.txt");
            gmesh->Print(outgmesh);
            
            //Solve H1 Problem
            config.exact.operator*().fSignConvention = 1;
            //Case1.exact.fSignConvention = 1;
            
            gmesh->ResetReference();
            TPZCompMesh *cmeshH1 = CompMeshH1(config);//CMeshPressure(Case1);
            
            SolveH1Problem(cmeshH1,config);
            
            //            TPZLinearAnalysis an(cmeshH1);
            //            an.SetExact(Case1.exact.ExactSolution());
            
            // Reconstruct Process
            TPZPostProcessError error(cmeshH1);
            error.SetAnalyticSolution(config.exact);
            
            TPZVec<STATE> estimatedelementerror;
            error.ComputeElementErrors(estimatedelementerror);
            //error.MultiPhysicsMesh()->LoadReferences();
            
            TPZFMatrix<STATE> true_elerror(cmeshH1->ElementSolution());
            TPZFMatrix<STATE> estimate_elerror(error.MultiPhysicsMesh()->ElementSolution());
            
            //true_elerror.Print("true error", std::cout);
            //estimate_elerror.Print("estimate error", std::cout);
            std::ofstream outTE("TrueErrorByElem.txt");
            true_elerror.Print(outTE);
            std::ofstream outEE("EstErrorByElem.txt");
            estimate_elerror.Print(outEE);
            
            //            std::ofstream outEE2("EstErrorByElem2.txt");
            //            estimate_elerror.Print(outEE2);
            
            std::ofstream outEE3("EstErrorByElem3.txt");
            estimatedelementerror.Print(outEE3);
            
            //std::cout << "true_elerror: " << true_elerror.GetVal(2, 3) << std::endl;
            
            PostProcessing(cmeshH1,true_elerror, estimate_elerror);
            
        }
    }
    
    //      SolvePoissonProblem(Case1);
    
    //    if(!SolvePoissonProblem(Case1)) {
    //        return 1;
    //    }
    
    return 0;
}

TPZGeoMesh *CreateGeoMesh() {
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    
    if (problemDimension == 2) {
        
        gmesh->SetDimension(2);
        
        // Creates matrix with quadrilateral node coordinates.
        const int quadNodeNumber = 4;
        REAL coordinates[quadNodeNumber][3] = {
            {-1., -1., 0.},
            {1., -1., 0.},
            {1., 1., 0.},
            {-1., 1., 0.}
        };
        
        // Inserts coordinates in the TPZGeoMesh object.
        for(int i = 0; i < quadNodeNumber; i++) {
            int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
            
            TPZVec<REAL> nodeCoord(3);
            nodeCoord[0] = coordinates[i][0];
            nodeCoord[1] = coordinates[i][1];
            nodeCoord[2] = coordinates[i][2];
            
            gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
        }
        
        // Creates quadrilateral element.
        int64_t index;
        TPZManVector<int64_t> nodeIDs(quadNodeNumber);
        
        for(int n = 0; n < quadNodeNumber; n++) {
            nodeIDs[n] = n;
        }
        gmesh->CreateGeoElement(EQuadrilateral, nodeIDs, matID, index);
        
        // Creates line elements where boundary conditions will be inserted.
        nodeIDs.Resize(2);
        for (int i = 0; i < 4; i++) {
            
            nodeIDs[0] = i % 4;
            nodeIDs[1] = (i + 1) % 4;
            
            gmesh->CreateGeoElement(EOned, nodeIDs, -1, index);
        }
        
        gmesh->BuildConnectivity();
        return gmesh;
    }
    else if (problemDimension == 3) {
        // Creates hexahedra element.
        gmesh->SetDimension(3);
        const int hexahedraNodes = 8;
        REAL coordinates[hexahedraNodes][3] = {
            {0., 0., 0.},
            {1., 0., 0.},
            {1., 1., 0.},
            {0., 1., 0.},
            {0., 0., 1.},
            {1., 0., 1.},
            {1., 1., 1.},
            {0., 1., 1.},
        };
        
        for(int n = 0; n < hexahedraNodes; n++) {
            int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
            TPZVec<REAL> coord(3);
            coord[0] = coordinates[n][0];
            coord[1] = coordinates[n][1];
            coord[2] = coordinates[n][2];
            gmesh->NodeVec()[nodeID] = TPZGeoNode(n, coord, *gmesh);
        }
        
        TPZVec<int64_t> nodeID(hexahedraNodes);
        for(int n = 0; n < hexahedraNodes; n++) {
            nodeID[n] = n;
        }
        
        // Inserts Dirichlet BC
        gmesh->BuildConnectivity();
        return gmesh;
    }
    else {
        DebugStop();
    }
}

TPZGeoMesh *CreateGeoMesh2() {
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    
    if (problemDimension == 2) {
        
        gmesh->SetDimension(2);
        
        // Creates matrix with quadrilateral node coordinates.
        const int quadNodeNumber = 4;
        REAL coordinates[quadNodeNumber][3] = {
            {0., 0., 0.},
            {1., 0., 0.},
            {1., 1., 0.},
            {0., 1., 0.}
        };
        
        // Inserts coordinates in the TPZGeoMesh object.
        for(int i = 0; i < quadNodeNumber; i++) {
            int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
            
            TPZVec<REAL> nodeCoord(3);
            nodeCoord[0] = coordinates[i][0];
            nodeCoord[1] = coordinates[i][1];
            nodeCoord[2] = coordinates[i][2];
            
            gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
        }
        
        // Creates quadrilateral element.
        int64_t index;
        TPZManVector<int64_t> nodeIDs(quadNodeNumber);
        
        for(int n = 0; n < quadNodeNumber; n++) {
            nodeIDs[n] = n;
        }
        gmesh->CreateGeoElement(EQuadrilateral, nodeIDs, matID, index);
        
        // Creates line elements where boundary conditions will be inserted.
        nodeIDs.Resize(2);
        for (int i = 0; i < 4; i++) {
            
            nodeIDs[0] = i % 4;
            nodeIDs[1] = (i + 1) % 4;
            
            gmesh->CreateGeoElement(EOned, nodeIDs, -1, index);
        }
        
        gmesh->BuildConnectivity();
        return gmesh;
    }
    else if (problemDimension == 3) {
        // Creates hexahedra element.
        gmesh->SetDimension(3);
        const int hexahedraNodes = 8;
        REAL coordinates[hexahedraNodes][3] = {
            {0., 0., 0.},
            {1., 0., 0.},
            {1., 1., 0.},
            {0., 1., 0.},
            {0., 0., 1.},
            {1., 0., 1.},
            {1., 1., 1.},
            {0., 1., 1.},
        };
        
        for(int n = 0; n < hexahedraNodes; n++) {
            int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
            TPZVec<REAL> coord(3);
            coord[0] = coordinates[n][0];
            coord[1] = coordinates[n][1];
            coord[2] = coordinates[n][2];
            gmesh->NodeVec()[nodeID] = TPZGeoNode(n, coord, *gmesh);
        }
        
        TPZVec<int64_t> nodeID(hexahedraNodes);
        for(int n = 0; n < hexahedraNodes; n++) {
            nodeID[n] = n;
        }
        
        // Inserts Dirichlet BC
        gmesh->BuildConnectivity();
        return gmesh;
    }
    else {
        DebugStop();
    }
}


bool SolvePoissonProblem(struct SimulationCase &sim_case) {
    
    // Creating the directory
    std::string command = "mkdir -p " + sim_case.dir_name;
    system(command.c_str());
    
    // Output files
    std::string file_name = sim_case.dir_name + "/" + "ErrorsHP_Poisson.txt";
    std::ofstream fileerrors(file_name, ios::app);   // To store all errors calculated by TPZAnalysis (PosProcess)
    
    // Initializing the auto adaptive process
    TPZVec<REAL> ervec, ErrorVec(100, 0.0);
    TPZVec<int64_t> NEquations(100, 0L);
    TPZVec<REAL> ervecbyel;
    TPZVec<REAL> gradervecbyel;
    
    /** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
    //scalnames.Push("POrder");
    scalnames.Push("Pressure");
    vecnames.Push("Derivative");
    
    fileerrors.flush();
    //le a malha geometrica
    TPZGeoMesh *gmesh = sim_case.gmesh;
    
    //refina a malha geometrica uniformemente
    {
        // Refines an element
        UniformRefinement(sim_case.numinitialrefine, gmesh);
        std::ofstream out("Gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        
    }
    
    // Creates computational mesh (approximation space and materials)
    TPZCompEl::SetgOrder(sim_case.porder);
    gmesh->SetName("Original GeoMesh");
    
    TPZManVector<TPZCompMesh *> meshvec(0);
    
    TPZCompMesh *pressuremesh = CMeshPressure(sim_case);
    pressuremesh->AdjustBoundaryElements();
    
    //    {
    //        std::ofstream out("CompMesh.vtk");
    //        TPZVTKGeoMesh::PrintCMeshVTK(pressuremesh, out);
    //        std::ofstream out2("CompMesh.txt");
    //        pressuremesh->Print(out2);
    //    }
    
    //define the model problem to be solve
    TLaplaceExample1 example;
    example.fExact = TLaplaceExample1::EX;//ESinSinDirNonHom;//ESinSin;//ESinSinDirNonHom;//ECosCos;
    example.fDimension = gmesh->Dimension();
    example.fSignConvention = -1;
    
    {
        for (auto it:pressuremesh->MaterialVec()) {
            TPZMaterialT<STATE> *mat = dynamic_cast<TPZMaterialT<STATE> *>(it.second);
            TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(mat);
            if (!bc) {
                mat->SetForcingFunction(example.ForceFunc(),3);
            }
            else {
                bc->SetForcingFunctionBC(example.ExactSolution(),3);
            }
        }
    }
    
    TPZLinearAnalysis an(pressuremesh, RenumType::EMetis);
    an.SetExact(example.ExactSolution());
    
    //    {
    //        std::stringstream sout;
    //      sout << sim_case.dir_name << "/" <<  "Poisson" << gmesh->Dimension() << "Porder" << sim_case.porder << ".vtk";
    ////        sout << sim_case.dir_name << "/" << "Poisson" << gmesh->Dimension() << "numref" << sim_case.numinitialrefine
    ////             << "Porder" << sim_case.porder << ".vtk";
    //        an.DefineGraphMesh(gmesh->Dimension(), scalnames, vecnames, sout.str());
    //    }
    
    pressuremesh->SetName("Adapted CompMesh");
    
    // Printing geometric and computational mesh
    //#ifdef ERRORESTIMATION_DEBUG
    //    {
    //        std::ofstream out("../PressureGeoMesh.txt");
    //        pressuremesh->Reference()->Print(out);
    //    }
    //#endif
    
#ifdef PZ_USING_MKL
    // Solves using a symmetric matrix then using Cholesky decomposition (direct method)
    //    TPZSymetricSpStructMatrix strmat(pressuremesh);
    TPZSkylineStructMatrix<STATE> strmat(pressuremesh);
    strmat.SetNumThreads(sim_case.nthreads);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(pressuremesh);
    strmat.SetNumThreads(sim_case.nthreads);
    strmat.SetDecomposeType(ECholesky);
#endif
    
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ECholesky);
    an.SetSolver(*direct);
    delete direct;
    
    an.Run();//resolveu o problema primal
    //reconstruction of flux and error estimation
    PostProcessProblem(an, gmesh, pressuremesh);
    
    return true;
}

bool PostProcessProblem(TPZAnalysis &an, TPZGeoMesh * gmesh, TPZCompMesh * pressuremesh) {
    // Post processing
    an.PostProcess(1, gmesh->Dimension());
    
    
    std::cout<<"Initializing reconstructed process"<<std::endl;
    TPZPostProcessError error(pressuremesh);
    
    TPZVec<STATE> estimatedelementerror, exactelementerror;
    
    error.ComputeElementErrors(estimatedelementerror);
    error.MultiPhysicsMesh()->LoadReferences();
    
    
    {
        int64_t nels = pressuremesh->ElementVec().NElements();
        pressuremesh->ElementSolution().Redim(nels, 6);
    }
    
    
    
    bool store_errors = true;
    an.PostProcessError(exactelementerror, store_errors);
    std::cout << "Exact error " << exactelementerror << std::endl;
    gmesh->ResetReference();
    error.MultiPhysicsMesh()->LoadReferences();
    
    //Compute the effectivity index
    std::cout<<"Computing effectivity index"<<std::endl;
    
    {
        TPZFMatrix<STATE> true_elerror(pressuremesh->ElementSolution());
        TPZFMatrix<STATE> estimate_elerror(error.MultiPhysicsMesh()->ElementSolution());
        int64_t nel = true_elerror.Rows();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = pressuremesh->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            TPZCompEl *mphys = gel->Reference();
            int64_t elindex2 = mphys->Index();
            true_elerror(el,0) = estimate_elerror(elindex2,2);
            true_elerror(el,1) = true_elerror(el,2);
            if (true_elerror(el,1) > 1.e-10) {
                true_elerror(el,2) = true_elerror(el,0)/true_elerror(el,1);
            }
        }
        pressuremesh->ElementSolution() = true_elerror;
    }
    
    {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("State");
        scalnames.Push("EstimatedError");
        scalnames.Push("TrueError");
        scalnames.Push("EffectivityIndex");
        an.DefineGraphMesh(pressuremesh->Dimension(), scalnames, vecnames, "ErrorEstimationH1.vtk");
        an.PostProcess(1);
    }
    
    if(gmesh) delete gmesh;
    
    return true;
}

void UniformRefinement(int nDiv, TPZGeoMesh *gmesh) {
    
    TPZManVector<TPZGeoEl*> children;
    for(int division = 0; division < nDiv; division++) {
        
        int64_t nels = gmesh->NElements();
        
        for(int64_t elem = 0; elem < nels; elem++) {
            
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            
            if(!gel || gel->HasSubElement()) continue;
            if(gel->Dimension() == 0) continue;
            
            if(0){ //if uniform refinament purely
                gel->Divide(children);
            }
            else{ // selected refinament
                if (division < nDiv-1){
                    gel->Divide(children);
                }
                else {
                    if (elem % 2 == 1 & gel->Dimension() != 1){
                        gel->Divide(children);
                    }
                }
            }
        }
    }

    {
        int64_t nel = gmesh->NElements();
        int dim = gmesh->Dimension();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->Dimension() == dim || gel->HasSubElement()) continue;
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while(neighbour != gelside) {
                if(neighbour.Element()->HasSubElement()) {
                    TPZStack<TPZGeoEl *> subs;
                    gel->Divide(subs);
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }
    
//    int nels = gmesh->NElements();
//    for(int64_t elem = 0; elem < nels; elem++) {
//
//        TPZGeoEl * gel = gmesh->ElementVec()[elem];
//
//        if(!gel || gel->HasSubElement()) continue;
//        if(gel->Dimension() != 1) continue;
//        TPZGeoElSide geoelside(gel);
//        TPZGeoElSide neig = geoelside.Neighbour();
//        if(neig.Element()->HasSubElement()){
//            gel->Divide(children);
//        }
//    }
}

TPZCompMesh *CMeshPressure(struct SimulationCase &sim_case) {
    
    int dim = sim_case.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;
    
    // Creates Poisson material
    TPZDarcyFlow *material = new TPZDarcyFlow(matID, dim);
    material->SetExactSol(sim_case.exact.ExactSolution(), 3);
    material->SetForcingFunction(sim_case.exact.ForceFunc(), 3);
    
    TPZCompMesh * cmesh = new TPZCompMesh(sim_case.gmesh);
    cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(material);
    
    
    for (auto matid : sim_case.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.);
        TPZManVector<REAL> val2(1, 0.);
        int bctype = 0;
        if (matid == -2) {
            bctype = 1;
        }
        TPZBndCondT<STATE> *bc = material->CreateBC(material, matid, bctype, val1, val2);
        bc->SetForcingFunctionBC(sim_case.exact.ExactSolution(),3);
        cmesh->InsertMaterialObject(bc);
    }
    
    
    cmesh->SetDefaultOrder(sim_case.porder);
    cmesh->SetAllCreateFunctionsContinuous();
    
    // Adjusts computational data structure
    cmesh->AutoBuild();
    return cmesh;
}


void SolveH1Problem(TPZCompMesh *cmeshH1, ProblemConfig &config){
    
    TPZLinearAnalysis an(cmeshH1);
    an.SetExact(config.exact.operator*().ExactSolution());
    
    
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
//    strmat.SetDecomposeType(ECholesky);
//    strmat.SetNumThreads(0);
    TPZSkylineStructMatrix<STATE> strmat(cmeshH1);
    strmat.SetNumThreads(0);
#endif
    
    std::set<int> matids;
    matids.insert(1);
    
    for(auto mat:config.bcmaterialids){
        matids.insert(mat);
    }
    
    strmat.SetMaterialIds(matids);
    an.SetStructuralMatrix(strmat);
    
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();//resolve o problema ate aqui
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Solution");
    vecnames.Push("Derivative");
    scalnames.Push("ExactSolution");
    
    int dim = cmeshH1->Reference()->Dimension();
    
    std::string plotname;
    {
        std::stringstream out;
        out << config.dir_name << "/" << "H1_Problem" << config.porder << "_" << dim
        << "D_" << "Ndiv_ " << config.ndivisions << config.problemname<<".vtk";
        plotname = out.str();
    }
    int resolution=3;
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(resolution,dim);
    
    TPZManVector<REAL> errorvec(10, 0.);
    int64_t nelem = cmeshH1->NElements();
    cmeshH1->LoadSolution(cmeshH1->Solution());
    cmeshH1->ExpandSolution();
    cmeshH1->ElementSolution().Redim(nelem, 10);
    
    an.PostProcessError(errorvec);//Error calculation with exact and approximate solution
    errorseminorm = errorvec[2];
    
    std::cout << "Computed errors " << errorvec << std::endl;
    
    
    //Erro
    
    ofstream myfile;
    myfile.open("ArquivosErros_exacto.txt", ios::app);
    myfile << "\n\n Error for H1 formulation ";
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << config.ndivisions << " Order = " << config.porder << "\n";
    myfile << "DOF Total = " << cmeshH1->NEquations() << "\n";
    myfile << "Energy norm = " << errorvec[0] << "\n";//norma energia
    myfile << "error norm L2 = " << errorvec[1] << "\n";//norma L2
    myfile << "Semi norm H1 = " << errorvec[2] << "\n";//norma L2
    myfile.close();
    
    
}


TPZCompMesh *CompMeshH1(ProblemConfig &problem) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZDarcyFlow *mat = 0;
    
    
    for (auto matid : problem.materialids) {
        TPZDarcyFlow *mix = new TPZDarcyFlow(matid, cmesh->Dimension());
        mix->SetExactSol(problem.exact.operator*().ExactSolution(),3);
        mix->SetForcingFunction(problem.exact.operator*().ForceFunc(),3);
        
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
        
    }
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.);
        TPZManVector<REAL> val2(1, 0.);
        int bctype = 0;
        val2.Fill(0.);
        TPZBndCondT<STATE> *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->SetForcingFunctionBC(problem.exact.operator*().ExactSolution(),3);
        
        cmesh->InsertMaterialObject(bc);
    }
    
    cmesh->SetDefaultOrder(problem.porder);//ordem
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    cmesh->AutoBuild();
    
    
    return cmesh;
}

TPZGeoMesh *GeometricMesh(int nel, TPZVec<int> &bcids) {
    
    TPZManVector<int> nx(2,nel);
    TPZManVector<REAL> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    TPZGenGrid2D gen(nx, x0, x1, 1, 0);
    //TPZGenGrid2D gen(nx,x0,x1,0);
    
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, bcids[0]);
    gen.SetBC(gmesh, 5, bcids[1]);
    gen.SetBC(gmesh, 6, bcids[2]);
    gen.SetBC(gmesh, 7, bcids[3]);
    
    gmesh->SetDimension(2);
    
    
    
    return gmesh;
}


bool PostProcessing(TPZCompMesh * pressuremesh, TPZFMatrix<STATE> true_elerror, TPZFMatrix<STATE> estimate_elerror) {
    
    TPZLinearAnalysis an(pressuremesh);
    
    int64_t nels = pressuremesh->ElementVec().NElements();
    pressuremesh->ElementSolution().Redim(nels, 6);
    
    
    //Compute the effectivity index
    std::cout<<"***** Computing effectivity index *****"<<std::endl;
    
    {
        REAL sum = 0.;
        int64_t nel = true_elerror.Rows();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = pressuremesh->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            TPZCompEl *mphys = gel->Reference();
            int64_t elindex2 = mphys->Index();
            REAL aux =  estimate_elerror(elindex2,2) + estimate_elerror(elindex2,3);
            sum += aux*aux; // To compute global effectivity index
            true_elerror(el,0) = estimate_elerror(elindex2,2);
            true_elerror(el,1) = true_elerror(el,2);
            if (true_elerror(el,1) > 1.e-10) {
                true_elerror(el,2) = true_elerror(el,0)/true_elerror(el,1);
            }
        }
        pressuremesh->ElementSolution() = true_elerror;
        
        TPZManVector<REAL> errorsum(5, 0.);
        pressuremesh->EvaluateError(false, errorsum);
        REAL globeffind = sqrt(sum)/errorsum[2];
        std::cout << "Global full error estimate: "<<sqrt(sum)<<std::endl;
        std::cout << "Global effectivity index: " << globeffind << std::endl;
        
    }
    
    {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Solution");
        scalnames.Push("ExactSolution");
        scalnames.Push("EstimatedError");
        scalnames.Push("TrueError");
        scalnames.Push("EffectivityIndex");
        vecnames.Push("Flux");
        vecnames.Push("ExactFlux");

        std::string plotname;
        {
            std::stringstream out;
            out <<"ErrorEstimationH1_" << pressuremesh->GetDefaultOrder() << "_" << pressuremesh->Reference()->Dimension()
            <<  "Ndofs " << pressuremesh->NEquations() << ".vtk";
            plotname = out.str();
        }
        
        an.DefineGraphMesh(pressuremesh->Dimension(), scalnames, vecnames, plotname);
        int resolution = 3;
        an.PostProcess(resolution);
    }

    return true;
}

TPZGeoMesh *CreateGeoCircleMesh() {
    TPZGeoMesh * gmesh = nullptr;
    if (readGMeshFromFile) {
        TPZGmshReader gmsh;
        std::string meshfilename = "LCircle.msh";
        
        gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        
        gmsh.PrintPartitionSummary(std::cout);
        
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmesh->SetDimension(2);
    } else {
        gmesh = CreateLCircleGeoMesh();
    }
    int initialRefinement = 1;
    UniformRefinement(initialRefinement, gmesh);
    
#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream out("OriginalGeometricMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
#endif
    return gmesh;
}


TPZGeoMesh *CreateLCircleGeoMesh() {
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    
    TPZVec<REAL> coord(3, 0.);
    
    // Inserts node at origin
    gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[0].Initialize(coord, *gmesh);
    
    // Inserts circumference nodes
    for (int64_t i = 0; i < 13; i++) {
        const REAL step = M_PI / 8;
        coord[0] = cos(i * step);
        coord[1] = sin(i * step);
        const int64_t newID = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
    }
    
    int matIdTriangle = 1, matIdArc = 2;
    
    // Inserts triangle elements
    TPZManVector<int64_t, 3> nodesIdVec(3);
    for (int64_t i = 0; i < 6; i++) {
        nodesIdVec[0] = 0;
        nodesIdVec[1] = 1 + 2 * i;
        nodesIdVec[2] = 3 + 2 * i;
        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodesIdVec, matIdTriangle, *gmesh);
    }
    // Inserts arc elements
    for (int64_t i = 0; i < 6; i++) {
        nodesIdVec[0] = 1 + 2 * i;
        nodesIdVec[1] = 3 + 2 * i;
        nodesIdVec[2] = 2 + 2 * i;
        new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
    }
    // Finally, inserts line elements to complete boundary
    nodesIdVec.Resize(2);
    nodesIdVec[0] = 0;
    nodesIdVec[1] = 1;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);
    
    nodesIdVec[0] = 0;
    nodesIdVec[1] = 13;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);
    
    gmesh->BuildConnectivity();
    
    return gmesh;
}


TPZGeoMesh* CreateLShapeMesh(int nel, TPZVec<int>& bcids){
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;
    
    // Creates matrix with quadrilateral node coordinates.
    const int NodeNumber = 8;
    REAL coordinates[NodeNumber][3] = {
        {0., 0., 0.},
        {1., 0., 0.},
        {1., 1., 0.},
        {0., 1., 0.},
        {-1.,1.,0.},
        {-1.,0.,0.},
        {-1.,-1.,0.},
        {0.,-1.,0.}
    };
    
    // Inserts coordinates in the TPZGeoMesh object.
    for(int i = 0; i < NodeNumber; i++) {
        int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
        
        TPZVec<REAL> nodeCoord(3);
        nodeCoord[0] = coordinates[i][0];
        nodeCoord[1] = coordinates[i][1];
        nodeCoord[2] = coordinates[i][2];
        
        gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
    }
    
    // Creates quadrilateral element.
    int64_t index =0;
    TPZManVector<int64_t> nodeIDs(3);
    //El 0
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 3;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    //El 1
    nodeIDs[0] = 2;
    nodeIDs[1] = 3;
    nodeIDs[2] = 1;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 2
    nodeIDs[0] = 3;
    nodeIDs[1] = 4;
    nodeIDs[2] = 0;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 3
    nodeIDs[0] = 5;
    nodeIDs[1] = 0;
    nodeIDs[2] = 4;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    //El 4
    nodeIDs[0] = 0;
    nodeIDs[1] = 5;
    nodeIDs[2] = 7;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 6
    nodeIDs[0] = 6;
    nodeIDs[1] = 7;
    nodeIDs[2] = 5;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    // Creates line elements where boundary conditions will be inserted.
    nodeIDs.Resize(2);
    
    for (int i = 0; i < NodeNumber-1; i++) {
        
        nodeIDs[0] = i;
        
        nodeIDs[1] = (i + 1);
        std::cout<<"xo "<<nodeIDs[0]<<" x1 "<<nodeIDs[1]<<" bcid "<<bcids[i]<< "\n";
        
        gmesh->CreateGeoElement(EOned, nodeIDs, bcids[i], index);
    }
    index ++;
    
    nodeIDs[0] = 7;
    nodeIDs[1] = 0;
    std::cout<<"xo "<<nodeIDs[0]<<" x1 "<<nodeIDs[1]<<" bcid "<<bcids[NodeNumber-1]<< "\n";
    
    gmesh->CreateGeoElement(EOned, nodeIDs, bcids[NodeNumber-1], index);
    
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}

//Origin centered square
TPZGeoMesh* CreateSquareShapeMesh(int nel, TPZVec<int>& bcids){
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;
    
    // Creates matrix with quadrilateral node coordinates.
    const int NodeNumber = 9;
    REAL coordinates[NodeNumber][3] = {
        {0., 0., 0.},
        {1., 0., 0.},
        {1., 1., 0.},
        {0., 1., 0.},
        {-1.,1.,0.},
        {-1.,0.,0.},
        {-1.,-1.,0.},
        {0.,-1.,0.},
        {1.,-1.,0.}
    };
    
    // Inserts coordinates in the TPZGeoMesh object.
    for(int i = 0; i < NodeNumber; i++) {
        int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
        
        TPZVec<REAL> nodeCoord(3);
        nodeCoord[0] = coordinates[i][0];
        nodeCoord[1] = coordinates[i][1];
        nodeCoord[2] = coordinates[i][2];
        
        gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
    }
    
    // Creates triangular element.
    int64_t index =0;
    TPZManVector<int64_t> nodeIDs(3);
    //El 0
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 3;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    //El 1
    nodeIDs[0] = 2;
    nodeIDs[1] = 3;
    nodeIDs[2] = 1;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 2
    nodeIDs[0] = 3;
    nodeIDs[1] = 4;
    nodeIDs[2] = 0;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 3
    nodeIDs[0] = 5;
    nodeIDs[1] = 0;
    nodeIDs[2] = 4;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    //El 4
    nodeIDs[0] = 0;
    nodeIDs[1] = 5;
    nodeIDs[2] = 7;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 5
    nodeIDs[0] = 6;
    nodeIDs[1] = 7;
    nodeIDs[2] = 5;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 6
    nodeIDs[0] = 7;
    nodeIDs[1] = 8;
    nodeIDs[2] = 0;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 7
    nodeIDs[0] = 1;
    nodeIDs[1] = 0;
    nodeIDs[2] = 8;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    // Creates line elements where boundary conditions will be inserted.
    nodeIDs.Resize(2);
    
    for (int i = 1; i < NodeNumber-1; i++) {
        
        nodeIDs[0] = i;
        
        nodeIDs[1] = (i + 1);
        std::cout<<"xo "<<nodeIDs[0]<<" x1 "<<nodeIDs[1]<<" bcid "<<bcids[i]<< "\n";
        
        gmesh->CreateGeoElement(EOned, nodeIDs, bcids[i], index);
    }
    index ++;
    
    nodeIDs[0] = 8;
    nodeIDs[1] = 1;
    std::cout<<"xo "<<nodeIDs[0]<<" x1 "<<nodeIDs[1]<<" bcid "<<bcids[NodeNumber-2]<< "\n";
    
    gmesh->CreateGeoElement(EOned, nodeIDs, bcids[NodeNumber-2], index);
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}

//Square in the first quadrant
TPZGeoMesh* CreateSquareShapeMesh2(int nel, TPZVec<int>& bcids){
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;
    
    // Creates matrix with quadrilateral node coordinates.
    const int NodeNumber = 9;
    REAL coordinates[NodeNumber][3] = {
        {0., 0., 0.},
        {0.5, 0., 0.},
        {1.,0.,0.},
        {1.,0.5,0.},
        {1.,1.,0.},
        {0.5,1.,0.},
        {0.,1,0.},
        {0., 0.5, 0.},
        {0.5, 0.5, 0.},
    };
    
    // Inserts coordinates in the TPZGeoMesh object.
    for(int i = 0; i < NodeNumber; i++) {
        int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
        
        TPZVec<REAL> nodeCoord(3);
        nodeCoord[0] = coordinates[i][0];
        nodeCoord[1] = coordinates[i][1];
        nodeCoord[2] = coordinates[i][2];
        
        gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
    }
    
    // Creates triangular element.
    int64_t index =0;
    TPZManVector<int64_t> nodeIDs(3);
    //El 0
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 7;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    //El 1
    nodeIDs[0] = 8;
    nodeIDs[1] = 7;
    nodeIDs[2] = 1;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 2
    nodeIDs[0] = 1;
    nodeIDs[1] = 2;
    nodeIDs[2] = 8;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 3
    nodeIDs[0] = 2;
    nodeIDs[1] = 3;
    nodeIDs[2] = 8;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    
    //El 4
    nodeIDs[0] = 3;
    nodeIDs[1] = 4;
    nodeIDs[2] = 5;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 5
    nodeIDs[0] = 8;
    nodeIDs[1] = 3;
    nodeIDs[2] = 5;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 6
    nodeIDs[0] = 8;
    nodeIDs[1] = 5;
    nodeIDs[2] = 6;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    //El 7
    nodeIDs[0] = 7;
    nodeIDs[1] = 8;
    nodeIDs[2] = 6;
    gmesh->CreateGeoElement(ETriangle, nodeIDs, matID, index);
    index++;
    // Creates line elements where boundary conditions will be inserted.
    nodeIDs.Resize(2);
    
    for (int i = 0; i < NodeNumber-2; i++) {
        
        nodeIDs[0] = i;
        
        nodeIDs[1] = (i + 1);
        std::cout<<"xo "<<nodeIDs[0]<<" x1 "<<nodeIDs[1]<<" bcid "<<bcids[i]<< "\n";
        
        gmesh->CreateGeoElement(EOned, nodeIDs, bcids[i], index);
    }
    index ++;
    
    nodeIDs[0] = 7;
    nodeIDs[1] = 0;
    std::cout<<"xo "<<nodeIDs[0]<<" x1 "<<nodeIDs[1]<<" bcid "<<bcids[NodeNumber-2]<< "\n";
    
    gmesh->CreateGeoElement(EOned, nodeIDs, bcids[NodeNumber-2], index);
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}

double Maximum(TPZVec<double> &vect){
    double max = vect[0];
    for(int i=1; i<3; i++){
        if(max < vect[i]){
            max = vect[i];
        }
    }
    return max;
}
