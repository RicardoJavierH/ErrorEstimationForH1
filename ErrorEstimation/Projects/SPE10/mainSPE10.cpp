//
// Created by Gustavo Batistela on 4/5/21.
//

#include "Tools.h"
#include <TPZGenGrid2D.h>
#include <TPZGenGrid3D.h>
#include <TPZMFSolutionTransfer.h>
#include <ToolsMHM.h>
#include <iostream>
#include <libInterpolate/Interpolate.hpp>
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include <memory>
#include <pzgmesh.h>

typedef _2D::BicubicInterpolator<REAL> Interpolator;

// Global variables
Interpolator interpolator;

// Function declarations
void ReadSPE10CellPermeabilities(TPZVec<REAL>*perm_vec, int layer);
TPZGeoMesh *CreateSPE10GeoMesh();
STATE PermeabilityFunction(const TPZVec<REAL> &x);
void InsertMaterials(TPZCompMesh *cmesh);
void CreateSPE10MHMCompMesh(TPZMHMixedMeshControl &mhm);
void SolveMHMProblem(TPZMHMixedMeshControl &mhm);
void EstimateError(TPZMHMixedMeshControl *mhm);

int main() {

    constexpr int layer = 36;
    constexpr int nx = 220;
    constexpr int ny = 60;
    constexpr int n_cells = nx * ny;

    auto perm_vec = TPZManVector<REAL, n_cells>(n_cells, 1);
    ReadSPE10CellPermeabilities(&perm_vec, layer);

    TPZGeoMesh *gmesh = CreateSPE10GeoMesh();
    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";
    Tools::PrintGeometry(gmesh, "SPE10GeoMesh", false, true);

    std::vector<REAL> x, y, perm;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            const int cell_id = ny * i + j;
            const double cell_perm = perm_vec[cell_id];
            x.push_back(0.5 + i);
            y.push_back(0.5 + j);
            perm.push_back(cell_perm);
        }
    }

    interpolator.setData(x.size(), x.data(), y.data(), perm.data());

    TPZMHMixedMeshControl mhm(gmesh);
    CreateSPE10MHMCompMesh(mhm);

    SolveMHMProblem(mhm);
    EstimateError(&mhm);

    return 0;
}

TPZGeoMesh *CreateSPE10GeoMesh() {
    std::cout << "Creating SPE10 initial grid...\n";

    const TPZManVector<REAL, 3> x0 = {0, 0, 0};
    const TPZManVector<REAL, 3> x1 = {220., 60., 0.};
    const TPZManVector<int, 3> ndiv = {20, 6, 0};

    TPZGenGrid2D gen(ndiv, x0, x1);

    gen.SetRefpatternElements(true);
    auto gmesh = new TPZGeoMesh;
    gen.Read(gmesh);

    gen.SetBC(gmesh, 5, -1);

    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    return gmesh;
}

void ReadSPE10CellPermeabilities(TPZVec<REAL> *perm_vec, const int layer) {

    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    int cell_id = 0;
    const int n_cells = perm_vec->size();
    const int start_line = 1 + n_cells * (layer - 1) / 6;

    int line_num = 0;
    int line_num2 = 0;
    while (perm_file) {
        line_num++;
        line_num2++;
        std::string line;
        std::getline(perm_file, line, '\n');

        if (line_num < start_line) continue;

        std::stringstream stream(line);
        for (int i = 0; i < 6; i++) {
            stream >> perm_vec->operator[](cell_id);
            cell_id++;
        }
        if (cell_id == n_cells) break;
    }
    std::cout << "Finished reading permeability data from input file!\n";
}

STATE PermeabilityFunction(const TPZVec<REAL> &x) {

    STATE perm;
    for (int i = 0; i < 2; i++) {
        perm = interpolator(x[0], x[1]);
        if (perm <= 1) {
            perm = 1;
        } else {
            perm += 1;
        }
    }
    std::cout << "[" << x[0] << ", " << x[1] << "]\n";
    std::cout << " perm = " << perm << std::endl;
    return perm;
    //std::cout << "[" << x[0] << ", " << x[1] << "]\n";
    //std::cout << "[" << res_mat(0, 0) << " " << res_mat(0, 1) << "\n";
    //std::cout        << res_mat(1, 0) << " " << res_mat(1, 1) << "\n";
    //std::cout        << res_mat(2, 0) << " " << res_mat(2, 1) << "\n";
    //std::cout        << res_mat(3, 0) << " " << res_mat(3, 1) << "]\n\n";
}

void CreateSPE10MHMCompMesh(TPZMHMixedMeshControl &mhm) {

    TPZGeoMesh *gmesh = mhm.GMesh().operator->();
    TPZManVector<int64_t, 22 * 6> coarse_indexes;
    ComputeCoarseIndices(gmesh, coarse_indexes);

    int nInternalRef = 3;
    Tools::UniformRefinement(nInternalRef, 2, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    mhm.DefinePartitionbyCoarseIndices(coarse_indexes);

    // Indicate material indices to the MHM control structure
    mhm.fMaterialIds = {1};
    mhm.fMaterialBCIds = {-1, -2, -3};

    // Insert the material objects in the multiphysics mesh
    TPZCompMesh *cmesh = mhm.CMesh().operator->();
    InsertMaterials(cmesh);

    // General approximation order settings
    mhm.SetInternalPOrder(1);
    mhm.SetSkeletonPOrder(1);
    mhm.SetHdivmaismaisPOrder(1);

    // Refine skeleton elements
    mhm.DivideSkeletonElements(0);
    mhm.DivideBoundarySkeletonElements();

    // Creates MHM mesh
    bool substructure = true;
    mhm.BuildComputationalMesh(substructure);
    //{
    //    std::string fileName = "CompMesh.txt";
    //    std::ofstream file(fileName);
    //    mhm.CMesh()->Print(file);
    //}
}

void SolveMHMProblem(TPZMHMixedMeshControl &mhm) {

    TPZAutoPointer<TPZCompMesh> cmesh = mhm.CMesh();

    bool should_renumber = true;
    TPZLinearAnalysis an(cmesh, should_renumber);

#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh.operator->());
    strmat.SetNumThreads(8);
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(8);
#endif

    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    std::cout << "Assembling\n";
    an.Assemble();

    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution();

    TPZMFSolutionTransfer transfer;
    transfer.BuildTransferData(cmesh.operator->());
    transfer.TransferFromMultiphysics();

    TPZStack<std::string> scalnames, vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }

    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    vecnames.Push("Flux");

    int resolution = 0;
    std::string plotname = "SPE10-Results.vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname);
    an.PostProcess(resolution, cmesh->Dimension());
}

void InsertMaterials(TPZCompMesh *cmesh) {

    typedef TPZMixedDarcyFlow TPZMixedPoisson;
    auto *mix = new TPZMixedPoisson(1, cmesh->Dimension());
    PermeabilityFunctionType a = PermeabilityFunction;
    mix->SetPermeabilityFunction(a);

    TPZFNMatrix<1, REAL> val1(1, 1, 0.);
    TPZManVector<REAL> val2(1, 0.);
    constexpr int dirichlet_bc = 0;

    // Pressure at reservoir boundary
    val2[0] = 1;
    TPZBndCond *pressure_left = mix->CreateBC(mix, -1, dirichlet_bc, val1, val2);

    cmesh->InsertMaterialObject(mix);
    cmesh->InsertMaterialObject(pressure_left);
}

void EstimateError(TPZMHMixedMeshControl *mhm) {

    std::cout << "\nError Estimation processing for MHM-Hdiv problem " << std::endl;

    // Error estimation
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!originalMesh) DebugStop();

    bool postProcWithHDiv = true;
    TPZMHMHDivErrorEstimator ErrorEstimator(*originalMesh, mhm, postProcWithHDiv);
    ErrorEstimator.PotentialReconstruction();

    std::string command = "mkdir SPE10";
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL> elementerrors;
    std::string outVTK = "SPE10-Errors.vtk";
    ErrorEstimator.ComputeErrors(errors, elementerrors, outVTK);
}
