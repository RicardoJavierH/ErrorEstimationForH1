//
//  TPZCreateHDivMesh.cpp
//  PZ
//
//  Created by Philippe Devloo on 6/21/16.
//
//

#include "TPZCreateMesh.h"
#include "pzbndcond.h"
#include "pzpoisson3d.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "mixedpoisson.h"
#include "problem.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "pzintel.h"
#include "pzcheckmesh.h"
#include "tpzcompmeshreferred.h"
#include "TPZMixedErrorEstimate.h"

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    //  TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    
    cmesh->SetDefaultOrder(pOrder);
    
    
    
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    return cmesh;
    
}

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim, bool disconnected, bool referred)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    material->NStateVariables();
    
    if(!referred)
    {
        //solucao exata
        TPZAutoPointer<TPZFunction<STATE> > solexata;
        int polorder = 4;
        TPZDummyFunction<STATE> *func = new TPZDummyFunction<STATE>(RightTermArcTangentBad);
        func->SetPolynomialOrder(polorder);
        solexata = func;
        material->SetForcingFunction(solexata);
    }
    
    //    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
    
    TPZCompMesh * cmesh;
    if (referred) {
        cmesh = new TPZCompMeshReferred(gmesh);
    }
    else
    {
        cmesh = new TPZCompMesh(gmesh);
    }
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    if (!referred)
    {
        TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
        TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
        TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
        TPZMaterial * BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
        TPZMaterial * BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);

        cmesh->InsertMaterialObject(BCond0);
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
        cmesh->InsertMaterialObject(BCond4);
        cmesh->InsertMaterialObject(BCond5);
    }
    cmesh->SetDefaultOrder(pOrder);
    //cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsContinuous();
    if (disconnected)
    {
        cmesh->ApproxSpace().CreateDisconnectedElements(disconnected);
    }
    if (referred)
    {
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuousReferred();
    }
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
        newnod.SetLagrangeMultiplier(1);
    }
    
    
#ifdef ERRORESTIMATION_DEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_2 pressure\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str());
    //	}
    //#endif
    return cmesh;
    
}


#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> &meshvec)
{
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim = gmesh->Dimension();
    
    
    //TPZMatPoissonD3 *material = new TPZMatPoissonD3(matId,dim); hdivantigo = true; // nesse material tem que ser true
    TPZMixedPoisson *material = new TPZMixedErrorEstimate<TPZMixedPoisson>(matId,dim);
    
    
    //incluindo os dados do problema
    TPZFNMatrix<9,REAL> PermTensor(3,3,0.);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3,0.);
    
    for (int i=0; i<3; i++) {
        PermTensor(i,i) = 1.;
        InvPermTensor(i,i) = 1.;
    }
    
    material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    int polorder = 4;
    TPZDummyFunction<STATE> *func = new TPZDummyFunction<STATE>(RightTermArcTangent);
    func->SetPolynomialOrder(polorder);
    solexata = func;
    material->SetForcingFunction(solexata);
    mphysics->SetDimModel(dim);
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    
    //Criando condicoes de contorno
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    
    
    
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    //    BCond1 = material->CreateBC(mat, bc1,neumann, val1, val2);
    
    BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    //    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2N);
    //    BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    //    BCond2->SetForcingFunction(FBCond2);
    
    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    //
    
    BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    //    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4N);
    //    BCond4 = material->CreateBC(mat, bc4,neumann, val1, val2);
    //    BCond4->SetForcingFunction(FBCond4);
    
    
    
    BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    mphysics->InsertMaterialObject(BCond5);
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    
    //        TPZMaterial * skeletonEl = material->CreateBC(mat, matskeleton, 3, val1, val2);
    //        mphysics->InsertMaterialObject(skeletonEl);
    
    //        TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(matskeleton, dim-1, 1);
    //        TPZMaterial * mat2(matskelet);
    //        mphysics->InsertMaterialObject(mat2);
    
    meshvec[0]->CleanUpUnconnectedNodes();
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    //------- Create and add group elements -------
    
    
    return mphysics;
    
}

/// adjust the polynomial orders of the hdiv elements such that the internal order is higher than the sideorders
void AdjustFluxPolynomialOrders(TPZCompMesh *fluxmesh, int hdivplusplus)
{
    int dim = fluxmesh->Dimension();
    /// loop over all the elements
    int64_t nel = fluxmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        // compute the maxorder
        int maxorder = -1;
        int ncon = intel->NConnects();
        for (int i=0; i<ncon-1; i++) {
            int conorder = intel->Connect(i).Order();
            maxorder = maxorder < conorder ? conorder : maxorder;
        }
        int nsides = gel->NSides();
        int nconside = intel->NSideConnects(nsides-1);
        // tive que tirar para rodar H1
//        if (nconside != 1 || maxorder == -1) {
//            DebugStop();
//        }
        int64_t cindex = intel->SideConnectIndex(nconside-1, nsides-1);
        TPZConnect &c = fluxmesh->ConnectVec()[cindex];
        if (c.NElConnected() != 1) {
            DebugStop();
        }
        if (c.Order()+hdivplusplus != maxorder) {
//            std::cout << "Changing the order of the central connect " << cindex << " from " << c.Order() << " to " << maxorder+hdivplusplus << std::endl;
            // change the internal connect order to be equal do maxorder
            intel->SetSideOrder(nsides-1, maxorder+hdivplusplus);
        }
    }
    fluxmesh->ExpandSolution();
}

/// set the pressure order acording to the order of internal connect of the elements of the fluxmesh
void SetPressureOrders(TPZCompMesh *fluxmesh, TPZCompMesh *pressuremesh)
{
    // build a vector with the required order of each element in the pressuremesh
    // if an element of the mesh dimension of the fluxmesh does not have a corresponding element in the pressuremesh DebugStop is called
    int meshdim = fluxmesh->Dimension();
    pressuremesh->Reference()->ResetReference();
    pressuremesh->LoadReferences();
    TPZManVector<int64_t> pressorder(pressuremesh->NElements(),-1);
    int64_t nel = fluxmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != meshdim) {
            continue;
        }
        int nsides = gel->NSides();
        int64_t cindex = intel->SideConnectIndex(0, nsides-1);
        TPZConnect &c = fluxmesh->ConnectVec()[cindex];
        int order = c.Order();
        TPZCompEl *pressureel = gel->Reference();
        TPZInterpolatedElement *pintel = dynamic_cast<TPZInterpolatedElement *>(pressureel);
        if (!pintel) {
            DebugStop();
        }
        pressorder[pintel->Index()] = order;
    }
    pressuremesh->Reference()->ResetReference();
    nel = pressorder.size();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        TPZInterpolatedElement *pintel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!pintel) {
            continue;
        }
        if (pressorder[el] == -1) {
            continue;
        }
        pintel->PRefine(pressorder[el]);
    }

    pressuremesh->ExpandSolution();
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

void CreateAllMeshes(TPZGeoMesh *gmesh, TPZVec<TPZCompMesh *> &meshvec, int porder, int dim, int hdivplusplus)
{
    meshvec.resize(5);
    // the HDiv flux mesh used to reconstruct the internal fluxes
    TPZCompMesh *fluxmesh = CMeshFlux(gmesh,porder,dim);
    meshvec[2] = fluxmesh;
    bool disconnected = false;
    bool referred = false;
    // this is the continuous pressure mesh for which we will test the error estimator
    TPZCompMesh *pressuremeshH1 = CMeshPressure(gmesh, porder, dim, disconnected,referred);
    meshvec[0] = pressuremeshH1;
    disconnected = true;
    referred = true;
    TPZCompMesh *pressuremeshHDiv = CMeshPressure(gmesh, porder, dim, disconnected,referred);
    meshvec[3] = pressuremeshHDiv;
    AdjustFluxPolynomialOrders(fluxmesh, hdivplusplus);
    SetPressureOrders(fluxmesh, pressuremeshHDiv);
    TPZManVector<TPZCompMesh *> meshvecloc(2);
    meshvecloc[0] = fluxmesh;
    meshvecloc[1] = pressuremeshHDiv;
    TPZCompMesh *mixed = CMeshMixed(gmesh, meshvecloc);
    meshvec[1] = mixed;
    
    disconnected = false;
    referred = true;
    TPZCompMesh *PartitionUnity = CMeshPressure(gmesh, 1, dim, disconnected, referred);
    TPZCompMeshReferred *PartitionRef = dynamic_cast<TPZCompMeshReferred *>(PartitionUnity);
    if (!PartitionRef) {
        DebugStop();
    }
    PartitionRef->LoadReferred(pressuremeshH1);
    TPZCompMeshReferred *refcomp = dynamic_cast<TPZCompMeshReferred *>(pressuremeshHDiv);
    refcomp->LoadReferred(PartitionUnity);
    // pressureHDiv -> PartitionUnity -> pressureH1
    meshvec[4] = PartitionUnity;
//    TPZCompMeshTools::GroupElements(mixed);
//    TPZCompMeshTools::CreatedCondensedElements(mixed, true);
}

/// Reconstruct the multi physics mesh from the fluxmesh
void ReconstructHDivMesh(TPZVec<TPZCompMesh *> &meshvec, int hdivplusplus)
{
    TPZCompMesh *HDivMesh = meshvec[1];
    TPZCompMesh *fluxmesh = meshvec[2];
    TPZCompMesh *pressuremesh = meshvec[3];
    TPZGeoMesh *gmesh = meshvec[0]->Reference();
    {
        std::map<int ,TPZMaterial * > HDivMat = HDivMesh->MaterialVec();
        HDivMesh->MaterialVec().clear();
        HDivMesh->CleanUp();
        HDivMesh->MaterialVec() = HDivMat;
        HDivMesh->CleanUpUnconnectedNodes();
    }
    {
        std::map<int ,TPZMaterial * > PressMat = meshvec[1]->MaterialVec();
        pressuremesh->MaterialVec().clear();
        gmesh->ResetReference();
        int64_t nel = meshvec[1]->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = pressuremesh->Element(el);
            if (cel) {
                delete cel;
            }
        }
        pressuremesh->CleanUp();
        pressuremesh->MaterialVec() = PressMat;
        pressuremesh->CleanUpUnconnectedNodes();
    }
    int meshdim = gmesh->Dimension();
    gmesh->ResetReference();
    fluxmesh->LoadReferences();
    fluxmesh->CleanUpUnconnectedNodes();
    fluxmesh->ExpandSolution();
    AdjustFluxPolynomialOrders(fluxmesh, hdivplusplus);
    {
        //        std::ofstream out("checkmesh.txt");
        gmesh->ResetReference();
        fluxmesh->LoadReferences();
        TPZCheckMesh check(fluxmesh,&std::cout);
        check.VerifyAllConnects();
    }
    
    // lookup which elements need to be created and with which order
    TPZVec<int64_t> loadels(gmesh->NElements(),-1);
    for (int64_t el=0; el<gmesh->NElements(); el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) {
            continue;
        }
        TPZCompEl *cel = gel->Reference();
        if (!cel) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            DebugStop();
        }
        int maxorder = intel->MaxOrder();
        loadels[el] = maxorder-1;
    }
    
    // create the pressure mesh
    gmesh->ResetReference();
    meshvec[1]->LoadReferences();
    for (int64_t el=0; el<gmesh->NElements(); el++) {
        if (loadels[el] != -1) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->Dimension() != gmesh->Dimension()) {
                continue;
            }
            meshvec[1]->SetDefaultOrder(loadels[el]);
            int64_t index;
            TPZCompEl *cel = pressuremesh->ApproxSpace().CreateCompEl(gmesh->Element(el), *meshvec[1], index);
            cel->Reference()->ResetReference();
        }
    }
    meshvec[1]->ExpandSolution();
    int64_t nconn = meshvec[1]->NConnects();
    for (int64_t ic = 0; ic<nconn; ic++) {
        pressuremesh->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    // create the multiphysics elements
    gmesh->ResetReference();
    HDivMesh->LoadReferences();
    for (int64_t el=0; el<gmesh->NElements(); el++) {
        if (loadels[el] != -1) {
            int64_t index;
            HDivMesh->ApproxSpace().CreateCompEl(gmesh->Element(el), *HDivMesh, index);
        }
    }
    
    TPZManVector<TPZCompMesh *,2> locmeshvec(2);
    locmeshvec[0] = meshvec[2];
    locmeshvec[1] = meshvec[3];
    
    TPZBuildMultiphysicsMesh::AddElements(locmeshvec, HDivMesh);
    
    TPZBuildMultiphysicsMesh::AddConnects(locmeshvec,HDivMesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(locmeshvec, HDivMesh);
    
    HDivMesh->ExpandSolution();
    TPZCompMeshTools::GroupElements(HDivMesh);
    HDivMesh->ComputeNodElCon();
    if(0)
    {
        std::ofstream out("ReconstructFlux.txt");
        meshvec[0]->Print(out);
    }
    if(0)
    {
        std::ofstream out("ReconstructPressure.txt");
        meshvec[1]->Print(out);
    }
    if(0)
    {
        std::ofstream out("Reconstruct.txt");
        HDivMesh->Print(out);
    }
    TPZCompMeshTools::CreatedCondensedElements(HDivMesh, true);
    
}

void UnitPressure(TPZCompMesh *PressMesh)
{
    int64_t nel=PressMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = PressMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int ncorner = gel->NCornerNodes();
        for (int ic=0; ic<ncorner; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = PressMesh->Block().Position(seqnum);
            PressMesh->Solution()(pos,0) = 1.;
        }
    }
}

