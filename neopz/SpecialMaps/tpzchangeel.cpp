/**
 * @file
 * @brief Contains the implementation of the TPZChangeEl methods. 
 */
#include "tpzchangeel.h"

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif
#include <iostream>
#include <cstdlib>

#include "tpzmathtools.h"

#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticpyramid.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticcube.h"
#include "tpzarc3d.h"
#include "TPZCylinderMap.h"
#include "TPZGeoElement.h"
#include "pzgeoelside.h"
#include "pzstack.h"
#include "tpzgeoelrefpattern.h"

#include <sstream>
using namespace std;
using namespace pzgeom;
using namespace pztopology;


TPZChangeEl::TPZChangeEl()
{
}
//------------------------------------------------------------------------------------------------------------

TPZChangeEl::~TPZChangeEl()
{
}
//------------------------------------------------------------------------------------------------------------

//#define verifyNeighbourhood
TPZGeoEl * TPZChangeEl::ChangeToQuadratic(TPZGeoMesh *Mesh, int64_t ElemIndex)
{
	TPZGeoEl * OldElem = Mesh->ElementVec()[ElemIndex];
    
    /////////////////////////
#ifdef verifyNeighbourhood
    std::ofstream before("before.txt");
    for(int s = 0; s < OldElem->NSides(); s++)
    {
        TPZGeoElSide oldSide(OldElem,s);
        TPZGeoElSide neighSide(oldSide.Neighbour());
        while(oldSide != neighSide)
        {
            before << s << "\t" << neighSide.Element()->Id() << "\t" << neighSide.Side() << "\n";
            neighSide = neighSide.Neighbour();
        }
    }
    before.close();
    TPZGeoEl * oldFather = OldElem->Father();
    int oldMePosition = -1;
    if(oldFather)
    {
        for(int s = 0; s < oldFather->NSubElements(); s++)
        {
            if(oldFather->SubElement(s) == OldElem)
            {
                oldMePosition = s;
                break;
            }
        }
    }
#endif
    /////////////////////////
    
#ifdef PZDEBUG
    if(!OldElem)
    {
        std::cout << "Null geoel on " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoEl * father = OldElem->Father();
    int which_subel = -1;
    if (father){
        which_subel = OldElem->WhichSubel();
    }
    
    int64_t midN;
    int nsides = OldElem->NSides();
    TPZVec<TPZGeoElSide> oldNeigh(nsides);
    StoreNeighbours(OldElem, oldNeigh);
    
    TPZVec<int64_t> NodesSequence(0);
    for(int s = 0; s < nsides; s++)
    {
        TPZGeoElSide mySide(OldElem,s);
        if(mySide.Dimension() == 0)
        {
            int64_t oldSz = NodesSequence.NElements();
            NodesSequence.resize(oldSz+1);
            NodesSequence[oldSz] = OldElem->NodeIndex(s);
        }
        if(CreateMiddleNodeAtEdge(Mesh, ElemIndex, s, midN))
        {
            int64_t oldSz = NodesSequence.NElements();
            NodesSequence.resize(oldSz+1);
            NodesSequence[oldSz] = midN;
        }
    }
    
    MElementType elType = OldElem->Type();
    int64_t oldId = OldElem->Id();
    int64_t oldMatId = OldElem->MaterialId();
    
	TPZGeoEl * NewElem = NULL;
    
    /** Deleting OldElem */
    Mesh->DeleteElement(OldElem);

    switch(elType) /** Inserting New Element in Mesh */
    {
        case(EOned) :
        {             
            NewElem = new TPZGeoElRefPattern<TPZQuadraticLine>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }  
        case(ETriangle) :
        {            
            NewElem = new TPZGeoElRefPattern<TPZQuadraticTrig>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(EQuadrilateral) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticQuad>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(ETetraedro) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticTetra>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(EPiramide) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticPyramid>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }    
        case(EPrisma) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticPrism>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(ECube) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticCube>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        default :
        {
            DebugStop();
            break;
        }
    }
    
#ifdef PZDEBUG
    if (!Mesh->Element(ElemIndex)) {
        DebugStop();
    }
#endif
    
    if(father)
    {
        NewElem->SetFather(father);
        father->SetSubElement(which_subel,NewElem);
    }


    RestoreNeighbours(NewElem, oldNeigh);
    
    if(NewElem->HasSubElement())
    {
        //Mudar subelementos para TPZGeoElMapped
    }
    
    /////////////////////////
#ifdef verifyNeighbourhood
    std::ofstream after("after.txt");
    for(int s = 0; s < NewElem->NSides(); s++)
    {
        TPZGeoElSide newSide(NewElem,s);
        TPZGeoElSide neighSide(newSide.Neighbour());
        while(newSide != neighSide)
        {
            after << s << "\t" << neighSide.Element()->Id() << "\t" << neighSide.Side() << "\n";
            neighSide = neighSide.Neighbour();
        }
    }
    after.close();
    TPZGeoEl * newFather = NewElem->Father();
    int newMePosition = -1;
    if(newFather)
    {
        for(int s = 0; s < newFather->NSubElements(); s++)
        {
            if(newFather->SubElement(s) == NewElem)
            {
                newMePosition = s;
                break;
            }
        }
    }
    if(oldFather != newFather || oldMePosition != newMePosition)
    {
        DebugStop();
    }
#endif
    /////////////////////////
    
	return NewElem;
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl* TPZChangeEl::ChangeToQuarterPoint(TPZGeoMesh *Mesh, int64_t ElemIndex, int targetSide)
{
    TPZGeoEl * gel = Mesh->ElementVec()[ElemIndex];
    
    #ifdef PZDEBUG
    if(!gel)
    {
        std::cout << "Null geoel on " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
    #endif
    
    if(gel->IsLinearMapping())
    {
        gel = TPZChangeEl::ChangeToQuadratic(Mesh, ElemIndex);
    }
    
    MElementType gelType = gel->Type();
    
    TPZStack<int> targetSideNodes, targetSideEdges;
    
    switch(gelType)
    {
        case (EPoint):
        {
            TPZPoint::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZPoint::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (EOned):
        {
            TPZLine::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZLine::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (ETriangle):
        {
            TPZTriangle::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZTriangle::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (EQuadrilateral):
        {
            TPZQuadrilateral::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZQuadrilateral::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (ETetraedro):
        {
            TPZTetrahedron::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZTetrahedron::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (EPiramide):
        {
            TPZPyramid::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZPyramid::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (EPrisma):
        {
            TPZPrism::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZPrism::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (ECube):
        {
            TPZCube::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZCube::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
        default:
        {
            std::cout << "Element type not found on " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
    }
    
    std::set<int> targetSideNodes_set, targetSideEdges_set, NOTtargetSideEdges_set;
    TPZGeoElSide targetSideEl(gel,targetSide);
    for(int nd = 0; nd < targetSideNodes.NElements(); nd++)
    {
        int tsNode = targetSideNodes[nd];
        targetSideNodes_set.insert(tsNode);
    }
    if(targetSideEl.Dimension() == 0)
    {
        targetSideNodes_set.insert(targetSide);
    }
    for(int ed = 0; ed < targetSideEdges.NElements(); ed++)
    {
        int tsEdge = targetSideEdges[ed];
        targetSideEdges_set.insert(tsEdge);
    }
    if(targetSideEl.Dimension() == 1)
    {
        targetSideEdges_set.insert(targetSide);
    }
    for(int sd = gel->NCornerNodes(); sd < gel->NSides(); sd++)
    {
        TPZGeoElSide gelSide(gel,sd);
        if(gelSide.Dimension() == 1 && targetSideEdges_set.find(sd) == targetSideEdges_set.end())
        {
            NOTtargetSideEdges_set.insert(sd);
        }
        else if(gelSide.Dimension() > 1)
        {
            break;
        }
    }
    
    const REAL dist = 0.25;
    
    std::set<int>::iterator it;
    for(it = NOTtargetSideEdges_set.begin(); it != NOTtargetSideEdges_set.end(); it++)
    {
        int edg = *it;
        TPZGeoElSide gelSide(gel,edg);
        
        int64_t initNode = gelSide.SideNodeLocIndex(0);
        int64_t finalNode = gelSide.SideNodeLocIndex(1); 
        int64_t meshInitNode = gelSide.SideNodeIndex(0);
        int64_t meshFinalNode = gelSide.SideNodeIndex(1);
        /**
         * Embora o elemento quadratico possua mais nohs (NNodes), a topologia segue igual ao
         * elemento linear no qual o elemento quadratico foi baseado.
         *
         * Isso quer dizer que Linear(Geom::NNodes) < Quadratic(Geom::NNodes), mas Linear(Geom::NSides) = Quadratic(Geom::NSides)
         *
         * Com isso a numeracao dos nohs nos meios das arestas coincide com a numeracao dos lados unidimensionais do elemento.
         */
        int64_t meshMiddleNode = gel->NodeIndex(edg);
        //********************
                                            
        double coordNear, coordFar;
        
        if(targetSideNodes_set.find(initNode) != targetSideNodes_set.end())//drag middle node to node 0 of this edge
        {
            for(int c = 0; c < 3; c++)
            {
                coordNear = Mesh->NodeVec()[meshInitNode].Coord(c);
                coordFar = Mesh->NodeVec()[meshFinalNode].Coord(c);
                Mesh->NodeVec()[meshMiddleNode].SetCoord(c, (1. - dist)*coordNear + (dist)*coordFar);
            }
        }
        else if(targetSideNodes_set.find(finalNode) != targetSideNodes_set.end())//drag middle node to node 1 of this edge
        {
            for(int c = 0; c < 3; c++)
            {
                coordNear = Mesh->NodeVec()[meshFinalNode].Coord(c);
                coordFar = Mesh->NodeVec()[meshInitNode].Coord(c);
                Mesh->NodeVec()[meshMiddleNode].SetCoord(c, (1. - dist)*coordNear + (dist)*coordFar);
            }
        }
    }				
    
    return gel;
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl * TPZChangeEl::ChangeToGeoBlend(TPZGeoMesh *Mesh, int64_t ElemIndex)
{
	TPZGeoEl * OldElem = Mesh->ElementVec()[ElemIndex];
	if(!OldElem)
    {
		PZError << "Error at " << __PRETTY_FUNCTION__ << " - NULL geometric element.\n";
		return NULL;
	}
    const MElementType oldType = OldElem->Type();
    int64_t oldId = OldElem->Id();
    const int64_t oldMatId = OldElem->MaterialId();
    const int nsides = OldElem->NSides();
    
    TPZVec<TPZGeoElSide> oldNeigh(nsides);
    StoreNeighbours(OldElem, oldNeigh);

    
	
    const int nnodes = OldElem->NCornerNodes();
    TPZManVector<int64_t> nodeindexes(nnodes);
    for(int i = 0; i < nnodes; i++)
    {
        nodeindexes[i] = OldElem->NodeIndex(i);
    }
    
    Mesh->DeleteElement(OldElem);
    
    TPZGeoEl * NewElem = Mesh->CreateGeoBlendElement(oldType, nodeindexes, oldMatId, oldId);

  RestoreNeighbours(NewElem, oldNeigh);
    
	NewElem->BuildBlendConnectivity();
	
	return NewElem;
}
//--------------------------------------------------------
TPZGeoEl * TPZChangeEl::ChangeToArc3D(TPZGeoMesh *mesh, const int64_t ElemIndex,
                                      const TPZVec<REAL> &xcenter, const REAL radius)
{

    auto CreateMidNode = [](const TPZVec<REAL> &x1, const TPZVec<REAL> &x2,
                            const REAL r,
                            const TPZVec<REAL> &xcenter){
        TPZVec<REAL> x3(3,0);
        const auto &xc = xcenter[0];
        const auto &yc = xcenter[1];
        const auto &zc = xcenter[2];
        
        //first we get its distance from xc
        x3[0] = (x1[0] + x2[0])/2 - xc;
        x3[1] = (x1[1] + x2[1])/2 - yc;
        x3[2] = (x1[2] + x2[2])/2 - zc;

        //norm of the vector
        const auto vecnorm = sqrt(x3[0]*x3[0] + x3[1]*x3[1] + x3[2]*x3[2]);

        //mid-arc coordinates
        x3[0] = xc + r * x3[0]/vecnorm;
        x3[1] = yc + r * x3[1]/vecnorm;
        x3[2] = zc + r * x3[2]/vecnorm;
        return x3;
    };
    
    TPZGeoEl * old_el = mesh->ElementVec()[ElemIndex];
    if(!old_el)
    {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - NULL geometric element.\n";
        return NULL;
    }
    const MElementType oldType = old_el->Type();
    if(oldType != EOned){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " geometric el is not 1d\n";
        return NULL;
    }
    const int64_t oldId = old_el->Id();
    const int64_t oldMatId = old_el->MaterialId();
    constexpr int nsides = 3;
    
    TPZVec<TPZGeoElSide> oldNeigh(nsides);
    StoreNeighbours(old_el, oldNeigh);

    //create new node
    TPZManVector<REAL,3> x1(3,0), x2(3,0), x3(3,0);
        
    old_el->Node(0).GetCoordinates(x1);
    old_el->Node(1).GetCoordinates(x2);

    x3 = CreateMidNode(x1, x2, radius, xcenter);
    const auto nodeidx = mesh->NodeVec().AllocateNewElement();
    mesh->NodeVec()[nodeidx].Initialize(x3,*mesh);

    TPZManVector<int64_t,3> nodeindexes =
        {old_el->NodeIndex(0), old_el->NodeIndex(1), nodeidx};
    
    mesh->DeleteElement(old_el);
    auto new_el =
        new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodeindexes, oldMatId, *mesh);

    RestoreNeighbours(new_el, oldNeigh);
    return new_el;
}

TPZGeoEl * TPZChangeEl::ChangeToCylinder(TPZGeoMesh *mesh, const int64_t ElemIndex,
                                         const TPZVec<REAL> &xcenter,
                                         const TPZFMatrix<REAL> &axis
                                         )
{
    
    auto SetCylData = [xcenter,axis,mesh](auto &cyl){
        cyl.SetOrigin(xcenter);
        cyl.SetCylinderAxis(axis);
        cyl.ComputeCornerCoordinates(*mesh);
    };

    TPZGeoEl * old_el = mesh->ElementVec()[ElemIndex];
    if(!old_el)
    {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - NULL geometric element.\n";
        return nullptr;
    }
    const MElementType oldType = old_el->Type();
    if (oldType == EPoint){
        return old_el;
    }
    
    const int64_t oldId = old_el->Id();
    const int oldMatId = old_el->MaterialId();
    const int nsides = old_el->NSides();
    const int nnodes = old_el->NCornerNodes();
    
    TPZVec<TPZGeoElSide> oldNeigh(nsides);
    StoreNeighbours(old_el, oldNeigh);
    TPZManVector<int64_t,4> nodeindexes(nnodes);
    for(int in = 0; in < nnodes; in++){
        nodeindexes[in] = old_el->NodeIndex(in);
    }
    
    mesh->DeleteElement(old_el);
    TPZGeoEl *new_el{nullptr};

    if (oldType == EOned){
        auto cyl =
        new TPZGeoElRefPattern<pzgeom::TPZCylinderMap<pzgeom::TPZGeoLinear>>(nodeindexes, oldMatId, *mesh);
        SetCylData(cyl->Geom());
        new_el = cyl;
    }
    else if(oldType == ETriangle){
        auto cyl =
            new TPZGeoElRefPattern<pzgeom::TPZCylinderMap<pzgeom::TPZGeoTriangle>>(nodeindexes, oldMatId, *mesh);
        SetCylData(cyl->Geom());
        new_el = cyl;
    }else if (oldType == EQuadrilateral){
        auto cyl =
            new TPZGeoElRefPattern<pzgeom::TPZCylinderMap<pzgeom::TPZGeoQuad>>(nodeindexes, oldMatId, *mesh);
        SetCylData(cyl->Geom());
        new_el = cyl;
    }else if (oldType == ETetraedro){
        auto cyl =
            new TPZGeoElRefPattern<pzgeom::TPZCylinderMap<pzgeom::TPZGeoTetrahedra>>(nodeindexes, oldMatId, *mesh);
        SetCylData(cyl->Geom());
        new_el = cyl;
    }else if (oldType == ECube){
        auto cyl =
            new TPZGeoElRefPattern<pzgeom::TPZCylinderMap<pzgeom::TPZGeoCube>>(nodeindexes, oldMatId, *mesh);
        SetCylData(cyl->Geom());
        new_el = cyl;
    }else if (oldType == EPrisma){
        auto cyl =
            new TPZGeoElRefPattern<pzgeom::TPZCylinderMap<pzgeom::TPZGeoPrism>>(nodeindexes, oldMatId, *mesh);
        SetCylData(cyl->Geom());
        new_el = cyl;
    }
    else if (oldType == EPiramide){
        auto cyl =
        new TPZGeoElRefPattern<pzgeom::TPZCylinderMap<pzgeom::TPZGeoPyramid>>(nodeindexes, oldMatId, *mesh);
        SetCylData(cyl->Geom());
        new_el = cyl;
    }
    else {
        DebugStop();
    }

    RestoreNeighbours(new_el, oldNeigh);
    return new_el;
}


//------------------------------------------------------------------------------------------------------------

bool TPZChangeEl::NearestNode(TPZGeoEl * gel, TPZVec<REAL> &x, int64_t &meshNode, double tol)
{    
	meshNode = -1;
	bool IsNearSomeNode = false;
	
	TPZVec<REAL> nodeCoord(3);
	int nnodes = gel->NNodes();
	
	for(int n = 0; n < nnodes; n++)
	{
		double dist = 0.;
		gel->NodePtr(n)->GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			dist += (x[c] - nodeCoord[c])*(x[c] - nodeCoord[c]);
		}
		dist = sqrt(dist);
		
		if(dist <= tol)
		{
			meshNode = gel->NodeIndex(n);
			IsNearSomeNode = true;
			break;
		}
	}
	
	return IsNearSomeNode;
}
//------------------------------------------------------------------------------------------------------------


int64_t TPZChangeEl::NearestNode(TPZGeoMesh * gmesh, TPZVec<REAL> &x, double tol)
{
	int meshNode = -1;
	
	TPZVec<REAL> nodeCoord(3);
	int64_t nnodes = gmesh->NNodes();
	
	for(int64_t n = 0; n < nnodes; n++)
	{
		double dist = 0.;
		gmesh->NodeVec()[n].GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			dist += (x[c] - nodeCoord[c])*(x[c] - nodeCoord[c]);
		}
		dist = sqrt(dist);
		
		if(dist <= tol)
		{
			meshNode = n;
			break;
		}
	}
	
	if(meshNode == -1)
	{
		std::cout << "Node not found for coordinates ( " << x[0] << " , " << x[1] << " , " << x[2] << " )" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		
		DebugStop();
	}
	
	return meshNode;
}
//------------------------------------------------------------------------------------------------------------

bool TPZChangeEl::CreateMiddleNodeAtEdge(TPZGeoMesh *Mesh, int64_t ElemIndex, int edge, int64_t &middleNodeId)
{
    TPZGeoEl * gel = Mesh->ElementVec()[ElemIndex];
    
    TPZGeoElSide myEdge(gel,edge);
    TPZGeoElSide neighEdge(myEdge.Neighbour());
    if(myEdge.Dimension() != 1)
    {
        return false;
    }
    
    TPZVec<REAL> n0Coord(3), n1Coord(3), middleCoordLocal(1), middleCoord(3);
    middleCoordLocal[0] = 0.;//middle node in edge (1D master element) is on qsi=0
    
    int64_t nearestN;
    myEdge.X(middleCoordLocal,middleCoord);
    if(NearestNode(gel, middleCoord, nearestN, 1.E-8))
    {
        middleNodeId = nearestN;//a malha jah contem um noh na coordenada desejada (middleCoord)
        return true;
    }
    //else...
    while(neighEdge != myEdge)
    {
        neighEdge.X(middleCoordLocal,middleCoord);
        if(NearestNode(neighEdge.Element(), middleCoord, nearestN, 1.E-8))
        {
            middleNodeId = nearestN;//a malha jah contem um noh na coordenada desejada (middleCoord)
            return true;
        }
        neighEdge = neighEdge.Neighbour();
    }

    const auto nodeidx = Mesh->NodeVec().AllocateNewElement();
    Mesh->NodeVec()[nodeidx].Initialize(middleCoord,*Mesh);
    middleNodeId = nodeidx;
    
    // //if not returned true...
    // TPZGeoNode midNode;
    // midNode.SetCoord(middleCoord);
    
    // /** Setting Midnodes Id's */
    // int64_t NewNodeId = Mesh->NNodes();
    // Mesh->SetNodeIdUsed(NewNodeId);
    // midNode.SetNodeId(NewNodeId);
    
    // /** Allocating Memory for MidNodes and Pushing Them */
    // middleNodeId = Mesh->NodeVec().AllocateNewElement();
    // Mesh->NodeVec()[middleNodeId] = midNode;

    return true;
}

void TPZChangeEl::StoreNeighbours(TPZGeoEl* gel, TPZVec<TPZGeoElSide> &neighs)
{
    const int nsides = gel->NSides();
    neighs.Resize(nsides);
    for(int s = 0; s < nsides; s++)
    {   
        TPZGeoElSide mySide(gel, s);
        TPZGeoElSide neigh = mySide.Neighbour();
        if(neigh.Exists() && neigh.Element() != gel){
            neighs[s] = neigh;
        }
    }
}
void TPZChangeEl::RestoreNeighbours(TPZGeoEl* gel, TPZVec<TPZGeoElSide> &neighs)
{
    const int nsides = gel->NSides();

    if(nsides != neighs.size()){
        PZError<<__PRETTY_FUNCTION__
                <<"\n neighbour vector should have size nsides. Aborting...\n";
        DebugStop();
    }

    for(int s = 0; s < nsides; s++)
    {
        TPZGeoElSide mygelside(gel,s);
        TPZGeoElSide &neigh = neighs[s];
        if(neigh.Exists()){
            neigh.SetConnectivity(mygelside);
        }else{
            gel->SetNeighbour(s,mygelside);
        }
    }
}
