//
//  TestTopology.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/6/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "pzmanvector.h"
#include "pztrnsform.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"
#include "tpzpyramid.h"
#include "pzgeotetrahedra.h"
#include "pzgeotriangle.h"
#include "TPZShapeHCurl.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeData.h"
#include "TPZMaterialData.h"
#include "pzelctemp.h"
#include "TPZShapeHDivKernel2D.h"

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapetetra.h"



#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

using namespace pztopology;
using namespace pzshape;
namespace topologytests{
    template <class top>
    void TestingSideProjections();
    template <class top>
    void TestingSideNodeProjections();
    template <class top,class TSHAPE>
    void TestingConstantDivergent();
    template <class top,class TSHAPE>
    void TestingConstantCurl();
    template <class top>
    void TestSideOrient();
}

TEST_CASE("verify faceorientation data structure")
{
    topologytests::TestSideOrient<TPZCube>();
    topologytests::TestSideOrient<TPZTetrahedron>();
    topologytests::TestSideOrient<TPZPrism>();
    topologytests::TestSideOrient<TPZQuadrilateral>();
    topologytests::TestSideOrient<TPZTriangle>();
}
TEST_CASE("projection_tests_1","[topology_tests]")
{
    topologytests::TestingSideProjections<TPZLine>();
    topologytests::TestingSideProjections<TPZTriangle>();
    topologytests::TestingSideProjections<TPZQuadrilateral>();
    topologytests::TestingSideProjections<TPZTetrahedron>();
    topologytests::TestingSideProjections<TPZCube>();
    topologytests::TestingSideProjections<TPZPrism>();
    topologytests::TestingSideProjections<TPZPyramid>();
}

TEST_CASE("projection_tests_2","[topology_tests]")
{
    topologytests::TestingSideNodeProjections<TPZTriangle>();
    topologytests::TestingSideNodeProjections<TPZQuadrilateral>();
    topologytests::TestingSideNodeProjections<TPZTetrahedron>();
    topologytests::TestingSideNodeProjections<TPZCube>();
    topologytests::TestingSideNodeProjections<TPZPrism>();
    topologytests::TestingSideNodeProjections<TPZPyramid>();
}

TEST_CASE("constant_divergent_test","[topology_tests]")
{
    topologytests::TestingConstantDivergent<TPZTriangle,TPZShapeTriang>();
    topologytests::TestingConstantDivergent<TPZQuadrilateral,TPZShapeQuad>();
    topologytests::TestingConstantDivergent<TPZTetrahedron,TPZShapeTetra>();
    topologytests::TestingConstantDivergent<TPZCube,TPZShapeCube>();
    topologytests::TestingConstantDivergent<TPZPrism,TPZShapePrism>();
}


TEST_CASE("constant_curl_test","[topology_tests]")
{
    topologytests::TestingConstantCurl<TPZTriangle,TPZShapeTriang>();
    topologytests::TestingConstantCurl<TPZQuadrilateral,TPZShapeQuad>();
    topologytests::TestingConstantCurl<TPZTetrahedron,TPZShapeTetra>();
    topologytests::TestingConstantCurl<TPZCube,TPZShapeCube>();
    topologytests::TestingConstantCurl<TPZPrism,TPZShapePrism>();
}

namespace topologytests{
    template <class top>
    void TestingSideProjections() {
        static std::string testName = __PRETTY_FUNCTION__;
        auto type = top::Type();
        int dim = top::Dimension;
        int nsides = top::NSides;
        for (auto iSide=0; iSide<nsides; iSide++) {
            TPZTransform<> tr1 = top::TransformSideToElement(iSide);
            TPZTransform<> tr2 = top::TransformElementToSide(iSide);
            TPZTransform<> tr3;
            tr3 = tr2.Multiply(tr1);
            TPZMatrix<REAL>& mat = tr3.Mult();
            const REAL tol = GetTolerance();
            const int nRows = mat.Rows();
            const int nCols = mat.Cols();
            for(int i = 0; i < nRows; i++){
                for(int j = 0; j < nCols; j++){
                    bool cond;
                    if (j == i) {
                        cond = mat(i,j)-1 < tol;
                    }else{
                        cond = std::abs(mat(i,j) ) < tol;
                    }
                    if(!cond){
                      std::cerr
                          << "\n" + testName + " failed" +
                                 "\ntopology: " + MElementType_Name(type) +
                                 "\n" + "side: " + std::to_string(iSide);
                    }
                    REQUIRE(cond);
                }
            }
        }
    }//TestingSideProjections

    template <class top>
    void TestingSideNodeProjections() {
        static std::string testName = __PRETTY_FUNCTION__;
        auto type = top::Type();
        int dim = top::Dimension;
        int nSides = top::NSides;
        int nNodes = top::NCornerNodes;
        TPZFNMatrix<24,REAL> cornerNodes(nNodes,dim);
        TPZManVector<REAL,3> node(dim);
        for (auto iNode=0; iNode<nNodes; iNode++) {
            top::ParametricDomainNodeCoord(iNode, node);
            for(auto x = 0; x < dim; x++) cornerNodes(iNode,x) = node[x];
        }
        const REAL tol = GetTolerance();
        for (auto iSide=nNodes; iSide<nSides - 1; iSide++) {
            const auto sideDim = top::SideDimension(iSide);
            TPZManVector<REAL,3> paramNode(sideDim);
            TPZTransform<> tr = top::TransformSideToElement(iSide);
            const auto nSideNodes = top::NSideNodes(iSide);
            for(auto iNode = 0; iNode < nSideNodes; iNode++){
                switch(nSideNodes){
                    case 2://edge
                        TPZLine::ParametricDomainNodeCoord(iNode,paramNode);
                        break;
                    case 3://triangle
                        TPZTriangle::ParametricDomainNodeCoord(iNode,paramNode);
                        break;
                    case 4://quadrilateral
                        TPZQuadrilateral::ParametricDomainNodeCoord(iNode,paramNode);
                        break;
                    default:
                        DebugStop();
                }
                tr.Apply(paramNode,node);
                const auto nodeId = top::SideNodeLocId(iSide,iNode);
                REAL diff = 0;
                for(auto x = 0; x < sideDim; x++) diff += (node[x] - cornerNodes(nodeId,x)) * (node[x] - cornerNodes(nodeId,x));
                bool cond = std::sqrt(diff) < tol;
                if(!cond){
                  std::cerr << "\n" + testName + " failed" +
                                   "\ntopology: " + MElementType_Name(type) +
                                   "\n" + "side: " + std::to_string(iSide);
                }
                REQUIRE(cond);
            }
        }
    }//TestingSideNodeProjections

    template <class top,class TSHAPE>
    void TestingConstantDivergent() {
        
        static std::string testName = __PRETTY_FUNCTION__;
        const  int64_t pOrderIntRule = 5;
        const auto nSides = top::NSides;
        const auto nCorner = top::NCornerNodes;
        const auto dim = top::Dimension;
        const auto nFaces = top::NFacets;
        auto type = top::Type();
        const REAL tol = 1.e-6;
        auto gel = pzgeom::TPZNodeRep<nCorner,top>();

        TPZIntPoints* intRule = gel.CreateSideIntegrationRule(nSides-1, pOrderIntRule);     
        int nEdges = 0;
        if (dim == 2){
            nEdges = nFaces;
        } else if (dim == 3){
            nEdges = nSides - 1 - nCorner - nFaces;
        }
        
        TPZVec<REAL> node(dim);
                
        TPZShapeData shapedata;
        TPZManVector<int64_t,nCorner> ids(nCorner,0);
        for(auto i=0; i<nCorner; i++) ids[i] = i;        
        auto &conOrders = shapedata.fHDivConnectOrders;
        constexpr auto nConnects = nFaces + 1;
        conOrders.Resize(nConnects,-1);
        for(auto i = 0; i < nConnects; i++) conOrders[i] = 1;
        TPZManVector<int,TSHAPE::NFacets> sideorient(TSHAPE::NFacets,0);
        for(int i=0; i<TSHAPE::NFacets; i++) sideorient[i] = 1;

        TPZShapeHDiv<TSHAPE>::Initialize(ids, conOrders, sideorient, shapedata);

        shapedata.fSideTransformationId.Resize(nSides-nCorner, 0);
        for (int iside = nCorner; iside< nSides ; iside++) {
            int pos = iside - nCorner;
            int trans_id = TSHAPE::GetTransformId(iside, ids); // Foi criado
            shapedata.fSideTransformationId[iside-nCorner] = trans_id;
        }

        const int npts = intRule->NPoints();
        auto nshape = shapedata.fSDVecShapeIndex.size();

        for (auto ipt = 0; ipt < npts; ipt++) {
            REAL w;
            intRule->Point(ipt, node, w);

            TPZFMatrix<REAL> phiHDiv(3,nshape);
            TPZFMatrix<REAL> divHDiv(3,nshape);
            TPZFMatrix<REAL> RT0Function(dim,nFaces);
            TPZVec<REAL> divRT0(nFaces);
            RT0Function.Zero();
            divRT0.Fill(0.);
            phiHDiv.Zero();
            divHDiv.Zero();

            TPZShapeHDiv<TSHAPE>::Shape(node,shapedata,phiHDiv,divHDiv);

            // std::cout << "phiHDiv = " << phiHDiv << std::endl;
            // std::cout << "divHDiv =  " << divHDiv << std::endl;
            
            //Compute the curl for each edge
            top::ComputeConstantHDiv(node,RT0Function,divRT0);

            // std::cout << "Constant phi = " << RT0Function << std::endl;
            // std::cout << "Constant div = " << divRT0 << std::endl;

            int nEdgesPerFace = 0;
            if ((top::Type() == ETriangle) || (top::Type() == EQuadrilateral)){
                nEdgesPerFace = 2;
            } else if (top::Type() == ETetraedro) {
                nEdgesPerFace = 3;
            } else if (top::Type() == ECube) {
                nEdgesPerFace = 4;
            }

            // Checks if all edges have the same curl value.
            bool condHdiv = true;
            bool condHdivRT0 = true;
             
            for (int i = 0; i < dim; i++){
                int aux = 0;
                for (int j = 0; j < nFaces; j++){
                    if (top::Type() == EPrisma){
                        if ( j==0 || j == nFaces-1){
                            nEdgesPerFace = 3;
                        } else {
                            nEdgesPerFace = 4;
                        }
                    }
                    REAL funVal = 0.;
                    for (int k = 0; k < (nEdgesPerFace); k++)
                    {
                        funVal += phiHDiv(i,aux+k) / (nEdgesPerFace);
                    }
                    if (fabs(funVal-RT0Function(i,j)) > tol){
                        condHdiv = false;
                    } 
                    aux += nEdgesPerFace;
                }
                
            }
            
            int aux = 0;
            for (int j = 0; j < nFaces; j++){
                if (top::Type() == EPrisma){
                    if ( j==0 || j == nFaces-1){
                        nEdgesPerFace = 3;
                    } else {
                        nEdgesPerFace = 4;
                    }
                }
                REAL divVal = 0.;
                for (int k = 0; k < (nEdgesPerFace); k++)
                {
                    divVal += divHDiv(aux+k) / (nEdgesPerFace);
                }
                if (fabs(divVal-divRT0[j]) > tol){
                    condHdivRT0 = false;
                } 
                aux += nEdgesPerFace;
            }
        
            if(!condHdiv || !condHdivRT0){
                std::cerr << "\n" + testName + " failed Hdiv" +
                                "\ntopology: " + MElementType_Name(type);
            }

            REQUIRE(condHdiv);
            REQUIRE(condHdivRT0);
        }
    }//Testing Constant Divergent


    template <class top,class TSHAPE>
    void TestingConstantCurl() {
        
        static std::string testName = __PRETTY_FUNCTION__;
        const  int64_t pOrderIntRule = 5;
        const auto nSides = top::NSides;
        const auto nCorner = top::NCornerNodes;
        const auto dim = top::Dimension;
        const auto nFaces = top::NFacets;
        auto type = top::Type();
        const REAL tol = 1.e-6;
        auto gel = pzgeom::TPZNodeRep<nCorner,top>();

        TPZIntPoints* intRule = gel.CreateSideIntegrationRule(nSides-1, pOrderIntRule);     
        int nEdges = 0;
        if (dim == 2){
            nEdges = nFaces;
        } else if (dim == 3){
            nEdges = nSides - 1 - nCorner - nFaces;
        }
        
        TPZVec<REAL> node(dim);
                
        TPZShapeData shapedata;
        TPZManVector<int64_t,nCorner> ids(nCorner,0);
        for(auto i=0; i<nCorner; i++) ids[i] = i;        
        auto &conOrders = shapedata.fHDivConnectOrders;
        constexpr auto nConnects = nSides - nCorner;
        conOrders.Resize(nConnects,-1);
        for(auto i = 0; i < nConnects; i++) conOrders[i] = 0;

        TPZShapeHCurl<TSHAPE>::Initialize(ids, conOrders, shapedata);

        const int npts = intRule->NPoints();
        auto nshape = shapedata.fSDVecShapeIndex.size();

        for (auto ipt = 0; ipt < npts; ipt++) {
            REAL w;
            intRule->Point(ipt, node, w);

            TPZFMatrix<REAL> phiHCurl(3,nshape);
            TPZFMatrix<REAL> curlHCurl(3,nshape);
            TPZFMatrix<REAL> N0Function(dim,nEdges);
            TPZFMatrix<REAL> curlN0(3,nEdges);
            N0Function.Zero();
            curlN0.Zero();
            phiHCurl.Zero();
            curlHCurl.Zero();

            TPZShapeHCurl<TSHAPE>::Shape(node,shapedata,phiHCurl,curlHCurl);

            // std::cout << "phiHCurl = " << phiHCurl << std::endl;
            // std::cout << "curlHCurl =  " << curlHCurl << std::endl;
            
            //Compute the curl for each edge
            top::ComputeConstantHCurl(node,N0Function,curlN0,shapedata.fSideTransformationId);

            // std::cout << "Constant phi = " << N0Function << std::endl;
            // std::cout << "Constant Curl = " << curlN0 << std::endl;

            for (int j = 0; j < nEdges; j++){
                TPZManVector<REAL,dim> phi_shape(dim,0.), phi_topol(dim,0.);
                for (int i = 0; i < dim; i++){
                    phi_shape[i] = phiHCurl(i,j);
                    phi_topol[i] = N0Function(i,j);
                }
                CAPTURE(phi_shape);
                CAPTURE(phi_topol);
                for (int i = 0; i < dim; i++){
                    REQUIRE(phiHCurl(i,j)-N0Function(i,j) == Catch::Approx(0).margin(tol));
                }
            }
            for (int j = 0; j < nEdges; j++){
                TPZManVector<REAL,3> curl_shape(dim,0.), curl_topol(dim,0.);
                for (int i = 0; i < dim; i++){
                    curl_shape[i] = curlHCurl(i,j);
                    curl_topol[i] = curlN0(i,j);
                }
                CAPTURE(curl_shape);
                CAPTURE(curl_topol);
                for (int i = 0; i < 3; i++){
                    REQUIRE(curlHCurl(i,j)-curlN0(i,j) == Catch::Approx(0).margin(tol));
                }
            }
        }
    }//Testing Constant Curl

#include "pzvec_extras.h"

static void Orthogonal(TPZFMatrix<REAL>&cornerco, TPZVec<REAL> &ortho){
    int nnodes = cornerco.Cols();
    if(nnodes == 2) {
        ortho.resize(2);
        ortho[1] = -cornerco(0,1)+cornerco(0,0);
        ortho[0] = -cornerco(1,0)+cornerco(1,1);
    } else if(nnodes == 3 || nnodes == 4) {
        TPZManVector<REAL,3> vec1(3),vec2(3);
        for (int i=0; i<3; i++) {
            vec1[i] = cornerco(i,1)-cornerco(i,0);
            vec2[i] = cornerco(i,nnodes-1)-cornerco(i,0);
        }
        ortho.resize(3);
        Cross(vec1, vec2, ortho);
    }
}
                                                                           
    template <class top>
    void TestSideOrient() {
        const int dim = top::Dimension;
        REAL cornerco[top::NCornerNodes][top::Dimension];
        for (int i=0; i<top::NCornerNodes; i++) {
            TPZManVector<REAL,3> co(dim,0.);
            top::CenterPoint(i,co);
            for(int c=0; c<top::Dimension; c++) cornerco[i][c] = co[c];
        }
        TPZManVector<REAL,3> centerelement(dim);
        top::CenterPoint(top::NSides-1,centerelement);
        int firstface = top::NSides-1-top::NumSides(dim-1);
        for(int side = firstface; side < top::NSides-1; side++) {
            TPZManVector<REAL,3> centerface(dim), dir(dim);
            top::CenterPoint(side,centerface);
            dir = centerface-centerelement;
            int nnodes = top::NSideNodes(side);
            TPZFNMatrix<12> faceco(dim,nnodes,0.);
            for (int in=0; in<nnodes; in++) {
                int locnod = top::SideNodeLocId(side,in);
                for(int c = 0; c<dim; c++) {
                    faceco(c,in) = cornerco[locnod][c];
                }
            }
            TPZManVector<REAL,3> normal(dim,0.);
            Orthogonal(faceco, normal);
            REAL inner = sdot(normal,dir);
            int faceorient = inner > 0 ? 1 : -1;
            int faceorientStored = top::GetFaceOrient(side-firstface);
            REQUIRE(faceorient == faceorientStored);
        }
    }


}//namespace
