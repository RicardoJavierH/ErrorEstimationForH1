/**
 * @file
 * @brief Contains the TPZTetrahedron class which defines the topology of the tetrahedron element. 
 */

#ifndef PZTOPOLOGYTPZTETRAHEDRON_H
#define PZTOPOLOGYTPZTETRAHEDRON_H


#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzeltype.h"
#include "pzaxestools.h"
#include "TPZTopologyUtils.h"
class TPZIntPoints;
class TPZIntTetra3D;
class TPZGraphElT3d;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

/// Groups all classes defining the structure of the master element
namespace pztopology {
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of the tetrahedron element. \ref topology "Topology"
	 * Sides 0 to 3 are vertices, sides 4 to 9 are lines, sides 10 to 13 are triangles 
	 * and side 14 is the tetrahedra.
	 */
	class TPZTetrahedron : public TPZSavable {
	public:

    friend void pztopology::GetPermutation<TPZTetrahedron>(const int permute, TPZVec<int> &permutation);
		/** @brief Topological characteristics */
    static constexpr int64_t NSides = 15;
		static constexpr int64_t NCornerNodes = 4;
    static constexpr int64_t Dimension = 3;
    static constexpr int64_t NFacets = 4;
    static constexpr int64_t NPermutations = 24;
      
		
    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;
                
		/** @brief Default constructor */
        TPZTetrahedron() : TPZRegisterClassId(&TPZTetrahedron::ClassId){
		}
		
		/** @brief Default destructor */
		virtual ~TPZTetrahedron() {
		}
		
		/** @name About sides of the topological element
		 * @{ */
		
		/** @brief Returns the dimension of the side */
		static int SideDimension(int side);
		
		/** @brief Get all sides with lower dimension on side */	
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
		/** @brief Get all sides with lower dimension but equal to DimTarget on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);
		
		/**
		 * @brief Returns all sides whose closure contains side
		 * @param side smaller dimension side
		 * @param high vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high);
		
		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NSideNodes(int side);
		/** @brief Returns the local node number of the node "node" along side "side" */
		static int SideNodeLocId(int side, int node);
		
		/** @brief Number of connects of the element (27) */
		static int NumSides();
		/** @brief Returns the number of connects for a set dimension */
		static int NumSides(int dimension);

		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NContainedSides(int side);
		/** @brief Returns the local connect number of the connect "c" along side "side" */
		static int ContainedSideLocId(int side, int c);


        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(const TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        /**
         * This method calculates the influence (a.k.a. the blend function) of the side side regarding an
         * interior point qsi. It is used by the TPZGeoBlend class.
         * @param side the index of the side
         * @param xi coordinates of the interior point
         * @param blendFactor influence (0 <= blendFactor <= 1)
         * * @param corrFactorDxi derivative of the blendFactor in respect to xi
         */
        template<class T>
        static void BlendFactorForSide(const int &side, const TPZVec<T> &xi, T &blendFactor,
                                      TPZVec<T> &corrFactorDxi);
		/** @} */
		
		/** @name About points at the parametric spaces
		 * @{ */

		/** @brief Returns the barycentric coordinates in the master element space of the original element */
		static void CenterPoint(int side, TPZVec<REAL> &center);
		
		/** @brief Verifies if the parametric point pt is in the element parametric domain */
		static bool IsInParametricDomain(const TPZVec<REAL> &pt, REAL tol = pztopology::gTolerance);

        /** @brief Verifies if the parametric point pt is in the element parametric domain (FAD version)*/
		static bool IsInParametricDomain(const TPZVec<Fad<REAL>> &pt, REAL tol = pztopology::gTolerance){
		    TPZVec<REAL> xi(pt.size());
		    for(int i = 0; i < pt.size(); i++) xi[i]= pt[i].val();
		    return IsInParametricDomain(xi,tol);
		}
        template<typename T,
                typename std::enable_if<std::is_same<T,Fad<REAL>>::value>::type* = nullptr>
        static bool IsInParametricDomain(const TPZVec<T> &pt, REAL tol){
            TPZVec<REAL> qsiReal(pt.size(),-1);
            for(int i = 0; i < qsiReal.size(); i++) qsiReal[i] = pt[i].val();
            return IsInParametricDomain(qsiReal,tol);
        }

        /** @brief Generates a random point in the master domain */
        static void RandomPoint(TPZVec<REAL> &pt);

        /**
         * This method will check if the projection to a certain side (MapToSide method) is regular,
         * i.e., if the interior point in the parametric domain is not too close to the projection's singularity.
         * @param side the index of the side upon which the interior point will be projected upon
         * @param xiInterior coordinates of the interior point
         * @return true if the interior point is far from the singularity
         */
        template<class T>
        static bool CheckProjectionForSingularity(const int &side, const TPZVec<T> &xiInterior);

        template<class T>
        static void MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide);
        
        static void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
		
		/** @} */
		
		/** @name About type of the topological element
		 * @{ */
		
		/** @brief Returns the type of the element as specified in file pzeltype.h */
		static constexpr MElementType Type(){return ETetraedro;}
		
		/** @brief Returns the type of the element side as specified in file pzeltype.h */
		static MElementType Type(int side) ;
		
		/** @} */
		
		/** @name About Transformations
		 * @{ */
		
		/**
		 * @brief Returns the transformation which takes a point from the side sidefrom to the side sideto
		 * @param sidefrom Side where the point resides
		 * @param sideto Side whose closure contains sidefrom
		 * @see the class TPZTransform
		 */
		static TPZTransform<> SideToSideTransform(int sidefrom, int sideto);

		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side Side from which the point will be tranformed (0<=side<=26)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformSideToElement(int side);
		/**
		 * @brief Returns the transformation which projects a point from the interior of the element to the side
		 * @param side Side to which the point will be tranformed (0<=side<=26)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformElementToSide(int side);
		
		/**
		 * @brief Method which identifies the transformation based on the IDs of the corner nodes
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */
		static int GetTransformId(const TPZVec<int64_t> &id);
		
		/**
		 * @brief Method which identifies the transformation of a side based on the IDs of the corner nodes
		 * @param side Index of side
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */	
		static int GetTransformId(const int side, const TPZVec<int64_t> &id);
		
		/** @} */
		
		/** @name Methods related over numeric integration
		 * @{ */
		
		/**
		 * @brief Create an integration rule over side
		 * @param side Side to create integration rule
		 * @param order Order of the integration rule to be created
		 */
		static TPZIntPoints * CreateSideIntegrationRule(int side, int order);
		
		/** @brief Typedef to numerical integration rule */
		typedef TPZIntTetra3D IntruleType;
		/** @brief Typedef to graphical element type */
		typedef TPZGraphElT3d GraphElType;
		
		/** @} */

		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side Side for which the permutation is needed
		 * @param id Ids of the corner nodes of the elements
		 * @param permgather Permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather);
		
		/** @brief Volume of the master element (measure) */
		static constexpr REAL RefElVolume() {
			return (1.0L/6.0L);
		}
        
        /* Given side and gradx the method returns directions needed for Hdiv space */
        static void ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors);
        static void GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao);
        static void GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao, TPZVec<int> &sidevectors);
        
        /// Compute the directions of the HDiv vectors
        template <class TVar>
        static void ComputeHDivDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions);

        /// Compute the directions of the HDiv vectors for constant divergent
        // template <class TVar>
        static void ComputeConstantHDiv(const TPZVec<REAL> &point, TPZFMatrix<REAL> &vecDiv, TPZVec<REAL> &div);
        static void ComputeConstantHDiv(const TPZVec<Fad<REAL>> &point, TPZFMatrix<Fad<REAL>> &vecDiv, TPZVec<Fad<REAL>> &div);
        template<class TVar>
        static void ComputeConstantHCurl(const TPZVec<TVar> &point, TPZFMatrix<TVar> &vecDiv, TPZFMatrix<TVar> &curl, const TPZVec<int> &transformationIds);
        static int GetFaceOrient(const int &face);

        /** Compute the directions of the HCurl vectors.
         * These vectors are combined with H1 shape functions to create the HCurl shape functions.
         * They *must be* computed in the following order:
         * - \f$v^{e,a}\f$: vector associated with edge \f$e\f$. It is normal to the edge \f$\hat{e}\f$ adjacent to \f$e\f$e by the vertex \f$a\f$a.
         * - \f$v^{e,T}\f$: vector associated with edge \f$e\f$. It is tangent to the edge \f$\hat{e}\f$.
         * - \f$v^{F,e}\f$: vector associated with face \f$F\f$. It is normal to the face \f$\hat{F}\f$ adjacent to \f$F\f$e by the edge \f$e\f$a.
         * - \f$v^{F,T}\f$: two orthornormal vectors associated with face \f$F\f$ and tangent to it.
         * - \f$v^{F,\perp}\f$: outward normal vector associated with face \f$F\f$ (3D only)
         * - \f$v^{K}\f$: set of orthonormal vectors associated with the volume of the element itself (3D only. In 2D \f$v^{F,T}\f$ does its job)
         * The side ordering should be respected. In the definition of the \f$v^{e,a}\f$ and the \f$v^{F,e}\f$ vectors, the subsides are ordered as the return of LowerDimensionSides.
         * @tparam TVar REAL or Fad<REAL>
         * @param gradx the gradient of the element mapping. if computing in normal element, gradx is the identity matrix.
         * @param directions computed directions
         * @param transformationIds transformation Ids associated with each side of dim > 0
         */
        template <class TVar>
        static void ComputeHCurlDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions, const TPZVec<int> &transformationIds);

        /**
         * Returns the number of bilinear sides to this shape. Needed to compute the number shapefunctions( NConnectShapeF )
         */
        static int NBilinearSides();

	protected:
		/** @name Data structure which defines the tetrahedral transformations */
		/** @{ */
		
		/** @brief Nodes over quadrilateral sides (2d - faces). */
    static constexpr int FaceNodes[4][3]  = { {0,1,2},{0,1,3},{1,2,3},{0,2,3} };

		/** @brief Nodes over lines sides (1d) */
    static constexpr int SideNodes[6][2]  = { {0,1},{1,2},{2,0},{0,3},{1,3},{2,3} };

		/** @brief Ids of the shape face */
    static constexpr int ShapeFaceId[4][3] = { {0,1,2},{0,1,3},{1,2,3},{0,2,3} };	

    /** @brief Valid permutations between nodes*/
    static constexpr int fPermutations[24][15] = {
      {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14},/*000*/
      {0,1,3,2,4,8,7,6,5,9,11,10,12,13,14},/*001*/
      {0,2,1,3,6,5,4,7,9,8,10,13,12,11,14},/*002*/
      {0,2,3,1,6,9,7,4,5,8,13,10,12,11,14},/*003*/
      {0,3,1,2,7,8,4,6,9,5,11,13,12,10,14},/*004*/
      {0,3,2,1,7,9,6,4,8,5,13,11,12,10,14},/*005*/
      {1,0,2,3,4,6,5,8,7,9,10,11,13,12,14},/*006*/
      {1,0,3,2,4,7,8,5,6,9,11,10,13,12,14},/*007*/
      {1,2,0,3,5,6,4,8,9,7,10,12,13,11,14},/*008*/
      {1,2,3,0,5,9,8,4,6,7,12,10,13,11,14},/*009*/
      {1,3,0,2,8,7,4,5,9,6,11,12,13,10,14},/*010*/
      {1,3,2,0,8,9,5,4,7,6,12,11,13,10,14},/*011*/
      {2,0,1,3,6,4,5,9,7,8,10,13,11,12,14},/*012*/
      {2,0,3,1,6,7,9,5,4,8,13,10,11,12,14},/*013*/
      {2,1,0,3,5,4,6,9,8,7,10,12,11,13,14},/*014*/
      {2,1,3,0,5,8,9,6,4,7,12,10,11,13,14},/*015*/
      {2,3,0,1,9,7,6,5,8,4,13,12,11,10,14},/*016*/
      {2,3,1,0,9,8,5,6,7,4,12,13,11,10,14},/*017*/
      {3,0,1,2,7,4,8,9,6,5,11,13,10,12,14},/*018*/
      {3,0,2,1,7,6,9,8,4,5,13,11,10,12,14},/*019*/
      {3,1,0,2,8,4,7,9,5,6,11,12,10,13,14},/*020*/
      {3,1,2,0,8,5,9,7,4,6,12,11,10,13,14},/*021*/
      {3,2,0,1,9,6,7,8,5,4,13,12,10,11,14},/*022*/
      {3,2,1,0,9,5,8,7,6,4,12,13,10,11,14}  /*023*/
    };
    /** @} */
  };
	
}

#endif
