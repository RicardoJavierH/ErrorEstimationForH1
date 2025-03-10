//
// Created by victor on 14/02/23.
//

#include "TPZHybridH1HdivFluxRecMaterial.h"

//
// Created by victor on 07/10/2020.
//

#include "TPZHybridH1ErrorEstimateMaterial.h"
#include "pzaxestools.h"
#include "TPZAnalyticSolution.h"
#include "TPZMaterialDataT.h"


#ifdef LOG4CXX
static LoggerPtr loggerF(Logger::getLogger("DebuggingF"));
#endif

TPZHybridH1HdivFluxRecMaterial::TPZHybridH1HdivFluxRecMaterial(int matid, int dim) : TPZMixedPoisson(matid, dim)
{
}

TPZHybridH1HdivFluxRecMaterial::TPZHybridH1HdivFluxRecMaterial() : TPZMixedPoisson()
{

}

TPZHybridH1HdivFluxRecMaterial::TPZHybridH1HdivFluxRecMaterial(const TPZHybridH1HdivFluxRecMaterial &copy) : TPZMixedPoisson(copy)
{
}

TPZHybridH1HdivFluxRecMaterial::TPZHybridH1HdivFluxRecMaterial(const TPZMixedPoisson &copy) : TPZMixedPoisson(copy)
{

}

TPZHybridH1HdivFluxRecMaterial::TPZHybridH1HdivFluxRecMaterial(TPZMatLaplacianHybrid &matlaplacian)
{
    this->SetId(matlaplacian.Id());
    this->SetDimension(matlaplacian.Dimension());

    STATE perm = matlaplacian.GetPermeability({0.,0.,0.});
    SetConstantPermeability(perm);

    if (matlaplacian.HasForcingFunction()) {
        TPZMatErrorCombinedSpaces<STATE> *me = this;
        TPZMatErrorCombinedSpaces<STATE> *lapl = &matlaplacian;
        auto funcpt = lapl->ExactSol();
        me->SetExactSol(funcpt,lapl->PolynomialOrderExact());
        this->SetForcingFunction(matlaplacian.ForcingFunction(),matlaplacian.ForcingFunctionPOrder());
    }
}

TPZHybridH1HdivFluxRecMaterial::~TPZHybridH1HdivFluxRecMaterial()
{

}

TPZHybridH1HdivFluxRecMaterial &TPZHybridH1HdivFluxRecMaterial::operator=(const TPZHybridH1HdivFluxRecMaterial &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}

void TPZHybridH1HdivFluxRecMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
     find (sigma_h,p_h,g_h) satisfying
        (invK.sigma_h,v_h) - (p_h,div v_h) = 0;
        -(div sigma_h,w_h) + (g_h,w_h) = (f,v);
        (p_h,t_h) = (u_avg,t_h),
    for all (v_h,w_h,t_h)
     **/

    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);
    int dim = datavec[0].axes.Rows();

    auto perm = GetPermeability(datavec[0].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);

    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    //TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &phip = datavec[1].phi;

    STATE force = 0.;
    if(fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction(datavec[1].x,res);
        force = res[0];
    }


    TPZFMatrix<REAL> &phiuk = datavec[0].phi;

    int phrq,phrp;
    phrq = datavec[0].fVecShapeIndex.NElements();
    phrp = phip.Rows();

    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
        }

        TPZFNMatrix<3,REAL> ivecZ(3,1,0.);
        TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
        for (int jq=0; jq<phrq; jq++)
        {

            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;

            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[0].fDeformedDirections(id,jvecind);
            }

            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<3; id++){
                for(int jd=0; jd<3; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }
            /// (invK.sigma_h,v)_K
            REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
            ek(iq,jq) += 1.*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
        }
    }

    // Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;

        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
            //ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
            //ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
        }
        TPZFNMatrix<3,REAL> axesvec(3,1,0.);
        datavec[0].axes.Multiply(ivec,axesvec);

        REAL divwq = 0.;
        divwq = datavec[0].divphi(ivecind);

        for (int jp=0; jp<phrp; jp++) {
            /// - (p_h,div v_h)
            REAL fact = (-1.)*weight*phip(jp,0)*divwq;
            // Matrix B
            ek(iq, phrq+jp) += fact;

            // Matrix B^T
            ek(phrq+jp,iq) += fact;

        }
    }
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        /// (f,v)
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }

    for(int ip=0; ip<phrp; ip++){
        /// (g_h,w_h)
        ek(phrq+ip,phrq+phrp) += phip(ip,0)*weight;
        ek(phrq+phrp,phrq+ip) += phip(ip,0)*weight;
    }
    ef(phrq+phrp,0) += weight*datavec[3].sol[0][0];
}

void TPZHybridH1HdivFluxRecMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const {

        datavec[3].SetAllRequirements(false);
        datavec[3].fNeedsSol = true;
        datavec[3].fNeedsNormal = true;

        datavec[4].SetAllRequirements(false);
        datavec[4].fNeedsSol = true;
        //datavec[4].fNeedsNormal = true;

}

void TPZHybridH1HdivFluxRecMaterial::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const{

        datavec[4].SetAllRequirements(false);
}

void TPZHybridH1HdivFluxRecMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    /**
      error[0] - error computed with exact pressure
      error[1] - error computed with reconstructed pressure
      error[2] - energy error computed with exact solution
     **/
    if(!ExactSol()) DebugStop();
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    auto x = data[fHDivConformPosition].x;

    /// Permeability data
    auto perm = GetPermeability(x);
    auto invPerm = 1./perm;

    /// Potential, velocity and source term from analytical data
    TPZVec<STATE> u(1);
    TPZFNMatrix<3,REAL> gradU(3,1),KGradU(3,1);
    ExactSol()(x, u, gradU);

    for (int i=0; i<3; i++) 
        KGradU(i,0) = perm*gradU(i,0);

    TPZVec<STATE> source(1);

    if(this->fForcingFunction){

        this->fForcingFunction(x,source);
    } 

    /// Broken H1 data from FEM simulation
    STATE u_h;
    TPZFMatrix<REAL> gradUhAxes(3,1),gradUh(3,1),KGradUh(3,1);

    u_h = data[fFEMbrokenH1Position].sol[0][0];

    gradUhAxes = data[fFEMbrokenH1Position].dsol[0];
    TPZAxesTools<REAL>::Axes2XYZ(gradUhAxes,gradUh,data[fFEMbrokenH1Position].axes);

    for (int i=0; i<3; i++) 
        KGradUh(i,0) = perm*gradUh(i,0);

    /// HDiv-conform data from reconstruction
    TPZFNMatrix<3,REAL> t_h(3,1);
    for(int ip = 0 ; ip < 3 ; ip++){
        t_h(ip,0) = data[0].sol[0][ip];
    }

    STATE divTh;
    divTh= data[fHDivConformPosition].divsol[0][0];

    REAL FEMsemiH1Error = 0.;
    REAL thFEMerrorEstimate = 0.;
    REAL residualFromHDivSpace = (source[0] - divTh)*(source[0] - divTh);

    for (int i=0; i<3; i++) {
            FEMsemiH1Error     += (KGradUh[i]-KGradU(i,0))*invPerm*(KGradUh[i]-KGradU(i,0));//Pq esta somando: o fluxo fem esta + e o exato -
            thFEMerrorEstimate += (KGradUh[i]+t_h[i])*invPerm*(KGradUh[i]+t_h[i]);
    }

    errors[0] = FEMsemiH1Error;          // ||grad(u_h)-grad(u)||
    errors[1] = thFEMerrorEstimate;      // ||grad(u_h)+sigma_h)||
    errors[2] = residualFromHDivSpace;   // ||f - Proj_divsigma||
}



int TPZHybridH1HdivFluxRecMaterial::VariableIndex(const std::string &name) const
{
    if(name == "minus_KGradUh") return 40;
    if(name == "minus_KGradU") return 42;
    if(name == "uh") return 43;
    if(name == "th") return 39;
    if(name == "POrder") return 46;

    if(name == "KGradUh_minus_KGradU") return 100;
    if(name == "th_plus_KGradUh") return 101;
    if(name == "residual") return 102;


    return -1;
}


int TPZHybridH1HdivFluxRecMaterial::NSolutionVariables(int var) const
{
    switch (var) {
    case 39:
    case 40:
    case 42:
        return 3;
        break;
    case 43:
    case 46:
    case 100:
    case 101:
    case 102:
        return 1;
        break;
    default:
        DebugStop();
        break;
    }
    return 0;
}

/**
 * @brief It return a solution to multiphysics simulation.
 * @param datavec [in] Data material vector
 * @param var [in] number of solution variables. See  NSolutionVariables() method
 * @param Solout [out] is the solution vector
 */

void TPZHybridH1HdivFluxRecMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout)
{

    /**
     datavec[0] HDiv reconstructed mesh
     datavec[1] L2 mesh,
     datavec[2] gMesh
     datavec[3] avg u_h from simulation
     datavec[4] H1 mesh from FEM
     **/

    int H1functionposition = 1;

    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);

    auto perm = GetPermeability(datavec[0].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);

    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);

    if(fExactSol)
    {
        this->fExactSol(datavec[H1functionposition].x, pressexact,gradu);

    }

    PermTensor.Multiply(gradu, fluxinv);

    int dim=this->fDim;
    switch (var)
    {
    case 40://FluxFem
    {
        TPZFMatrix<REAL> &dsolaxes = datavec[4].dsol[0];
        TPZFNMatrix<9, REAL> dsol(3, 0);
        TPZFNMatrix<9, REAL> KGradsol(3, 0);
        TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[4].axes);

        PermTensor.Multiply(dsol, KGradsol);

        for (int i = 0; i < 3; i++) Solout[i]  = -KGradsol(i,0);
    }
    break;
    //Flux reconstrucion
    case 39: // sigma_h
        for(int id=0 ; id<3; id++) {
            Solout[id] = datavec[0].sol[0][id];
        }
    break;
    case 42://flux exact
        for(int i=0; i<dim; i++) Solout[i] = -fluxinv(i);
        break;

    case 43://Pressure fem
        Solout[0] = datavec[fFEMbrokenH1Position].sol[0][0];
        break;
    case 46://order p
        Solout[0] = datavec[1].p;
        break;

    default:
        //TPZMixedPoisson::Solution(datavec,var, Solout);
        DebugStop();
    }
}

