//
//  TPZMixedHDivErrorEstimate.cpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 20/04/18.
//

#include "TPZMixedHdivErrorEstimate.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "pzaxestools.h"
#include "TPZMaterialDataT.h"


template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat>::TPZMixedHDivErrorEstimate() : MixedMat()
{
    
}

template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat>::TPZMixedHDivErrorEstimate(int matid, int dim) : MixedMat(matid,dim)
{
    
}

template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat>::TPZMixedHDivErrorEstimate(const MixedMat &copy) : MixedMat(copy)
{
    
}

template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat>::~TPZMixedHDivErrorEstimate()
{
    
}

template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat>::TPZMixedHDivErrorEstimate(const TPZMixedHDivErrorEstimate &cp) : MixedMat(cp)
{
    
}

template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat> &TPZMixedHDivErrorEstimate<MixedMat>::operator=(const TPZMixedHDivErrorEstimate &copy)
{
    MixedMat *me = this;
    me->operator=(copy);
    return *this;
}

template<class MixedMat>
void TPZMixedHDivErrorEstimate<MixedMat>::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const
{
    MixedMat::FillDataRequirements(datavec);
    {
        int i = 1;
        datavec[i].SetAllRequirements(false);
        //datavec[0].fNeedsSol = true;
        datavec[i].fNeedsSol = true;
    }
    
}

template<class MixedMat>
int TPZMixedHDivErrorEstimate<MixedMat>::VariableIndex(const std::string &name) const
{
    if(name == "FluxFem") return 40;
    if(name == "FluxReconstructed") return 41;
    if(name == "FluxExact") return 42;
    if(name == "PressureFEM") return 43;
    if(name == "PressureReconstructed") return 44;
    if(name == "PressureExact") return 45;
    if(name == "PressureErrorExact") return 100;
    if(name == "PressureErrorEstimate") return 101;
    if(name == "EnergyErrorExact") return 102;
    if(name == "EnergyErrorEstimate") return 103;
    if(name == "ResidualError") return 104;
    if(name == "PressureEffectivityIndex") return 105;
    if(name == "EnergyEffectivityIndex") return 106;
    if(name == "POrder") return 46;
   
    return -1;
}

template<class MixedMat>
int TPZMixedHDivErrorEstimate<MixedMat>::NSolutionVariables(int var) const
{
    switch (var) {
        case 40:
        case 41:
        case 42:
            return 3;
            break;
        case 43:
        case 44:
        case 45:
        case 46:
        case 100:
        case 101:
        case 102:
        case 103:
        case 104:
        case 105:
        case 106:
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
template<class MixedMat>
void TPZMixedHDivErrorEstimate<MixedMat>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout)
{
    /**
     datavec[0]= Hdiv Resconstructed
     datavec[1]= Pressure Resconstructed
     datavec[2]= Hdiv FEM
     datavec[3]= Pressure FEM
    
     **/
    
    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);
    
    STATE perm = MixedMat::GetPermeability(datavec[1].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);
    int dim = MixedMat::fDim;


    STATE pressureexact = 0.;
    TPZManVector<STATE,2> pressvec(1,0.);
    TPZFNMatrix<9, STATE> gradu(dim, 1, 0.), fluxinv(dim, 1);

    if(MixedMat::fExactSol)
    {

        MixedMat::fExactSol(datavec[0].x, pressvec,gradu);
        gradu.Resize(3, 1);
        //gradu(2,0) = 0.;
    }
    
    PermTensor.Multiply(gradu, fluxinv);
    pressureexact = pressvec[0];
    switch (var)
    {
        case 40://FluxFem
            for(int i=0; i<dim; i++) Solout[i] = datavec[2].sol[0][i];
            break;
        case 41://FluxReconstructed
            for (int i=0; i<dim; i++) Solout[i] = datavec[0].sol[0][i];
            break;
        case 42:
            for(int i=0; i<dim; i++) Solout[i] = -fluxinv(i);
            break;
        case 43://PressureFem
            Solout[0] = datavec[3].sol[0][0];
            break;
        case 44://PressureReconstructed
            Solout[0] = datavec[1].sol[0][0];
            break;
        case 45:
            Solout[0] = pressureexact;
            break;
        case 46:
            Solout[0] = datavec[1].p;
            break;
        default:
            DebugStop();
    }
}


/// make a contribution to the error computation
template<class MixedMat>
void TPZMixedHDivErrorEstimate<MixedMat>::Errors(const TPZVec<TPZMaterialDataT<STATE> > &data, TPZVec<REAL> &errors)
{
    
    /**
     datavec[0]= Hdiv Reconstructed
     datavec[1]= Pressure Reconstructed
     datavec[2]= Hdiv FEM
     datavec[3]= Pressure FEM
     **/
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    const int dim = MixedMat::fDim;

    TPZManVector<STATE,3> fluxfem(3), fluxreconstructed(3), pressurefem(1), pressurereconstructed(1);

    
    fluxreconstructed = data[0].sol[0];
    fluxfem = data[2].sol[0];
    STATE divsigmafem=data[2].divsol[0][0];
    
    
    TPZVec<STATE> divsigma(1,0.);

    TPZManVector<STATE, 1> u_exact(1);
    TPZFNMatrix<9, STATE> du_exact(3, 3);
    if(this->HasExactSol()){
        this->ExactSol()(data[0].x,u_exact,du_exact);
    }
    if(this->ForcingFunction())
    {
        this->ForcingFunction()(data[0].x,divsigma);

    }

    
    REAL residual = 0.;
//    std::cout<<" divsigma "<<divsigma[0] << " divsigmafem "<<divsigmafem<<"\n";
    
    residual = (divsigma[0] - divsigmafem)*(divsigma[0] - divsigmafem);
    
    
    pressurereconstructed[0] = data[1].sol[0][0];
    pressurefem[0] = data[3].sol[0][0];
    
    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);
    
    REAL perm = MixedMat::GetPermeability(data[1].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);
    
    
    TPZFNMatrix<3,REAL> fluxexactneg;

    //sigma=-K grad(u)
    {
        TPZFNMatrix<9,REAL> gradpressure(dim,1);
        for (int i=0; i<dim; i++) {
            gradpressure(i,0) = du_exact[i];
        }
        PermTensor.Multiply(gradpressure,fluxexactneg);
    }
    
    
    REAL innerexact = 0.;
    REAL innerestimate = 0.;
    for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++) {
            innerexact += (fluxfem[i]+fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]+fluxexactneg(j,0));//Pq esta somando: o fluxo fem esta + e o exato -
            innerestimate += (fluxfem[i]-fluxreconstructed[i])*InvPermTensor(i,j)*(fluxfem[j]-fluxreconstructed[j]);
        }
    }
    errors[0] = (pressurefem[0]-u_exact[0])*(pressurefem[0]-u_exact[0]);//exact error pressure
    errors[1] = (pressurefem[0]-pressurereconstructed[0])*(pressurefem[0]-pressurereconstructed[0]);//error pressure reconstructed
    errors[2] = innerexact;//error flux exact
    errors[3] = innerestimate;//error flux reconstructed
    errors[4] = residual;
    for(int i = 0; i<5; i++)
    {
        if(std::isnan(errors[i])) DebugStop();
    }
}





template class TPZMixedHDivErrorEstimate<TPZMixedDarcyFlow>;
