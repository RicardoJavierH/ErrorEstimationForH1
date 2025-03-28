//
//  TPZHDivErrorEstimateMaterial.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//
/**
 This material implement the Ainsworth proposal
 **/

#ifndef TPZHDivErrorEstimateMaterial_hpp
#define TPZHDivErrorEstimateMaterial_hpp

#include <stdio.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"

typedef TPZMixedDarcyFlow TPZMixedPoisson;

class TPZHDivErrorEstimateMaterial : public TPZMixedPoisson {

public:
    TPZHDivErrorEstimateMaterial(int matid, int dim);

    TPZHDivErrorEstimateMaterial();
    
    TPZHDivErrorEstimateMaterial(const TPZHDivErrorEstimateMaterial &copy);
    
    TPZHDivErrorEstimateMaterial(const TPZMixedPoisson &copy);
    
    virtual ~TPZHDivErrorEstimateMaterial();
    
    TPZHDivErrorEstimateMaterial &operator=(const TPZHDivErrorEstimateMaterial &copy);

    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    virtual void FillBoundaryConditionDataRequirements(int type,  TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    bool fNeumannLocalProblem = false;

    void SetNeumannProblem(bool neumannProblem) {
        fNeumannLocalProblem = neumannProblem;
    }

    virtual int NEvalErrors() const override { return 5; } // erro de oscilacao de dados tbem

    /// Compute the error and error estimate
    // error[0] - error computed with exact pressure
    // error[1] - error computed with reconstructed pressure
    // error[2] - energy error computed with exact solution
    // error[3] - energy error computed with reconstructed solution
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;

    void ErrorsBC(TPZVec<TPZMaterialDataT<STATE>> &data,
                  TPZVec<REAL> &errors, TPZBndCondT<STATE> &bc);

    virtual int VariableIndex(const std::string &name) const override;
    virtual int NSolutionVariables(int var) const override;
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    // Returns the first non-null approximation space index, which will be the
    // H1 reconstruction space
    int FirstNonNullApproxSpaceIndex(const TPZVec<TPZMaterialDataT<STATE>> &datavec) const;
};

#endif /* TPZHDivErrorEstimateMaterial_hpp */
