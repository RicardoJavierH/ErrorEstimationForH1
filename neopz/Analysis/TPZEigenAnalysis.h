//
// Created by Francisco Teixeira Orlandini on 11/17/17.
//

#ifndef PZ_TPZEIGENANALYSIS_H
#define PZ_TPZEIGENANALYSIS_H
#include "TPZAnalysis.h"     //For TPZAnalysis
#include "tpzautopointer.h" //For TPZAutoPointer
#include "pzmatrix.h"       //For TPZFMatrix

template<typename TVar>
class TPZEigenSolver;

/**
 * @ingroup Analysis
 * @brief Performs the Finite Element Analysis of a (standard or generalised) eigenvalue problem.
 */
class TPZEigenAnalysis : public TPZAnalysis{
public:
    enum class Mat{A,B};//< Which matrix should be assembled
    /** @name Constructors*/
    /** @{ */
    /** @brief Create an empty TPZEigenAnalysis object*/
    TPZEigenAnalysis();
    
#ifdef PZ_USING_METIS
    /** @brief Create an TPZEigenAnalysis object from one mesh pointer */
    TPZEigenAnalysis(TPZCompMesh *mesh, const RenumType& renumtype = RenumType::EMetis, std::ostream &out = std::cout);
    /** @brief Create an TPZEigenAnalysis object from one mesh auto pointer object */
    TPZEigenAnalysis(TPZAutoPointer<TPZCompMesh> mesh, const RenumType& renumtype = RenumType::EMetis, std::ostream &out = std::cout);
#else
  /** @brief Create an TPZEigenAnalysis object from one mesh pointer */
  TPZEigenAnalysis(TPZCompMesh *mesh, const RenumType& renumtype = RenumType::ESloan, std::ostream &out = std::cout);
  /** @brief Create an TPZEigenAnalysis object from one mesh auto pointer object */
  TPZEigenAnalysis(TPZAutoPointer<TPZCompMesh> mesh, const RenumType& renumtype = RenumType::ESloan, std::ostream &out = std::cout);

#endif

    /** @} */
    /** @name FEM*/
    /** @{ */
    /** @brief Gets the eigensolver */
    template<class TVar>
    TPZEigenSolver<TVar> &EigenSolver();
    /** @brief Set the solver
      @note In this function it will be checked if the solver is a TPZEigenSolver*/
    void SetSolver(const TPZSolver &solver) override;
    /** @brief Assemble the matrices associated with the EVP*/
    void Assemble() override;
    /** @brief Assemble one of the matrices associated with the EVP*/
    void AssembleMat(const TPZEigenAnalysis::Mat mat);
    /** @brief Solve the EVP problem*/
    void Solve() override;
    /** @brief Gets the eigenvectors calculated by the Solve method*/
    TPZFMatrix<CSTATE> &GetEigenvectors()
    {return fEigenvectors;}
    /** @brief Gets the eigenvalues  by the Solve method*/
    TPZVec<CSTATE> &GetEigenvalues()
    {return fEigenvalues;}
    /** @brief Gets the eigenvectors calculated by the Solve method*/
    const TPZFMatrix<CSTATE> &GetEigenvectors() const
    {return fEigenvectors;}
    /** @brief Gets the eigenvalues  by the Solve method*/
    const TPZVec<CSTATE> &GetEigenvalues() const
    {return fEigenvalues;}
    /** @brief Sets the eigenvectors*/
    void SetEigenvectors(const TPZFMatrix<CSTATE> &ev){fEigenvectors = ev;}
    /** @brief Sets the eigenvalues*/
    void SetEigenvalues(const TPZVec<CSTATE> &ev){fEigenvalues = ev;}
    /** @brief Set to compute eigenvectors or just eigenvalues*/
    inline void SetComputeEigenvectors(const bool opt);
    /** @brief Whether to compute eigenvectors or just eigenvalues*/
    inline bool ComputeEigenvectors() const;
    /** @} */
    /** @name ReadWrite*/
    /**  @{ */
    //! Class unique identifier
    int ClassId() const override;
    //! Write attributes to TPZStream
    void Write(TPZStream &buf, int withclassid) const override;
    //! Read attributes from TPZStream to replicate instance from file
    void Read(TPZStream &buf, void *context) override;
    /** @} */
protected:
    /**
    * @brief Stores the computed eigenvalues
    */
    TPZManVector<CSTATE,10> fEigenvalues;
    /**
     * @brief Stores the computed eigenvectors
     */
    TPZFMatrix<CSTATE> fEigenvectors;
    //! Whether to compute eigenvectors
    bool fCalcVectors{true};
private:
    /*the following methods are set to private until they are properly implemented*/
    using TPZAnalysis::PrePostProcessTable;
    using TPZAnalysis::PostProcessTable;
    using TPZAnalysis::PostProcessError;
    using TPZAnalysis::PostProcessErrorSerial;
    using TPZAnalysis::PostProcessErrorParallel;
    template<class TVar, Mat MAT>
    void AssembleT();
    template<class TVar>
    void SolveT();
    template<class TVar>
    void ConfigAssemble();
    bool fIsSetUp{false};
};


void TPZEigenAnalysis::SetComputeEigenvectors(const bool opt)
{
    fCalcVectors = opt;
}

bool TPZEigenAnalysis::ComputeEigenvectors() const
{
    return fCalcVectors;
}

#define INSTANTIATE_TEMPLATES(TVar)                                            \
  extern template TPZEigenSolver<TVar> &TPZEigenAnalysis::EigenSolver<TVar>();

INSTANTIATE_TEMPLATES(STATE)
INSTANTIATE_TEMPLATES(CSTATE)
#undef INSTANTIATE_TEMPLATES

#endif //PZ_TPZEIGENANALYSIS_H
