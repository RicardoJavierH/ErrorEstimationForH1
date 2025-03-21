//
// Created by victor on 29/08/2022.
//

#include "TPZStructMatrixOMPorTBB.h"
/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixOMPorTBB methods.
 */

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include "TPZLinearAnalysis.h"
#include "pzsfulmat.h"
#include "TPZParallelUtils.h"
#include "pzgnode.h"
#include "TPZTimer.h"
#include "TPZElementMatrixT.h"
#include "TPZSYSMPMatrix.h"
#include "TPZThreadPool.h"

#include "TPZLapack.h"

#ifdef USING_TBB
#include <tbb/parallel_for.h>

#include <tbb/info.h>
#include <tbb/task_arena.h>
#include <tbb/global_control.h>
#endif

#ifdef USING_OMP
#include "omp.h"
#endif
#include "pzcheckconsistency.h"
#include "TPZMaterial.h"
#include "TPZStrMatParInterface.h"

using namespace std;

#include "pzlog.h"

#include "run_stats_table.h"
#include <thread>

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.TPZStructMatrixOMPorTBB");
static TPZLogger loggerel("pz.strmatrix.element");
static TPZLogger loggerel2("pz.strmatrix.elementinterface");
static TPZLogger loggerelmat("pz.strmatrix.elementmat");
static TPZLogger loggerCheck("pz.strmatrix.checkconsistency");
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif

RunStatsTable stat_ass_graph_LCC("-ass_graph_ot", "Run statistics table for the graph creation and coloring TPZStructMatrixOMPorTBB.");

template<class TVar>
TPZStructMatrixOMPorTBB<TVar>::TPZStructMatrixOMPorTBB(): TPZStrMatParInterface()
{
this->SetNumThreads(TPZThreadPool::globalInstance().threadCount());
}

template<class TVar>
TPZStructMatrixOMPorTBB<TVar>::TPZStructMatrixOMPorTBB(const TPZStructMatrixOMPorTBB &copy) : TPZStrMatParInterface(copy) {
    fElSequenceColor = copy.fElSequenceColor;
    fElBlocked = copy.fElBlocked;
    fElementsComputed = copy.fElementsComputed;
    fElementCompleted = copy.fElementCompleted;
    fSomeoneIsSleeping = copy.fSomeoneIsSleeping;
    fShouldColor = copy.fShouldColor;
    fUsingTBB = copy.fUsingTBB;
}


static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
#ifdef USING_MKL
    mkl_domain_set_num_threads(1, MKL_DOMAIN_BLAS);
    //mkl_set_num_threads_local(1);
#endif
    ass_stiff.start();

    const auto &equationFilter =
            (dynamic_cast<TPZStructMatrix *>(this))->EquationFilter();
    ass_stiff.start();

    if (equationFilter.IsActive()) {
        int64_t neqcondense = equationFilter.NActiveEquations();
#ifdef PZDEBUG
        if (mat.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<STATE> rhsloc(neqcondense, rhs.Cols(), 0.);
        if (this->fNumThreads) {
            this->MultiThread_Assemble(mat, rhsloc, guiInterface);
        } else {
            this->Serial_Assemble(mat, rhsloc, guiInterface);
        }

        equationFilter.Scatter(rhsloc, rhs);
    } else {
        if (this->fNumThreads) {
            this->MultiThread_Assemble(mat, rhs, guiInterface);
        } else {
            this->Serial_Assemble(mat, rhs, guiInterface);
        }
    }
    ass_stiff.stop();
#ifdef USING_MKL
    mkl_set_num_threads_local(0);
#endif
}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_rhs.start();
    const auto &equationFilter =
            (dynamic_cast<TPZStructMatrix *>(this))->EquationFilter();
    if(equationFilter.IsActive())
    {
        int64_t neqcondense = equationFilter.NActiveEquations();
        int64_t neqexpand = equationFilter.NEqExpand();

        TPZFMatrix<TVar> &source = dynamic_cast<TPZFMatrix<TVar>&>(rhs);
        if(rhs.Rows() != neqexpand || Norm(source) != 0.)
        {
            DebugStop();
        }
        TPZFMatrix<STATE> rhsloc(neqcondense,1,0.);
        if(this->fNumThreads)
        {
#ifdef HUGEDEBUG
            TPZFMatrix<STATE> rhsserial(rhsloc);
            this->Serial_Assemble(rhsserial, guiInterface);
#endif
            this->MultiThread_Assemble(rhsloc,guiInterface);
#ifdef HUGEDEBUG
            rhsserial -= rhsloc;
            REAL norm = Norm(rhsserial);
            std::cout << "difference between serial and parallel " << norm << std::endl;
#endif
        }
        else
        {
            this->Serial_Assemble(rhsloc,guiInterface);
        }
        equationFilter.Scatter(rhsloc,rhs);
    }
    else
    {
        if(this->fNumThreads){
#ifdef HUGEDEBUG
            TPZFMatrix<STATE> rhsserial(rhs);
            this->Serial_Assemble(rhsserial, guiInterface);
#endif
            this->MultiThread_Assemble(rhs,guiInterface);
#ifdef HUGEDEBUG
            REAL normrhs = Norm(rhs);
            REAL normrhsserial = Norm(rhsserial);
            std::cout << "normrhs = " << normrhs << " normrhsserial " << normrhsserial << std::endl;
            rhsserial -= rhs;
            REAL norm = Norm(rhsserial);
            std::cout << "difference between serial and parallel " << norm << std::endl;
#endif
        }
        else{
            this->Serial_Assemble(rhs,guiInterface);
        }
    }
    ass_rhs.stop();
}


template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::Serial_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    //mkl_set_num_threads_local(1);
    //mkl_domain_set_num_threads(1, MKL_DOMAIN_BLAS);

    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();

    TPZMatrix<TVar> &stiffness = dynamic_cast<TPZMatrix<TVar>&>(mat);
    TPZFMatrix<TVar> &source = dynamic_cast<TPZFMatrix<TVar>&>(rhs);

    const auto &equationFilter =
            (dynamic_cast<TPZStructMatrix *>(this))->EquationFilter();

    int64_t nelem = cmesh->NElements();
    std::cout << "LCC_Serial_Assemble\n";
    for (int64_t iel = 0; iel < nelem; iel++)
    {
        TPZCompEl *el = cmesh->Element(iel);
        if (!el) continue;


        TPZElementMatrixT<TVar> ek(cmesh, TPZElementMatrix::EK), ef(cmesh, TPZElementMatrix::EF);

        el->CalcStiff(ek, ef);

        if(!el->HasDependency()) {
            ek.ComputeDestinationIndices();
            equationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            //            TPZSFMatrix<STATE> test(stiffness);
            //            TPZFMatrix<STATE> test2(stiffness.Rows(),stiffness.Cols(),0.);
            //            stiffness.Print("before assembly",std::cout,EMathematicaInput);
            stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            source.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //            stiffness.Print("stiffness after assembly STK = ",std::cout,EMathematicaInput);
            //            rhs.Print("rhs after assembly Rhs = ",std::cout,EMathematicaInput);
            //            test2.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //            test -= stiffness;
            //            test.Print("matriz de rigidez diference",std::cout);
            //            test2.Print("matriz de rigidez interface",std::cout);


        } else {
            // the element has dependent nodes
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            equationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);

            stiffness.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            source.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);

        }
    }
#ifdef PZDEBUG
    VerifyStiffnessSum(mat);
#endif
#ifdef USING_MKL
    mkl_set_num_threads_local(0);
#endif
}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::VerifyStiffnessSum(TPZBaseMatrix & mat){

    TPZMatrix<TVar> &stiffness = dynamic_cast<TPZMatrix<TVar>&>(mat);

    REAL totalSum = 0;
    for(int irow =0; irow < mat.Rows();irow++){
        for (int icol = 0; icol < mat.Cols(); icol++)
            totalSum += abs(stiffness.GetVal( irow, icol));
    }
    std::cout << "totalSum stisffnes =" << totalSum <<std::endl;
}
template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::Serial_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){


}

template<class TVar>
int TPZStructMatrixOMPorTBB<TVar>::GetNumberColors(){
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();
    int64_t nelem = cmesh->NElements();
    if (fElVecColor.size() != nelem) DebugStop();
    int ncolor = -1;
    for (int iel=0; iel<nelem; iel++){
        if (ncolor < fElVecColor[iel])
            ncolor = fElVecColor[iel];
    }
    ncolor++;
    return ncolor;
}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::MultiThread_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    //mkl_domain_set_num_threads(1, MKL_DOMAIN_BLAS);
    //mkl_set_num_threads_local(1);
    if (fShouldColor){
        if (fUsingTBB){
            AssemblingUsingTBBandColoring(mat,rhs);

        }else{
            AssemblingUsingOMPandColoring(mat,rhs);
        }

    }else{
        if (fUsingTBB){
            AssemblingUsingTBBbutNotColoring(mat,rhs);

        }else{
            AssemblingUsingOMPbutNotColoring(mat,rhs);
        }
    }

#ifdef PZDEBUG
    VerifyStiffnessSum(mat);
#endif
}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::AssemblingUsingTBBandColoring(TPZBaseMatrix & mat, TPZBaseMatrix & rhs ){
#ifndef USING_TBB
    DebugStop();
#endif
#ifdef USING_TBB
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();

    int64_t nelem = cmesh->NElements();
    const int nthread = this->fNumThreads;

    TPZStructMatrixOMPorTBB::OrderElements();
    int ncolor = GetNumberColors();

    for (int icol=0; icol<ncolor; icol++){

        //tbb::task_scheduler_init init(nthread); //dont work in computer of LABMEC
        tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, nthread);
        tbb::parallel_for( tbb::blocked_range<int64_t>(0,nelem),
                          [&](tbb::blocked_range<int64_t> r){
        for (int64_t iel = r.begin(); iel < r.end(); iel++)
            {
                if (icol != fElVecColor[iel]) continue;
                TPZCompEl *el = cmesh->Element(iel);
                if (!el) continue;

                ComputingCalcstiffAndAssembling(mat,rhs,el);

            }
        });
    }
#endif
}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::AssemblingUsingOMPandColoring(TPZBaseMatrix & mat,TPZBaseMatrix & rhs ){
#ifdef USING_OMP

    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();

    int64_t nelem = cmesh->NElements();
    const int nthread = this->fNumThreads;

    TPZStructMatrixOMPorTBB::OrderElements();
    int ncolor = GetNumberColors();

    for (int icol=0; icol<ncolor; icol++){

        omp_set_num_threads(nthread);
        #pragma omp parallel for schedule(dynamic,1)
        for (int64_t iel = 0; iel < nelem; iel++){
        if (icol != fElVecColor[iel]) continue;
        TPZCompEl *el = cmesh->Element(iel);
        if (!el) continue;

            ComputingCalcstiffAndAssembling(mat,rhs,el);

        }
    }
#else
    DebugStop();
#endif
}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::AssemblingUsingTBBbutNotColoring(TPZBaseMatrix & mat, TPZBaseMatrix & rhs ){
#ifdef USING_TBB
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();
    int64_t nelem = cmesh->NElements();
    const int nthread = this->fNumThreads;


        //tbb::task_scheduler_init init(nthread); //dont work in computer of LABMEC
        tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, nthread);
        tbb::parallel_for( tbb::blocked_range<int64_t>(0,nelem),
                          [&](tbb::blocked_range<int64_t> r){
        for (int64_t iel = r.begin(); iel < r.end(); iel++)
        {
            TPZCompEl *el = cmesh->Element(iel);
            if (!el) continue;

            ComputingCalcstiffAndAssembling(mat,rhs,el);

        }
        });
#else
    DebugStop();
#endif

}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::AssemblingUsingOMPbutNotColoring(TPZBaseMatrix & mat, TPZBaseMatrix & rhs ){
#ifdef USING_OMP

    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();

    int64_t nelem = cmesh->NElements();
    const int nthread = this->fNumThreads;

    omp_set_num_threads(nthread);
    #pragma omp parallel for schedule(dynamic,1)
    for (int64_t iel = 0; iel < nelem; iel++){
        {
            TPZCompEl *el = cmesh->Element(iel);
            if (!el) continue;

            ComputingCalcstiffAndAssembling(mat,rhs,el);
        }
    }
#else
    DebugStop();
#endif
}



template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::MultiThread_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    std::cout <<"\n\n\nPlease implement TPZStructMatrixOMPorTBB::MultiThread_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)\n\n\n";
    DebugStop();
}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::ComputingCalcstiffAndAssembling(TPZBaseMatrix & mat,TPZBaseMatrix & rhs,TPZCompEl *el){

    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();

    TPZMatrix<TVar> &stiffness = dynamic_cast<TPZMatrix<TVar>&>(mat);
    TPZFMatrix<TVar> &source = dynamic_cast<TPZFMatrix<TVar>&>(rhs);

    TPZElementMatrixT<TVar> ek(cmesh, TPZElementMatrix::EK), ef(cmesh, TPZElementMatrix::EF);
    el->CalcStiff(ek, ef);

    const auto &equationFilter =
            (dynamic_cast<TPZStructMatrix *>(this))->EquationFilter();

    if(!el->HasDependency()) {
        ek.ComputeDestinationIndices();
        equationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);

        if (fShouldColor){
            stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            source.AddFelNonAtomic(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        }else{
            stiffness.AddKelAtomic(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            source.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        }

    } else {
        // the element has dependent nodes
        ek.ApplyConstraints();
        ef.ApplyConstraints();
        ek.ComputeDestinationIndices();
        equationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
        if (fShouldColor){
            stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            source.AddFelNonAtomic(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        }else{
            stiffness.AddKelAtomic(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            source.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        }

    }
}



static void AssembleColor(int64_t el, TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute) {
    for (auto connect : connectlist) {
        elContribute[connect] = el;
    }
}

//Coloring computational elements 2022
template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::OrderElements(){
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();
    int64_t nconnects = cmesh->NConnects();
    int64_t nelement = cmesh->ElementVec().NElements();
    fElVecColor.Resize(nelement);
    fElVecColor.Fill(-1);
    TPZVec<int> used(nconnects,0);
    bool haswork = true;
    int color = 0;

    while (haswork){
        haswork = false;
        for (int64_t iel = 0; iel < nelement; iel++){
            TPZCompEl *compel = cmesh->Element(iel);
            if (!compel) continue;
            if (fElVecColor[iel] != -1) continue;

            TPZStack<int64_t> connectlist;
            compel->BuildConnectList(connectlist);

            bool canColor = true;
            for (auto ic : connectlist){
                if (used[ic]){
                    canColor = false;
                    haswork = true;
                }
            }
            if (canColor){
                for (auto ic : connectlist) used[ic] = 1;
                fElVecColor[iel] = color;
            }
        }
        if (haswork) {
            color++;
            used.Fill(0);
        }
    }
}


template<class TVar>
int TPZStructMatrixOMPorTBB<TVar>::ClassId() const{
    return Hash("TPZStructMatrixOMPorTBB") ^
        ClassIdOrHash<TVar>() << 1 ^
        TPZStrMatParInterface::ClassId();
}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::Read(TPZStream& buf, void* context) {
    TPZStrMatParInterface::Read(buf, context);

    buf.Read(fElBlocked);
    buf.Read(fElSequenceColor);
    buf.Read(&fElementCompleted);
    buf.Read(&fSomeoneIsSleeping);
    int64_t fCurrentIndexLong;
    buf.Read(&fCurrentIndexLong);
    fCurrentIndex = fCurrentIndexLong;
}

template<class TVar>
void TPZStructMatrixOMPorTBB<TVar>::Write(TPZStream& buf, int withclassid) const {
    TPZStrMatParInterface::Write(buf, withclassid);

    buf.Write(fElBlocked);
    buf.Write(fElSequenceColor);
    buf.Write(&fElementCompleted);
    buf.Write(&fSomeoneIsSleeping);
    int64_t fCurrentIndexLong;
    fCurrentIndexLong = fCurrentIndex;
    buf.Write(&fCurrentIndexLong);
}

template class TPZStructMatrixOMPorTBB<STATE>;
template class TPZRestoreClass<TPZStructMatrixOMPorTBB<STATE>>;
template class TPZStructMatrixOMPorTBB<CSTATE>;
template class TPZRestoreClass<TPZStructMatrixOMPorTBB<CSTATE>>;
