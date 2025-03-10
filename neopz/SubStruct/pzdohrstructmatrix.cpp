/**
 * @file
 * @brief Contains the implementation of the TPZDohrStructMatrix methods.
 */

#include "pzdohrstructmatrix.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrsubstructCondense.h"
#include "tpzdohrprecond.h"
#include "tpznodesetcompute.h"
#include "TPZRenumbering.h"
#include "pzmetis.h"

#include "pzskylstrmatrix.h"
#include "pzmatred.h"
#include "tpzmatredstructmatrix.h"
#include "tpzpairstructmatrix.h"
#include "pzfstrmatrix.h"

#include "pzsubcmesh.h"
#include "pzintel.h"

#include "TPZBoostGraph.h"
#include "pzsloan.h"
#include "pzvisualmatrix.h"
#include "TPZRefPatternTools.h"
#include "tpzverysparsematrix.h"

#include <sstream>
#include <map>
#include "pzlog.h"

#include "TPZSimpleTimer.h"
#include "TPZTimeTemp.h"
#include "TPZVTKGeoMesh.h"
#include <stdlib.h>

#ifdef PZ_LOG
static TPZLogger logger("structmatrix.dohrstructmatrix");
static TPZLogger loggerasm("structmatrix.dohrstructmatrix.asm");
#endif

#include "clock_timer.h"
#include "timing_analysis.h"
#include "arglib.h"
#include "run_stats_table.h"

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
using namespace tbb;
#endif

#include <thread>

using namespace std;

#include <rcm.h>
#include <TPZBFileStream.h>

#ifdef USING_PAPI
#include <papi.h>

static float stiff_sum = 0;
#endif

/** @brief Return the number of submeshes */
static int64_t NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh);

/** @brief return a pointer to the isub submesh */
static TPZSubCompMesh *SubMesh(TPZAutoPointer<TPZCompMesh> compmesh, int isub);

template<class TVar>
static void AssembleMatrices(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense<TVar> > substruct, TPZAutoPointer<TPZDohrAssembly<TVar> > dohrassembly,
                             mutex* TestThread);
template<class TVar>
static void DecomposeBig(TPZAutoPointer<TPZDohrSubstructCondense<TVar> > substruct, int numa_node);
template<class TVar>

static void DecomposeInternal(TPZAutoPointer<TPZDohrSubstructCondense<TVar> > substruct, int numa_node);

#define NOMETIS \
    PZError<<"TPZDohrStructMatrix requires Metis library\n";\
    PZError<<"Please reconfigure NeoPZ library with:\n";\
    PZError<<"PZ_USING_METIS=ON"<<std::endl;\
    DebugStop();
template<class TVar, class TPar>
 TPZDohrStructMatrix<TVar,TPar>::TPZDohrStructMatrix() :
TPZStructMatrixT<TVar>(), fDohrAssembly(0), fDohrPrecond(0), fAccessElement()
{
#ifndef PZ_USING_METIS
    NOMETIS
#endif
}

template<class TVar, class TPar>
 TPZDohrStructMatrix<TVar,TPar>::TPZDohrStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh) :
TPZStructMatrixT<TVar>(cmesh), fDohrAssembly(0),
fDohrPrecond(0), fAccessElement()
{
#ifndef PZ_USING_METIS
    NOMETIS
#endif
}

template<class TVar, class TPar>
 TPZDohrStructMatrix<TVar,TPar>::TPZDohrStructMatrix(const TPZDohrStructMatrix &copy) :
TPZStructMatrixT<TVar>(copy), fDohrAssembly(copy.fDohrAssembly), fDohrPrecond(copy.fDohrPrecond), fAccessElement()
{
#ifndef PZ_USING_METIS
    NOMETIS
#endif
}

template<class TVar, class TPar>
 TPZDohrStructMatrix<TVar,TPar>::~TPZDohrStructMatrix()
{
}


// this will create a DohrMatrix
template<class TVar, class TPar>
TPZMatrix<TVar> * TPZDohrStructMatrix<TVar,TPar>::Create()
{
#ifndef PZ_USING_METIS
    NOMETIS
#endif    
	TPZSimpleTimer timeforcopute; // init of timer for compute
	this->fMesh->ComputeNodElCon();
	TPZAutoPointer<TPZDohrAssembly<TVar> > assembly = new TPZDohrAssembly<TVar>;
	fDohrAssembly = assembly;
	
	this->fMesh->InitializeBlock();
	{
		TPZVec<int64_t> perm,iperm;
		TPZStack<int64_t> elgraph,elgraphindex;
		
		
		int nindep = this->fMesh->NIndependentConnects();
		this->fMesh->ComputeElGraph(elgraph,elgraphindex);
		int nel = elgraphindex.NElements()-1;
#ifdef USING_BOOST
		TPZBoostGraph boost(nel,nindep);
		boost.setGType(TPZBoostGraph::KMC);
		boost.SetElementGraph(elgraph, elgraphindex);
		boost.CompressedResequence(perm, iperm);
#else
		TPZSloan sloan(nel,nindep);
		sloan.SetElementGraph(elgraph, elgraphindex);
		sloan.Resequence(perm, iperm);
#endif
		this->fMesh->Permute(perm);
	}
	int nsub = NSubMesh(this->fCompMesh);
	int isub;
	
	for(isub=0; isub<nsub; isub++)
	{
		TPZSubCompMesh *submesh = SubMesh(this->fCompMesh, isub);
#ifdef PZDEBUG
		std::cout << '.'; std::cout.flush();
#endif
		if(!submesh) 
		{
			continue;
		}
		TPZVec<int64_t> perm,iperm;
		TPZStack<int64_t> elgraph,elgraphindex;
		int64_t nindep = submesh->NIndependentConnects();
		submesh->ComputeElGraph(elgraph,elgraphindex);
		int64_t nel = elgraphindex.NElements()-1;
#ifdef USING_BOOST
		TPZBoostGraph boost(nel,nindep);
		boost.setGType(TPZBoostGraph::KMC);
		boost.SetElementGraph(elgraph, elgraphindex);
		boost.CompressedResequence(perm, iperm);
#else
		TPZSloan sloan(nel,nindep);
		sloan.SetElementGraph(elgraph, elgraphindex);
		sloan.Resequence(perm, iperm);
#endif
		
		submesh->Permute(perm);
#ifdef PZDEBUG 
		std::stringstream filename;
		filename << "SubMatrix" << submesh->Index() << ".vtk";
		TPZFMatrix<REAL> fillin(50,50);
		submesh->ComputeFillIn(50, fillin);
		VisualMatrix(fillin, filename.str().c_str());
#endif
	}		
	
	tempo.ft1comput = timeforcopute.ReturnTimeDouble(); //end of time for compute
#ifdef PZDEBUG
	std::cout << tempo.ft1comput << std::endl;
	std::cout << "Identifying corner nodes\n";
	TPZSimpleTimer timefornodes; // init of timer
#endif
	
	IdentifyCornerNodes();
    
#ifdef PZDEBUG
	tempo.ft4identcorner = timefornodes.ReturnTimeDouble();
	std::cout << "Total for Identifying Corner Nodes: " << tempo.ft4identcorner << std::endl; // end of timer
#endif
    
	TPZDohrMatrix<TVar,TPZDohrSubstructCondense<TVar> > *dohr = new TPZDohrMatrix<TVar,TPZDohrSubstructCondense<TVar> >(assembly);

	int64_t neq = this->fMesh->NEquations();
	dohr->Resize(neq,neq);
	// fCornerEqs was initialized during the mesh generation process
	dohr->SetNumCornerEqs(this->fCornerEqs.size());
	
	assembly->fFineEqs.Resize(nsub);
	assembly->fCoarseEqs.Resize(nsub);
	for(isub=0; isub<nsub; isub++)
	{
		TPZSubCompMesh *submesh = SubMesh(this->fCompMesh, isub);
		if(!submesh) 
		{
			continue;
		}
		TPZAutoPointer<TPZDohrSubstructCondense<TVar> > substruct = new TPZDohrSubstructCondense<TVar>();
		submesh->ComputeNodElCon();
		int64_t neq = ((TPZCompMesh *)submesh)->NEquations();
		//    int neq = substruct->fStiffness->Rows();
		
		substruct->fNEquations = neq;
		
		
		// identify the equation numbers of the submesh
		std::map<int,int> globinv,globaleqs;
		// initialize the globaleqs data structure
		// the global eqs are ordered in the sequence the connects appear
		IdentifyEqNumbers(submesh, globaleqs ,globinv);
		int next = globaleqs.size();
		substruct->fNumExternalEquations = next;
        substruct->fNumInternalEquations = submesh->NumInternalEquations();
		assembly->fFineEqs[isub].Resize(next);
		std::map<int,int>::iterator it;
		int count = 0;
		for(it=globaleqs.begin(); it!=globaleqs.end(); it++)
		{
			assembly->fFineEqs[isub][count++] = it->second; 
		}
        
		
		// initialize the permutations from the mesh enumeration to the external enumeration
		typedef typename TPZDohrSubstructCondense<TVar>::ENumbering ENumbering;
		typedef std::pair<ENumbering,ENumbering> Numberingpair;
		ENumbering tsub,text,tint;
		tsub = TPZDohrSubstructCondense<TVar>::Submesh;
		text = TPZDohrSubstructCondense<TVar>::ExternalFirst;
		tint = TPZDohrSubstructCondense<TVar>::InternalFirst;
		
		TPZVec<int> &toexternal = substruct->fPermutationsScatter[Numberingpair(tsub,text)];
		TPZVec<int> &fromexternal = substruct->fPermutationsScatter[Numberingpair(text,tsub)];
		toexternal.Resize(neq,-1);
		fromexternal.Resize(neq,-1);
		int nel = globaleqs.size();
		count = 0;
		for(it=globaleqs.begin(); it!=globaleqs.end(); it++)
		{
			toexternal[it->first] = count++;
		}
		count = nel++;
		int ieq;
		for(ieq=0; ieq<neq; ieq++)
		{
			if(toexternal[ieq] == -1) toexternal[ieq] = count++;
		}
		for(ieq=0; ieq<neq; ieq++)
		{
			fromexternal[toexternal[ieq]] = ieq;
		}
		
		ComputeInternalEquationPermutation(submesh, substruct->fPermutationsScatter[Numberingpair(tsub,tint)], substruct->fPermutationsScatter[Numberingpair(tint,tsub)]);
		//		IdentifyEqNumbers(submesh, substruct->fGlobalIndex,globinv);
		
		// initialize the fC matrix
		// associate each column of the fC matrix with a coarse index
		IdentifySubCornerEqs(globinv,substruct->fCoarseNodes,assembly->fCoarseEqs[isub]);
        
        
		//		int ncoarse = substruct->fCoarseNodes.NElements();
		
		// reorder by internal nodes
		// the fInternalEqs data structure will not be filled if the connects are made internal
		
		// this permutes the nodes of the submesh
		// This is a lengthy process which should run on the remote processor
		dohr->AddSubstruct(substruct);
	}
	return dohr;
}


template<class TVar>
class parallel_assemble_task_t
{
private:
    
    /** We divide the assemble procedure into N work items, which will
     be executed by one or several tasks. The TBB parallel_for
     construct automatically divide the work items in subsets and
     "creates" tasks to execute the work in each subset. Each task
     invokes the operator(blocked_range subset), which will be
     responsible for executing the work items in the subset. */
    template<class TTVar>
    struct work_item_t
    {
        work_item_t (unsigned submesh_idx, const TPZAutoPointer<TPZDohrSubstructCondense<TTVar> >& substruct) :
        fSubMeshIndex(submesh_idx), fSubstruct(substruct) {}
        
        unsigned fSubMeshIndex;
        TPZAutoPointer<TPZDohrSubstructCondense<TTVar> > fSubstruct;
    };
    
    /** Array of work items. */
    std::vector<work_item_t<TVar> > work_items;
    // TODO: Try the cache_aligned_allocator for improved performance.
    //std::vector<work_item_t<TVar>,cache_aligned_allocator<work_item_t<TVar> > > work_items;
    
    /* Pointers to shared data structures. */
    TPZAutoPointer<TPZDohrAssembly<TVar> > fAssembly;
    TPZAutoPointer<TPZCompMesh> fMesh;
    
public:
    
    parallel_assemble_task_t(TPZAutoPointer<TPZDohrAssembly<TVar> > assembly,
                             TPZAutoPointer<TPZCompMesh> mesh) :
    fAssembly(assembly), fMesh(mesh) {}
    
    /** Add a new work item to be list. */
    void push_work_item(unsigned submesh_idx, const TPZAutoPointer<TPZDohrSubstructCondense<TVar> >& substruct)
    {
        work_items.push_back(work_item_t<TVar>(submesh_idx, substruct));
    }
    
    /** Execute work items serially. */
    void run_serial()
    {
        typename std::vector<work_item_t<TVar> >::iterator it = work_items.begin();
        typename std::vector<work_item_t<TVar> >::iterator end = work_items.end();
        
        for (;it != end; it++)
        {
            work_item_t<TVar>& wi = *it;
            TPZSubCompMesh* submesh = SubMesh(fMesh, wi.fSubMeshIndex);
            ::AssembleMatrices(submesh, wi.fSubstruct, fAssembly,NULL);
            ::DecomposeBig(wi.fSubstruct, -2 /* Do not realloc */);
            ::DecomposeInternal(wi.fSubstruct, -2 /* Do not realloc */);
        }
    }
    
#ifdef USING_TBB
    /** Computing operator for the parallel for. */
    void operator()(const blocked_range<size_t>& range) const
    {
        for(size_t i=range.begin(); i!=range.end(); ++i ) {
            const work_item_t<TVar>& wi = work_items[i];
            TPZSubCompMesh* submesh = SubMesh(fMesh, wi.fSubMeshIndex);
            ::AssembleMatrices(submesh, wi.fSubstruct, fAssembly,NULL);
            ::DecomposeBig(wi.fSubstruct,-2 /* Do not realloc */);
            ::DecomposeInternal(wi.fSubstruct,-2 /* Do not realloc */);
        }
    }
    
    /** Execute work items in parallel. */
    void run_parallel_for()
    {
        /* TBB Parallel for. It will split the range into N sub-ranges and
         invoke the operator() for each sub-range.*/
        parallel_for(blocked_range<size_t>(0,work_items.size(), 1 /*IdealGrainSize*/), *this);
    }
#endif
    
};

template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::AssembleTBB(TPZBaseMatrix & mat, TPZBaseMatrix & rhs_base,
                                      TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    if(!dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" incompatible type. Aborting...\n";
        DebugStop();
    }
    auto& rhs = dynamic_cast<TPZFMatrix<TVar>&>(rhs_base);
    TPZBaseMatrix *dohrgeneric = &mat;
    TPZDohrMatrix<TVar,TPZDohrSubstructCondense<TVar> > *dohr =
    dynamic_cast<TPZDohrMatrix<TVar,TPZDohrSubstructCondense<TVar> > *> (dohrgeneric);
    
    const std::list<TPZAutoPointer<TPZDohrSubstructCondense<TVar> > > &sublist = dohr->SubStructures();
    unsigned isub;
    unsigned nsub = NSubMesh(this->fMesh);
    auto it = sublist.begin();
    parallel_assemble_task_t<TVar> parallel_tasks(fDohrAssembly, this->fMesh);
    
    /* Initialize work items. */
    std::cout << "Assembling " << nsub << " submeshes" << std::endl;
    for (isub=0; isub<nsub ; isub++) {
        TPZSubCompMesh *submesh = SubMesh(this->fMesh, isub);
        if(!submesh) continue;
        parallel_tasks.push_work_item(isub, *it);
        it++;
    }
    
    /* Run assemble and decompose. */
#ifdef USING_TBB
    parallel_tasks.run_parallel_for();
#else
    parallel_tasks.run_serial();
#endif
    
    /* Post processing. */
    for (isub=0, it=sublist.begin(); it != sublist.end(); it++, isub++) {
        TPZFMatrix<TVar> rhsloc((*it)->fNumExternalEquations,1,0.);
        (*it)->ContributeRhs(rhsloc);
        fDohrAssembly->Assemble(isub,rhsloc,rhs);
    }
    
    dohr->Initialize();
    
    TPZDohrPrecond<TVar,TPZDohrSubstructCondense<TVar> > *precond = new TPZDohrPrecond<TVar,TPZDohrSubstructCondense<TVar> > (*dohr,fDohrAssembly);
    
    precond->Initialize();
    
    fDohrPrecond = precond;
    
    return; // dohrgeneric;
}

template<class T>
struct ThreadDohrmanAssemblyList_ThreadArgs_t
{
    ThreadDohrmanAssemblyList_ThreadArgs_t() : thread_idx(-1), list(NULL) {}
    
    /* Thread index. */
    unsigned thread_idx;
    /* Thread descriptor. */
    std::thread thread;
    /* List of items to be assembled. */
    ThreadDohrmanAssemblyList<T>* list;
};

/**
 * @brief Assemble the global system of equations into the matrix which has already been created
 */

    

/* Run statistics. */
/** Jorge comments this code because is missing a file ARGLIB.CPP. 
 This file is in PERFUTIL directory and must to be added to solve linking problems.
*/
RunStatsTable dohr_ass   ("-tpz_dohr_ass", "Raw data table statistics for TPZDohrStructMatrix::Assemble assemble (first)");
RunStatsTable dohr_dec   ("-tpz_dohr_dec", "Raw data table statistics for TPZDohrStructMatrix::Assemble decompose (second)");


template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs_base,
                                   TPZAutoPointer<TPZGuiInterface> guiInterface,
                                   unsigned numthreads_assemble, unsigned numthreads_decompose)
{
  if (!dynamic_cast<TPZFMatrix<TVar> *>(&rhs_base)) {
    PZError << __PRETTY_FUNCTION__;
    PZError << " incompatible type. Aborting...\n";
    DebugStop();
  }
  auto &rhs = dynamic_cast<TPZFMatrix<TVar> &>(rhs_base);
#ifdef PERF_ANALYSIS
    ClockTimer timer;
    TimingAnalysis ta;
#endif
    
    TPZBaseMatrix *dohrgeneric = &mat;
    TPZDohrMatrix<TVar,TPZDohrSubstructCondense<TVar> > *dohr = dynamic_cast<TPZDohrMatrix<TVar,TPZDohrSubstructCondense<TVar> > *> (dohrgeneric);
    const std::list<TPZAutoPointer<TPZDohrSubstructCondense<TVar> > > &sublist = dohr->SubStructures();
    
    int nsub = NSubMesh(this->fCompMesh); // mod this->fMesh
    auto it = sublist.begin();
    
    /* Create a list of items to assemble. */
    ThreadDohrmanAssemblyList<TVar> worklist;
    for (int isub=0; isub<nsub ; isub++) {
        TPZSubCompMesh *submesh = SubMesh(this->fCompMesh, isub); // mod this->fMesh
        if(!submesh) continue;
        ThreadDohrmanAssembly<TVar> *work = new ThreadDohrmanAssembly<TVar>(this->fCompMesh,isub,*it,fDohrAssembly); //mod this->fMesh
        worklist.Append(work);
        it++;
    }
    
    if(guiInterface){
        if(guiInterface->AmIKilled()){
            return ;//0;
        }
    }
    
    // First pass : assembling the matrices
    ThreadDohrmanAssemblyList<TVar> worklistAssemble(worklist);
    auto itwork =
    worklistAssemble.fList.begin();
    while (itwork != worklistAssemble.fList.end()) {
        (*itwork)->fTask = ThreadDohrmanAssembly<TVar>::EComputeMatrix;
        itwork++;
    }
    
    
#ifdef USING_PAPI
    float rtime, ptime, mflops;
    int64_t flpops;
    PAPI_flops ( &rtime, &ptime, &flpops, &mflops );
#endif
    
    dohr_ass.start();
    if (numthreads_assemble == 0) {
        /* Put the main thread to work on all items. */
        ThreadDohrmanAssemblyList_ThreadArgs_t<TVar> targ;
        targ.thread_idx=0;
        targ.list = &worklistAssemble;
        ThreadDohrmanAssemblyList<TVar>::ThreadWork(&targ);
    }
    else {
        /* Threads arguments. */
        std::vector<ThreadDohrmanAssemblyList_ThreadArgs_t<TVar> > args(numthreads_assemble);
        
        /* Assemble multi-threaded */
        for(unsigned itr=0; itr<numthreads_assemble; itr++)
        {
            ThreadDohrmanAssemblyList_ThreadArgs_t<TVar>* targ = &(args[itr]);
            targ->thread_idx=itr;
            targ->list = &worklistAssemble;
            targ->thread = thread(ThreadDohrmanAssemblyList<TVar>::ThreadWork, targ);
        }
        /* Sync. */
        for(unsigned itr=0; itr<numthreads_assemble; itr++)
        {
            args[itr].thread.join();
        }
    }
    dohr_ass.stop();
    
#ifdef USING_PAPI
    float ltime;
    PAPI_flops ( &ltime, &ptime, &flpops, &mflops );
    
    printf("Assemble Time: %.2f \t", ltime-rtime);
    printf("Assemble Stiffness : %.2f seconds\n", stiff_sum);
    
#endif
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        int isub = 0;
        for (it=sublist.begin(); it!=sublist.end(); it++) {
            std::stringstream sout;
            sout << "Substructure number " << isub <<std::endl;
            isub++;
           // TPZDohrSubstructCondense<TVar> *ptr = (*it).operator->();
            (*it)->fMatRed->Print("Matred",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
    }
#endif
    
    // Second  pass : decomposing
    // Philippe: this may be easier to adapt the code for NUMA.
    // Edson: TODO: Measure time again.
    ThreadDohrmanAssemblyList<TVar> worklistDecompose;
    itwork = worklist.fList.begin();
    while (itwork != worklist.fList.end()) {
        TPZAutoPointer<ThreadDohrmanAssembly<TVar> > pt1 = new ThreadDohrmanAssembly<TVar> (*itwork);
        pt1->fTask = ThreadDohrmanAssembly<TVar>::EDecomposeBig;
        worklistDecompose.Append(pt1);
        TPZAutoPointer<ThreadDohrmanAssembly<TVar> > pt2 = new ThreadDohrmanAssembly<TVar>(*itwork);
        pt2->fTask = ThreadDohrmanAssembly<TVar>::EDecomposeInternal;
        worklistDecompose.Append(pt2);
        itwork++;
    }
    
    dohr_dec.start();
    if (numthreads_decompose == 0) {
        /* Compute it sequentialy */
        ThreadDohrmanAssemblyList_ThreadArgs_t<TVar> targ;
        targ.thread_idx = 0;
        targ.list = &worklistDecompose;
        ThreadDohrmanAssemblyList<TVar>::ThreadWork(&targ);
    }
    else {
        /* Threads arguments. */
        std::vector<ThreadDohrmanAssemblyList_ThreadArgs_t<TVar> >
        args(numthreads_decompose);
        for(unsigned itr=0; itr<numthreads_decompose; itr++)
        {
            ThreadDohrmanAssemblyList_ThreadArgs_t<TVar>& targ = args[itr];
            targ.thread_idx=itr;
            targ.list = &worklistDecompose;
            targ.thread = thread(ThreadDohrmanAssemblyList<TVar>::ThreadWork,
                                 &targ);
        }
        for(unsigned itr=0; itr<numthreads_decompose; itr++)
        {
            args[itr].thread.join();
        }
    }
    dohr_dec.stop();
    
    // Post processing (TODO: check whethe it is time consuming
    int isub;
    for (it=sublist.begin(), isub=0; it != sublist.end(); it++,isub++) {
        TPZFMatrix<TVar> rhsloc((*it)->fNumExternalEquations,1,0.);
        (*it)->ContributeRhs(rhsloc);
        fDohrAssembly->Assemble(isub,rhsloc,rhs);
    }
    
    dohr->SetNumThreads(this->fNumThreads);

    dohr->Initialize();
    TPZDohrPrecond<TVar,TPZDohrSubstructCondense<TVar> > *precond = new TPZDohrPrecond<TVar,TPZDohrSubstructCondense<TVar> > (*dohr,fDohrAssembly);
    precond->Initialize();
    fDohrPrecond = precond;
    
    return; // dohrgeneric;
}

/**
 * @brief Assemble the global right hand side
 */
template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::Assemble(TPZBaseMatrix & rhs_base, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    if(!dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" incompatible type. Aborting...\n";
        DebugStop();
    }
    auto& rhs = dynamic_cast<TPZFMatrix<TVar>&>(rhs_base);
    int rows = this->fMesh->NEquations();
    rhs.Redim(rows,1);
    TPZDohrPrecond<TVar,TPZDohrSubstructCondense<TVar> > *precond = dynamic_cast<TPZDohrPrecond<TVar,TPZDohrSubstructCondense<TVar> > *>(fDohrPrecond.operator->());
    const std::list<TPZAutoPointer<TPZDohrSubstructCondense<TVar> > > &sublist = precond->Global();
    
    int nsub = NSubMesh(this->fMesh);
    auto it = sublist.begin();
    
    
    int isub;
    for (isub=0; isub<nsub ; isub++) {
        TPZSubCompMesh *submesh = SubMesh(this->fCompMesh, isub);
        if(!submesh)
        {
            DebugStop();
            continue;
        }
        TPZFStructMatrix<TVar> fullstr(submesh);
        (*it)->fLocalLoad.Zero();
        fullstr.Assemble((*it)->fLocalLoad,guiInterface);
        it++;
    }
    for (it=sublist.begin(), isub=0; it != sublist.end(); it++,isub++) {
        
        // const std::list<TPZAutoPointer<TPZDohrSubstructCondense> > &sublist
        // *it represents the substructure
        TPZFMatrix<TVar> rhsloc((*it)->fNumExternalEquations,1,0.);
        (*it)->ContributeRhs(rhsloc);
        fDohrAssembly->Assemble(isub,rhsloc,rhs);
    }
    
    
}


// identify cornernodes
template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::IdentifyCornerNodes()
{
    fCornerEqs.clear();
    TPZStack<int64_t> elementgraph,elementgraphindex;
    TPZStack<int64_t> expelementgraph,expelementgraphindex;
    std::set<int> subelindexes;
    int nelem = this->fMesh->NElements();
    int iel;
    for (iel=0; iel<nelem ; iel++) {
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *> (this->fMesh->ElementVec()[iel]);
        if (sub) {
            subelindexes.insert(iel);
        }
    }
    // Determine the eligible connect sequence numbers
    std::set<int64_t> cornerconnind;
	std::set<int> cornerconnseq;
    this->fMesh->BuildCornerConnectList(cornerconnind);
    std::set<int64_t>::iterator it;
    for (it=cornerconnind.begin(); it!=cornerconnind.end(); it++) {
        TPZConnect &c = this->fMesh->ConnectVec()[*it];
        int seqnum = c.SequenceNumber();
        cornerconnseq.insert(seqnum);
    }
    
    //    this->fCompMesh->ComputeElGraph(elementgraph,elementgraphindex);
    int nindep = this->fMesh->NIndependentConnects();
    //  int neq = fCMesh->NEquations();
    this->fMesh->ComputeElGraph(elementgraph,elementgraphindex);
    int nel = elementgraphindex.NElements()-1;
    // expand the element graph to include a ficticious internal node to all elements
    expelementgraphindex.Push(0);
    int nelprev = nel;
    
    
    int count = 0;
    for (iel=0; iel<nel; iel++) {
        int nc = elementgraphindex[iel+1]-elementgraphindex[iel];
        if (nc) {
            int index = elementgraphindex[iel];
            int ic;
            for (ic=0; ic<nc; ic++) {
                expelementgraph.Push(0);
                expelementgraph[count++] = elementgraph[index++];
            }
            expelementgraph.Push(0);
            expelementgraph[count++] = nindep;
            nindep++;
        }
        expelementgraphindex.Push(count);
    }
    
    
    int next = fExternalConnectIndexes.NElements();
    
    
    if(next)
        //	if(0)
    {
        TPZManVector<int> externalconnect(nindep,0);
        // add the external connects
        int iext;
        for (iext=0; iext<next; iext++) {
            int extindex = fExternalConnectIndexes[iext];
            int seqnum = this->fMesh->ConnectVec()[extindex].SequenceNumber();
            if (seqnum >= 0) {
                externalconnect[seqnum] = 1;
            }
        }
        nel = expelementgraphindex.NElements()-1;
        for (iel=0; iel<nel; iel++) {
            bool hasext = false;
            int firstnode = expelementgraphindex[iel];
            int lastnode = expelementgraphindex[iel+1];
            int nodeindex;
            for (nodeindex= firstnode; nodeindex < lastnode; nodeindex++) {
                int node = expelementgraph[nodeindex];
                if (externalconnect[node] ==1) {
                    hasext = true;
                    break;
                }
            }
            if (hasext) {
                for (nodeindex= firstnode; nodeindex < lastnode; nodeindex++) {
                    int node = expelementgraph[nodeindex];
                    if (externalconnect[node] ==1) {
                        expelementgraph.Push(node);
                    }
                    expelementgraph.Push(nindep++);
                }
                expelementgraphindex.Push(expelementgraph.NElements());
            }
        }
    }
    
    
    
    
    // Put a global external element on top of everything
    //	if (next) {
    if (0) {
        count = expelementgraph.NElements();
        int iext;
        for (iext=0; iext<next; iext++) {
            int extindex = fExternalConnectIndexes[iext];
            int seqnum = this->fMesh->ConnectVec()[extindex].SequenceNumber();
            if (seqnum >= 0) {
                expelementgraph.Push(0);
                expelementgraph[count++] = seqnum;
            }
        }
        expelementgraphindex.Push(count);
    }
    nel = expelementgraphindex.NElements()-1;
    //	nel = elementgraphindex.NElements()-1;
    TPZRenumbering renum(nel,nindep);
    renum.SetElementGraph(expelementgraph, expelementgraphindex);
#ifdef PZ_LOG
    {
        std::stringstream sout;
        renum.Print(expelementgraph, expelementgraphindex,"Expanded graph",sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str())
		}
    }
#endif
    //	renum.SetElementGraph(elementgraph, elementgraphindex);
    std::set<int> othercornereqs;
    renum.CornerEqs(3,nelprev,cornerconnseq,othercornereqs);
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream str;
        int nelem = this->fMesh->NElements();
        int iel;
        int sub = 0;
        for (iel=0; iel<nelem; iel++) {
            TPZCompEl *cel = this->fMesh->ElementVec()[iel];
            if (!cel) {
                continue;
            }
            str << "SubCMesh " << sub << std::endl;
            int nc = cel->NConnects();
            int ic;
            for (ic=0; ic<nc; ic++) {
                TPZConnect &con = cel->Connect(ic);
                int seqnum = con.SequenceNumber();
                if (othercornereqs.find(seqnum) != othercornereqs.end()) {
                    str << seqnum << " ";
                }
            }
            str << std::endl;
            sub++;
        }
        LOGPZ_DEBUG(logger,str.str());
    }
#endif
#ifdef PZDEBUG
    std::set<int> cornerseqnums;
#endif
    int nnodes = this->fMesh->Block().NBlocks();
    int in;
    for (in=0; in<nnodes; in++) {
        if (othercornereqs.find(in) != othercornereqs.end()) {
#ifdef PZDEBUG
            cornerseqnums.insert(in);
#endif
            int pos = this->fMesh->Block().Position(in);
            int size = this->fMesh->Block().Size(in);
            int ieq;
            for(ieq=0; ieq<size; ieq++)
            {
                this->fCornerEqs.insert(pos+ieq);
            }
            
        }
    }
#ifdef PZDEBUG
    std::cout << "Number cornereqs " << fCornerEqs.size() << std::endl;

    cornerseqnums = othercornereqs;
    std::set<int> connectindices;
    TPZStack<int> geonodeindices;
    int ncon = this->fMesh->ConnectVec().NElements();
    int ic;
    for (ic=0; ic<ncon; ic++) {
        if (cornerseqnums.find(this->fMesh->ConnectVec()[ic].SequenceNumber()) != cornerseqnums.end()) {
            connectindices.insert(ic);
        }
    }
    int el;
    int numcel = this->fMesh->NElements();
    for (el=0; el<numcel; el++) {
        TPZCompEl *cel = this->fMesh->ElementVec()[el];
        if(!cel) continue;
        TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (cel);
        if(!submesh) continue;
        int elsub;
        int nelsub = submesh->NElements();
        for (elsub=0; elsub<nelsub; elsub++) {
            TPZCompEl *cel = submesh->ElementVec()[elsub];
            if (!cel) {
                continue;
            }
            int ic;
            int nc = cel->NConnects();
            for (ic=0; ic<nc ; ic++) {
                int connectindex = cel->ConnectIndex(ic);
                int fatherindex = submesh->NodeIndex(connectindex,this->fMesh);
                if(fatherindex != -1)
                {
                    if (connectindices.find(fatherindex) != connectindices.end())
                    {
                        // good one
                        TPZGeoEl *gel = cel->Reference();
                        int ncornernodes = gel->NCornerNodes();
                        if(ic<ncornernodes)
                        {
                            int nodeindex = gel->NodeIndex(ic);
                            geonodeindices.Push(nodeindex);
                        }
                        connectindices.erase(fatherindex);
                    }
                }
            }
        }
    }
    TPZAutoPointer<TPZGeoMesh> pointgmesh = new TPZGeoMesh;
    pointgmesh->NodeVec() = this->fMesh->Reference()->NodeVec();
    TPZManVector<int64_t> nodeindices(1,0);
    int ngeo = geonodeindices.NElements();
    int igeo;
    for (igeo=0; igeo<ngeo; igeo++) {
        nodeindices[0] = geonodeindices[igeo];
        int64_t index;
        pointgmesh->CreateGeoElement(EPoint,nodeindices,1,index);
    }
    pointgmesh->BuildConnectivity();
    std::ofstream arquivo("PointMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(pointgmesh.operator->(),arquivo,true);
#endif
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream str;
        str << "number of corner equations " << fCornerEqs.size() << std::endl;
        int count = 0;
        str << " corner equations ";
        std::set<int>::iterator it;
        for(it=fCornerEqs.begin(); it!=fCornerEqs.end(); it++)
        {
            str << *it << " ";
            count++;
            if (!(count%100)) {
                str << std::endl;
            }
        }
        str << std::endl;
        
        count = 0;
        
        str << "\nnumber of corner block indices after " << othercornereqs.size() << std::endl;
        for(it=othercornereqs.begin(); it!=othercornereqs.end(); it++)
        {
            str << *it << " ";
            count++;
            if (!(count%100)) {
                str << std::endl;
            }
            
        }
        LOGPZ_DEBUG(logger,str.str());
    }
#endif
}

// get the global equation numbers of a substructure (and their inverse)
template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::IdentifyEqNumbers(TPZSubCompMesh *sub, std::map<int,int> &global, std::map<int,int> &globinv)
{
    int64_t ncon = sub->ConnectVec().NElements();
    // ncon is the number of connects of the subcompmesh
    TPZCompMesh *super = this->fMesh;
    int64_t ic;
#ifdef PZ_LOG_STOP
    std::stringstream sout;
    sout << "total submesh connects/glob/loc ";
#endif
    for(ic=0; ic<ncon; ic++)
    {
        int64_t glob = sub->NodeIndex(ic,super);
        // continue is the connect is internal
        if(glob == -1) continue;
        int64_t locseq = sub->ConnectVec()[ic].SequenceNumber();
        int64_t globseq = super->ConnectVec()[glob].SequenceNumber();
        int64_t locpos = sub->Block().Position(locseq);
        int64_t globpos = super->Block().Position(globseq);
        int locsize = sub->Block().Size(locseq);
        //    int globsize = super->Block().Size(globseq);
        int ieq;
        for(ieq =0; ieq<locsize; ieq++)
        {
#ifdef PZ_LOG_STOP
            sout << ic << "/" << globpos+ieq << "/" << locpos+ieq << " ";
#endif
            global[locpos+ieq] = globpos+ieq;
            globinv[globpos+ieq] = locpos+ieq;
        }
    }
#ifdef PZ_LOG_STOP
	if (logger.isDebugEnabled())
	{
		LOGPZ_DEBUG(logger, sout.str())
	}
#endif
}

// return the number of submeshes
int64_t NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh)
{
    int64_t nel = compmesh->NElements();
    TPZCompEl *cel;
    int64_t iel, count = 0;
    for(iel=0; iel<nel; iel++)
    {
        cel = compmesh->ElementVec()[iel];
        if(!cel) continue;
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub) count++;
    }
    return count;
}

// return a pointer to the isub submesh
TPZSubCompMesh *SubMesh(TPZAutoPointer<TPZCompMesh> compmesh, int isub)
{
    int64_t nel = compmesh->NElements();
    TPZCompEl *cel;
    int64_t iel, count = 0;
    for(iel=0; iel<nel; iel++)
    {
        cel = compmesh->ElementVec()[iel];
        if(!cel) continue;
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub && isub == count) return sub;
        if(sub) count++;
    }
    return NULL;
}

// computes the permutation vectors from the subcompmesh ordening to the "internal first" ordering
// the mesh is modified during this method but is returned to its original state at the end of execution
template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::ComputeInternalEquationPermutation(TPZSubCompMesh *sub,
                                                             TPZVec<int> &scatterpermute, TPZVec<int> &gatherpermute)
{
    // This permutation vector is with respect to the blocks of the mesh
    TPZVec<int64_t> scatterpermuteblock;
    sub->ComputePermutationInternalFirst(scatterpermuteblock);
    TPZBlock destblock = sub->Block();
    TPZBlock &origblock = sub->Block();
    int64_t nblocks = origblock.NBlocks();
    if(scatterpermuteblock.NElements() != origblock.NBlocks())
    {
        std::cout << __PRETTY_FUNCTION__ << " something seriously wrong!!!\n";
    }
    int64_t ib;
    for(ib=0; ib<nblocks; ib++)
    {
        destblock.Set(scatterpermuteblock[ib],origblock.Size(ib));
    }
    destblock.Resequence();
    
    int64_t neq = ((TPZCompMesh *)sub)->NEquations();
    scatterpermute.Resize(neq);
    gatherpermute.Resize(neq);
    scatterpermute.Fill(-1);
    gatherpermute.Fill(-1);
    int64_t ncon = sub->ConnectVec().NElements();
#ifdef PZ_LOG_STOP
    std::stringstream sout;
    sout << "internal submesh connects/glob/loc ";
#endif
    int64_t ic;
    for(ic=0; ic<ncon; ic++)
    {
        // skip dependent connects
        TPZConnect &con = sub->ConnectVec()[ic];
        if(con.HasDependency() || con.IsCondensed() ) continue;
        int64_t locseq = sub->ConnectVec()[ic].SequenceNumber();
        // skip unused connects
        if(locseq < 0) continue;
        int destseq = scatterpermuteblock[locseq];
        int64_t locpos = origblock.Position(locseq);
        int64_t destpos = destblock.Position(destseq);
        int size = origblock.Size(locseq);
        //    int globsize = super->Block().Size(globseq);
        int ieq;
        for(ieq =0; ieq<size; ieq++)
        {
#ifdef PZ_LOG_STOP
            sout << ic << "/" << locpos+ieq << "/" << destpos+ieq << " ";
#endif
            scatterpermute[locpos+ieq] = destpos+ieq;
        }
    }
    int64_t ieq;
    for(ieq = 0; ieq < neq; ieq++)
    {
        gatherpermute[scatterpermute[ieq]] = ieq;
    }
#ifdef PZ_LOG_STOP
	if (logger.isDebugEnabled())
	{
		LOGPZ_DEBUG(logger, sout.str())
	}
#endif
    
}

// Identify the corner equations associated with a substructure
template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::IdentifySubCornerEqs(std::map<int,int> &globaltolocal, TPZVec<int> &cornereqs,
                                               TPZVec<int> &coarseindex)
{
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Input data for IdentifySubCornerEqs \nglobaltolocal";
        std::map<int,int>::iterator mapit;
        for(mapit = globaltolocal.begin(); mapit != globaltolocal.end(); mapit++)
        {
            sout << " [" << mapit->first << " , " << mapit->second << "] ";
        }
        sout << "\nCorner equations stored in the GenSubStructure data ";
        std::set<int>::iterator setit;
        for(setit = fCornerEqs.begin(); setit != fCornerEqs.end(); setit++)
        {
            sout << *setit << " , ";
        }
        sout << "\ncornereqs " << cornereqs;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    
    cornereqs.Resize(fCornerEqs.size());
    coarseindex.Resize(fCornerEqs.size());
    std::set<int>::iterator it;
    int64_t count = 0;
    int64_t localcount = 0;
    for(it = fCornerEqs.begin(); it!= fCornerEqs.end(); it++,count++)
    {
        if(globaltolocal.find(*it) != globaltolocal.end())
        {
            cornereqs[localcount] = globaltolocal[*it];
            coarseindex[localcount] = count;
            localcount++;
        }
    }
    cornereqs.Resize(localcount);
    coarseindex.Resize(localcount);
}


// partition the mesh in submeshes
template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::SubStructure(int nsub )
{
    
    int64_t nel = this->fMesh->NElements();
    int meshdim = this->fMesh->Dimension();
    int64_t nnodes = this->fMesh->NIndependentConnects();
    
    TPZMetis metis(nel,nnodes);
    TPZStack<int64_t> elgraph,elgraphindex;
    this->fMesh->ComputeElGraph(elgraph,elgraphindex);
    metis.SetElementGraph(elgraph, elgraphindex);
    TPZManVector<int> domain_index(nel,-1);
    metis.Subdivide(nsub, domain_index);
    CorrectNeighbourDomainIndex(this->fMesh, domain_index);
#ifdef PZDEBUG
    {
        TPZGeoMesh *gmesh = this->fMesh->Reference();
        int64_t nelgeo = gmesh->NElements();
        TPZVec<int> domaincolor(nelgeo,-999);
        int64_t cel;
        for (cel=0; cel<nel; cel++) {
            TPZCompEl *compel = this->fMesh->ElementVec()[cel];
            if(!compel) continue;
            TPZGeoEl *gel = compel->Reference();
            if (!gel) {
                continue;
            }
            domaincolor[gel->Index()] = domain_index[cel];
        }
        std::ofstream vtkfile("partitionbefore.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, domaincolor);
    }
#endif
    if(meshdim == 3)
    {
        int nsubnew = 0;
        while (nsubnew != nsub)
        {
            nsubnew = SeparateUnconnected(domain_index,nsub,meshdim-1);
            nsub = nsubnew;
        }
        nsub = ClusterIslands(domain_index,nsub,meshdim-1);
    }
    CorrectNeighbourDomainIndex(this->fMesh, domain_index);
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Geometric mesh and domain indices\n";
        this->fMesh->Reference()->Print(sout);
        sout << "Domain indices : \n";
        int64_t nel = this->fMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            sout << "el " << el << " domain " << domain_index[el] << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    
#ifdef PZDEBUG
    {
        TPZGeoMesh *gmesh = this->fMesh->Reference();
        int64_t nelgeo = gmesh->NElements();
        TPZVec<int> domaincolor(nelgeo,-999);
        int64_t cel;
        for (cel=0; cel<nel; cel++) {
            TPZCompEl *compel = this->fMesh->ElementVec()[cel];
            if(!compel) continue;
            TPZGeoEl *gel = compel->Reference();
            if (!gel) {
                continue;
            }
            domaincolor[gel->Index()] = domain_index[cel];
        }
        std::ofstream vtkfile("partition.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, domaincolor);
    }
#endif
    int isub;
    TPZManVector<TPZSubCompMesh *> submeshes(nsub,0);
    for (isub=0; isub<nsub; isub++) {
        int64_t index;
#ifdef PZDEBUG
        std::cout << '^'; std::cout.flush();
#endif
        submeshes[isub] = new TPZSubCompMesh(*(this->fMesh));
        index = submeshes[isub]->Index();
        if (index < domain_index.NElements()) {
            domain_index[index] = -1;
        }
    }
    int64_t iel;
    for (iel=0; iel<nel; iel++) {
        int domindex = domain_index[iel];
        if (domindex >= 0) {
            TPZCompEl *cel = this->fMesh->ElementVec()[iel];
            if (!cel) {
                continue;
            }
            submeshes[domindex]->TransferElement(this->fMesh,iel);
        }
    }
    for (isub = 0; isub<nsub; isub++) {
        int64_t nel = submeshes[isub]->NElements();
        if (nel == 0) {
            delete submeshes[isub];
            submeshes[isub] = 0;
        }
    }
    this->fMesh->ComputeNodElCon();
    for (isub=0; isub<nsub; isub++) {
        if (submeshes[isub])
        {
            submeshes[isub]->MakeAllInternal();
            submeshes[isub]->PermuteExternalConnects();
#ifdef PZDEBUG
            std::cout << '*'; std::cout.flush();
#endif
        }
    }
    
    this->fMesh->ComputeNodElCon();
    this->fMesh->CleanUpUnconnectedNodes();
}

// This is a lengthy process which should run on the remote processor assembling all
template<class TVar>
void AssembleMatrices(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense<TVar> > substruct, TPZAutoPointer<TPZDohrAssembly<TVar> > dohrassembly,
                      mutex* TestThread)
{
    //	static std::set<int> subindexes;
    //	int index = submesh->Index();
    //	if (subindexes.find(index) != subindexes.end()) {
    //		DebugStop();
    //	}
    //	subindexes.insert(index);
    
    
    {
        typedef typename TPZDohrSubstructCondense<TVar>::ENumbering ENumbering;
        typedef std::pair<ENumbering,ENumbering> pairnumbering;
        pairnumbering fromsub(TPZDohrSubstructCondense<TVar>::Submesh,TPZDohrSubstructCondense<TVar>::InternalFirst);
        TPZVec<int> &permutescatter = substruct->fPermutationsScatter[fromsub];
        // create a skyline matrix based on the current numbering of the mesh
        // put the stiffness matrix in a TPZMatRed object to facilitate the computation of phi and zi
        TPZSkylineStructMatrix<TVar> skylstr(submesh);
        skylstr.EquationFilter().Reset();
        
        
        TPZAutoPointer<TPZMatrix<TVar> > Stiffness = skylstr.Create();
        
        
        TPZMatRed<TVar, TPZFMatrix<TVar> > *matredbig = new TPZMatRed<TVar,TPZFMatrix<TVar> >(Stiffness->Rows()+substruct->fCoarseNodes.NElements(),Stiffness->Rows());
        
        
        matredbig->SetK00(Stiffness);
        substruct->fMatRedComplete = matredbig;
        
        
        
        TPZVec<int64_t> permuteconnectscatter;
        
        substruct->fNumInternalEquations = submesh->NumInternalEquations();
        
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "SubMesh Index = " << submesh->Index() << " Before permutation sequence numbers ";
            int64_t i;
            int64_t ncon = submesh->ConnectVec().NElements();
            for (i=0; i<ncon; i++) {
                sout << i << '|' << submesh->ConnectVec()[i].SequenceNumber() << " ";
            }
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        // change the sequencing of the connects of the mesh, putting the internal connects first
        submesh->PermuteInternalFirst(permuteconnectscatter);
        

#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "SubMesh Index = " << submesh->Index() << " After permutation sequence numbers ";
            int64_t i;
            int64_t ncon = submesh->ConnectVec().NElements();
            for (i=0; i<ncon; i++) {
                sout << i << '|' << submesh->ConnectVec()[i].SequenceNumber() << " ";
            }
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "SubMesh Index = " << submesh->Index() << "\nComputed scatter vector ";
            sout << permuteconnectscatter;
            sout << "\nStored scatter vector " << permutescatter;
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        
        // create a "substructure matrix" based on the submesh using a skyline matrix structure as the internal matrix
        TPZMatRedStructMatrix<TPZSkylineStructMatrix<TVar>,TPZVerySparseMatrix<TVar> > redstruct(submesh);
        TPZMatRed<TVar, TPZVerySparseMatrix<TVar> > *matredptr = dynamic_cast<TPZMatRed<TVar, TPZVerySparseMatrix<TVar> > *>(redstruct.Create());
        //TPZAutoPointer<TPZMatRed<TPZVerySparseMatrix> > matred = matredptr;
        
        // create a structural matrix which will assemble both stiffnesses simultaneously
        // permutescatter will reorder the equations to internal first
        TPZPairStructMatrix pairstructmatrix(submesh,permutescatter);
        
        // reorder the sequence numbering of the connects to reflect the original ordering
        TPZVec<int64_t> invpermuteconnectscatter(permuteconnectscatter.NElements());
        int64_t iel;
        for (iel=0; iel < permuteconnectscatter.NElements(); iel++) {
            invpermuteconnectscatter[permuteconnectscatter[iel]] = iel;
        }
        TPZAutoPointer<TPZMatrix<TVar> > InternalStiffness = matredptr->K00();
        
#ifdef PZDEBUG
        std::stringstream filename;
        filename << "SubMatrixInternal" << submesh->Index() << ".vtk";
        TPZFMatrix<REAL> fillin(50,50);
        submesh->ComputeFillIn(50, fillin);
        VisualMatrix(fillin, filename.str().c_str());
#endif
        
        // put the equation back in the optimized ordering for all equations (original ordering)
        submesh->Permute(invpermuteconnectscatter);
        
        
        

        // compute both stiffness matrices simultaneously
        substruct->fLocalLoad.Redim(Stiffness->Rows(),1);
#ifdef USING_PAPI
        float rtime, ptime, mflops, ltime;
        int64_t flpops;
        
        PAPI_flops ( &rtime, &ptime, &flpops, &mflops );
#endif
        pairstructmatrix.Assemble(Stiffness.operator->(), matredptr, substruct->fLocalLoad);
#ifdef USING_PAPI
        PAPI_flops ( &ltime, &ptime, &flpops, &mflops );
        //printf("Stiff: %.2f \t", ltime-rtime);
        
        stiff_sum += ltime-rtime;
#endif
        // fLocalLoad is in the original ordering of the submesh
        matredbig->SimetrizeMatRed();
        matredptr->SimetrizeMatRed();
        
        substruct->fWeights.Resize(Stiffness->Rows());
        int64_t i;
        for(i=0; i<substruct->fWeights.NElements(); i++)
        {
            substruct->fWeights[i] = Stiffness->GetVal(i,i);
        }
        // Desingularize the matrix without affecting the solution
        int64_t ncoarse = substruct->fCoarseNodes.NElements(), ic;
        int64_t neq = Stiffness->Rows();
        for(ic=0; ic<ncoarse; ic++)
        {
            int coarse = substruct->fCoarseNodes[ic];
            Stiffness->operator()(coarse,coarse) += 10.;
            //Philippe 7/6/2012
            //matredbig->operator()(coarse,coarse) += 10.;
            matredbig->operator()(neq+ic,coarse) = 1.;
            matredbig->operator()(coarse,neq+ic) = 1.;
        }
        //substruct->fStiffness = Stiffness;
        TPZStepSolver<TVar> *InvertedStiffness = new TPZStepSolver<TVar>(Stiffness);
        InvertedStiffness->SetMatrix(Stiffness);
        
        //EBORIN: Uncomment the following line to replace Cholesky by LDLt decomposition
        //#ifdef USE_LDLT_DECOMPOSITION
        
#ifdef USE_LDLT_DECOMPOSITION
        InvertedStiffness->SetDirect(ELDLt);
#else
        InvertedStiffness->SetDirect(ECholesky);
#endif
        matredbig->SetSolver(InvertedStiffness);
        
        
        TPZStepSolver<TVar> *InvertedInternalStiffness = new TPZStepSolver<TVar>(InternalStiffness);
        InvertedInternalStiffness->SetMatrix(InternalStiffness);
#ifdef DUMP_LDLT_MATRICES
        InvertedInternalStiffness->SetDirect(ELDLt);
#else
        InvertedInternalStiffness->SetDirect(ECholesky);
#endif
        matredptr->SetSolver(InvertedInternalStiffness);
        matredptr->SetReduced();
        TPZMatRed<TVar,TPZFMatrix<TVar> > *matredfull = new TPZMatRed<TVar,TPZFMatrix<TVar> >(*matredptr);
        substruct->fMatRed = matredfull;
        
        
    }
}

#ifdef DUMP_LDLT_MATRICES

#include <TPZBFileStream.h>
std::mutex dump_matrix_mutex;
unsigned matrix_unique_id = 0;

void dump_matrix(TPZAutoPointer<TPZMatrix<TVar> > Stiffness)
{
    std::scoped_lock<std::mutex> lock(dump_matrix_mutex);
    std::cout << "Dump stiffness matrix at DecomposeBig..." << std::endl;
    std::stringstream fname;
    fname << "matrix_" << matrix_unique_id++ << ".bin";
    TPZBFileStream fs;
    fs.OpenWrite(fname.str());
    Stiffness->Write(fs, 0);
    std::cout << "Dump stiffness matrix at DecomposeBig... [Done]" << std::endl;
}

#endif

//EBORIN: consumes tasks from the ThreadDohrmanAssemblyList list. The tasks
//        are ThreadDohrmanAssembly::AssembleMatrices

#ifdef USING_LIBNUMA
#include<numa.h>
class NUMA_mng_t {
    
public:
    NUMA_mng_t() {
        max_node_id = numa_max_node();
        next_rr_node = 0;
    }
    /** Return the number of nodes on the system.
     *  Nodes are identified from 0 to num_nodes-1. */
    unsigned get_num_nodes() {return (max_node_id+1);}
    /** Return the next node on a round-robin fashion. */
    unsigned get_rr_node_id() {return (next_rr_node++);}
    
private:
    
    unsigned max_node_id;
    /** Next round-robin node. */
    unsigned next_rr_node;
};

NUMA_mng_t NUMA;
#endif

clarg::argBool naa("-naDALora", "NUMA aware Dohrman Assembly List thread work objects re-allocation.", false);
clarg::argInt  naat("-naDALorat", "NUMA aware Dohrman Assembly List thread work objects re-allocation threshold.", 0);

#ifdef USING_LIBNUMA
clarg::argBool nats("-naDALtws", "NUMA aware (node round-robin) Dohrman Assembly List thread work scheduling.", false);
#endif

template<class TVar>
void DecomposeBig(TPZAutoPointer<TPZDohrSubstructCondense<TVar> > substruct, int numa_node)
{
    TPZAutoPointer<TPZMatRed<TVar,TPZFMatrix<TVar> > > matredbig = substruct->fMatRedComplete;
    TPZAutoPointer<TPZMatrix<TVar> > Stiffness = matredbig->K00();

    if (Stiffness->MemoryFootprint() > naat.get_value()) {
      Stiffness.ReallocForNuma(numa_node);
    }
    
#ifdef USE_LDLT_DECOMPOSITION
    Stiffness->Decompose(ELDLt);
#else
    Stiffness->Decompose(ECholesky);
#endif
    
    substruct->Initialize();
}

template<class TVar>
void DecomposeInternal(TPZAutoPointer<TPZDohrSubstructCondense<TVar> > substruct, int numa_node)
{
    TPZAutoPointer<TPZMatRed<TVar,TPZFMatrix<TVar> > > matred = substruct->fMatRed;
    TPZAutoPointer<TPZMatrix<TVar> > InternalStiffness = matred->K00();
    
    if (InternalStiffness->MemoryFootprint() > naat.get_value()) {
      InternalStiffness.ReallocForNuma(numa_node);
    }
    
#ifdef USE_LDLT_DECOMPOSITION
    InternalStiffness->Decompose(ELDLt);
#else
    InternalStiffness->Decompose(ECholesky);
#endif
}

//EComputeMatrix, EDecomposeInternal, EDecomposeBig
template<class TVar>
void ThreadDohrmanAssembly<TVar>::AssembleMatrices(mutex &threadtest, int numa_node)
{
    ThreadDohrmanAssembly *threadData = this;
    TPZSubCompMesh *submesh = SubMesh(threadData->fMesh,threadData->fSubMeshIndex);
    switch (fTask) {
        case EComputeMatrix:
            ::AssembleMatrices(submesh,threadData->fSubstruct,threadData->fAssembly,&threadtest);
            break;
        case EDecomposeInternal:
            DecomposeInternal(threadData->fSubstruct, numa_node);
            break;
        case EDecomposeBig:
            DecomposeBig(threadData->fSubstruct, numa_node);
            break;
        default:
            DebugStop();
            break;
    }
#ifdef PZ_LOG
    if (fTask == EComputeMatrix)
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            /*      sout << "Submesh for element " << iel << std::endl;
             submesh->Print(sout);*/
            sout << "Substructure for submesh " << fSubMeshIndex << std::endl;
            fSubstruct->Print(sout);
            LOGPZ_DEBUG(loggerasm,sout.str())
        }
#endif
    
}

template<class TVar>
ThreadDohrmanAssemblyList<TVar>::ThreadDohrmanAssemblyList() : fAccessElement(), fTestThreads()
{
}

template<class TVar>
ThreadDohrmanAssemblyList<TVar>::ThreadDohrmanAssemblyList(ThreadDohrmanAssemblyList<TVar> &cpy) : fList(cpy.fList), fAccessElement(), fTestThreads()
{
}

template<class TVar>
ThreadDohrmanAssemblyList<TVar>::~ThreadDohrmanAssemblyList()
{
}

template<class TVar>
void ThreadDohrmanAssemblyList<TVar>::Append(TPZAutoPointer<ThreadDohrmanAssembly<TVar> > object)
{
    unique_lock<mutex> lock(fAccessElement);
    fList.push_back(object);
}

template<class TVar>
TPZAutoPointer<ThreadDohrmanAssembly<TVar> > ThreadDohrmanAssemblyList<TVar>::NextObject()
{
    TPZAutoPointer<ThreadDohrmanAssembly<TVar> > result;
    unique_lock<mutex> lock(fAccessElement);
    if (fList.begin() != fList.end()) {
        result = *fList.begin();
        fList.pop_front();
    }
    return result;
}

template<class TVar>
void *ThreadDohrmanAssemblyList<TVar>::ThreadWork(void *voidptr)
{
    ThreadDohrmanAssemblyList_ThreadArgs_t<TVar>* args =
    (ThreadDohrmanAssemblyList_ThreadArgs_t<TVar>*) (voidptr);
    
    /* bind thread and newlly allocated memory to node if -naDALtws is set. */
    int node_id = -2 /* Do not realloc */;
    
#ifdef USING_LIBNUMA
    if (nats.was_set()) {
        struct bitmask* nodemask = numa_allocate_nodemask();
        numa_bitmask_clearall(nodemask);
        numa_bitmask_setbit(nodemask,args->thread_idx%NUMA.get_num_nodes());
        numa_bind(nodemask);
        numa_free_nodemask(nodemask);
    }
    if (naa.was_set()) {
        node_id = args->thread_idx%NUMA.get_num_nodes();
    }
#else
    if (naa.was_set()) {
        node_id = -1; /* Realloc */
    }
#endif
    
    TPZAutoPointer<ThreadDohrmanAssembly<TVar> > runner = args->list->NextObject();
    
    while (runner) {
        runner->AssembleMatrices(args->list->fTestThreads,node_id);
        runner = args->list->NextObject();
    }
    
    return 0;
}

// Identify the external connects
template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::IdentifyExternalConnectIndexes()
{
    // for each computational element
    std::set<int64_t> connectindexes;
    int64_t iel;
    int64_t nel = this->fMesh->NElements();
    for (iel=0; iel<nel; iel++) {
        // if it has a neighbour along its interior, skip
        TPZCompEl *cel = this->fMesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        int is;
        int ns = gel->NSides();
        int dim = gel->Dimension();
        TPZStack<TPZCompElSide> compneigh;
        
        // if there is a neighbour along the side of dimension dim skip
        TPZGeoElSide gelside(gel,ns-1);
        gelside.ConnectedCompElementList(compneigh,0,0);
        if (compneigh.NElements()) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) {
            continue;
        }
        // loop over the sides of dimension dim-1
        for (is=0; is<ns; is++)
        {
            // if there is a neighbour of dimension >= dim skip
            // the side connects are external
            TPZGeoElSide gelside(gel,is);
            if (gelside.Dimension() != dim-1) {
                continue;
            }
            compneigh.Resize(0);
            gelside.ConnectedCompElementList(compneigh, 0, 0);
            int64_t ncomp = compneigh.NElements();
            int64_t ic;
            for (ic=0; ic<ncomp; ic++) {
                TPZCompElSide celside = compneigh[ic];
                TPZGeoElSide gelside = celside.Reference();
                if (gelside.Element()->Dimension() == dim) {
                    break;
                }
            }
            // if no neighbour has dimension dim
            if (ic == ncomp) {
                int nsconnect = intel->NSideConnects(is);
                int isc;
                for (isc=0; isc<nsconnect; isc++) {
                    int64_t ind = intel->SideConnectIndex(isc,is);
                    connectindexes.insert(ind);
                }
            }
        }
    }
    std::set<int64_t>::iterator it;
    fExternalConnectIndexes.Resize(connectindexes.size());
    int64_t i = 0;
    for (it=connectindexes.begin(); it != connectindexes.end(); it++,i++) {
        fExternalConnectIndexes[i] = *it;
    }
}

// Verifies if the subdomains are connected by sides of connectdimension and separate them if not
// returns the new number of subdomains
template<class TVar, class TPar>
int  TPZDohrStructMatrix<TVar,TPar>::SeparateUnconnected(TPZVec<int> &domain_index, int nsub, int connectdimension)
{
    std::map<int,int> domain_index_count;
    int64_t iel;
    int64_t nel = this->fMesh->NElements();
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = this->fMesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        int mydomainindex = domain_index[cel->Index()];
        domain_index_count[mydomainindex]++;
    }
    std::set<int> domain_check;
    
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = this->fMesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        int mydomainindex = domain_index[cel->Index()];
        if (domain_check.find(mydomainindex) != domain_check.end()) {
            continue;
        }
        domain_check.insert(mydomainindex);
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        
        TPZStack<TPZGeoEl *> gelstack;
        gelstack.Push(gel);
        std::set<TPZCompEl *> gelcluster;
        while (gelstack.NElements())
        {
            TPZGeoEl *gel = gelstack.Pop();
            if (gelcluster.find(gel->Reference()) != gelcluster.end()) {
                continue;
            }
            int beforesize = gelcluster.size();
            gelcluster.insert(gel->Reference());
            int checksize = gelcluster.size();
            if (checksize == beforesize) {
                DebugStop();
            }
            
            int nsides = gel->NSides();
            int is;
            for (is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != connectdimension) {
                    continue;
                }
                TPZGeoElSide gelside(gel,is);
                TPZStack<TPZCompElSide> elsidevec;
                gelside.ConnectedCompElementList(elsidevec, 0, 0);
                int64_t nneigh = elsidevec.NElements();
                int64_t neigh;
                for (neigh = 0; neigh <nneigh; neigh++) {
                    TPZCompElSide celside = elsidevec[neigh];
                    TPZCompEl *celloc = celside.Element();
                    if (domain_index[celloc->Index()] != mydomainindex) {
                        continue;
                    }
                    if (gelcluster.find(celloc) == gelcluster.end()) {
                        gelstack.Push(celloc->Reference());
                    }
                }
            }
        }
        
        if (gelcluster.size() != (std::set<TPZCompEl *>::size_type)domain_index_count[mydomainindex]) {
            if (gelcluster.size() > (std::set<TPZCompEl *>::size_type)domain_index_count[mydomainindex]) {
                DebugStop();
            }
            domain_index_count[mydomainindex] -= gelcluster.size();
            domain_index_count[nsub] = gelcluster.size();
            std::set<TPZCompEl *>::iterator it;
            domain_check.erase(mydomainindex);
            domain_check.insert(nsub);
            for (it=gelcluster.begin(); it!=gelcluster.end(); it++) {
                domain_index[(*it)->Index()]=nsub;
            }
            nsub++;
        }
    }
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Number of elements per domain ";
        std::map<int,int>::iterator it;
        int64_t count = 0;
        for (it=domain_index_count.begin(); it != domain_index_count.end(); it++) {
            if (! (count++ %40)) {
                sout << std::endl;
            }
            sout << it->first << " " << it->second << " " << std::endl;
        }
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    return nsub;
}

// Eliminate subdomains who are embedded in other subdomains
// returns the number of subdomains
template<class TVar, class TPar>
int  TPZDohrStructMatrix<TVar,TPar>::ClusterIslands(TPZVec<int> &domain_index,int nsub,int connectdimension)
{
    int meshdim = this->fMesh->Dimension();
    int64_t nel = this->fMesh->NElements();
    int64_t mincount = nel/nsub/20;
    // contains for each subdomain the set of neighbouring domains
    TPZVec<std::set<int> > domain_neighbours(nsub);
    // contains for each domain the number of cells within that domain
    std::map<int,int> domain_index_count;
    int64_t iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = this->fMesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        int mydomainindex = domain_index[cel->Index()];
        //        if (mydomainindex == 0) {
        //            std::stringstream sout;
        //            cel->Print(sout);
        //            TPZGeoEl *gel = cel->Reference();
        //            if (gel) {
        //                gel->Print(sout);
        //            }
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        domain_index_count[mydomainindex]++;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        int nsides = gel->NSides();
        int geldim = gel->Dimension();
        int is;
        for (is=0; is<nsides; is++) {
            int sidedim = gel->SideDimension(is);
            if (sidedim != connectdimension && geldim>=connectdimension) {
                continue;
            }
            if (geldim < connectdimension && is != nsides-1)
            {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> elsidevec;
            gelside.ConnectedCompElementList(elsidevec, 0, 0);
            int64_t nneigh = elsidevec.NElements();
            int64_t neigh;
            int64_t nneighvalid = 0;
            for (neigh = 0; neigh <nneigh; neigh++) {
                TPZCompElSide celside = elsidevec[neigh];
                TPZCompEl *celloc = celside.Element();
                TPZGeoEl *gelloc = celloc->Reference();
                if (gelloc->Dimension() != meshdim) {
                    continue;
                }
                nneighvalid++;
                int celdomain = domain_index[celloc->Index()];
                if (celdomain != mydomainindex)
                {
                    domain_neighbours[mydomainindex].insert(celdomain);
                }
            }
            if (nneighvalid == 0)
            {
                // include the boundary as a ficticious neighbour index
                domain_neighbours[mydomainindex].insert(-1);
            }
        }
    }
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        for (int64_t i=0; i<domain_neighbours.size(); i++) {
            std::set<int>::const_iterator it;
            sout << "Domain index " << i << " neighbours ";
            for (it=domain_neighbours[i].begin(); it != domain_neighbours[i].end(); it++) {
                sout << *it << " ";
            }
            sout << std::endl;
        }
        std::map<int,int>::const_iterator it = domain_index_count.begin();
        while (it != domain_index_count.end()) {
            sout << "Domain index " << it->first << " number of elements " << it->second << std::endl;
            it++;
        }
        LOGPZ_DEBUG(logger, sout.str())
        
    }
#endif
    // compute a destination domain index for each domain (used for clustering domains)
    int isub;
    TPZManVector<int> domain_dest(nsub,-1);
    int64_t count = 0;
    for (isub=0; isub < nsub; isub++)
    {
        // if the subdomain is neighbour to only one subdomain
        // this means that the subdomain is isolated (only boundaries as neighbours) (not treated)
        // or that the domain is embedded in another domain
        if (domain_neighbours[isub].size() == 1 )
        {
            // merge both subdomains
            int target = *(domain_neighbours[isub].begin());
            // target == -1 is not treated here
            if (target == -1) {
                continue;
            }
            if (domain_dest[target] == -1 && domain_dest[isub] == -1)
            {
                domain_dest[isub] = count;
                domain_dest[target] = count;
                count++;
            }
            else if (domain_dest[target] == -1)
            {
                domain_dest[target] = domain_dest[isub];
            }
            else
            {
                domain_dest[isub] = domain_dest[target];
            }
            
        }
        else if(domain_dest[isub] == -1 && domain_index_count[isub] < mincount)
        {
            // the domain has very little elements
            // the domain has at least two neighbouring domains (may include the ficticious -1 domain)
            std::map<int,int> sizeDomain;
            std::set<int>::iterator it;
            for (it = domain_neighbours[isub].begin(); it != domain_neighbours[isub].end(); it++) {
                if (*it != -1) {
                    sizeDomain[domain_index_count[isub]] = *it;
                }
            }
            int domaintargetindex = sizeDomain.rbegin()->second;
            int destdomainindexcount = domain_index_count[domaintargetindex];
            int domainshrinkcount = domain_index_count[isub];
            domain_index_count[domaintargetindex] = destdomainindexcount+domainshrinkcount;
            domain_index_count[isub] = 0;
            if(domain_dest[domaintargetindex] == -1)
            {
                domain_dest[domaintargetindex] = count;
                domain_dest[isub] = count;
                count++;
            }
            else {
                domain_dest[isub] = domain_dest[domaintargetindex];
            }
            
        }
        else if (domain_dest[isub] == -1)
        {
            domain_dest[isub] = count++;
        }
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        int isub;
        for (isub=0; isub < nsub; isub++) {
            sout << "isub = " << isub << " number of neighbours " << domain_neighbours[isub].size() << " domains ";
            std::set<int>::iterator it;
            for (it = domain_neighbours[isub].begin(); it != domain_neighbours[isub].end(); it++) {
                sout << *it << " ";
            }
            sout << std::endl;
        }
        sout << "Destination domain " << domain_dest << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    int domsize = domain_index.NElements();
    int d;
    for (d=0; d<domsize; d++) {
        domain_index[d] = domain_dest[domain_index[d]];
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Number of elements per domain ";
        std::map<int,int>::iterator it;
        int64_t count = 0;
        for (it=domain_index_count.begin(); it != domain_index_count.end(); it++) {
            if (! (count++ %40)) {
                sout << std::endl;
            }
            sout << it->first << " " << it->second << " ";
        }
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    return count;
}


template<class TVar, class TPar>
int  TPZDohrStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZDohrStructMatrix") ^ TPZStructMatrix::ClassId() << 1;
}

template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::Write( TPZStream &str, int withclassid ) const
{
    TPZPersistenceManager::WritePointer(this->fMesh, &str);
    int hasdohrassembly = 0;
    if (fDohrAssembly) {
        hasdohrassembly = 1;
    }
    str.Write(&hasdohrassembly);
    if (hasdohrassembly) {
        TPZPersistenceManager::WritePointer(fDohrAssembly.operator ->(), &str);
    }
    str.Write( fExternalConnectIndexes);
    str.Write(fCornerEqs);
}

template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::Read(TPZStream &str, void *context )
{
    this->SetMesh(TPZAutoPointerDynamicCast<TPZCompMesh>(TPZPersistenceManager::GetAutoPointer(&str)));
    int hasdohrassembly;
    str.Read(&hasdohrassembly);
    if (hasdohrassembly) {
        fDohrAssembly = TPZAutoPointerDynamicCast<TPZDohrAssembly<TVar>>(TPZPersistenceManager::GetAutoPointer(&str));
    }
    str.Read( fExternalConnectIndexes);
    str.Read( fCornerEqs);
}

/** @brief Set the domain index of the lower dimension elements equal to the domain index of their neighbour */
template<class TVar, class TPar>
void  TPZDohrStructMatrix<TVar,TPar>::CorrectNeighbourDomainIndex(TPZCompMesh *cmesh, TPZVec<int> &domainindex)
{
    int64_t nel = cmesh->NElements();
    TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh->ElementVec();
    bool changed = true;
    while(changed)
    {
        changed = false;
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = elvec[el];
            if (! cel) {
                continue;
            }
            TPZGeoEl *gel = cel->Reference();
            if (!gel) {
                continue;
            }
            int nsides = gel->NSides();
            TPZGeoElSide neighbour = gel->Neighbour(nsides-1);
            if (neighbour.Element() != gel) {
                TPZCompEl *neighcel = neighbour.Element()->Reference();
                if (! neighcel) {
                    continue;
                }
                int64_t neighindex = neighcel->Index();
                if (domainindex[el] != domainindex[neighindex]) {
                    domainindex[el] = domainindex[neighindex];
                    changed = true;
                }
            }
        }
    }
}


#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZDohrStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZDohrStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZDohrStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;

#undef NOMETIS
