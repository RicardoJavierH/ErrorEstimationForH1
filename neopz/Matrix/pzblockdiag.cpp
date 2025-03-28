/**
 * @file
 * @brief Contains the implementation of the TPZBlockDiagonal methods.
 */

#include "pzfmatrix.h"
#include "pzblockdiag.h"

#include <random>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "pzlog.h"
#include <sstream>
#include "pzstack.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template<class TVar>
void TPZBlockDiagonal<TVar>::AddBlock(int64_t i, TPZFMatrix<TVar> &block){

	int64_t firstpos = fBlockPos[i];
	int64_t bsize = fBlockSize[i];
	
	int64_t r,c;
	for(r=0; r<bsize; r++) {
		for(c=0; c<bsize; c++) {
			fStorage[firstpos+r+bsize*c] += block(r,c);
		}
	}
}

template<class TVar>
void TPZBlockDiagonal<TVar>::SetBlock(int64_t i, TPZFMatrix<TVar> &block){
	int64_t firstpos = fBlockPos[i];
	int64_t bsize = fBlockSize[i];

#ifdef PZDEBUG
	if(block.Rows() != bsize || block.Cols() != bsize){
		DebugStop();
	}
#endif
	int64_t r,c;
	for(r=0; r<bsize; r++) {
		for(c=0; c<bsize; c++) {
			fStorage[firstpos+r+bsize*c] = block(r,c);
		}
	}
}

template<class TVar>
void TPZBlockDiagonal<TVar>::GetBlock(int64_t i, TPZFMatrix<TVar> &block){
    int64_t firstpos = fBlockPos[i];
    int64_t bsize = fBlockSize[i];
    block.Redim(bsize,bsize);
    int64_t r,c;
    for(r=0; r<bsize; r++) {
        for(c=0; c<bsize; c++) {
            block(r,c) = fStorage[firstpos+r+bsize*c];
        }
    }
}

template<class TVar>
TPZAutoPointer<TPZFMatrix<TVar>>
TPZBlockDiagonal<TVar>::GetBlockPtr(int64_t i){
	if(i > -1 && i < this->NumberofBlocks()){
		return fBlockMats[i];
	}else{
		return nullptr;
	}
}


template<class TVar>
void TPZBlockDiagonal<TVar>::Initialize(const TPZVec<int> &blocksize){
	int64_t nblock = blocksize.NElements();
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Number of blocks \t" << nblock;
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	fBlockSize = blocksize;
	fBlockPos.Resize(nblock+1,0);
	fBlockMats.Resize(nblock, nullptr);
	int64_t b;
	int64_t ndata = 0;
	int64_t neq = 0;
	int bsize;
	for(b=0; b<nblock; b++) {
		bsize = blocksize[b];
		fBlockPos[b+1] = fBlockPos[b]+bsize*bsize;
		//position doesnt matter if blocksize is zero,
		//and this will avoid runtime errors if zero sized blocks are at
		//the end of block list
		if(bsize==0){fBlockPos[b] = 0;}
		ndata += bsize*bsize;
		neq += bsize;
	}
#ifdef PZ_LOG
	if(ndata > 10000000)
	{
		std::stringstream sout;
		sout << "Calling fStorage.Resize(ndata,0.) with ndata = " << ndata;
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	
	fStorage.Fill(0.,0);
	fStorage.Resize(ndata,0.);
	this->fDecomposed = ENoDecompose;
	this->fRow = neq;
	this->fCol = neq;
	
	for(b=0; b<nblock; b++) {
		const auto bsize = blocksize[b];
		const auto pos = fBlockPos[b];
		if(bsize){
			fBlockMats[b] = new TPZFMatrix<TVar>(bsize,bsize, &fStorage[pos], bsize*bsize);
		}else{
			fBlockMats[b] = nullptr;
		}
	}
}

template<class TVar>
void TPZBlockDiagonal<TVar>::BuildFromMatrix(TPZMatrix<TVar> &mat) {
	if(mat.Rows() != this->Rows()) {
		cout << "TPZBlockDiagonal::BuildFromMatrix wrong data structure\n";
		return;
	}
	int64_t nblock = fBlockSize.NElements();
	int64_t b,eq=0;
	for(b=0; b<nblock; b++) {
		int r,c,bsize = fBlockSize[b];
		int64_t pos = fBlockPos[b];
		for(r=0; r<bsize; r++){
			for(c=0; c<bsize; c++) {
				fStorage[pos+r+c*bsize] = mat.GetVal(eq+r,eq+c);
			}
		}
		eq += bsize;
	}
}

/*******************/
/*** Constructor ***/

template<class TVar>
TPZBlockDiagonal<TVar>::TPZBlockDiagonal()
: TPZRegisterClassId(&TPZBlockDiagonal::ClassId),
TPZMatrix<TVar>(), fStorage(), fBlockPos(1,0),fBlockSize()
{
	
}

template<class TVar>
TPZBlockDiagonal<TVar>::TPZBlockDiagonal(const TPZVec<int> &blocksize)
: TPZRegisterClassId(&TPZBlockDiagonal::ClassId),
TPZMatrix<TVar>(), fStorage(), fBlockPos(1,0),fBlockSize()
{
	Initialize(blocksize);
}

/********************/
/*** Constructors ***/

template<class TVar>
TPZBlockDiagonal<TVar>::TPZBlockDiagonal(const TPZVec<int> &blocksizes, const TPZFMatrix<TVar> &glob)
: TPZRegisterClassId(&TPZBlockDiagonal::ClassId),
TPZMatrix<TVar>(), fBlockSize(blocksizes)
{
	Initialize(blocksizes);
	const auto nblock = blocksizes.NElements();
	int eq = 0;
	for(auto b=0; b<nblock; b++) {
		const auto bsize = fBlockSize[b];
		const auto pos = fBlockPos[b];
		for(auto r=0; r<bsize; r++) {
			for(auto c=0; c<bsize; c++) {
				fStorage[pos+r+bsize*c]= glob.GetVal(eq+r,eq+c);
			}
		}
		eq += bsize;
	}
	
}

/*********************************/
/*** Constructor( TPZBlockDiagonal& ) ***/
template<class TVar>
TPZBlockDiagonal<TVar>::TPZBlockDiagonal (const TPZBlockDiagonal<TVar> & A)
: TPZRegisterClassId(&TPZBlockDiagonal::ClassId),
TPZMatrix<TVar>( A.Dim(), A.Dim() ), fStorage(A.fStorage),
	fBlockPos(A.fBlockPos), fBlockSize(A.fBlockSize)
{
	const int nblock = fBlockSize.NElements();
	fBlockMats.Resize(nblock);
	for(int b=0; b<nblock; b++) {
		const auto bsize = fBlockSize[b];
		const auto pos = fBlockPos[b];
		fBlockMats[b] = new TPZFMatrix<TVar>(bsize,bsize, &fStorage[pos], bsize*bsize);
	}
}

/******************/
/*** Destructor ***/
template<class TVar>
TPZBlockDiagonal<TVar>::~TPZBlockDiagonal ()
{
}

/***********/
/*** Put ***/
template<class TVar>
int TPZBlockDiagonal<TVar>::Put(const int64_t row,const int64_t col,const TVar& value )
{
	//  cout << "TPZBlockDiagonal.Put should not be called\n";
	if ( (row >= Dim()) || (col >= Dim()) || row<0 || col<0)
    {
		cout << "TPZBlockDiagonal::Put: " << row << "," << col << "," << Dim();
		cout << "\n";
		return( 0 );
    }
	
	return( PutVal( row, col, value ) );
}

/***********/
/*** PutVal ***/
template<class TVar>
int TPZBlockDiagonal<TVar>::PutVal(const int64_t row,const int64_t col,const TVar& value )
{
	int64_t b = 0;
	int64_t nb = fBlockSize.NElements();
	if(nb==0) {
		cout << "TPZBlockDiagonal::PutVal called with parameters out of range\n";
		return -1;
	}
	int64_t eq=0;
	int64_t bsize = fBlockSize[b];
	while(eq+bsize <= row && b < nb) {
		eq+=bsize;
		b++;
		bsize = fBlockSize[b];
	}
	if(b==nb) {
    	cout << "TPZBlockDiagonal::PutVal wrong data structure\n";
		return -1;
	}
	if(col < eq || col >= eq+bsize) {
		if(value != TVar(0.)) {
    		cout << "TPZBlockDiagonal::PutVal, indices row col out of range\n";
			return -1;
		} else {
			return 0;
		}
	}
	fStorage[fBlockPos[b]+row-eq+bsize*(col-eq)] = value;
	return 0;
}



/***********/
/*** Get ***/

template<class TVar>
const TVar
TPZBlockDiagonal<TVar>::Get(const int64_t row,const int64_t col ) const
{
	if ( (row >= Dim()) || (col >= Dim()) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZBlockDiagonal::Get <indices out of band matrix range>" );
	
	return( GetVal( row, col ) );
}

template<class TVar>
TVar &
TPZBlockDiagonal<TVar>::operator()(const int64_t row, const int64_t col) {
	
	int64_t b = 0;
	int64_t nb = fBlockSize.NElements();
	if(nb==0) {
		cout << "TPZBlockDiagonal::operator() called with parameters out of range\n";
		static TVar zero = 0.;
		return zero;
	}
	int64_t eq=0;
	int64_t bsize = fBlockSize[b];
	while(eq+bsize <= row && b < nb) {
		eq+=bsize;
		b++;
		bsize = fBlockSize[b];
	}
	if(b==nb) {
		cout << "TPZBlockDiagonal::operator() wrong data structure\n";
		static TVar zero = 0.;
		return zero;
	}
	if(col < eq || col >= eq+bsize) {
		cout << "TPZBlockDiagonal::operator(), indices row col out of range\n";
		static TVar zero = 0.;
		return zero;
	}
	return fStorage[fBlockPos[b]+row-eq+bsize*(col-eq)];
}

/***********/
/*** GetVal ***/
template<class TVar>
const TVar TPZBlockDiagonal<TVar>::GetVal(const int64_t row,const int64_t col ) const
{	
	int64_t b = 0;
	int64_t nb = fBlockSize.NElements();
	if(nb==0) {
		cout << "TPZBlockDiagonal::GetVal called with parameters out of range\n";
	}
	int64_t eq=0;
	int64_t bsize = fBlockSize[b];
	while(eq+bsize <= row && b < nb) {
		eq+=bsize;
		b++;
		bsize = fBlockSize[b];
	}
	if(b==nb) {
		cout << "TPZBlockDiagonal::GetVal wrong data structure\n";
	}
	if(col < eq || col >= eq+bsize) {
		return (TVar)0;
	}
	return fStorage[fBlockPos[b]+row-eq+bsize*(col-eq)];
}


/******** Operacoes com MATRIZES GENERICAS ********/

/*******************/
/*** MultiplyAdd ***/
//
//  perform a multiply add operation to be used by iterative solvers
//

template<class TVar>
void TPZBlockDiagonal<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
							   const TVar alpha,const TVar beta ,const int opt) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	
	if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZBlockDiagonal::MultAdd <matrixs with incompatible dimensions>" );
	if(beta != TVar(0) && (x.Cols() != y.Cols() || x.Rows() != y.Rows())) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZBlockDiagonal::MultAdd incompatible dimensions\n");
	}
	
	this->PrepareZ(y,z,beta,opt);
	int64_t xcols = x.Cols();
	int64_t nb= fBlockSize.NElements();
	int64_t ic, b, eq=0;
	int bsize, r, c;
	if(opt == 0) {
		for (ic = 0; ic < xcols; ic++) {
			eq=0;
			for(b=0; b<nb; b++) {
				bsize = fBlockSize[b];
				int64_t pos = fBlockPos[b];
				for(r=0; r<bsize; r++) {
					for(c=0; c<bsize; c++) {
						z(eq+r,ic) += alpha*fStorage[pos+r+bsize*c]*x.GetVal((eq+c),ic);
					}
				}
				eq += bsize;
			}
		}
	} else {
		for (ic = 0; ic < xcols; ic++) {
			eq=0;
			for(b=0; b<nb; b++) {
				bsize = fBlockSize[b];
				int64_t pos = fBlockPos[b];
				for(r=0; r<bsize; r++) {
					for(c=0; c<bsize; c++) {
						z(eq+r,ic) += alpha*fStorage[pos+c+bsize*r]*x.GetVal((eq+c),ic);
					}
				}
				eq+=bsize;
			}
		}
	}
}


/***************/
/**** Zero ****/
template<class TVar>
int TPZBlockDiagonal<TVar>::Zero()
{
	
	fStorage.Fill(0.,0);
	this->fDecomposed = ENoDecompose;
	
	return( 1 );
}



/********************/
/*** Transpose () ***/
template<class TVar>
void TPZBlockDiagonal<TVar>::Transpose (TPZMatrix<TVar> *const T) const
{
	T->Resize( Dim(), Dim() );
	
	int64_t b, eq = 0, pos;
	int bsize, r, c;
	int64_t nb = fBlockSize.NElements();
	for ( b=0; b<nb; b++) {
		pos= fBlockPos[b];
		bsize = fBlockSize[b];
		for(r=0; r<bsize; r++) {
			for(c=0; c<bsize; c++) {
				T->PutVal(eq+r,eq+c,fStorage[pos+c+r*bsize]);
			}
		}
		eq += bsize;
	}
}


/*****************/
/*** Decompose_LU ***/
template<class TVar>
int TPZBlockDiagonal<TVar>::Decompose_LU()
{
#ifdef PZ_LOG
    if(logger.isDebugEnabled())	LOGPZ_DEBUG(logger, "TPZBlockDiagonal::Decompose_LU");
#endif
	
	if (  this->fDecomposed && this->fDecomposed == ELU) {
		return ELU;
	} else if(this->fDecomposed) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZBlockDiagonal::Decompose_LU is already decomposed with other scheme");
	}
	
	int64_t b,nb,pos;
	int bsize;
	nb = fBlockSize.NElements();
	for(b=0;b<nb; b++) {
		pos = fBlockPos[b];
		bsize = fBlockSize[b];
		if(!bsize) continue;
		
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream mess;
            mess << "TPZBlockDiagonal::Decompose_LU() - bsize = " << bsize << ", bsize*bsize = " << bsize*bsize;
            LOGPZ_DEBUG(logger,mess.str());
        }
#endif


		fBlockMats[b]->Decompose_LU();
	}
	this->fDecomposed = ELU;
	return 1;
}

template<class TVar>
int
TPZBlockDiagonal<TVar>::Substitution( TPZFMatrix<TVar> *B) const
{
	if(this->fDecomposed != ELU) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZBlockDiagonal::Decompose_LU is decomposed with other scheme");
	}
	
	int64_t b,nb,pos,bsize,eq=0;
	nb = fBlockSize.NElements();
	int64_t c, nc = B->Cols();
	for(c=0; c<nc; c++) {
		eq = 0;
		for(b=0;b<nb; b++) {
			bsize = fBlockSize[b];
			if(!bsize) continue;
			TPZFMatrix<TVar> BTemp(bsize,1,&(B->operator()(eq,c)),bsize);
			fBlockMats[b]->Substitution(&BTemp);
			eq+= bsize;
		}
	}
	return 1;
}

/************************** Private **************************/

/*************/
/*** Clear ***/
template<class TVar>
int TPZBlockDiagonal<TVar>::Clear()
{
	fStorage.Resize(0);
	fBlockPos.Resize(0);
	fBlockSize.Resize(0);
	fBlockMats.Resize(0);
	this->fRow = 0;
	this->fCol = 0;
	this->fDecomposed = ENoDecompose;
	return( 1 );
}

template<class TVar>
int TPZBlockDiagonal<TVar>::main() {
	
	cout << "Entering the main program\n";
	cout.flush();
	TPZFMatrix<TVar> ref(7,7,0.);
	int r,c;
	for(r=0; r<7; r++) {
		for(c=0; c<7; c++) {
			ref(r,c) = ((TVar)(float)(5+r*c+3*r));
		}
		ref(r,r) += (TVar)1000;
	}
	TPZVec<int> blocksize(3);
	blocksize[0] = 2;
	blocksize[1] = 4;
	blocksize[2] = 1;
	TPZBlockDiagonal bd1(blocksize,ref);
	TPZBlockDiagonal bd2(bd1);
	ref.Print("original matrix",std::cout);
	bd1.Solve_LU(&ref);
	bd1.Solve_LU(&ref);
	ref.Print("after inverting the diagonal",std::cout);
	TPZFMatrix<TVar> ref2;
	bd2.Multiply(ref,ref2);
	bd2.Multiply(ref2,ref);
	ref.Print("restoring the original matrix",std::cout);
	return 1;
	
}

template<class TVar>
void TPZBlockDiagonal<TVar>::Print(const char *msg, std::ostream &out, const MatrixOutputFormat format) const {
	
	if(format != EFormatted)
	{
		TPZMatrix<TVar>::Print(msg,out,format);
		return;
	}
	out << "TPZBlockDiagonal matrix ";
	if(msg) out << msg;
	out  << std::endl;
	
	int64_t nblock = fBlockSize.NElements();
	out << "Number of blocks " << nblock << std::endl; 
	int64_t b,bsize,pos;
	for(b=0; b<nblock; b++) {
		bsize = fBlockSize[b];
		out << "block number " << b << " size : " << bsize << std::endl;
		int64_t r,c;
		pos = fBlockPos[b];
        for(r=0; r<bsize ; r++) {
            for(c=0; c<bsize; c++) {
				out << fStorage[pos+r+bsize*c] << ' ';
			}
			out << std::endl;
		}
	}
}

/**
 * Updates the values of the matrix based on the values of the matrix
 */
template<class TVar>
void TPZBlockDiagonal<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)
{
  if(!mat)
  {
    cout << "TPZBlockDiagonal::UpdateFrom" << " called with zero argument\n";
    return;
  }
  TPZMatrix<TVar>& matref = mat;
  UpdateFrom(matref);
}

template<class TVar>
void TPZBlockDiagonal<TVar>::UpdateFrom(TPZMatrix<TVar>& matref)
{
  this->fDecomposed = ENoDecompose;
  int64_t nblock = fBlockSize.NElements();
  int64_t b,bsize,pos,firsteq = 0;
  for(b=0; b<nblock; b++) {
    bsize = fBlockSize[b];
    if(bsize){
      //    int r,c;
      pos = fBlockPos[b];
      TPZFMatrix<TVar> block(bsize,bsize,&fStorage[pos],bsize*bsize);
      matref.GetSub(firsteq,firsteq,bsize,bsize,block);
      firsteq += bsize;
    }
  }
}

/** Fill the matrix with random values (non singular matrix) */
template<class TVar>
void TPZBlockDiagonal<TVar>::AutoFill(int64_t neq, int64_t jeq, SymProp sp) {
    if (neq != jeq) {
        DebugStop();
    }
    TPZStack<int> blsizes;
    int64_t totalsize = 0;
    while (totalsize < neq) {
        int64_t blsize = (neq*rand())/RAND_MAX;
        blsize = blsize < neq-totalsize ? blsize : neq-totalsize;
        blsizes.Push(blsize);
        totalsize += blsize;
    }
    Initialize(blsizes);
    // Initialize the blocksizes!!
	int64_t b, bsize, eq = 0, pos;
	int64_t nb = fBlockSize.NElements(), r, c;

	const bool must_conj = is_complex<TVar>::value && sp == SymProp::Herm;
	for ( b=0; b<nb; b++) {
		pos= fBlockPos[b];
		bsize = fBlockSize[b];
		for(c=0; c<bsize; c++) {
            RTVar sum = 0.;
            r=0;
						if (must_conj){
							if constexpr (is_complex<TVar>::value){
								for (r=0; r<c; r++) {
									fStorage[pos+c+r*bsize]=std::conj(fStorage[pos+r+c*bsize]);
									sum += fabs(fStorage[pos+r+c*bsize]);
								}
							}else{
								DebugStop();//how the hell did we get here?
							}
						}else{
							for (r=0; r<c; r++) {
								fStorage[pos+c+r*bsize]=fStorage[pos+r+c*bsize];
								sum += fabs(fStorage[pos+r+c*bsize]);
							}
						}
			for(; r<bsize; r++) {
				auto val = this->GetRandomVal();
				fStorage[pos+c+r*bsize] = (val);
                if(c!= r) sum += fabs(val);
			}
            //totototo
//            if(r==4)
//            {
//                std::cout << "sum " << sum << std::endl;
//            }
            if (fabs(fStorage[pos+c+c*bsize]) < sum) {
                fStorage[pos+c+c*bsize] = (TVar)(sum + (RTVar)1.);
            }
		}
		eq += bsize;
	}
}

template<class TVar>
int TPZBlockDiagonal<TVar>::ClassId() const{
    return Hash("TPZBlockDiagonal") ^ TPZMatrix<TVar>::ClassId() << 1;
}

template class TPZBlockDiagonal<float>;
template class TPZBlockDiagonal<double>;
template class TPZBlockDiagonal<long double>;

template class TPZBlockDiagonal<std::complex<float> >;
template class TPZBlockDiagonal<std::complex<double> >;
template class TPZBlockDiagonal<std::complex<long double> >;
