
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

  //
  // MultXM
  //

  template <class T, class T1> void RowMajorMultXM(
      const T1 alpha, const UpperTriMatrixView<T>& A)
  {
    TMVAssert(A.isrm());
    TMVAssert(!A.isunit());
    TMVAssert(alpha != T1(1));
    TMVAssert(A.size() > 0);
    TMVAssert(A.ct() == NonConj);

    const int ds = A.stepi()+1;
    T* Aptr = A.ptr();
    for(size_t i=0,len=A.rowsize();len>0;++i,--len,Aptr+=ds) 
      // A.row(i,i,N) *= alpha;
      DoMultXV(alpha,VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
	  len);
  }

  template <class T, class T1> void ColMajorMultXM(
      const T1 alpha, const UpperTriMatrixView<T>& A)
  {
    TMVAssert(A.iscm());
    TMVAssert(!A.isunit());
    TMVAssert(alpha != T1(1));
    TMVAssert(A.size() > 0);
    TMVAssert(A.ct() == NonConj);

    T* Aptr = A.ptr();
    for(size_t j=0;j<A.rowsize();++j,Aptr+=A.stepj()) 
      // A.col(j,0,j+1) *= alpha;
      DoMultXV(alpha,VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
	  j+1);
  }

  template <class T> void MultXM(const T alpha, const UpperTriMatrixView<T>& A)
    // A = alpha * A
  {
#ifdef XDEBUG
    Matrix<T> A2 = alpha * Matrix<T>(A);
    Matrix<T> A0 = A;
#endif
    TMVAssert(!A.isunit());

    if (A.size() > 0 && alpha != T(1)) {
      if (alpha == T(0)) A.Zero();
      else {
	if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
	else if (A.isrm())
	  if (IMAG(alpha) == RealType(T)(0))
	    RowMajorMultXM(REAL(alpha),A);
	  else
	    RowMajorMultXM(alpha,A);
	else if (A.iscm())
	  if (IMAG(alpha) == RealType(T)(0))
	    ColMajorMultXM(REAL(alpha),A);
	  else
	    ColMajorMultXM(alpha,A);
	else 
	  for(size_t i=0;i<A.colsize();++i) 
	    A.row(i,i,A.rowsize()) *= alpha;
      }
    }
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*Norm(A)) {
      cerr<<"TriMultXM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriMatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv

