
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"

namespace tmv {

  //
  // AddMM
  //

  template <class T1, class T2, class T3, class T4> inline void RowMajorAddMM(
      const T1 alpha, const GenMatrix<T2>& A, 
      const T3 beta, const MatrixView<T4>& B)
  {
    TMVAssert(A.isrm());
    TMVAssert(B.isrm());
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(alpha != T1(0));
    TMVAssert(beta != T3(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);

    const T2* Aptr = A.cptr();
    T4* Bptr = B.ptr();

    for(size_t i=0;i<A.colsize();++i,Aptr+=A.stepi(),Bptr+=B.stepi()) {
      if (beta != T3(1)) 
	DoMultXV(beta,VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),
	    B.rowsize());
      if (A.isconj())
	DoAddVV(alpha,CVIt<T2,Unit,Conj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),
	    B.rowsize());
      else
	DoAddVV(alpha,CVIt<T2,Unit,NonConj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),
	    B.rowsize());
    }
  }

  template <class T, class Ta> inline void RowAddMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const T beta, const MatrixView<T>& B)
  {
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(beta != T(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);

    MultXM(beta,B);
    for(size_t i=0;i<A.colsize();++i) AddVV(alpha,A.row(i),B.row(i));
  }


  template <class T1, class T2, class T3, class T4> inline void ColMajorAddMM(
      const T1 alpha, const GenMatrix<T2>& A, 
      const T3 beta, const MatrixView<T4>& B)
  {
    TMVAssert(A.iscm());
    TMVAssert(B.iscm());
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(alpha != T1(0));
    TMVAssert(beta != T3(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);

    const T2* Aptr = A.cptr();
    T4* Bptr = B.ptr();
    for(size_t j=0;j<A.rowsize();++j,Aptr+=A.stepj(),Bptr+=B.stepj()) {
      if (beta != T3(1)) 
	DoMultXV(beta,VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),
	    B.colsize());
      if (A.isconj())
	DoAddVV(alpha,CVIt<T2,Unit,Conj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),
	    B.colsize());
      else
	DoAddVV(alpha,CVIt<T2,Unit,NonConj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),
	    B.colsize());
    }
  }

  template <class T, class Ta> inline void ColAddMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const T beta, const MatrixView<T>& B)
  {
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(beta != T(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);

    MultXM(beta,B);
    for(size_t j=0;j<A.rowsize();++j) AddVV(alpha,A.col(j),B.col(j));
  }

  template <class T, class Ta> inline void DoAddMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const T beta, const MatrixView<T>& B)
  { 
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(beta != T(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);

    if (B.isconj()) DoAddMM(CONJ(alpha),A.QuickConjugate(),CONJ(beta),
	B.QuickConjugate());
    else if (A.stepi()==B.stepi() && A.stepj()==B.stepj() && 
	((B.stepj()==int(1) && B.stepi()==int(B.rowsize())) ||
	 (B.stepi()==int(1) && B.stepj()==int(B.colsize())) ||
	 (B.stepj()<=B.stepi() && B.stepi()==B.stepj()*int(A.rowsize())) ||
	 (B.stepi()<=B.stepj() && B.stepj()==B.stepi()*int(A.colsize())))) {
      ConstVectorView<Ta> AA = A.LinearView();
      VectorView<T> BB = B.LinearView();
      MultXV(beta,BB);
      AddVV(alpha,AA,BB);
    }
    else if (A.isrm() && B.isrm()) 
      if (IMAG(alpha) == RealType(T)(0))
	if (IMAG(beta) == RealType(T)(0))
	  RowMajorAddMM(REAL(alpha),A,REAL(beta),B); 
	else
	  RowMajorAddMM(REAL(alpha),A,beta,B); 
      else
	if (IMAG(beta) == RealType(T)(0))
	  RowMajorAddMM(alpha,A,REAL(beta),B); 
	else
	  RowMajorAddMM(alpha,A,beta,B); 
    else if (A.iscm() && B.iscm()) 
      if (IMAG(alpha) == RealType(T)(0))
	if (IMAG(beta) == RealType(T)(0))
	  ColMajorAddMM(REAL(alpha),A,REAL(beta),B); 
	else
	  ColMajorAddMM(REAL(alpha),A,beta,B); 
      else
	if (IMAG(beta) == RealType(T)(0))
	  ColMajorAddMM(alpha,A,REAL(beta),B); 
	else
	  ColMajorAddMM(alpha,A,beta,B); 
    else if (A.colsize() < A.rowsize()) RowAddMM(alpha,A,beta,B);
    else ColAddMM(alpha,A,beta,B);
  }

  template <class T, class Ta> inline void AddMM(const T alpha,
      const GenMatrix<Ta>& A, const T beta, const MatrixView<T>& B)
    // B = alpha * A + beta * B
  {
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    if (B.colsize() > 0 && B.rowsize() > 0) {
      if (alpha == T(0)) {
	if (beta != T(1)) MultXM(beta,B);
      } else if (beta == T(0)) {
	B = A;
	if (alpha != T(1)) MultXM(alpha,B);
      }
      else {
	if (B.SameStorageAs(A)) {
	  if (B.SameAs(A)) { 
	    T ab = alpha+beta;
	    if (ab != T(1)) MultXM(ab,B);
	  } else {
	    if (B.isrm()) {
	      Matrix<T,RowMajor> A2 = A;
	      DoAddMM(alpha,A2,beta,B);
	    } else {
	      Matrix<T,ColMajor> A2 = A;
	      DoAddMM(alpha,A2,beta,B);
	    }
	  }
	} 
	else DoAddMM(alpha,A,beta,B);
      }
    }
  }

#define InstFile "TMV_MatrixArith_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


