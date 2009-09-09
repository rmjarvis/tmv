
#include "TMV.h"

namespace tmv {

  //
  // AddMM
  //

  template <bool rm, bool a1, bool ca, class T, class T1, class T2> 
    inline void DoRowAddMM(const T1 alpha, const GenMatrix<T2>& A, 
	const MatrixView<T>& B)
    {
      TMVAssert(A.colsize() == B.colsize());
      TMVAssert(A.rowsize() == B.rowsize());
      TMVAssert(alpha != T1(0));
      TMVAssert(B.colsize() > 0);
      TMVAssert(B.rowsize() > 0);
      TMVAssert(B.ct() == NonConj);
      TMVAssert(!B.SameStorageAs(A));
      TMVAssert(rm == (A.isrm() && B.isrm()));
      TMVAssert(a1 == (alpha == T1(1)));
      TMVAssert(ca == A.isconj());

      const T2* Arowi = A.cptr();
      T* Browi = B.ptr();
      const size_t M = A.colsize();
      const size_t N = A.rowsize();
      const int Asi = A.stepi();
      const int Asj = (rm ? 1 : A.stepj());
      const int Bsi = B.stepi();
      const int Bsj = (rm ? 1 : B.stepj());

      for(size_t i=M;i>0;--i,Arowi+=Asi,Browi+=Bsi) {
	const T2* Aij = Arowi;
	T* Bij = Browi;
	for(size_t j=N;j>0;--j,(rm?++Aij:Aij+=Asj),(rm?++Bij:Bij+=Bsj)) {
	  if (a1) *Bij += (ca ? CONJ(*Aij) : *Aij);
	  else *Bij += alpha * (ca ? CONJ(*Aij) : *Aij);
	}
      }
    }

  template <bool rm, class T, class Ta> inline void RowAddMM(
      const T alpha, const GenMatrix<Ta>& A, const MatrixView<T>& B)
  { 
    if (IMAG(alpha) == RealType(T)(0))
      if (REAL(alpha) == RealType(T)(1))
	if (A.isconj()) DoRowAddMM<rm,true,true>(REAL(alpha),A,B);
	else DoRowAddMM<rm,true,false>(REAL(alpha),A,B);
      else
	if (A.isconj()) DoRowAddMM<rm,false,true>(REAL(alpha),A,B);
	else DoRowAddMM<rm,false,false>(REAL(alpha),A,B);
    else
      if (A.isconj()) DoRowAddMM<rm,false,true>(alpha,A,B);
      else DoRowAddMM<rm,false,false>(alpha,A,B);
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
    TMVAssert(B.ct() == NonConj);
    TMVAssert(!B.SameStorageAs(A));

    if (A.stor() == B.stor() && A.CanLinearize() && B.CanLinearize()) {
      TMVAssert(A.stepi() == B.stepi() && A.stepj() == B.stepj());
      AddVV(alpha,A.ConstLinearView(),beta,B.LinearView());
    } else {
      if (beta != T(1)) B *= beta;

      if (A.isrm() && B.isrm())
	RowAddMM<true>(alpha,A,B); 
      else if (A.iscm() && B.iscm())
	RowAddMM<true>(alpha,A.Transpose(),B.Transpose()); 
      else if (A.rowsize() > A.colsize())
	RowAddMM<false>(alpha,A,B); 
      else
	RowAddMM<false>(alpha,A.Transpose(),B.Transpose()); 
    }
  }

  template <class T, class Ta> void AddMM(const T alpha,
      const GenMatrix<Ta>& A, const T beta, const MatrixView<T>& B)
    // B = alpha * A + beta * B
  {
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    if (B.colsize() > 0 && B.rowsize() > 0) {
      if (B.isconj()) 
	AddMM(CONJ(alpha),A.Conjugate(),CONJ(beta),B.Conjugate());
      else if (alpha == T(0)) {
	B *= beta;
      } else if (beta == T(0)) {
	B = A;
	B *= alpha;
      } else {
	if (B.SameStorageAs(A)) {
	  if (B.SameAs(A)) { 
	    B *= alpha+beta;
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


