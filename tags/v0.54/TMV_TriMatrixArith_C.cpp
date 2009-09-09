
#include "TMV.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

  //
  // AddMM
  //

  template <bool rm, bool a1, bool ca, class T1, class T2, class T3> 
    inline void DoRowAddMM(
	const T1 alpha, const GenUpperTriMatrix<T2>& A, 
	const UpperTriMatrixView<T3>& B)
    {
      TMVAssert(!A.isunit());
      TMVAssert(!B.isunit());
      TMVAssert(A.size() == B.size());
      TMVAssert(A.size() > 0);
      TMVAssert(alpha != T1(0));
      TMVAssert(B.ct() == NonConj);

      const size_t N = A.size();
      const int Astepj = rm ? 1 : A.stepj();
      const int Ads = A.stepi() + Astepj;
      const int Bstepj = rm ? 1 : B.stepj();
      const int Bds = B.stepi() + Bstepj;
      const T2* Aii = A.cptr();
      T3* Bii = B.ptr();

      for(size_t len=N;len>0;--len,Aii+=Ads,Bii+=Bds) {
	const T2* Aij = Aii;
	T3* Bij = Bii;
	for(size_t j=len;j>0;--j,(rm?++Aij:Aij+=Astepj),
	    (rm?++Bij:Bij+=Bstepj)) {
	  if (a1) *Bij += (ca ? CONJ(*Aij) : *Aij);
	  else *Bij += alpha * (ca ? CONJ(*Aij) : *Aij);
	}
      }
    }

  template <bool rm, class T, class Ta> inline void RowAddMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A, 
      const UpperTriMatrixView<T>& B)
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

  template <bool cm, bool a1, bool ca, class T1, class T2, class T3> 
    inline void DoColAddMM(
	const T1 alpha, const GenUpperTriMatrix<T2>& A, 
	const UpperTriMatrixView<T3>& B)
    {
      TMVAssert(!A.isunit());
      TMVAssert(!B.isunit());
      TMVAssert(A.size() == B.size());
      TMVAssert(A.size() > 0);
      TMVAssert(alpha != T1(0));
      TMVAssert(B.ct() == NonConj);

      const size_t N = A.size();
      const int Astepi = (cm ? 1 : A.stepi());
      const int Astepj = A.stepj();
      const int Bstepi = (cm ? 1 : B.stepi());
      const int Bstepj = B.stepj();
      const T2* Acolj = A.cptr()+(N-1)*Astepj;
      T3* Bcolj = B.ptr()+(N-1)*Bstepj;

      for(size_t j=N;j>0;--j,Acolj-=Astepj,Bcolj-=Bstepj) {
	const T2* Aij = Acolj;
	T3* Bij = Bcolj;
	for(size_t i=j;i>0;--i,(cm?++Aij:Aij+=Astepi),(cm?++Bij:Bij+=Bstepi)) {
	  if (a1) *Bij += (ca ? CONJ(*Aij) : *Aij);
	  else *Bij += alpha * (ca ? CONJ(*Aij) : *Aij);
	}
      }
    }

  template <bool cm, class T, class Ta> inline void ColAddMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A, 
      const UpperTriMatrixView<T>& B)
  {
    if (IMAG(alpha) == RealType(T)(0)) 
      if (REAL(alpha) == RealType(T)(1))
	if (A.isconj()) DoColAddMM<cm,true,true>(REAL(alpha),A,B); 
	else DoColAddMM<cm,true,false>(REAL(alpha),A,B);
      else
	if (A.isconj()) DoColAddMM<cm,false,true>(REAL(alpha),A,B); 
	else DoColAddMM<cm,false,false>(REAL(alpha),A,B);
    else
      if (A.isconj()) DoColAddMM<cm,false,true>(alpha,A,B); 
      else DoColAddMM<cm,false,false>(alpha,A,B);
  }

  template <class T, class Ta> inline void DoAddMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  { 
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);

    if (A.isunit()) {
      B.diag().AddToAll(alpha);
      DoAddMM(alpha,A.OffDiag(),B.OffDiag());
    } else {
      if (A.isrm() && B.isrm()) RowAddMM<true>(alpha,A,B); 
      else if (A.iscm() && B.iscm()) ColAddMM<true>(alpha,A,B); 
      else RowAddMM<false>(alpha,A,B); 
    }
  }

  template <class T, class Ta> void AddMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, 
      const T beta, const UpperTriMatrixView<T>& B)
    // B = alpha * A + beta * B
  {
#ifdef XDEBUG
    Matrix<T> B2 = alpha*Matrix<T>(A) + beta*Matrix<T>(B);
    Matrix<T> B0 = B;
#endif

    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    if (A.size() > 0) {
      if (B.isconj()) AddMM(CONJ(alpha),A.Conjugate(),
	  CONJ(beta),B.Conjugate());
      else if (alpha == T(0)) {
	B *= beta;
      } else if (beta == T(0)) {
	B = A;
	B *= alpha;
      }
      else {
	if (B.SameStorageAs(A)) {
	  if (B.SameAs(A)) { 
	    B *= (alpha + beta);
	  } else {
	    if (B.isrm()) {
	      UpperTriMatrix<Ta,NonUnitDiag,RowMajor> A2 = A;
	      B *= beta;
	      DoAddMM(alpha,A2,B);
	    } else {
	      UpperTriMatrix<Ta,NonUnitDiag,ColMajor> A2 = A;
	      B *= beta;
	      DoAddMM(alpha,A2,B);
	    }
	  }
	} 
	else {
	  B *= beta;
	  DoAddMM(alpha,A,B);
	}
      }
    }
#ifdef XDEBUG
    if (Norm(B-B2) > 0.001*Norm(B)) {
      cerr<<"TriAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> void AddMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, 
      const T beta, const GenUpperTriMatrix<Tb>& B, const MatrixView<T>& C)
    // C = alpha * A + beta * B
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C2 = alpha*A0 + beta*B0;
#endif

    TMVAssert(A.size() == B.size());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());

    if (C.isconj()) AddMM(CONJ(alpha),A.Conjugate(),CONJ(beta),
	B.Conjugate(),C.Conjugate());
    else if (A.size() > 0) {
      if (C.SameStorageAs(A) && C.SameStorageAs(B)) {
	if (A.isrm() || B.isrm()) {
	  Matrix<T,RowMajor> temp(C.colsize(),C.rowsize());
	  AddMM(alpha,A,beta,B,temp.View());
	  C = temp;
	} else {
	  Matrix<T,ColMajor> temp(C.colsize(),C.rowsize());
	  AddMM(alpha,A,beta,B,temp.View());
	  C = temp;
	}
      } else {
	LowerTriMatrixView<T> LC = LowerTriMatrixViewOf(C,UnitDiag);
	UpperTriMatrixView<T> UC = UpperTriMatrixViewOf(C,UnitDiag);
	if (C.SameStorageAs(A)) {
	  LC.OffDiag() = alpha * A.OffDiag();
	  UC.OffDiag() = beta * B.OffDiag();
	} else {
	  UC.OffDiag() = beta * B.OffDiag();
	  LC.OffDiag() = alpha * A.OffDiag();
	}
	if (A.isunit()) {
	  if (B.isunit()) {
	    C.diag().SetAllTo(alpha+beta);
	  } else {
	    C.diag() = B.diag();
	    C.diag() *= beta;
	    C.diag().AddToAll(alpha);
	  }
	} else {
	  if (B.isunit()) {
	    C.diag() = A.diag();
	    C.diag() *= alpha;
	    C.diag().AddToAll(beta);
	  } else if (C.SameStorageAs(A)) {
	    C.diag() = alpha * A.diag();
	    C.diag() += beta * B.diag();
	  } else {
	    C.diag() = beta * B.diag();
	    C.diag() += alpha * A.diag();
	  }
	}
      }
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(abs(alpha)*Norm(A0)+abs(beta)*Norm(B0))) {
      cerr<<"TriAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriMatrixArith_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


