
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

  //
  // AddMM
  //

  template <class T1, class T2, class T3> inline void RowMajorAddMM(
      const T1 alpha, const GenUpperTriMatrix<T2>& A, 
      const UpperTriMatrixView<T3>& B)
  {
    TMVAssert(A.isrm());
    TMVAssert(B.isrm());
    TMVAssert(!A.isunit());
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T1(0));
    TMVAssert(B.ct() == NonConj);

    const size_t N = A.size();
    const int Ads = A.stepi()+1;
    const int Bds = B.stepi()+1;
    const T2* Aptr = A.cptr();
    T3* Bptr = B.ptr();

    // B.row(i,i,N) += alpha * A.row(i,i,N);
    if (A.isconj())
      for(size_t i=0,len=N;len>0;++i,--len,Aptr+=Ads,Bptr+=Bds) 
	DoAddVV(alpha,CVIt<T2,Unit,Conj>(Aptr,1),
	    VIt<T3,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),len);
    else
      for(size_t i=0,len=N;len>0;++i,--len,Aptr+=Ads,Bptr+=Bds) 
	DoAddVV(alpha,CVIt<T2,Unit,NonConj>(Aptr,1),
	    VIt<T3,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),len);
  }

  template <class T, class T2> inline void DoRowMajorAddMM(
      const T alpha, const GenUpperTriMatrix<T2>& A, 
      const UpperTriMatrixView<T>& B)
  {
    if (IMAG(alpha)==RealType(T)(0)) RowMajorAddMM(REAL(alpha),A,B); 
    else RowMajorAddMM(alpha,A,B); 
  }

  template <class T, class Ta> inline void RowAddMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, 
      const UpperTriMatrixView<T>& B)
  {
    TMVAssert(!A.isunit());
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);

    const size_t N = A.size();

    for(size_t i=0;i<N;++i) B.row(i,i,N) += alpha * A.row(i,i,N);
  }

  template <class T1, class T2, class T3> inline void ColMajorAddMM(
      const T1 alpha, const GenUpperTriMatrix<T2>& A, 
      const UpperTriMatrixView<T3>& B)
  {
    TMVAssert(A.iscm());
    TMVAssert(B.iscm());
    TMVAssert(!A.isunit());
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T1(0));
    TMVAssert(B.ct() == NonConj);

    const size_t N = A.size();
    const T2* Aptr = A.cptr();
    T3* Bptr = B.ptr();

    // B.col(j,0,j+1) += alpha * A.col(j,0,j+1);
    if (A.isconj())
      for(size_t j=0;j<N;++j,Aptr+=A.stepj(),Bptr+=B.stepj()) 
	DoAddVV(alpha,CVIt<T2,Unit,Conj>(Aptr,1),
	    VIt<T3,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),j+1);
    else
      for(size_t j=0;j<N;++j,Aptr+=A.stepj(),Bptr+=B.stepj()) 
	DoAddVV(alpha,CVIt<T2,Unit,NonConj>(Aptr,1),
	    VIt<T3,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),j+1);
  }

  template <class T, class T2> inline void DoColMajorAddMM(
      const T alpha, const GenUpperTriMatrix<T2>& A, 
      const UpperTriMatrixView<T>& B)
  {
    if (IMAG(alpha)==RealType(T)(0)) ColMajorAddMM(REAL(alpha),A,B); 
    else ColMajorAddMM(alpha,A,B); 
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
    }
    else {
      if (A.isrm() && B.isrm()) DoRowMajorAddMM(alpha,A,B); 
      else if (A.iscm() && B.iscm()) DoColMajorAddMM(alpha,A,B); 
      else RowAddMM(alpha,A,B); 
    }
  }

  template <class T, class Ta> inline void AddMM(const T alpha,
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
	    B *= alpha + beta;
	  } else {
	    if (B.isrm()) {
	      UpperTriMatrix<Ta,NonUnitDiag,RowMajor> A2 = A;
	      B *= beta;
	      DoAddMM(alpha,A2,B);
	    } else if (B.iscm()) {
	      UpperTriMatrix<Ta,NonUnitDiag,ColMajor> A2 = A;
	      B *= beta;
	      DoAddMM(alpha,A2,B);
	    } else {
	      UpperTriMatrix<Ta,NonUnitDiag,DiagMajor> A2 = A;
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

  template <class T, class Ta, class Tb> inline void AddMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, 
      const T beta, const GenUpperTriMatrix<Tb>& B, const MatrixView<T>& C)
    // C = alpha * A + beta * B
  {
#ifdef XDEBUG
    Matrix<T> C2 = alpha*Matrix<T>(A) + beta*Matrix<T>(B);
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
	LowerTriMatrixView<T> LC(C,UnitDiag);
	UpperTriMatrixView<T> UC(C,UnitDiag);
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
    if (Norm(C-C2) > 0.001*Norm(C)) {
      cerr<<"TriAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
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


