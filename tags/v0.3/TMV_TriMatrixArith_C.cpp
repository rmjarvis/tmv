
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Tri.h"

namespace tmv {

  //
  // AddMM
  //

  template <class T1, class T2, class T3, class T4> inline void RowMajorAddMM(
      const T1 alpha, const GenUpperTriMatrix<T2>& A, 
      const T3 beta, const UpperTriMatrixView<T4>& B)
  {
    TMVAssert(A.isrm());
    TMVAssert(B.isrm());
    TMVAssert(!A.isunit());
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T1(0));
    TMVAssert(beta != T3(0));
    TMVAssert(B.ct() == NonConj);

    const size_t N = A.size();
    const int Ads = A.stepi()+1;
    const int Bds = B.stepi()+1;
    const T2* Aptr = A.cptr();
    T4* Bptr = B.ptr();

    if (beta != T3(1)) B *= beta;
    if (A.isconj())
      for(size_t i=0,len=N;len>0;++i,--len,Aptr+=Ads,Bptr+=Bds) 
	DoAddVV(alpha,CVIt<T2,Unit,Conj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1),len);
    else
      for(size_t i=0,len=N;len>0;++i,--len,Aptr+=Ads,Bptr+=Bds) 
	DoAddVV(alpha,CVIt<T2,Unit,NonConj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1),len);
  }

  template <class T, class Ta> inline void RowAddMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, 
      const T beta, const UpperTriMatrixView<T>& B)
  {
    TMVAssert(!A.isunit());
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(beta != T(0));
    TMVAssert(B.ct() == NonConj);

    const size_t N = A.size();

    if (beta != T(1)) B *= beta;
    if (alpha == T(1)) {
      for(size_t i=0;i<N;++i) B.row(i,i,N) += A.row(i,i,N);
    } else {
      for(size_t i=0;i<N;++i) B.row(i,i,N) += alpha * A.row(i,i,N);
    }
  }

  template <class T1, class T2, class T3, class T4> inline void ColMajorAddMM(
      const T1 alpha, const GenUpperTriMatrix<T2>& A, 
      const T3 beta, const UpperTriMatrixView<T4>& B)
  {
    TMVAssert(A.iscm());
    TMVAssert(B.iscm());
    TMVAssert(!A.isunit());
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T1(0));
    TMVAssert(beta != T3(0));
    TMVAssert(B.ct() == NonConj);

    const size_t N = A.size();
    const T2* Aptr = A.cptr();
    T4* Bptr = B.ptr();

    if (beta != T3(1)) B *= beta;
    if (A.isconj())
      for(size_t j=0;j<N;++j,Aptr+=A.stepj(),Bptr+=B.stepj()) 
	DoAddVV(alpha,CVIt<T2,Unit,Conj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1),j+1);
    else
      for(size_t j=0;j<N;++j,Aptr+=A.stepj(),Bptr+=B.stepj()) 
	DoAddVV(alpha,CVIt<T2,Unit,NonConj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1),j+1);
  }

  template <class T, class Ta> inline void DoAddMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, 
      const T beta, const UpperTriMatrixView<T>& B)
  { 
    TMVAssert(!A.isunit());
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(beta != T(0));
    TMVAssert(B.ct() == NonConj);

    if (A.isrm() && B.isrm()) 
      if (IMAG(alpha)==RealType(T)(0))
	if (IMAG(beta)==RealType(T)(0))
	  RowMajorAddMM(REAL(alpha),A,REAL(beta),B); 
	else
	  RowMajorAddMM(REAL(alpha),A,beta,B); 
      else
	if (IMAG(beta)==RealType(T)(0))
	  RowMajorAddMM(alpha,A,REAL(beta),B); 
	else
	  RowMajorAddMM(alpha,A,beta,B); 
    else if (A.iscm() && B.iscm())
      if (IMAG(alpha)==RealType(T)(0))
	if (IMAG(beta)==RealType(T)(0))
	  ColMajorAddMM(REAL(alpha),A,REAL(beta),B); 
	else
	  ColMajorAddMM(REAL(alpha),A,beta,B); 
      else
	if (IMAG(beta)==RealType(T)(0))
	  ColMajorAddMM(alpha,A,REAL(beta),B); 
	else
	  ColMajorAddMM(alpha,A,beta,B); 
    else RowAddMM(alpha,A,beta,B); 
  }

  // MJ: SameStorage checks
  template <class T, class Ta> inline void AddMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, 
      const T beta, const UpperTriMatrixView<T>& B)
    // B = alpha * A + beta * B
  {
    //cerr<<"TriAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    if (A.size() > 0) {
      if (B.isconj()) AddMM(CONJ(alpha),A.QuickConjugate(),
	  CONJ(beta),B.QuickConjugate());
      else if (alpha == T(0)) {
	if (beta != T(1)) MultXM(beta,B);
      } else if (beta == T(0)) {
	B = A;
	if (alpha != T(1)) MultXM(alpha,B);
      }
      else
	if (A.isunit()) {
	  if (beta != T(1)) B.diag() *= beta;
	  B.diag().AddToAll(alpha);
	  DoAddMM(alpha,A.OffDiag(),beta,B.OffDiag());
	}
	else DoAddMM(alpha,A,beta,B); 
    }
    //cerr<<"->B = "<<B<<endl;
  }

  template <class T, class Ta, class Tb> inline void AddMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, 
      const T beta, const GenUpperTriMatrix<Tb>& B, const MatrixView<T>& C)
    // C = alpha * A + beta * B
  {
    //cerr<<"TriAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    if (C.isconj()) AddMM(CONJ(alpha),A.QuickConjugate(),CONJ(beta),
	B.QuickConjugate(),C.QuickConjugate());
    else if (C.SameStorageAs(A) && C.SameStorageAs(B)) {
      if (A.isrm() || B.isrm()) {
	Matrix<T,RowMajor> temp(C.colsize(),C.rowsize());
	AddMM(alpha,A,beta,B,temp.QuickView());
	C = temp;
      } else {
	Matrix<T,ColMajor> temp(C.colsize(),C.rowsize());
	AddMM(alpha,A,beta,B,temp.QuickView());
	C = temp;
      }
    } else {
      //cerr<<"Normal:\n";
      LowerTriMatrixView<T> LC(C,UnitDiag);
      UpperTriMatrixView<T> UC(C,UnitDiag);
      if (C.SameStorageAs(A)) {
	LC.OffDiag() = alpha * A.OffDiag();
	UC.OffDiag() = beta * B.OffDiag();
      } else {
	UC.OffDiag() = beta * B.OffDiag();
	//cerr<<"B.Off = "<<B.OffDiag()<<endl;
	//cerr<<"UC.Off = "<<UC.OffDiag()<<endl;
	//cerr<<"UC = "<<UC<<endl;
	LC.OffDiag() = alpha * A.OffDiag();
	//cerr<<"B.Off = "<<A.OffDiag()<<endl;
	//cerr<<"LC.Off = "<<LC.OffDiag()<<endl;
	//cerr<<"LC = "<<LC<<endl;
      }
      if (A.isunit()) {
	if (B.isunit()) {
	  C.diag().SetAllTo(alpha+beta);
	} else {
	  C.diag() = beta*B.diag();
	  C.diag().AddToAll(alpha);
	}
      } else {
	if (B.isunit()) {
	  C.diag() = alpha*A.diag();
	  C.diag().AddToAll(beta);
	} else {
	  C.diag() = alpha*A.diag();
	  C.diag() +=+ beta*B.diag();
	}
      }
    }
    //cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
  }

#define InstFile "TMV_TriMatrixArith_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


