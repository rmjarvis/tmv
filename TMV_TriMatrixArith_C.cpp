
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Tri.h"

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

    for(size_t i=0;i<N;++i) AddVV(alpha,A.row(i,i,N),B.row(i,i,N));
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
    //cerr<<"TriAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    TMVAssert(!B.isunit());
    TMVAssert(A.size() == B.size());
    if (A.size() > 0) {
      if (B.isconj()) AddMM(CONJ(alpha),A.QuickConjugate(),
	  CONJ(beta),B.QuickConjugate());
      else if (alpha == T(0)) {
	MultXM(beta,B);
      } else if (beta == T(0)) {
	B = A;
	MultXM(alpha,B);
      }
      else {
	if (B.SameStorageAs(A)) {
	  if (B.SameAs(A)) { 
	    MultXM(alpha+beta,B);
	  } else {
	    if (B.isrm()) {
	      UpperTriMatrix<T,NonUnitDiag,RowMajor> A2 = A;
	      MultXM(beta,A2.View());
	      DoAddMM(alpha,A2,B);
	    } else if (B.iscm()) {
	      UpperTriMatrix<T,NonUnitDiag,ColMajor> A2 = A;
	      MultXM(beta,A2.View());
	      DoAddMM(alpha,A2,B);
	    } else {
	      UpperTriMatrix<T,NonUnitDiag,DiagMajor> A2 = A;
	      MultXM(beta,A2.View());
	      DoAddMM(alpha,A2,B);
	    }
	  }
	} 
	else {
	  MultXM(beta,B);
	  DoAddMM(alpha,A,B);
	}
      }
    }
    //cerr<<"->B = "<<B<<endl;
  }

  template <class T, class Ta, class Tb> inline void AddMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, 
      const T beta, const GenUpperTriMatrix<Tb>& B, const MatrixView<T>& C)
    // C = alpha * A + beta * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());
    //cerr<<"TriAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    if (C.isconj()) AddMM(CONJ(alpha),A.QuickConjugate(),CONJ(beta),
	B.QuickConjugate(),C.QuickConjugate());
    else if (A.size() > 0) {
      if (C.SameStorageAs(A) && C.SameStorageAs(B)) {
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
	LowerTriMatrixView<T> LC(C,UnitDiag);
	UpperTriMatrixView<T> UC(C,UnitDiag);
	if (C.SameStorageAs(A)) {
	  LC.OffDiag() = A.OffDiag();
	  MultXM(alpha,LC.OffDiag());
	  UC.OffDiag() = B.OffDiag();
	  MultXM(beta,UC.OffDiag());
	} else {
	  UC.OffDiag() = B.OffDiag();
	  MultXM(beta,UC.OffDiag());
	  LC.OffDiag() = A.OffDiag();
	  MultXM(alpha,LC.OffDiag());
	}
	if (A.isunit()) {
	  if (B.isunit()) {
	    C.diag().SetAllTo(alpha+beta);
	  } else {
	    C.diag()=B.diag();
	    MultXV(beta,C.diag());
	    C.diag().AddToAll(alpha);
	  }
	} else {
	  if (B.isunit()) {
	    C.diag()=A.diag();
	    MultXV(alpha,C.diag());
	    C.diag().AddToAll(beta);
	  } else {
	    C.diag()=A.diag();
	    MultXV(alpha,C.diag());
	    AddVV(beta,B.diag(),C.diag());
	  }
	}
      }
    }
    //cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
  }

#define InstFile "TMV_TriMatrixArith_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


