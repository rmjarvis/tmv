
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Band.h"

//#define XDEBUG

namespace tmv {

  //
  // AddMM
  //

  template <class T1, class T2, class T3, class T4> inline void RowMajorAddMM(
      const T1 alpha, const GenBandMatrix<T2>& A, 
      const T3 beta, const BandMatrixView<T4>& B)
  {
    //cerr<<"RowMajor Band AddMM\n";
    TMVAssert(A.isrm());
    TMVAssert(B.isrm());
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(A.nlo() == B.nlo());
    TMVAssert(A.nhi() == B.nhi());
    TMVAssert(alpha != T1(0));
    TMVAssert(beta != T1(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.ct() == NonConj);

    size_t j1=0;
    size_t j2=B.nhi()+1;
    size_t k=B.nlo();
    size_t len = B.nhi()+1;
    const int Ads = A.stepi()+1;
    const int Bds = B.stepi()+1;
    const T2* Aptr = A.cptr();
    T4* Bptr = B.ptr();

    for(size_t i=0;i<B.colsize();++i) {
      if (beta!=T3(1)) 
	DoMultXV(beta,VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),
	    len);
      if (A.isconj())
	DoAddVV(alpha,CVIt<T2,Unit,Conj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),len);
      else
	DoAddVV(alpha,CVIt<T2,Unit,NonConj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),len);
      if (k>0) {--k; ++len; Aptr+=A.stepi(); Bptr+=B.stepi(); } 
      else { ++j1; Aptr+=Ads; Bptr+=Bds; }
      if (j2<B.rowsize()) ++j2;
      else if (j1==B.rowsize()) break;
      else --len;
    }
  }

  template <class T, class T2> inline void DoRowMajorAddMM(
      const T alpha, const GenBandMatrix<T2>& A, 
      const T beta, const BandMatrixView<T>& B)
  {
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
  }

  template <class T1, class T2, class T3, class T4> inline void ColMajorAddMM(
      const T1 alpha, const GenBandMatrix<T2>& A, 
      const T3 beta, const BandMatrixView<T4>& B)
  {
    //cerr<<"ColMajor Band AddMM\n";
    TMVAssert(A.iscm());
    TMVAssert(B.iscm());
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(A.nlo() == B.nlo());
    TMVAssert(A.nhi() == B.nhi());
    TMVAssert(alpha != T1(0));
    TMVAssert(beta != T1(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.ct() == NonConj);

    size_t i1=0;
    size_t i2=B.nlo()+1;
    size_t k=B.nhi();
    size_t len=B.nlo()+1;
    const int Ads=A.stepj()+1;
    const int Bds=B.stepj()+1;
    const T2* Aptr = A.cptr();
    T4* Bptr = B.ptr();

    for(size_t j=0;j<A.rowsize();++j) {
      if (beta!=T3(1)) 
	DoMultXV(beta,VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),
	    len);
      if (A.isconj())
	DoAddVV(alpha,CVIt<T2,Unit,Conj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),len);
      else
	DoAddVV(alpha,CVIt<T2,Unit,NonConj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),len);
      if (k>0) { --k; ++len; Aptr+=A.stepj(); Bptr+=B.stepj(); }  
      else { ++i1; Aptr+=Ads; Bptr+=Bds; }
      if (i2<B.colsize()) ++i2;
      else if (i1==B.colsize()) break;
      else --len;
    }
  }

  template <class T, class T2> inline void DoColMajorAddMM(
      const T alpha, const GenBandMatrix<T2>& A, 
      const T beta, const BandMatrixView<T>& B)
  {
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
  }

  template <class T1, class T2, class T3, class T4> inline void DiagMajorAddMM(
      const T1 alpha, const GenBandMatrix<T2>& A, 
      const T3 beta, const BandMatrixView<T4>& B)
  {
    //cerr<<"DiagMajor Band AddMM\n";
    TMVAssert(A.isdm());
    TMVAssert(B.isdm());
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(A.nlo() == B.nlo());
    TMVAssert(A.nhi() == B.nhi());
    TMVAssert(alpha != T1(0));
    TMVAssert(beta != T1(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.ct() == NonConj);

    size_t len=min(B.colsize()-B.nlo(),B.rowsize());
    size_t j2 = len;
    const T2* Aptr = A.cptr()+B.nlo()*A.stepi();
    T4* Bptr = B.ptr()+B.nlo()*B.stepi();

    for(int i=-B.nlo();i<=B.nhi();++i) {
      if (beta!=T3(1)) 
	DoMultXV(beta,VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),
	    len);
      if (A.isconj())
	DoAddVV(alpha,CVIt<T2,Unit,Conj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),len);
      else
	DoAddVV(alpha,CVIt<T2,Unit,NonConj>(Aptr,1),
	    VIt<T4,Unit,NonConj>(Bptr,1 FIRSTLAST1(B.first,B.last) ),len);
      if (i<0) { 
	Aptr -= A.stepi();
	Bptr -= B.stepi();
	if (j2<A.rowsize()) { ++len; ++j2; } 
      }
      else { 
	Aptr += A.stepj();
	Bptr += B.stepj();
	if (j2<A.rowsize()) ++j2; else --len; 
      }
    }
  }

  template <class T, class T2> inline void DoDiagMajorAddMM(
      const T alpha, const GenBandMatrix<T2>& A, 
      const T beta, const BandMatrixView<T>& B)
  {
    if (IMAG(alpha)==RealType(T)(0))
      if (IMAG(beta)==RealType(T)(0))
	DiagMajorAddMM(REAL(alpha),A,REAL(beta),B); 
      else
	DiagMajorAddMM(REAL(alpha),A,beta,B); 
    else
      if (IMAG(beta)==RealType(T)(0))
	DiagMajorAddMM(alpha,A,REAL(beta),B); 
      else
	DiagMajorAddMM(alpha,A,beta,B); 
  }

  template <class T, class Ta> inline void DiagAddMM(
      const T alpha, const GenBandMatrix<Ta>& A, 
      const T beta, const BandMatrixView<T>& B)
  {
    //cerr<<"Diag Band AddMM\n";
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(A.nlo() == B.nlo());
    TMVAssert(A.nhi() == B.nhi());
    TMVAssert(alpha != T(0));
    TMVAssert(beta != T(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.ct() == NonConj);

    for(int i=-B.nlo();i<=B.nhi();++i) {
      MultXV(beta,B.diag(i));
      AddVV(alpha,A.diag(i),B.diag(i));
    }
  }

  template <class T, class Ta> inline void DoAddMM(
      const T alpha, const GenBandMatrix<Ta>& A, 
      const T beta, const BandMatrixView<T>& B)
  { 
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(B.nlo() >= A.nlo());
    TMVAssert(B.nhi() >= A.nhi());
    TMVAssert(alpha != T(0));
    TMVAssert(beta != T(0));
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.rowsize() > 0);

    if (B.isconj()) DoAddMM(CONJ(alpha),A.QuickConjugate(),CONJ(beta),
	B.QuickConjugate());
    else {
      BandMatrixView<T> BB = BandMatrixViewOf(B,A.nlo(),A.nhi());

      if (A.isrm() && B.isrm()) DoRowMajorAddMM(alpha,A,beta,BB);
      else if (A.iscm() && B.iscm()) DoColMajorAddMM(alpha,A,beta,BB);
      else if (A.isdm() && B.isdm()) DoDiagMajorAddMM(alpha,A,beta,BB);
      else DiagAddMM(alpha,A,beta,BB);

      if (beta != T(1)) {
	for(int i=-B.nlo();i<-A.nlo();++i) MultXV(beta,B.diag(i));
	for(int i=A.nhi()+1;i<=B.nhi();++i) MultXV(beta,B.diag(i));
      }
    }
  }

  template <class T, class Ta> inline void AddMM(const T alpha,
      const GenBandMatrix<Ta>& A, 
      const T beta, const BandMatrixView<T>& B)
    // B = alpha * A + beta * B
  {
    //cerr<<"Band AddMM: alpha = "<<alpha<<"  beta = "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
#ifdef XDEBUG
    Matrix<T> B2 = alpha*Matrix<T>(A)+beta*Matrix<T>(B);
#endif
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(B.nlo() >= A.nlo());
    TMVAssert(B.nhi() >= A.nhi());

    if (B.colsize() > 0 && B.rowsize() > 0) {
      if (alpha == T(0)) {
	if (beta != T(1)) MultXM(beta,B);
      }
      else if (beta == T(0)) {
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
	      BandMatrix<T,RowMajor> A2 = A;
	      DoAddMM(alpha,A2,beta,B);
	    } else if (B.iscm()) {
	      BandMatrix<T,ColMajor> A2 = A;
	      DoAddMM(alpha,A2,beta,B);
	    } else {
	      BandMatrix<T,DiagMajor> A2 = A;
	      DoAddMM(alpha,A2,beta,B);
	    }
	  }
	} 
	else DoAddMM(alpha,A,beta,B);
      }
    }
#ifdef XDEBUG
    if (Norm(B2-B) > 0.0001) {
      cerr<<"Band AddMM\n";
      cerr<<"alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"   "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"   "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_BandMatrixArith_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


