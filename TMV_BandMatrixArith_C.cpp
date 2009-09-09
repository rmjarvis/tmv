
#include "TMV.h"
#include "TMV_Band.h"

//#define XDEBUG

namespace tmv {

  //
  // AddMM
  //

  template <class T, class Ta> void DoAddMM(
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
    TMVAssert(B.ct() == NonConj);
    TMVAssert(!B.SameStorageAs(A));

    if (A.stor() == B.stor() && A.nlo() == B.nlo() && A.nhi() == B.nhi() && 
	A.CanLinearize() && B.CanLinearize()) {
      TMVAssert(A.stepi() == B.stepi() && A.stepj() == B.stepj());
      AddVV(alpha,A.ConstLinearView(),beta,B.LinearView());
    } else {
      for(int i=-A.nlo();i<=A.nhi();++i) 
	AddVV(alpha,A.diag(i),beta,B.diag(i));
      if (beta != T(1)) {
	for(int i=-B.nlo();i<-A.nlo();++i) B.diag(i) *= beta;
	for(int i=A.nhi()+1;i<=B.nhi();++i) B.diag(i) *= beta;
      }
    }
  }

  template <class T, class Ta> void AddMM(const T alpha,
      const GenBandMatrix<Ta>& A, 
      const T beta, const BandMatrixView<T>& B)
    // B = alpha * A + beta * B
  {
#ifdef XDEBUG
    //cerr<<"Band AddMM: alpha = "<<alpha<<"  beta = "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    Matrix<T> B2 = alpha*Matrix<T>(A)+beta*Matrix<T>(B);
#endif
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(B.nlo() >= A.nlo());
    TMVAssert(B.nhi() >= A.nhi());

    if (B.colsize() > 0 && B.rowsize() > 0) {
      if (B.isconj())
	AddMM(CONJ(alpha),A.Conjugate(),CONJ(beta),B.Conjugate());
      else if (alpha == T(0)) 
	B *= beta;
      else if (beta == T(0)) {
	B = A; 
	B *= alpha;
      }
      else {
	if (B.SameStorageAs(A)) {
	  if (B.SameAs(A)) {
	    B *= alpha+beta;
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
    if (Norm(B2-B) > 0.001*Norm(B)) {
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


