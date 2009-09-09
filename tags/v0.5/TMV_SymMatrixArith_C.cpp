
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  //
  // AddMM
  //

  template <class T, class Ta> void AddMM(const T alpha,
      const GenSymMatrix<Ta>& A, 
      const T beta, const SymMatrixView<T>& B)
    // B = alpha * A + beta * B
  {
#ifdef XDEBUG
    Matrix<T> B2 = alpha*Matrix<T>(A) + beta*Matrix<T>(B);
    Matrix<T> B0 = B;
#endif

    AddMM(alpha,A.UpperTri(),beta,B.UpperTri()); 

#ifdef XDEBUG
    if (Norm(B-B2) > 0.001*Norm(B)) {
      cerr<<"SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> void AddMM(const T alpha,
      const GenSymMatrix<Ta>& A, 
      const T beta, const MatrixView<T>& B)
    // B = alpha * A + beta * B
  {
#ifdef XDEBUG
    Matrix<T> B2 = alpha*Matrix<T>(A) + beta*Matrix<T>(B);
    Matrix<T> B0 = B;
#endif

    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());

    AddMM(alpha,A.UpperTri(),beta,UpperTriMatrixViewOf(B));
    if (A.size() > 0)
      AddMM(alpha,A.LowerTri().OffDiag(),beta,
	  LowerTriMatrixViewOf(B).OffDiag());

#ifdef XDEBUG
    if (Norm(B-B2) > 0.001*Norm(B)) {
      cerr<<"SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymMatrixArith_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


