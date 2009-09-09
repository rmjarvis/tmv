
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  //
  // MultXM
  //

  template <class T> void MultXM(const T alpha, const SymMatrixView<T>& A)
    // A = alpha * A
  {
#ifdef XDEBUG
    Matrix<T> A2 = alpha * Matrix<T>(A);
    Matrix<T> A0 = A;
#endif

    TMVAssert(!A.isherm() || IMAG(alpha) == RealType(T)(0));

    A.UpperTri() *= alpha;

#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*max(RealType(T)(1),Norm(A))) {
      cerr<<"SymMultXM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymMatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv

