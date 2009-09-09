
#include "TMV.h"
#include "TMV_Diag.h"

namespace tmv {

  //
  // LDiv
  //

  template <class T, class T1, class T2> void SV_LDiv(
      const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
      const GenMatrix<T1>& V, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T>& x)
  {
    // A x = m
    // U S V x = m
    // x = Vt S^-1 Ut m
    TMVAssert(m.colsize() == U.colsize()); // = M
    TMVAssert(x.colsize() == V.rowsize()); // = N
    TMVAssert(x.rowsize() == m.rowsize()); // = R
    TMVAssert(kmax <= V.rowsize()); // = K
    TMVAssert(kmax <= U.colsize());
    Matrix<T> m2 = U.Adjoint().Rows(0,kmax) * m; // KxR
    m2 /= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = V.Adjoint().Cols(0,kmax) * m2; // NxR
  }

  //
  // RDiv
  //

  template <class T, class T1, class T2> void SV_RDiv(
      const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
      const GenMatrix<T1>& V, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T>& x) 
  {
    // x A = m
    // x U S V = m
    // x = m Vt S^-1 Ut
    TMVAssert(m.rowsize() == V.rowsize()); // = N
    TMVAssert(x.rowsize() == U.colsize()); // = M
    TMVAssert(x.colsize() == m.colsize()); // = R
    TMVAssert(kmax <= U.colsize()); // = K
    TMVAssert(kmax <= V.rowsize());

    Matrix<T,ColMajor> m2 = m * V.Adjoint().Cols(0,kmax); // = RxK
    m2 %= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = m2 * U.Adjoint().Rows(0,kmax); // = RxM
  }

#define InstFile "TMV_SVDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


