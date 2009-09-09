
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

  template <class T, class T1, class T2> void SV_LDiv(
      const GenMatrix<T1>& UV,
      const GenVector<T1>& Ubeta, const GenVector<T1>& Vbeta, 
      const GenMatrix<T1>* UV2, const GenVector<T1>* Qbeta,
      const GenMatrix<RealType(T1)>& U1, const GenMatrix<RealType(T1)>& V1,
      const GenVector<RealType(T1)>& S, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T>& x)
  {
    // A x = m
    // U U1 S V1 V x = m
    // x = Vt V1t S^-1 U1t Ut m
    TMVAssert(m.colsize() == UV.colsize()); // = M
    TMVAssert(x.colsize() == UV.rowsize()); // = N
    TMVAssert(x.rowsize() == m.rowsize()); // = R
    TMVAssert(kmax <= UV.rowsize());
    TMVAssert(UV.rowsize() <= UV.colsize());
    
    const size_t N = UV.rowsize();

    Matrix<T,ColMajor> m2 = m; // MxR
    if (UV2) {
      TMVAssert(Qbeta);
      Q_LDivEq(UV,*Qbeta,m2.View());
      Q_LDivEq(*UV2,Ubeta,m2.Rows(0,N));
    } else {
      Q_LDivEq(UV,Ubeta,m2.View());
    }
    Matrix<T,ColMajor> m3 = U1.Adjoint().Rows(0,kmax) * m2.Rows(0,N); 
      // = KxR
    m3 /= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = V1.Adjoint().Cols(0,kmax) * m3; // NxR
    if (UV2) {
      Q_RDivEq(UV2->SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  x.Rows(1,N).Transpose());
    } else {
      Q_RDivEq(UV.SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  x.Rows(1,N).Transpose());
    }
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

  template <class T, class T1, class T2> void SV_RDiv(
      const GenMatrix<T1>& UV,
      const GenVector<T1>& Ubeta, const GenVector<T1>& Vbeta, 
      const GenMatrix<T1>* UV2, const GenVector<T1>* Qbeta,
      const GenMatrix<RealType(T1)>& U1, const GenMatrix<RealType(T1)>& V1,
      const GenVector<RealType(T1)>& S, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T>& x)
  {
    // x A = m
    // x U U1 S V1 V = m
    // x = m Vt V1t S^-1 U1t Ut
    TMVAssert(m.rowsize() == UV.rowsize()); // = N
    TMVAssert(x.rowsize() == UV.colsize()); // = M
    TMVAssert(x.colsize() == m.colsize()); // = R
    TMVAssert(kmax <= UV.rowsize()); // = K
    TMVAssert(UV.rowsize() <= UV.colsize());
    
    const size_t N = UV.rowsize();

    Matrix<T,ColMajor> m2 = m; // = RxN
    if (UV2) {
      Q_LDivEq(UV2->SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  m2.Cols(1,N).Transpose());
    } else {
      Q_LDivEq(UV.SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  m2.Cols(1,N).Transpose());
    }
    Matrix<T,ColMajor> m3 = m2 * V1.Adjoint().Cols(0,kmax); // = RxK
    m3 %= DiagMatrixViewOf(S.SubVector(0,kmax));
    x.Cols(0,N) = m3 * U1.Adjoint().Rows(0,kmax); // = RxN
    x.Cols(N,x.rowsize()).Zero(); // x = RxM
    if (UV2) {
      TMVAssert(Qbeta);
      Q_RDivEq(*UV2,Ubeta,x.Cols(0,N));
      Q_RDivEq(UV,*Qbeta,x);
    } else {
      Q_RDivEq(UV,Ubeta,x);
    }
  }

#define InstFile "TMV_SVDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


