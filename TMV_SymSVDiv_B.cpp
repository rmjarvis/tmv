
#include "TMV_Sym.h"
#include "TMV.h"
#include "TMV_Givens.h"
#include "TMV_Diag.h"
#include "TMV_Band.h"
#include "TMV_Householder.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t SYM_TRIDIAG_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t SYM_TRIDIAG_BLOCKSIZE = 64;
#endif

  template <class T> void HermSVDiv<T>::Inverse(
      const SymMatrixView<T>& sinv) const
  {
    TMVAssert(sinv.size() == S.size());
    TMVAssert(sinv.isherm());
    // MJ: This is a placeholder.  Obviously there is a faster version.
    Matrix<T> minv(S.size(),S.size());
    Inverse(minv.View());
    sinv = SymMatrixViewOf(minv,sinv.uplo());
  }

  template <class T> void HermSVDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  { 
    TMVAssert(U);
    // A^-1 = U S^-1 Ut
    Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = U->Cols(0,kmax) * SinvUt;
  }

  template <class T> void HermSVDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(U);
    // A = U S Ut
    // At = U S Ut
    // AtA = U S^2 Ut
    // (AtA)^-1 = U S^-2 Ut
    //
    Matrix<T,RowMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = SinvUt.Adjoint() * SinvUt;
  }

  template <class T> void SymSVDiv<T>::Inverse(
      const SymMatrixView<T>& sinv) const
  {
    TMVAssert(sinv.size() == S.size());
    TMVAssert(sinv.issym());
    // MJ: This is a placeholder.  Obviously there is a faster version.
    Matrix<T> minv(S.size(),S.size());
    Inverse(minv.View());
    sinv = HermMatrixViewOf(minv,sinv.uplo());
  }

  template <class T> void SymSVDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  { 
    TMVAssert(U&&V);
    // A = U S V
    // A^-1 = Vt S^-1 Ut
    Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = V->Adjoint().Cols(0,kmax) * SinvUt;
  }

  template <class T> void SymSVDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(U&&V);
    // A = U S V
    // At = Vt S Ut
    // AtA = Vt S^2 V
    // (AtA)^-1 = Vt S^-2 V
    //
    Matrix<T,RowMajor> SinvV = V->Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = SinvV.Adjoint() * SinvV;
  }

#define InstFile "TMV_SymSVDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


