
#include "TMV_Sym.h"
#include "TMV.h"
#include "TMV_Diag.h"

namespace tmv {

  template <class T> void HermSVDiv<T>::Inverse(
      const SymMatrixView<T>& sinv) const
  {
    TMVAssert(sinv.size() == S.size());
    TMVAssert(sinv.isherm());
    Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    SymMultMM(T(1),U->Cols(0,kmax),SinvUt,0,sinv);
  }

  template <class T> void HermSVDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  { 
    TMVAssert(U.get());
    // A^-1 = U S^-1 Ut
    Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = U->Cols(0,kmax) * SinvUt;
  }

  template <class T> void HermSVDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(U.get());
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
    Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    SymMultMM(T(1),V->Adjoint().Cols(0,kmax),SinvUt,0,sinv);
  }

  template <class T> void SymSVDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  { 
    TMVAssert(U.get() && V.get());
    // A = U S V
    // A^-1 = Vt S^-1 Ut
    Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = V->Adjoint().Cols(0,kmax) * SinvUt;
  }

  template <class T> void SymSVDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(U.get() && V.get());
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


