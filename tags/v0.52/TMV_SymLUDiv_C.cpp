
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

  template <class T> void HermLUDiv<T>::Inverse(
      const SymMatrixView<T>& sinv) const
  {
    TMVAssert(sinv.size() == LLx.size());
    TMVAssert(sinv.isherm());
    // MJ: This is a placeholder.  Obviously there is a faster version.
    Matrix<T> minv(LLx.size(),LLx.size());
    Inverse(minv.View());
    sinv = HermMatrixViewOf(minv,sinv.uplo());
  }

  template <class T> void HermLUDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LLx.size());
    TMVAssert(minv.rowsize() == LLx.size());

    minv.SetToIdentity();
    LDivEq(minv.View());
  }

  template <class T> void HermLUDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LLx.size());
    TMVAssert(minv.rowsize() == LLx.size());

    Matrix<T,ColMajor> temp(LLx.size(),LLx.size());
    temp.SetToIdentity();
    LDivEq(temp.View());
    minv =  temp*temp.Adjoint();
  }

  template <class T> void SymLUDiv<T>::Inverse(
      const SymMatrixView<T>& sinv) const
  {
    TMVAssert(sinv.size() == LLx.size());
    TMVAssert(sinv.issym());
    // MJ: This is a placeholder.  Obviously there is a faster version.
    Matrix<T> minv(LLx.size(),LLx.size());
    Inverse(minv.View());
    sinv = SymMatrixViewOf(minv,sinv.uplo());
  }

  template <class T> void SymLUDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LLx.size());
    TMVAssert(minv.rowsize() == LLx.size());

    minv.SetToIdentity();
    LDivEq(minv.View());
  }

  template <class T> void SymLUDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LLx.size());
    TMVAssert(minv.rowsize() == LLx.size());

    Matrix<T,ColMajor> temp(LLx.size(),LLx.size());
    temp.SetToIdentity();
    LDivEq(temp.View());
    minv =  temp*temp.Adjoint();
  }

#define InstFile "TMV_SymLUDiv_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


