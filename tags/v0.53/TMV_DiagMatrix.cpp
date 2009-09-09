
#include "TMV_Diag.h"

namespace tmv {

  template <class T> void GenDiagMatrix<T>::Inverse(
      const DiagMatrixView<T>& minv) const
  {
    this->SetDiv();
    const DiagDiv<T>* ddiv =
      dynamic_cast<const DiagDiv<T>*>(this->GetDiv());
    TMVAssert(ddiv);
    ddiv->Inverse(minv);
  }

  template <class T> void GenDiagMatrix<T>::InverseATA(
      const DiagMatrixView<T>& minv) const
  {
    this->SetDiv();
    const DiagDiv<T>* ddiv =
      dynamic_cast<const DiagDiv<T>*>(this->GetDiv());
    TMVAssert(ddiv);
    ddiv->InverseATA(minv);
  }

  template <class T> template <class T1> void GenDiagMatrix<T>::DivEq(
      const DiagMatrixView<T1>& m) const
  {
    this->SetDiv();
    const DiagDiv<T>* ddiv =
      dynamic_cast<const DiagDiv<T>*>(this->GetDiv());
    TMVAssert(ddiv);
    ddiv->DivEq(m);
  }

  template <class T> template <class T1, class T0> void GenDiagMatrix<T>::Div(
      const GenDiagMatrix<T1>& m1, const DiagMatrixView<T0>& m0) const
  {
    this->SetDiv();
    const DiagDiv<T>* ddiv =
      dynamic_cast<const DiagDiv<T>*>(this->GetDiv());
    TMVAssert(ddiv);
    ddiv->Div(m1,m0);
  }

  template <class T> void GenDiagMatrix<T>::NewDivider() const
  { this->itsdiv = new DiagDiv<T>(*this); }

  template <class T, IndexStyle I> istream& operator>>(istream& is,
      DiagMatrix<T,I>*& m)
  {
    char d;
    is >> d;
    if (!is || d != 'D') 
      tmv_error("reading D in DiagMatrix::Read");
    size_t size;
    is >> size;
    if (!is) 
      tmv_error("reading size in DiagMatrix::Read");
    m = new DiagMatrix<T,I>(size);
    m->diag().Read(is);
    return is;
  }

  template <class T> istream& operator>>(istream& is,
      const DiagMatrixView<T>& m)
  {
    char d;
    is >> d;
    if (!is || d != 'D') 
      tmv_error("reading D in DiagMatrix::Read");
    size_t size;
    is >> size;
    if (!is) 
      tmv_error("reading size in DiagMatrix::Read");
    if (m.size() != size)
      tmv_error("size does not match in DiagMatrix::Read");
    TMVAssert(m.colsize() == size);
    m.diag().Read(is);
    return is;
  }



#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


