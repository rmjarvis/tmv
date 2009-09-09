
#include "TMV_Diag.h"

namespace tmv {

  template <class T> void GenDiagMatrix<T>::DInverse(
      const DiagMatrixView<T>& minv) const
  {
    this->SetDiv();
    const DiagDiv<T>* ddiv =
      dynamic_cast<const DiagDiv<T>*>(this->GetDiv());
    TMVAssert(ddiv);
    ddiv->DInverse(minv);
  }

  template <class T> void GenDiagMatrix<T>::DInverseATA(
      const DiagMatrixView<T>& minv) const
  {
    this->SetDiv();
    const DiagDiv<T>* ddiv =
      dynamic_cast<const DiagDiv<T>*>(this->GetDiv());
    TMVAssert(ddiv);
    ddiv->DInverseATA(minv);
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

  template <class T> istream& operator>>(istream& fin, DiagMatrix<T>* m)
  {
    size_t s;
    fin >> s;
    if (!fin) tmv_error("reading size in DiagMatrix::Read");
    m = new DiagMatrix<T>(s);
    m->diag().Read(fin);
    return fin;
  }

  template <class T> istream& operator>>(istream& fin,
      const DiagMatrixView<T>& m)
  {
    size_t s;
    fin >> s;
    if (!fin) tmv_error("reading size in DiagMatrix::Read");
    TMVAssert(m.colsize() == s);
    m.diag().Read(fin);
    return fin;
  }



#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


