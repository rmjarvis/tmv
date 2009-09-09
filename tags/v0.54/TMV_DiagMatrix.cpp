
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
  {
    this->itsdiv.reset(new DiagDiv<T>(*this)); 
  }

  template <class T> void DiagMatrixReadError<T>::Write(
      ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for DiagMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (m.get() && s != m->size()) {
      os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
    }
    if (!is) {
      if (iseof) {
	os<<"Input stream reached end-of-file prematurely.\n";
      } else if (isbad) {
	os<<"Input stream is corrupted.\n";
      } else {
	os<<"Input stream cannot read next character.\n";
      }
    }
    if (m.get()) {
      os<<"The portion of the DiagMatrix which was successfully read is: \n";
      ConstDiagMatrixView<T> mm = m->View();
      os<<"( ";
      for(size_t ii=0;ii<i;++ii)
	os<<' '<<mm(ii,ii)<<' ';
      os<<" )\n";
    }
  }

  template <class T, IndexStyle I> istream& operator>>(istream& is,
      auto_ptr<DiagMatrix<T,I> >& m)
  {
    char d;
    is >> d;
    if (!is || d != 'D') 
      throw DiagMatrixReadError<T>(is,'D',d);
    size_t size;
    is >> size;
    if (!is) 
      throw DiagMatrixReadError<T>(is);
    m.reset(new DiagMatrix<T,I>(size));
    try {
      m->diag().Read(is);
    }
    catch (VectorReadError<T>& ve) {
      throw DiagMatrixReadError<T>(ve.i,*m,ve.exp,ve.got,ve.s,
	  ve.is,ve.iseof,ve.isbad);
    }
    return is;
  }

  template <class T> istream& operator>>(istream& is,
      const DiagMatrixView<T>& m)
  {
    char d;
    is >> d;
    if (!is || d != 'D') 
      throw DiagMatrixReadError<T>(is,'D',d);
    size_t s;
    is >> s;
    if (!is) 
      throw DiagMatrixReadError<T>(is);
    if (m.size() != s)
      throw DiagMatrixReadError<T>(m,is,s);
    TMVAssert(m.size() == s);
    try {
      m.diag().Read(is);
    }
    catch (VectorReadError<T>& ve) {
      throw DiagMatrixReadError<T>(ve.i,m,ve.exp,ve.got,ve.s,
	  ve.is,ve.iseof,ve.isbad);
    }
    return is;
  }

#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


