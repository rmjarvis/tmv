
#include "TMV_Tri.h"
#include "TMV_Diag.h"

namespace tmv {

  template <class T> bool UpperTriDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenUpperTriMatrix<T>* tm = dynamic_cast<const GenUpperTriMatrix<T>*>(&m);
    TMVAssert(tm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*tm)<<"  "<<*tm<<endl;
      *fout << "T = "<<*itsm<<endl;
    }
    RealType(T) nm = Norm(*itsm-*tm);
    nm /= Norm(*itsm);
    if (fout) {
      *fout << "Norm(M-T)/Norm(T) = "<<nm<<endl;
    }
    return nm < Epsilon<T>();
  }

  template <class T> bool LowerTriDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenLowerTriMatrix<T>* tm = dynamic_cast<const GenLowerTriMatrix<T>*>(&m);
    TMVAssert(tm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*tm)<<"  "<<*tm<<endl;
      *fout << "T = "<<*itsm<<endl;
    }
    RealType(T) nm = Norm(*itsm-*tm);
    nm /= Norm(*itsm);
    if (fout) {
      *fout << "Norm(M-T)/Norm(T) = "<<nm<<endl;
    }
    return nm < Epsilon<T>();
  }

  template <class T> T UpperTriDiv<T>::Det() const
  {
    if (!donedet) {
      if (!itsm->isunit()) det *= DiagMatrixViewOf(itsm->diag()).Det();
      donedet = true;
    }
    return det;  
  }                  

  template <class T> T LowerTriDiv<T>::Det() const
  {
    if (!donedet) {
      if (!itsm->isunit()) det *= DiagMatrixViewOf(itsm->diag()).Det();
      donedet = true;
    }
    return det;  
  }                  

#define InstFile "TMV_TriDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


