
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

namespace tmv {

#define APTR1 inplace ? 0 : new T[A.colsize()*A.rowsize()]
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()
#define LUX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) : \
  inplace ? A.NonConst().View() : \
  MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)

  template <class T> LUDiv<T>::LUDiv(const GenMatrix<T>& A, bool _inplace) :
    istrans(A.isrm()), inplace(_inplace && (A.iscm() || A.isrm())),
    Aptr1(APTR1), Aptr(APTR), LUx(LUX), P(new size_t[A.colsize()]),
    det(T(1)), donedet(false)
  {
    TMVAssert(A.IsSquare());
    if (istrans) {
      if (inplace) TMVAssert(A.Transpose() == LUx);
      else LUx = A.Transpose();
    }
    else {
      if (inplace) TMVAssert(A == LUx);
      else LUx = A;
    }
    LU_Decompose(LUx,P.get(),det);
  }
#undef LUX
#undef APTR

  template <class T> bool LUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenMatrix<T>* mm = dynamic_cast<const GenMatrix<T>*>(&m);
    TMVAssert(mm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*mm)<<"  ";
      *fout << (istrans ? mm->Transpose() : mm->View()) <<endl;
      *fout << "L = "<<GetL()<<endl;
      *fout << "U = "<<GetU()<<endl;
    }
    Matrix<T> lu = GetL()*GetU();
    lu.ReversePermuteRows(GetP());
    RealType(T) nm = Norm(lu- (istrans ? mm->Transpose() : mm->View()) );
    nm /= Norm(GetL())*Norm(GetU());
    if (fout) {
      *fout << "PLU = "<<lu<<endl;
      *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<endl;
    }
    Matrix<T> m2 = *mm;
    m2.DivideUsing(SVS);
    m2.SetDiv();
    return nm < m2.Condition()*m2.colsize()*Epsilon<T>();
  }

  template <class T> T LUDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(LUx.diag()).Det();
      donedet = true;
    }
    return det;
  }

#define InstFile "TMV_LUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


