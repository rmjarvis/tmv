
#include "TMV.h"
#include "TMV_Diag.h"

namespace tmv {

#define APTR1 inplace ? 0 : new T[A.colsize()*A.rowsize()]     
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()
#define UX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) : \
  inplace ? A.NonConst().View() : \
  MatrixViewOf(Aptr,A.colsize(),A.rowsize(), \
      A.IsSquare() ? BaseStorOf(A) : ColMajor)

  template <class T> SVDiv<T>::SVDiv(const GenMatrix<T>& A,
      bool _inplace, bool StoreU, bool StoreV) :
    istrans(A.colsize() < A.rowsize()), 
    inplace(_inplace && (A.isrm() || A.iscm())), Aptr1(APTR1), Aptr(APTR),
    U(new MatrixView<T>(UX)), S(U->rowsize()), V(0), det(T(1))
  {
    if (istrans) {
      if (inplace) TMVAssert(A.Transpose() == *U); 
      else *U = A.Transpose();
    } else {
      if (inplace) TMVAssert(A == *U); 
      else *U = A;
    }

    if (istrans) swap(StoreU,StoreV);

    if (StoreV) {
      V.reset(new Matrix<T,ColMajor>(U->rowsize(),U->rowsize()));
      SV_Decompose(*U,S.View(),V->View(),det,StoreU);
    }
    else SV_Decompose(*U,S.View(),det,StoreU);

    if (!StoreU && !inplace) { Aptr1.reset(0); }
    if (!StoreU) { U.reset(0); }

    // Set kmax for actual 0 elements (to within machine precision).
    // Any further cut in the number of singular values to use
    // should be done by the user.
    Thresh(Epsilon<T>());
  }

#undef UX
#undef APTR

  template <class T> bool SVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenMatrix<T>* mm = dynamic_cast<const GenMatrix<T>*>(&m);
    TMVAssert(mm);
    if (fout) {
      *fout << "SVDiv CheckDecomp\n";
      *fout << "M = "<<tmv::Type(*mm)<<"  "<<*mm<<endl;
      *fout << "U = "<<GetU()<<endl;
      *fout << "S = "<<GetS()<<endl;
      *fout << "V = "<<GetV()<<endl;
    }
    Matrix<T> usv = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(usv-*mm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    RealType(T) cond = S(0) / (S(kmax-1));
    if (fout) {
      *fout << "USV = "<<usv<<endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
      *fout <<"  "<<cond<<" * "<<Epsilon<T>()<<endl;
    }
    return nm < cond*Epsilon<T>();
  }

  template <class T> void SVDiv<T>::Thresh(RealType(T) toler,
      ostream* debugout) const
  {
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = S(0)*toler;
    for(kmax=S.size(); kmax>0 && S(kmax-1)<=thresh; --kmax);
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"Smax = "<<S(0)<<", thresh = "<<thresh<<endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

  template <class T> void SVDiv<T>::Top(size_t neigen,
      ostream* debugout) const
  {
    TMVAssert(neigen > 0);
    TMVAssert(neigen <= S.size());
    kmax = neigen;
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

#define InstFile "TMV_SVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


