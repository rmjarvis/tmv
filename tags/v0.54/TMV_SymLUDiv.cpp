
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Diag.h"

namespace tmv {

#define APTR1 inplace ? 0 : new T[A.size()*A.size()]
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()
#define LLX \
  inplace ? \
    A.uplo()==Upper ? \
      (A.isherm() ? A.NonConst().Adjoint() : A.NonConst().Transpose()) : \
      A.NonConst() : \
    A.isherm() ? \
      HermMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri())) : \
      SymMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri()))

  template <class T> SymLUDiv<T>::SymLUDiv(const GenSymMatrix<T>& A,
      bool _inplace) :
    inplace(_inplace && (A.isrm() || A.iscm())), 
    Aptr1(APTR1), Aptr(APTR), LLx(LLX), 
    xD(Vector<T>(A.size()-1)),
    P(new size_t[A.colsize()]), det(T(1))
  {
    if (inplace) TMVAssert(A == LLx); 
    else LLx = A;

    SymLU_Decompose(LLx,xD.View(),P.get(),det);
    //cerr<<"After LUDecompose: \n";
    //cerr<<"xD = "<<xD<<endl;
  }

#undef APTR
#undef LLX
#undef XDPTR
#undef XD

  template <class T> bool SymLUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    //cerr<<"Start CheckDecomp: xD = "<<xD<<endl;
    const GenSymMatrix<T>* sm = dynamic_cast<const GenSymMatrix<T>*>(&m);
    TMVAssert(sm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*sm)<<"  "<<*sm<<endl;
      *fout << "L = "<<GetL()<<endl;
      *fout << "D = "<<GetD()<<endl;
      *fout << "P = ";
      for(size_t i=0;i<sm->size();i++) *fout<<P.get()[i]<<" ";
      *fout <<endl;
    }
    Matrix<T> lu = GetL()*GetD()*
      (LLx.isherm()?GetL().Adjoint():GetL().Transpose());
    lu.ReversePermuteRows(P.get());
    lu.ReversePermuteCols(P.get());
    RealType(T) nm = Norm(lu-*sm);
    nm /= SQR(Norm(GetL()))*Norm(GetD());
    if (fout) {
      *fout << "LDLt = "<<lu<<endl;
      *fout << "M-LDLt = "<<(*sm-lu)<<endl;
      *fout << "Norm(M-LDLt)/Norm(LDLt) = "<<nm<<endl;
      *fout << "Norm(m) = "<<Norm(*sm)<<endl;
      *fout << "Norm(LDLt) = "<<SQR(Norm(GetL()))*Norm(GetD())<<endl;
    }
    RealType(T) kappa;
    if (LLx.isherm()) {
      HermMatrix<T> sm2 = *sm;
      sm2.DivideUsing(SVS);
      kappa = sm2.Condition();
    } else {
      SymMatrix<T> sm2 = *sm;
      sm2.DivideUsing(SVS);
      kappa = sm2.Condition();
    }
    //cerr<<"Done CheckDecomp: xD = "<<xD<<endl;
    return nm < kappa*sm->size()*Epsilon<T>();
  }

  template <class T> const Matrix<T> SymLUDiv<T>::GetD() const 
  { 
    Matrix<T> temp(LLx.size(),LLx.size(),T(0));
    temp.diag() = LLx.diag();
    temp.diag(-1) = xD;
    temp.diag(1) = LLx.isherm() ? xD.Conjugate() : xD.View();
    return temp;
  }

#define InstFile "TMV_SymLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


