
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#define APTR (inplace ? A.NonConst().ptr() : new T[A.size()*A.size()])
#define LLX \
  (inplace ? A.uplo()==Upper ? A.NonConst().Adjoint() : A.NonConst() : \
  HermMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri())))
#define XDPTR (inplace ? new T[A.size()-1] : Aptr+LLx.stepj())
#define XD VectorViewOf(xDptr,A.size()-1,inplace ? 1 : int(A.size()+1))

  template <class T> HermLUDiv<T>::HermLUDiv(const GenSymMatrix<T>& A,
      bool _inplace) :
    inplace(_inplace && (A.isrm() || A.iscm())), Aptr(APTR), LLx(LLX), 
    xDptr(XDPTR), xD(XD),
    P(new size_t[A.colsize()]), det(T(1))
  {
    TMVAssert(IsReal(T()) || A.isherm());
    if (inplace) { TMVAssert(A == LLx); }
    else LLx = A;

    SymLU_Decompose(LLx,xD,P,det);
  }

#undef LLX
#define LLX \
  inplace ? A.uplo()==Upper ? A.NonConst().Transpose() : A.NonConst() : \
  SymMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri()))

  template <class T> SymLUDiv<T>::SymLUDiv(const GenSymMatrix<T>& A,
      bool _inplace) :
    inplace(_inplace && (A.isrm() || A.iscm())), Aptr(APTR), LLx(LLX), 
    xDptr(XDPTR), xD(XD),
    P(new size_t[A.colsize()]), det(T(1))
  {
    TMVAssert(IsComplex(T()) && !A.isherm());
    if (inplace) { TMVAssert(A == LLx); }
    else LLx = A;

    SymLU_Decompose(LLx,xD,P,det);
  }

#undef APTR
#undef LLX
#undef XDPTR
#undef XD

  template <class T> HermLUDiv<T>::~HermLUDiv() 
  { delete[] P; if (!inplace) delete[] Aptr; else delete[] xDptr; }

  template <class T> SymLUDiv<T>::~SymLUDiv() 
  { delete[] P; if (!inplace) delete[] Aptr; else delete[] xDptr; }

  template <class T> bool SymLUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenSymMatrix<T>* sm = dynamic_cast<const GenSymMatrix<T>*>(&m);
    TMVAssert(sm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*sm)<<"  "<<*sm<<endl;
      *fout << "L = "<<GetL()<<endl;
      *fout << "D = "<<GetD()<<endl;
    }
    Matrix<T> lu = GetL()*GetD()*GetL().Transpose();
    lu.ReversePermuteRows(P);
    lu.ReversePermuteCols(P);
    RealType(T) nm = Norm(lu-*sm);
    nm /= SQR(Norm(GetL()))*Norm(GetD());
    if (fout) {
      *fout << "LDLT = "<<lu<<endl;
      *fout << "Norm(M-LDLT)/Norm(LDLT) = "<<nm<<endl;
      *fout << "M-LDLT = "<<(*sm-lu)<<endl;
      *fout << "Norm(m) = "<<Norm(*sm)<<endl;
      *fout << "Norm(LDLT) = "<<SQR(Norm(GetL()))*Norm(GetD())<<endl;
    }
    SymMatrix<T> sm2 = *sm;
    sm2.DivideUsing(SVS);
    sm2.SetDiv();
    return nm < sm2.Condition()*sm2.size()*Epsilon<T>();
  }

  template <class T> bool HermLUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenSymMatrix<T>* sm = dynamic_cast<const GenSymMatrix<T>*>(&m);
    TMVAssert(sm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*sm)<<"  "<<*sm<<endl;
      *fout << "L = "<<GetL()<<endl;
      *fout << "D = "<<GetD()<<endl;
      *fout << "P = ";
      for(size_t i=0;i<sm->size();i++) *fout<<P[i]<<" ";
      *fout <<endl;
    }
    Matrix<T> lu = GetL()*GetD()*GetL().Adjoint();
    lu.ReversePermuteRows(P);
    lu.ReversePermuteCols(P);
    RealType(T) nm = Norm(lu-*sm);
    nm /= SQR(Norm(GetL()))*Norm(GetD());
    if (fout) {
      *fout << "LDLt = "<<lu<<endl;
      *fout << "M-LDLt = "<<(*sm-lu)<<endl;
      *fout << "Norm(M-LDLt)/Norm(LDLt) = "<<nm<<endl;
    }
    HermMatrix<T> sm2 = *sm;
    sm2.DivideUsing(SVS);
    sm2.SetDiv();
    return nm < sm2.Condition()*sm2.size()*Epsilon<T>();
  }

  template <class T> const Matrix<T> HermLUDiv<T>::GetD() const 
  { 
    Matrix<T> temp(LLx.size(),LLx.size(),T(0));
    temp.diag() = LLx.diag();
    temp.diag(1) = xD.Conjugate();
    temp.diag(-1) = xD;
    return temp;
  }

  template <class T> const Matrix<T> SymLUDiv<T>::GetD() const 
  { 
    Matrix<T> temp(LLx.size(),LLx.size(),T(0));
    temp.diag() = LLx.diag();
    temp.diag(1) = temp.diag(-1) = xD;
    return temp;
  }

#define InstFile "TMV_SymLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


