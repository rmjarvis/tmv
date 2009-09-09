
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"
#include "TMV_Sym.h"

namespace tmv {

#define RecursiveCH

#ifdef TMV_BLOCKSIZE
  const size_t CH_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t CH_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t CH_BLOCKSIZE = 64;
  const size_t CH_BLOCKSIZE2 = 2;
#endif

#define APTR inplace ? A.NonConst().ptr() : new T[A.size()*A.size()]
#define LLX \
  inplace ? A.uplo()==Upper ? A.NonConst().Adjoint() : A.NonConst() : \
  HermMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri()))

  template <class T> HermCHDiv<T>::HermCHDiv(const GenSymMatrix<T>& A,
      bool _inplace) :
    inplace(_inplace && (A.iscm() || A.isrm())), Aptr(APTR), LLx(LLX), 
    det(T(0)),donedet(false)
  {
    TMVAssert(IsReal(T()) || A.isherm());
    if (inplace) { TMVAssert(A==LLx); }
    else LLx = A;
    HermCH_Decompose(LLx);
  }

  template <class T> HermCHDiv<T>::~HermCHDiv() 
  { if (!inplace) delete[] Aptr; }

#undef APTR
#undef LLX

  template <class T> bool HermCHDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenSymMatrix<T>* sm = dynamic_cast<const GenSymMatrix<T>*>(&m);
    TMVAssert(sm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*sm)<<"  "<<*sm<<endl;
      *fout << "L = "<<GetL()<<endl;
    }
    Matrix<T> lu = GetL()*GetL().Adjoint();
    RealType(T) nm = Norm(lu-*sm);
    nm /= SQR(Norm(GetL()));
    if (fout) {
      *fout << "LLt = "<<lu<<endl;
      *fout << "Norm(M-LLt)/Norm(LLt) = "<<nm<<endl;
    }
    HermMatrix<T> sm2 = *sm;
    sm2.DivideUsing(SVS);
    sm2.SetDiv();
    return nm < sm2.Condition()*sm2.size()*Epsilon<T>();
  }

  //
  // LDivEq
  //

  template <class T1, class T2> void HermCH_LDivEq(
      const GenSymMatrix<T1>& LL, const MatrixView<T2>& m)
  {
    TMVAssert(LL.size() == m.colsize());
    // m = (LLt)^-1 m
    //   = Lt^-1 L^-1 m
    m /= LL.LowerTri();
    m /= LL.UpperTri();
  }

  //
  // RDivEq Matrix
  //

  template <class T1, class T2> void HermCH_RDivEq(
      const GenSymMatrix<T1>& LL, const MatrixView<T2>& m)
  {
    TMVAssert(LL.size() == m.rowsize());
    // m = m (LLt)^-1 
    //   = m Lt^-1 L^-1
    m %= LL.UpperTri();
    m %= LL.LowerTri();
  }

  template <class T> T HermCHDiv<T>::Det() const
  {
    if (!donedet) {
      det = DiagMatrixViewOf(LLx.diag()).Det();
      det *= det;
      donedet = true;
    }
    return det;
  }

#define InstFile "TMV_SymCHDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


