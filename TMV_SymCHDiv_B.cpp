
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

#define RecursiveCH

#ifdef TMV_BLOCKSIZE
  const size_t CH_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t CH_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t CH_BLOCKSIZE = 64;
  const size_t CH_BLOCKSIZE2 = 2;
#endif

  template <class T> void NonLapCHInverse(const SymMatrixView<T>& sinv)
  {
    TMVAssert(sinv.isherm());
    //cerr<<"LLt = "<<sinv<<endl;
    // inv = (L Lt)^-1 = Lt^-1 L^-1
    LowerTriMatrixView<T> L = sinv.LowerTri();
    //cerr<<"L = "<<L<<endl;
    L = L.Inverse();
    //cerr<<"Linv = "<<L<<endl;
    //cerr<<"PUL = "<<Type(L.Adjoint()*L)<<endl;
    sinv = L.Adjoint() * L;
    //cerr<<"sinv = "<<sinv<<endl;
  }

#ifdef LAP
  template <class T> void LapCHInverse(const SymMatrixView<T>& sinv)
  { NonLapCHInverse(sinv); }
  template <> void LapCHInverse(const SymMatrixView<double>& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    char uplo = (sinv.uplo() == Upper ? 'U' : 'L');
    int n = sinv.size();
    int lda = sinv.stepj();
    int info;
    dpotri(&uplo,&n,sinv.ptr(),&lda,&info);
    if (info < 0) tmv_error("dpotri returned info < 0");
  }
  template <> void LapCHInverse(const SymMatrixView<complex<double> >& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    char uplo = (sinv.uplo() == Upper ? 'U' : 'L');
    int n = sinv.size();
    int lda = sinv.stepj();
    int info;
    zpotri(&uplo,&n,LAP_Complex(sinv.ptr()),&lda,&info);
    if (info < 0) tmv_error("zpotri returned info < 0");
  }
#ifndef NOFLOAT
  template <> void LapCHInverse(const SymMatrixView<float>& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    char uplo = (sinv.uplo() == Upper ? 'U' : 'L');
    int n = sinv.size();
    int lda = sinv.stepj();
    int info;
    spotri(&uplo,&n,sinv.ptr(),&lda,&info);
    if (info < 0) tmv_error("spotri returned info < 0");
  }
  template <> void LapCHInverse(const SymMatrixView<complex<float> >& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    char uplo = (sinv.uplo() == Upper ? 'U' : 'L');
    int n = sinv.size();
    int lda = sinv.stepj();
    int info;
    cpotri(&uplo,&n,LAP_Complex(sinv.ptr()),&lda,&info);
    if (info < 0) tmv_error("cpotri returned info < 0");
  }
#endif // FLOAT
#endif // LAP

  template <class T> void DoCHInverse(const GenSymMatrix<T>& LLx,
      const SymMatrixView<T>& sinv) 
  {
    if (sinv.size() > 0) {
      if (!(sinv.iscm() || sinv.isrm())) {
	HermMatrix<T,Lower,ColMajor> temp(sinv.size());
	DoCHInverse(LLx,temp.View());
	sinv = temp;
#ifdef LAP
      } else if (sinv.isrm()) {
	DoCHInverse(LLx.Transpose(),sinv.Transpose());
#endif
      } else {
	sinv = LLx;
#ifdef LAP
	LapCHInverse(sinv);
#else
	NonLapCHInverse(sinv);
#endif
      }
    }
  }

  template <class T> void HermCHDiv<T>::Inverse(
      const SymMatrixView<T>& sinv) const
  {
    TMVAssert(sinv.size() == LLx.size());
    TMVAssert(sinv.isherm());

    DoCHInverse(LLx,sinv);
  }

  template <class T> void HermCHDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LLx.size());
    TMVAssert(minv.rowsize() == LLx.size());

    if (IsComplex(T())) minv.diag().Imag().Zero();
    Inverse(HermMatrixViewOf(minv,Lower));
    UpperTriMatrixViewOf(minv).OffDiag() =
      LowerTriMatrixViewOf(minv).OffDiag().Adjoint();
  }

  template <class T> void Square(const MatrixView<T>& A)
  {
    const size_t N = A.colsize();
    if (N == 1) {
      const T A00 = *A.ptr();
      *A.ptr() = NORM(REAL(A00));
    } else {
      const size_t K = N/2;
      MatrixView<T> A00 = A.SubMatrix(0,K,0,K);
      MatrixView<T> A10 = A.SubMatrix(K,N,0,K);
      MatrixView<T> A01 = A.SubMatrix(0,K,K,N);
      MatrixView<T> A11 = A.SubMatrix(K,N,K,N);

      // [ A00 A10t ] [ A00 A10t ] 
      // [ A10 A11  ] [ A10 A11  ]
      // = [ A00^2 + A10t A10    A00 A10t + A10t A11 ]
      //   [ A10 A00 + A11 A10   A10 A10t + A11^2    ]

      // A10 stores the actual data for A10
      // We can therefore write to A01, sort of as a temp matrix.
      A01 = A00 * A10.Adjoint();
      A01 += A10.Adjoint() * A11;

      Square(A00);
      A00 += A10.Adjoint()*A10;
      Square(A11);
      A11 += A10*A10.Adjoint();

      A10 = A01.Adjoint();
    }
  }

  
  template <class T> void HermCHDiv<T>::InverseATA(
      const MatrixView<T>& ata) const
  {
    TMVAssert(ata.colsize() == LLx.size());
    TMVAssert(ata.rowsize() == LLx.size());
    // ata = (At A)^-1 = A^-1 (A^-1)t
    //     = A^-1 A^-1
    SymMatrixView<T> hermata = HermMatrixViewOf(ata,Lower);
    Inverse(hermata);
    Square(ata);
  }

#define InstFile "TMV_SymCHDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


