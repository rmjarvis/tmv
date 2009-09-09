
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_SymCHDiv_A.h"

//#define XDEBUG

namespace tmv {

  template <class T> inline void NonLapCHInverse(const SymMatrixView<T>& sinv)
  {
    TMVAssert(sinv.isherm());
    // inv = (L Lt)^-1 = Lt^-1 L^-1
    LowerTriMatrixView<T> L = sinv.LowerTri();
    L = L.Inverse();
    sinv = L.Adjoint() * L;
  }

#ifdef ALAP
  template <class T> inline void LapCHInverse(const SymMatrixView<T>& sinv)
  { NonLapCHInverse(sinv); }
  template <> inline void LapCHInverse(const SymMatrixView<double>& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    int n = sinv.size();
    int lda = sinv.stepj();
    LAPNAME(dpotri) (LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
	LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("dpotri");
  }
  template <> inline void LapCHInverse(
      const SymMatrixView<complex<double> >& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    int n = sinv.size();
    int lda = sinv.stepj();
    LAPNAME(zpotri) (LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
	LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("zpotri");
  }
#ifndef NOFLOAT
  template <> inline void LapCHInverse(const SymMatrixView<float>& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    int n = sinv.size();
    int lda = sinv.stepj();
    LAPNAME(spotri) (LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
	LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("spotri");
  }
  template <> inline void LapCHInverse(
      const SymMatrixView<complex<float> >& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    int n = sinv.size();
    int lda = sinv.stepj();
    LAPNAME(cpotri) (LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
	LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("cpotri");
  }
#endif // FLOAT
#endif // LAP

  template <class T> inline void DoCHInverse(const GenSymMatrix<T>& LLx,
      const SymMatrixView<T>& sinv) 
  {
    if (sinv.size() > 0) {
      if (!(sinv.iscm() || sinv.isrm())) {
	HermMatrix<T,Lower,ColMajor> temp(sinv.size());
	DoCHInverse(LLx,temp.View());
	sinv = temp;
#ifdef ALAP
      } else if (sinv.isrm()) {
	DoCHInverse(LLx.Transpose(),sinv.Transpose());
#endif
      } else {
	sinv = LLx;
#ifdef ALAP
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
#ifdef XDEBUG
    Matrix<T> A = LLx.LowerTri() * LLx.UpperTri();
#endif

    DoCHInverse(LLx,sinv);

#ifdef XDEBUG
    Matrix<T> eye = A * sinv;
    RealType(T) kappa = Norm(A) * Norm(sinv);
    if (Norm(eye-T(1)) > 0.0001*kappa*sinv.size()) {
      cerr<<"A = "<<A<<endl;
      cerr<<"sinv = "<<sinv<<endl;
      cerr<<"A*sinv = "<<A*sinv<<endl;
      cerr<<"sinv*A = "<<sinv*A<<endl;
      cerr<<"Norm(A*sinv-1) = "<<Norm(A*sinv-T(1))<<endl;
      cerr<<"kappa = "<<kappa<<endl;
      abort();
    }
#endif
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

  template <bool herm, class T> void SymATASquare(const MatrixView<T>& A)
  {
    const size_t N = A.colsize();
    if (N == 1) {
      const T A00 = *A.ptr();
      if (herm)
	*A.ptr() = NORM(REAL(A00));
      else 
	*A.ptr() = SQR(A00);
    } else {
      const size_t K = N/2;
      MatrixView<T> A00 = A.SubMatrix(0,K,0,K);
      MatrixView<T> A10 = A.SubMatrix(K,N,0,K);
      MatrixView<T> A01 = A.SubMatrix(0,K,K,N);
      MatrixView<T> A11 = A.SubMatrix(K,N,K,N);
      MatrixView<T> A10t = herm ? A10.Adjoint() : A10.Transpose();

      // [ A00 A10t ] [ A00 A10t ] 
      // [ A10 A11  ] [ A10 A11  ]
      // = [ A00^2 + A10t A10    A00 A10t + A10t A11 ]
      //   [ A10 A00 + A11 A10   A10 A10t + A11^2    ]

      // A10 stores the actual data for A10
      // We can therefore write to A01, sort of as a temp matrix.
      
      A01 = A00 * A10t;
      A01 += A10t * A11;

      SymATASquare<herm>(A00);
      A00 += A10t*A10;
      SymATASquare<herm>(A11);
      A11 += A10*A10t;

      A10t = A01;
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
    SymATASquare<true>(ata);
  }

#define InstFile "TMV_SymCHDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


