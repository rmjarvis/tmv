
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

  template <class T> inline void NonLapLUInverse(const MatrixView<T>& minv)
  {
    // m = P L U
    // m^-1 = U^-1 L^-1 Pt
    UpperTriMatrixView<T> U = UpperTriMatrixViewOf(minv);
    LowerTriMatrixView<T> L = LowerTriMatrixViewOf(minv,UnitDiag);
    U = U.Inverse();
    L = L.Inverse();
    minv = U*L;
    // Do Pt back in DoLUInverse
  }

#ifdef ALAP
  template <class T> inline void LapLUInverse(const MatrixView<T>& minv)
  { NonLapLUInverse(minv); }
  template <> inline void LapLUInverse(const MatrixView<double>& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    LAPNAME(dgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("dgetri");
  }
  template <> inline void LapLUInverse(const MatrixView<complex<double> >& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
#endif
    LAPNAME(zgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("zgetri");
  }
#ifndef NOFLOAT
  template <> inline void LapLUInverse(const MatrixView<float>& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    LAPNAME(sgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("sgetri");
  }
  template <> inline void LapLUInverse(const MatrixView<complex<float> >& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
#endif
    LAPNAME(cgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("cgetri");
  }
#endif // FLOAT
#endif // ALAP

  template <class T> inline void DoLUInverse(const GenMatrix<T>& LUx,
      const MatrixView<T>& minv, size_t* P)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.colsize() == LUx.colsize());
#ifdef XDEBUG
    Matrix<T> m = LowerTriMatrixViewOf(LUx,UnitDiag) *
      UpperTriMatrixViewOf(LUx);
    m.ReversePermuteRows(P);
#ifdef ALAP
    Matrix<T,ColMajor> minv2(minv.colsize(),minv.colsize());
    minv2 = LUx;
    NonLapLUInverse(minv2.View());
    minv2.ReversePermuteCols(P);
#endif
#endif

    if (minv.colsize() > 0) {
      if (!(minv.iscm()
#ifndef ALAP
	    || minv.isrm()
#endif
	   )) {
	Matrix<T,ColMajor> temp(minv.colsize(),minv.colsize());
	DoLUInverse(LUx,temp.View(),P);
	minv = temp;
      } else {
	minv = LUx;
#ifdef ALAP
	LapLUInverse(minv);
#else
	NonLapLUInverse(minv);
#endif
	minv.ReversePermuteCols(P);
      }
    }

#ifdef XDEBUG
    RealType(T) normdiff = Norm(m*minv - T(1));
    RealType(T) kappa = Norm(m)*Norm(minv);
    if (normdiff > 0.001*kappa*minv.colsize()) {
      cerr<<"LUInverse:\n";
      cerr<<"m = "<<m<<endl;
      cerr<<"LUx = "<<LUx<<endl;
      cerr<<"minv = "<<minv<<endl;
      cerr<<"m*minv = "<<m*minv<<endl;
      cerr<<"minv*m = "<<minv*m<<endl;
#ifdef ALAP
      cerr<<"Non-lap inverse = "<<minv2<<endl;
      cerr<<"m*minv2 = "<<m*minv2<<endl;
      cerr<<"minv2*m = "<<minv2*m<<endl;
#endif
      cerr<<"Norm(m*minv - 1) = "<<normdiff<<endl;
      cerr<<"kappa = "<<kappa<<endl;
      abort();
    }
#endif
  }

  template <class T> void LUDiv<T>::Inverse(const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LUx.colsize());
    TMVAssert(minv.rowsize() == LUx.colsize());
    // m = P L U
    // m^-1 = U^-1 L^-1 Pt
    if (istrans) DoLUInverse(LUx,minv.Transpose(),P.get());
    else DoLUInverse(LUx,minv,P.get());
  }

  template <class T> void LUDiv<T>::InverseATA(const MatrixView<T>& ata) const
  {
    TMVAssert(ata.colsize() == LUx.colsize());
    TMVAssert(ata.rowsize() == LUx.colsize());
    // (At A)^-1 = A^-1 (A^-1)t
    // = (U^-1 L^-1 Pt) (P L^-1t U^-1t)
    // = U^-1 L^-1 L^-1t U^-1t
    //
    // if PLU is really AT, then
    // A^-1 = P L^-1T U^-1T
    // (At A)^-1 = P L^-1T U^-1T U^-1* L^-1* Pt

    LowerTriMatrixView<T> L = LowerTriMatrixViewOf(LUx,UnitDiag);
    UpperTriMatrixView<T> U = UpperTriMatrixViewOf(LUx);

    if (istrans) {
      UpperTriMatrixView<T> uinv = UpperTriMatrixViewOf(ata);
      uinv = U.Inverse();
      ata = uinv.Transpose() * uinv.Conjugate();
      ata /= L.Transpose();
      ata %= L.Conjugate();
      ata.ReversePermuteCols(P.get());
      ata.ReversePermuteRows(P.get());
    } else {
      LowerTriMatrixView<T> linv = LowerTriMatrixViewOf(ata,UnitDiag);
      linv = L.Inverse();
      ata = linv * linv.Adjoint();
      ata /= U;
      ata %= U.Adjoint();
    }
  }

#define InstFile "TMV_LUDiv_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


