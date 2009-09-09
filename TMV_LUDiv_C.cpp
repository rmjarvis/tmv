
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

  template <class T> void NonLapLUInverse(const MatrixView<T>& minv)
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

#ifdef LAP
  template <class T> void LapLUInverse(const MatrixView<T>& minv)
  { NonLapLUInverse(minv); }
  template <> void LapLUInverse(const MatrixView<double>& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    dgetri(&n,minv.ptr(),&lda,LAP_IPiv(n),work,&lwork,&info);
    if (info < 0) tmv_error("dgetri returned info < 0");
  }
  template <> void LapLUInverse(const MatrixView<complex<double> >& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    zgetri(&n,LAP_Complex(minv.ptr()),&lda,LAP_IPiv(n),LAP_Complex(work),
	&lwork,&info);
    if (info < 0) tmv_error("zgetri returned info < 0");
  }
#ifndef NOFLOAT
  template <> void LapLUInverse(const MatrixView<float>& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    sgetri(&n,minv.ptr(),&lda,LAP_IPiv(n),work,&lwork,&info);
    if (info < 0) tmv_error("sgetri returned info < 0");
  }
  template <> void LapLUInverse(const MatrixView<complex<float> >& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    cgetri(&n,LAP_Complex(minv.ptr()),&lda,LAP_IPiv(n),LAP_Complex(work),
	&lwork,&info);
    if (info < 0) tmv_error("cgetri returned info < 0");
  }
#endif // FLOAT
#endif // LAP

  template <class T> void DoLUInverse(const GenMatrix<T>& LUx,
      const MatrixView<T>& minv, size_t* P)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.colsize() == LUx.colsize());
#ifdef XDEBUG
    Matrix<T> m = LowerTriMatrixViewOf(LUx,UnitDiag) *
      UpperTriMatrixViewOf(LUx);
    m.ReversePermuteRows(P);
#ifdef LAP
    Matrix<T,ColMajor> minv2(minv.colsize(),minv.colsize());
    minv2 = LUx;
    NonLapLUInverse(minv2.View());
    minv2.ReversePermuteCols(P);
#endif
#endif

    if (minv.colsize() > 0) {
#ifdef LAP
      if (!(minv.iscm())) 
#else
      if (!(minv.iscm() || minv.isrm())) 
#endif
      {
	Matrix<T,ColMajor> temp(minv.colsize(),minv.colsize());
	DoLUInverse(LUx,temp.View(),P);
	minv = temp;
      } else {
	minv = LUx;
#ifdef LAP
	LapLUInverse(minv);
#else
	NonLapLUInverse(minv);
#endif
	minv.ReversePermuteCols(P);
      }
    }

#ifdef XDEBUG
    RealType(T) normdiff = Norm(m*minv - RealType(T)(1));
    if (normdiff > 0.001*Norm(m)*Norm(minv)) {
      cerr<<"LUInverse:\n";
      cerr<<"m = "<<m<<endl;
      cerr<<"LUx = "<<LUx<<endl;
      cerr<<"minv = "<<minv<<endl;
      cerr<<"m*minv = "<<m*minv<<endl;
      cerr<<"minv*m = "<<minv*m<<endl;
#ifdef LAP
      cerr<<"Non-lap inverse = "<<minv2<<endl;
      cerr<<"m*minv2 = "<<m*minv2<<endl;
      cerr<<"minv2*m = "<<minv2*m<<endl;
#endif
      cerr<<"Norm(m*minv - 1) = "<<normdiff<<endl;
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
    if (istrans) DoLUInverse(LUx,minv.Transpose(),P);
    else DoLUInverse(LUx,minv,P);
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
      ata.ReversePermuteCols(P);
      ata.ReversePermuteRows(P);
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


