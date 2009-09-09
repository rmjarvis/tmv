
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t TRI_DIV_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t TRI_DIV_BLOCKSIZE = 64;
#endif
 
  template <bool unit, class T> void RecursiveInverse(
      const UpperTriMatrixView<T>& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    TMVAssert(unit == m.isunit());

    const size_t N = m.size();
    const size_t nb = TRI_DIV_BLOCKSIZE;

    if (N == 1) {
      if (!unit) {
	T*const mptr = m.ptr();
	if (*mptr == T(0)) 
	  tmv_error("Singular Matrix found in Tri Inverse");
	*mptr = RealType(T)(1) / (*mptr);
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      UpperTriMatrixView<T> m00 = m.SubTriMatrix(0,k);
      MatrixView<T> m01 = m.SubMatrix(0,k,k,N);
      UpperTriMatrixView<T> m11 = m.SubTriMatrix(k,N);

      // m00 m01' + m01 m11' = 0
      // m00 m01' = -m01 m11'
      // m01' = -m00' m01 m11'

      RecursiveInverse<unit>(m00);
      RecursiveInverse<unit>(m11);
      m01 = -m00 * m01;
      m01 *= m11;
    }
  }

  template <class T> void NonLapInverse(const UpperTriMatrixView<T>& m)
  {
    if (m.isunit()) RecursiveInverse<true>(m);
    else RecursiveInverse<false>(m);
  }

#ifdef ALAP
  template <class T> void LapInverse(const UpperTriMatrixView<T>& m)
  { NonLapInverse(m); }
  template <> void LapInverse(const UpperTriMatrixView<double>& m)
  {
#ifdef CLAP
    clapack_dtrtri(m.iscm()?CblasColMajor:CblasRowMajor,
	CblasUpper, m.isunit()?CblasUnit:CblasNonUnit,
	m.size(),m.ptr(),m.iscm()?m.stepj():m.stepi());
#else
    char uplo = m.iscm() ? 'U' : 'L';
    char diag = m.isunit() ? 'U' : 'N';
    int n = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
    int info;
    dtrtri(&uplo,&diag,&n,m.ptr(),&lda,&info);
    if (info < 0) tmv_error("dtrtri returned info < 0");
#endif
  }
  template <> void LapInverse(const UpperTriMatrixView<complex<double> >& m)
  {
#ifdef CLAP
    clapack_ztrtri(m.iscm()?CblasColMajor:CblasRowMajor,
	CblasUpper, m.isunit()?CblasUnit:CblasNonUnit,
	m.size(),m.ptr(),m.iscm()?m.stepj():m.stepi());
#else
    char uplo = m.iscm() ? 'U' : 'L';
    char diag = m.isunit() ? 'U' : 'N';
    int n = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
    int info;
    ztrtri(&uplo,&diag,&n,LAP_Complex(m.ptr()),&lda,&info);
    if (info < 0) tmv_error("ztrtri returned info < 0");
#endif
  }
#ifndef NOFLOAT
  template <> void LapInverse(const UpperTriMatrixView<float>& m)
  {
#ifdef CLAP
    clapack_strtri(m.iscm()?CblasColMajor:CblasRowMajor,
	CblasUpper, m.isunit()?CblasUnit:CblasNonUnit,
	m.size(),m.ptr(),m.iscm()?m.stepj():m.stepi());
#else
    char uplo = m.iscm() ? 'U' : 'L';
    char diag = m.isunit() ? 'U' : 'N';
    int n = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
    int info;
    strtri(&uplo,&diag,&n,m.ptr(),&lda,&info);
    if (info < 0) tmv_error("strtri returned info < 0");
#endif
  }
  template <> void LapInverse(const UpperTriMatrixView<complex<float> >& m)
  {
#ifdef CLAP
    clapack_ctrtri(m.iscm()?CblasColMajor:CblasRowMajor,
	CblasUpper, m.isunit()?CblasUnit:CblasNonUnit,
	m.size(),m.ptr(),m.iscm()?m.stepj():m.stepi());
#else
    char uplo = m.iscm() ? 'U' : 'L';
    char diag = m.isunit() ? 'U' : 'N';
    int n = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
    int info;
    ctrtri(&uplo,&diag,&n,LAP_Complex(m.ptr()),&lda,&info);
    if (info < 0) tmv_error("ctrtri returned info < 0");
#endif
  }
#endif // FLOAT
#endif // ALAP

  template <class T> void DoInverse(const GenUpperTriMatrix<T>& m,
      const UpperTriMatrixView<T>& minv)
  {
    TMVAssert(minv.size() == m.size());
    TMVAssert(!minv.isunit() || m.isunit());

    if (m.size() > 0) {
      if (!(minv.iscm() || minv.isrm())) {
	UpperTriMatrix<T> temp(m.size());
	DoInverse(m,temp.View());
	minv = temp;
      } else {
	minv = m;
#ifdef ALAP
	LapInverse(minv);
#else
	NonLapInverse(minv);
#endif
      }
    }
  }

  template <class T> void UpperTriDiv<T>::Inverse(
      const UpperTriMatrixView<T>& minv) const
  {
    TMVAssert(minv.size() == itsm->size());
    TMVAssert(!minv.isunit() || itsm->isunit());
#ifdef XDEBUG
    UpperTriMatrix<T> m0 = *itsm;
#endif

    DoInverse(*itsm,minv);

#ifdef XDEBUG
    Matrix<T> eye = minv*m0;
    if (Norm(eye-T(1)) > 0.0001*(Norm(m0)+Norm(minv))) {
      cerr<<"UpperTriMatrix Inverse:\n";
      cerr<<"m = "<<tmv::Type(*itsm)<<"  "<<m0<<endl;
      cerr<<"minv = "<<minv<<endl;
      cerr<<"minv*m = "<<minv*m0<<endl;
      cerr<<"m*minv = "<<m0*minv<<endl;
      abort();
    }
#endif
  }

  template <class T> void LowerTriDiv<T>::Inverse(
      const LowerTriMatrixView<T>& minv) const
  {
    TMVAssert(minv.size() == itsm->size());
    TMVAssert(!minv.isunit() || itsm->isunit());
#ifdef XDEBUG
    LowerTriMatrix<T> m0 = *itsm;
#endif

    DoInverse(itsm->Transpose(),minv.Transpose());

#ifdef XDEBUG
    Matrix<T> eye = minv*m0;
    if (Norm(eye-T(1)) > 0.0001*(Norm(m0)+Norm(minv))) {
      cerr<<"LowerTriMatrix Inverse:\n";
      cerr<<"m = "<<tmv::Type(*itsm)<<"  "<<m0<<endl;
      cerr<<"minv = "<<minv<<endl;
      cerr<<"minv*m = "<<minv*m0<<endl;
      cerr<<"m*minv = "<<m0*minv<<endl;
      abort();
    }
#endif
  }

  template <class T> void UpperTriDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == itsm->size());
    TMVAssert(minv.rowsize() == itsm->size());

    UpperTriMatrixView<T> U = UpperTriMatrixViewOf(minv);
    DoInverse(*itsm,U);
    minv = U * U.Adjoint();

    /*
    UpperTriMatrix<T,NonUnitDiag,ColMajor> temp(itsm->size());
    Inverse(temp.View());
    minv = temp*Adjoint(temp);
    */
  }

  template <class T> void LowerTriDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == itsm->size());
    TMVAssert(minv.rowsize() == itsm->size());

    LowerTriMatrixView<T> L = LowerTriMatrixViewOf(minv);
    DoInverse(itsm->Transpose(),L.Transpose());
    minv = L * L.Adjoint();

    /*
    LowerTriMatrix<T,NonUnitDiag,ColMajor> temp(itsm->size());
    Inverse(temp.View());
    minv = temp*Adjoint(temp);
    */
  }

#define InstFile "TMV_TriDiv_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


