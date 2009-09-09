
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define QR_BLOCKSIZE TMV_BLOCKSIZE
#else
#define QR_BLOCKSIZE 64
#endif

  //
  // QR Downdate
  //

  template <class T> inline bool NonBlockQR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);
    // Given that A0 = Q0 R0
    // Given that [ A0 ] = Q1 R1
    //            [ A  ] 
    // Find R0 so that A0 = Q0 R0
    // Input R is R1, output is R0

    const size_t N = A.rowsize();

    T* Rdiag = R.ptr();
    const size_t ds = R.stepi()+R.stepj();

    for(size_t j=0;j<N;++j,Rdiag+=ds) {
      // Apply the Householder Reflection for this column
      const VectorView<T> v = A.col(j);
      T beta;
      if (!(Householder_UnReflect(*Rdiag,v,beta))) return false;
      TMVAssert(beta != T(1));
      VectorView<T> m0 = R.row(j,j+1,N);
      MatrixView<T> mx = A.Cols(j+1,N);

      // m0' = m0 - beta m0 - beta vtmx
      // m0 = (m0' + beta btmx)/(1-beta)
      Vector<T> bvtmx = beta*v.Conjugate()*mx;
      m0 += bvtmx;
      m0 /= T(1)-beta;

      // mx' = mx - beta v (m0 + vtmx)
      bvtmx += beta*m0;
      mx -= v ^ bvtmx;
    }
    return true;
  }

  template <class T> inline bool RecursiveQR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A,
      const UpperTriMatrixView<T>& Z, bool makeZ)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const size_t N = A.rowsize();

    if (N==1) {
      T b;
      if (!(Householder_UnReflect(R(0,0),A.col(0),b))) return false;
      Z(0,0) = CONJ(b);
    } else if (N==2) {
      T b0;
      if (!(Householder_UnReflect(R(0,0),A.col(0),b0))) return false;
      Z(0,0) = CONJ(b0);
      if (b0 != T(0)) {
	TMVAssert(b0 != T(1));
	T vtmx = A.col(0).Conjugate() * A.col(1);
	R(0,1) = (R(0,1) + b0*vtmx)/(T(1)-b0);
	A.col(1) -= b0*(vtmx + R(0,1)) * A.col(0);
      }

      T b1;
      if (!(Householder_UnReflect(R(1,1),A.col(1),b1))) return false;
      Z(1,1) = CONJ(b1);

      if (makeZ) {
	T vtmx = A.col(0).Conjugate() * A.col(1);
	Z(0,1) = -CONJ(b0*b1)*vtmx;
      }
    } else {
      size_t j1 = N/2;

      UpperTriMatrixView<T> R1 = R.SubTriMatrix(0,j1);
      MatrixView<T> Rx = R.SubMatrix(0,j1,j1,N);
      UpperTriMatrixView<T> R2 = R.SubTriMatrix(j1,N);

      MatrixView<T> A1 = A.Cols(0,j1);
      MatrixView<T> A2 = A.Cols(j1,N);

      UpperTriMatrixView<T> Z1 = Z.SubTriMatrix(0,j1);
      MatrixView<T> Zx = Z.SubMatrix(0,j1,j1,N);
      UpperTriMatrixView<T> Z2 = Z.SubTriMatrix(j1,N);

      if (!(RecursiveQR_Downdate(R1,A1,Z1,true))) return false;

      Zx = A1.Adjoint() * A2;
      Zx = Z1.Adjoint() * Zx;

      Rx += Zx;
      LowerTriMatrix<T> ImZt = T(1)-Z1.Adjoint();
      Rx /= ImZt;

      Zx += Z1.Adjoint() * Rx;
      A2 -= A1 * Zx;

      if (!(RecursiveQR_Downdate(R2,A2,Z2,makeZ))) return false;

      if (makeZ) {
	Zx = A1.Adjoint() * A2; 
	Zx = -Z1*Zx;
	Zx *= Z2;
      }
    }
    return true;
  }

  template <class T> inline bool BlockQR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const size_t N = A.rowsize();

    /*
    const size_t ds = R.stepi() + R.stepj();
    T* Rdiag = R.ptr();
    */

    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(min(size_t(QR_BLOCKSIZE),N));
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QR_BLOCKSIZE);
      MatrixView<T> A1 = A.Cols(j1,j2);
      UpperTriMatrixView<T> R1 = R.SubTriMatrix(j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,j2-j1);

      if (!(RecursiveQR_Downdate(R1,A1,Z,j2<N))) return false;

      /*
      for(size_t j=j1,jj=0;j<j2;++j,++jj,Rdiag+=ds) {
	const VectorView<T> v = A.col(j);
	T b;
	if (!(Householder_UnReflect(*Rdiag,v,b))) return false;
	if (b != T(0)) {
	  TMVAssert(b != T(1));
	  VectorView<T> m0 = R.row(j,j+1,j2);
	  MatrixView<T> mx = A.Cols(j+1,j2);

	  // m0' = m0 - beta m0 - beta vtmx
	  // m0 = (m0' + beta btmx)/(1-beta)
	  Vector<T> bvtmx = b*v.Conjugate()*mx;
	  m0 += bvtmx;
	  m0 /= T(1)-b;

	  // mx' = mx - beta v (m0 + vtmx)
	  bvtmx += b*m0;
	  mx -= v ^ bvtmx;

	  if (j2 < N) {
	    if (jj > 0) {
	      VectorView<T> z = Z.col(jj,0,jj);
	      z = A.Cols(j1,j).Adjoint() * v;
	      z = -b * Z.SubTriMatrix(0,jj) * z;
	    }
	    Z(jj,jj) = CONJ(b);
	  }
	} else {
	  Z.col(jj,0,jj+1).Zero();
	}
      }
      */

      if (j2 < N) {
	// m0' = m0 - Zt(Ytmx+m0)
	// m0' + ZtYtmx = (I-Zt) m0;
	MatrixView<T> m0 = R.SubMatrix(j1,j2,j2,N);
	MatrixView<T> Y = A.Cols(j1,j2);
	MatrixView<T> mx = A.Cols(j2,N);

	Matrix<T,ColMajor> ZtYtm = Y.Adjoint() * mx;
	ZtYtm = Z.Adjoint() * ZtYtm;

	m0 += ZtYtm;
	LowerTriMatrix<T> ImZt = T(1)-Z.Adjoint();
	m0 /= ImZt;

	ZtYtm += Z.Adjoint() * m0;
	mx -= Y * ZtYtm;
      }
      j1 = j2;
    }
    return true;
  }

  template <class T> bool QR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

#ifdef XDEBUG
    Matrix<T> R0 = R;
    Matrix<T> A0 = A;
    UpperTriMatrix<T> R2 = R;
    Matrix<T> A2 = A;
    bool ret2 = NonBlockQR_Downdate(R2.View(),A2.View());
#endif

    bool ret;
    if (A.rowsize() > 0) {
      if (A.rowsize() > QR_BLOCKSIZE)
	return BlockQR_Downdate(R,A);
      else {

	UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
	ret = RecursiveQR_Downdate(R,A,Z.View(),false);

	/*
	ret = NonBlockQR_Downdate(R,A);
	*/
      }
    }
    else ret = true;

#ifdef XDEBUG
    if (ret && (!ret2 || Norm(R2-R) > 1.e-5*Norm(A0)*Norm(R0)) ) {
      cerr<<"QR_Downdate\n";
      cerr<<"Succeed? = "<<ret<<endl;
      cerr<<"R0 = "<<Type(R)<<"  "<<R0<<endl;
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"R -> "<<R<<endl;
      cerr<<"NonBlock R -> "<<R2<<endl;
      abort();
    }
#endif

    return ret;
  }

#define InstFile "TMV_QRDiv_F.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


