
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
  // QR Update
  //

  template <class T> inline void NonBlockQR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);
    // Given that A0 = Q0 R0
    // Find R1, so that [ A0 ] = Q1 R1
    //                  [ A  ] 
    // Input R is R0, output is R1

    const size_t N = A.rowsize();

    T* Rdiag = R.ptr();
    const size_t ds = R.stepi()+R.stepj();
    T det(0);

    for(size_t j=0;j<N;++j,Rdiag+=ds) {
      // Apply the Householder Reflection for this column
      const VectorView<T> v = A.col(j);
      T beta = Householder_Reflect(*Rdiag,v,det);
      if (beta != T(0))
	Householder_LMult(v,beta,R.row(j,j+1,N),A.Cols(j+1,N));
    }
  }

  template <class T> inline void RecursiveQR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A,
      const UpperTriMatrixView<T>& Z, bool makeZ)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const size_t N = A.rowsize();
    T det(0);

    if (N==1) {
      T b = Householder_Reflect(R(0,0),A.col(0),det);
      Z(0,0) = CONJ(b);
    } else if (N==2) {
      T b0 = Householder_Reflect(R(0,0),A.col(0),det);
      if (b0 != T(0)) {
	T temp = b0*(A.col(0).Conjugate()*A.col(1) + R(0,1));
	R(0,1) -= temp;
	A.col(1) -= temp * A.col(0);
      }
      Z(0,0) = CONJ(b0);
      T b1 = Householder_Reflect(R(1,1),A.col(1),det);
      Z(1,1) = CONJ(b1);

      if (makeZ) {
	T temp = A.col(0).Conjugate()*A.col(1);
	Z(0,1) = -Z(0,0)*Z(1,1)*temp;
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

      RecursiveQR_Update(R1,A1,Z1,true);

      // Zx is a temporary here - it happens to be the right shape.
      Zx = A1.Adjoint() * A2; 
      Zx += Rx;
      Zx = Z1.Adjoint()*Zx;
      Rx -= Zx;
      A2 -= A1 * Zx;

      RecursiveQR_Update(R2,A2,Z2,makeZ);

      if (makeZ) {
	Zx = A1.Adjoint() * A2; 
	Zx = -Z1*Zx;
	Zx *= Z2;
      }
    }
  }

  template <class T> inline void BlockQR_Update(
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

      RecursiveQR_Update(R1,A1,Z,j2<N);

      /*
      T det(0);
      for(size_t j=j1,jj=0;j<j2;++j,++jj,Rdiag+=ds) {
	const VectorView<T> v = A.col(j);
	T b = Householder_Reflect(*Rdiag,v,det);
	if (b != T(0)) {
	  Householder_LMult(v,b,R.row(j,j+1,j2),A.Cols(j+1,j2));
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
	Matrix<T,ColMajor> ZtYtm = A.Cols(j1,j2).Adjoint() * A.Cols(j2,N);
	ZtYtm += R.SubMatrix(j1,j2,j2,N);
	ZtYtm = Z.Adjoint() * ZtYtm;
	R.SubMatrix(j1,j2,j2,N) -= ZtYtm;
	A.Cols(j2,N) -= A.Cols(j1,j2) * ZtYtm;
      }
      j1 = j2;
    }
  }

  template <class T> void QR_Update(
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
    NonBlockQR_Update(R2.View(),A2.View());
#endif
    if (A.rowsize() > 0)
      if (A.rowsize() > QR_BLOCKSIZE)
	BlockQR_Update(R,A);
      else {

	UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
	RecursiveQR_Update(R,A,Z.View(),false);

	/*
	NonBlockQR_Update(R,A);
	*/
      }
#ifdef XDEBUG
    if (Norm(R2-R) > 1.e-5*Norm(R0)*Norm(A0)) {
      cerr<<"QR_Update\n";
      cerr<<"R0 = "<<Type(R)<<"  "<<R0<<endl;
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"R -> "<<R<<endl;
      cerr<<"NonBlock R -> "<<R2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_QRDiv_E.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


