///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_Matrix.h"
#include "TMV_Givens.h"
#include "TMV_DiagMatrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_Householder.h"
#include "TMV_QRDiv.h"
#include "TMV_SVDiv.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"

#define XDEBUG

#ifdef XDEBUG
#include "TMV_DiagMatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define HESS_BLOCKSIZE TMV_BLOCKSIZE
#else
#define HESS_BLOCKSIZE 16
#endif

  template <class T> inline const MatrixView<T>* ZMV() 
  { return (const MatrixView<T>*)(0); }

  //
  // Reduce Matrix to Hessenberg Form (upper tri with one lower sub-diag)
  //
  
  template <class T> inline void NonBlockHessenberg(
      const MatrixView<T>& A, const VectorView<T>& Ubeta)
  {
#ifdef XDEBUG
    cerr<<"Start NonBlock Hessenberg Reduction: A = "<<A<<endl;
    Matrix<T> A0(A);
#endif
    // Decompose A into U H Ut
    // H is a Hessenberg Matrix
    // U is a Unitary Matrix
    // On output, H is stored in the upper-Hessenberg part of A
    // U is stored in compact form in the rest of A along with 
    // the vector Ubeta.
    const size_t N = A.rowsize();

    TMVAssert(A.colsize() == A.rowsize());
    TMVAssert(N > 0);
    TMVAssert(Ubeta.size() == N-1);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(!Ubeta.isconj());
    TMVAssert(Ubeta.step()==1);

    // We use Householder reflections to reduce A to the Hessenberg form:
    T* Uj = Ubeta.ptr();
    T det = 0; // Ignore Householder det calculations
    for(size_t j=0;j<N-1;++j,++Uj) {
      *Uj = Householder_Reflect(A.SubMatrix(j+1,N,j,N),det);
      if (*Uj != T(0))
	Householder_LMult(A.col(j+2,N),*Uj,A.SubMatrix(0,N,j+1,N).Adjoint());
    }

#ifdef XDEBUG
    Matrix<T> U(N,N,T(0));
    U.SubMatrix(1,N,1,N) = A.SubMatrix(1,N,0,N-1);
    U.UpperTri().Zero();
    Vector<T> Ubeta2(N);
    Ubeta2.SubVector(1,N) = Ubeta;
    Ubeta2(0) = T(0);
    GetQFromQR(U.View(),Ubeta2);
    Matrix<T> H = A;
    if (N>2) LowerTriMatrixViewOf(H).OffDiag().OffDiag().Zero();
    Matrix<T> AA = U*H*U.Adjoint();
    if (Norm(A0-AA) > 0.001*Norm(A0)) {
      cerr<<"NonBlock Hessenberg: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"A = "<<A<<endl;
      cerr<<"Ubeta = "<<Ubeta<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"H = "<<H<<endl;
      cerr<<"UHUt = "<<AA<<endl;
      abort();
    }
#endif
  }

  template <class T> inline void BlockHessenberg(
      const MatrixView<T>& A, const VectorView<T>& Ubeta)
  {
    // Much like the block version of Bidiagonalize, we try to maintain
    // the operation of several successive Householder matrices in
    // a block form, where the net Block Householder is I - YZYt.
    //
    // But as with the bidiagonlization algorithm (and unlike a simple
    // block QR decomposition), we update the matrix from both the left 
    // and the right, so we also need to keep track of the product
    // ZYtm in addition.
    //
    // The block update at the end of the block loop is
    // m' = (I-YZYt) m (I-YZtYt)
    //
    // The Y matrix is stored in the first K columns of m,
    // and the Hessenberg portion of these columns is updated as we go.
    // For the right-hand-side update, m -= mYZtYt, the m on the right
    // needs to be the full original matrix m, including the original
    // versions of these K columns.  Therefore, we can't wait until 
    // the end for this calculation.  
    //
    // Instead, we keep track of mYZt as we progress, so the final update
    // is:
    //
    // m' = (I-YZYt) (m - mYZt Y)
    //
    // We also need to do this same calculation for each column as we
    // progress through the block.
    //
    const size_t N = A.rowsize();

#ifdef XDEBUG
    Matrix<T> A0(A);
#endif

    TMVAssert(A.rowsize() == A.colsize());
    TMVAssert(N > 0);
    TMVAssert(Ubeta.size() == N-1);
    TMVAssert(!Ubeta.isconj());
    TMVAssert(Ubeta.step()==1);

    size_t ncolmax = std::min(size_t(HESS_BLOCKSIZE),N-1);
    Matrix<T,RowMajor> mYZt_full(N,ncolmax);
    UpperTriMatrix<T,NonUnitDiag,ColMajor> Z_full(ncolmax);

    T det(0); // Ignore Householder Determinant calculations
    T* Uj = Ubeta.ptr();
    for(size_t j1=0;j1<N-1;) {
      size_t j2 = std::min(N-1,j1+HESS_BLOCKSIZE);
      size_t ncols = j2-j1;
      MatrixView<T> mYZt = mYZt_full.SubMatrix(0,N-j1,0,ncols);
      UpperTriMatrixView<T> Z = Z_full.SubTriMatrix(0,ncols);

      for(size_t j=j1,jj=0;j<j2;++j,++jj,++Uj) { // jj = j-j1

	// Update current column of A
	//
	// m' = (I - YZYt) (m - mYZt Yt)
	// A(0:N,j)' = A(0:N,j) - mYZt(0:N,0:j) Y(j,0:j)t
	A.col(j,j1+1,N) -= mYZt.Cols(0,j) * A.row(j,0,j).Conjugate();
	//
	// A(0:N,j)'' = A(0:N,j) - Y Z Yt A(0:N,j)'
	// 
	// Let Y = (L)     where L is unit-diagonal, lower-triangular,
	//         (M)     and M is rectangular
	//
	LowerTriMatrixView<T> L = 
	  LowerTriMatrixViewOf(A.SubMatrix(j1+1,j+1,j1,j),UnitDiag);
	MatrixView<T> M = A.SubMatrix(j+1,N,j1,j);
	// Use the last column of Z as temporary storage for Yt A(0:N,j)'
	VectorView<T> YtAj = Z.col(jj,0,jj);
	YtAj = L.Adjoint() * A.col(j,j1+1,j+1);
	YtAj += M.Adjoint() * A.col(j,j+1,N);
	YtAj = Z.SubTriMatrix(0,jj) * YtAj;
	A.col(j,j1+1,j+1) -= L * YtAj;
	A.col(j,j+1,N) -= M * YtAj;

	// Do the Householder reflection 
	VectorView<T> u = A.col(j,j+1,N);
	T bu = Householder_Reflect(u,det);
	*Uj = bu;

	// Save the top of the u vector, which isn't actually part of u
	T& Atemp = *u.cptr();
	TMVAssert(IMAG(Atemp) == RealType(T)(0));
	RealType(T) Aorig = REAL(Atemp);
	Atemp = RealType(T)(1);

	// Update Z
	VectorView<T> Zj = Z.col(jj,0,jj);
	Zj = -bu * M.Adjoint() * u;
	Zj = Z * Zj;
	Z(jj,jj) = -bu;

	// Update mYtZt:
	//
	// mYZt(0:N,j) = m(0:N,0:N) Y(0:N,0:j) Zt(0:j,j)
	//             = m(0:N,j+1:N) Y(j+1:N,j) Zt(j,j)
	//             = bu* m(0:N,j+1:N) u 
	//
	mYZt.col(jj) = CONJ(bu) * A.SubMatrix(j1,N,j+1,N) * u;

	// Restore Aorig, which is actually part of the Hessenberg matrix.
	Atemp = Aorig;
      }

      // Update the rest of the matrix:
      // A(j2,j2-1) needs to be temporarily changed to 1 for use in Y
      T& Atemp = *(A.ptr() + j2*A.stepi() + (j2-1)*A.stepj());
      TMVAssert(IMAG(Atemp) == RealType(T)(0));
      RealType(T) Aorig = Atemp;
      Atemp = RealType(T)(1);

      // m' = (I-YZYt) (m - mYZt Y)
      MatrixView<T> m = A.SubMatrix(j1,N,j2,N);
      ConstMatrixView<T> Y = A.SubMatrix(j2+1,N,j1,j2);
      m -= mYZt * Y.Adjoint();
      BlockHouseholder_LMult(Y,Z,m);

      // Restore A(j2,j2-1)
      Atemp = Aorig;
      j1 = j2;
    }

#ifdef XDEBUG
    Matrix<T> U(N,N,T(0));
    U.SubMatrix(1,N,1,N) = A.SubMatrix(1,N,0,N-1);
    U.UpperTri().Zero();
    U(0,0) = T(1);
    Vector<T> Ubeta2(N);
    Ubeta2.SubVector(1,N) = Ubeta;
    Ubeta2(0) = T(0);
    GetQFromQR(U.View(),Ubeta2);
    Matrix<T> H = A;
    if (N>2) LowerTriMatrixViewOf(H).OffDiag().OffDiag().Zero();
    Matrix<T> AA = U*H*U.Adjoint();
    if (Norm(A0-AA) > 0.001*Norm(A0)) {
      cerr<<"NonBlock Hessenberg: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"A = "<<A<<endl;
      cerr<<"Ubeta = "<<Ubeta<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"H = "<<H<<endl;
      cerr<<"UHUt = "<<AA<<endl;
      Matrix<T,ColMajor> A2 = A0;
      Vector<T> Ub2(Ubeta.size());
      NonBlockHessenberg(A2.View(),Ub2.View());
      cerr<<"cf NonBlock: A -> "<<A2<<endl;
      cerr<<"Ubeta = "<<Ub2<<endl;
      abort();
    }
#endif
  }

  template <class T> inline void NonLapHessenberg(
      const MatrixView<T>& A, const VectorView<T>& Ubeta)
  {
    TMVAssert(A.rowsize() == A.colsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(Ubeta.size() == A.rowsize()-1);

#if 0
    if (A.rowsize() > HESS_BLOCKSIZE)
      BlockHessenberg(A,Ubeta,Vbeta,D,E,det);
    else
#endif
      NonBlockHessenberg(A,Ubeta);
  }

#ifdef LAP
  template <class T> inline void LapHessenberg(
      const MatrixView<T>& A, const VectorView<T>& Ubeta)
  { NonLapHessenberg(A,Ubeta,Vbeta,D,E,det); }
#ifdef INST_DOUBLE
  template <> inline void LapHessenberg(
      const MatrixView<double>& A, const VectorView<double>& Ubeta)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.colsize() == A.rowsize());
    TMVAssert(Ubeta.size() == A.rowsize()-1);
    TMVAssert(A.ct()==NonConj);

    int n = A.rowsize();
    int ilo = 1;
    int ihi = n;
    int lda = A.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    LAPNAME(dgehrd) (LAPCM LAPV(n),LAPV(ilo),LAPV(ihi),
	LAPP(A.ptr()),LAPV(lda),LAPP(Ubeta.ptr())
	LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
    LAP_Results("dgehrd");
#else
    LAP_Results(int(work[0]),m,n,lwork,"dgehrd");
#endif
  }
  template <> inline void LapHessenberg(
      const MatrixView<std::complex<double> >& A, 
      const VectorView<std::complex<double> >& Ubeta)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.colsize() == A.rowsize());
    TMVAssert(Ubeta.size() == A.rowsize()-1);
    TMVAssert(A.ct()==NonConj);

    int n = A.rowsize();
    int ilo = 1;
    int ihi = n;
    int lda = A.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    std::complex<double>* work = LAP_ZWork(lwork);
#endif
    LAPNAME(zgehrd) (LAPCM LAPV(n),LAPV(ilo),LAPV(ihi),
	LAPP(A.ptr()),LAPV(lda),LAPP(Ubeta.ptr())
	LAPWK(work) LAPVWK(lwork) LAPINFO);
    Ubeta.ConjugateSelf();
#ifdef LAPNOWORK
    LAP_Results("zgehrd");
#else
    LAP_Results(int(REAL(work[0])),m,n,lwork,"zgehrd");
#endif
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapHessenberg(
      const MatrixView<float>& A, const VectorView<float>& Ubeta)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.colsize() == A.rowsize());
    TMVAssert(Ubeta.size() == A.rowsize()-1);
    TMVAssert(A.ct()==NonConj);

    int n = A.rowsize();
    int ilo = 1;
    int ihi = n;
    int lda = A.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    LAPNAME(sgebrd) (LAPCM LAPV(n),LAPV(ilo),LAPV(ihi),
	LAPP(A.ptr()),LAPV(lda),LAPP(Ubeta.ptr())
	LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
    LAP_Results("sgehrd");
#else
    LAP_Results(int(work[0]),m,n,lwork,"sgehrd");
#endif
  }
  template <> inline void LapHessenberg(
      const MatrixView<std::complex<float> >& A, 
      const VectorView<std::complex<float> >& Ubeta) 
  {
    TMVAssert(A.iscm());
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Ubeta.size() == A.rowsize());
    TMVAssert(A.ct()==NonConj);

    int n = A.rowsize();
    int ilo = 1;
    int ihi = n;
    int lda = A.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    std::complex<float>* work = LAP_CWork(lwork);
#endif
    LAPNAME(cgehrd) (LAPCM LAPV(n),LAPV(ilo),LAPV(ihi),
	LAPP(A.ptr()),LAPV(lda),LAPP(Ubeta.ptr())
	LAPWK(work) LAPVWK(lwork) LAPINFO);
    Ubeta.ConjugateSelf();
#ifdef LAPNOWORK
    LAP_Results("cgehrd");
#else
    LAP_Results(int(REAL(work[0])),m,n,lwork,"cgehrd");
#endif
  }
#endif 
#endif // LAP

  template <class T> inline void Hessenberg(
      const MatrixView<T>& A, const VectorView<T>& Ubeta)
  {
    TMVAssert(A.colsize() == A.rowsize());
    TMVAssert(Ubeta.size() == A.rowsize()-1);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(Ubeta.step() == 1);
    
    if (A.rowsize() > 0) {
#ifdef LAP
      if (A.iscm()) 
	LapHessenberg(A,Ubeta);
      else 
#endif
	NonLapHessenberg(A,Ubeta);
    }
  }

  /*
  //
  // Schur_From_Hessenberg: QR method
  //
  
  // MJ: Stopped here
  template <class T> inline void BidiagonalChopSmallElements(
      const VectorView<T>& D, const VectorView<T>& E)
  {
    // This routines sets to 0 any elements in E or D which
    // are essentially 0, given the machine precision:
    // if |E(i)| < Epsilon() * (|D(i)| + |D(i+1)|) then E(i) <- 0
    // if |D(i)| < Epsilon() * |B| then D(i) <- 0
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    RealType(T) eps = Epsilon<T>();
    RealType(T) dthresh = eps * SQRT(NormSq(D) + NormSq(E));

    T* Di = D.ptr();
    T* Ei = E.ptr();
    RealType(T) absDim1 = ABS(*Di);
    if (absDim1 < dthresh) *Di = T(0);
    ++Di;
    for(size_t k=E.size();k>0;--k,++Di,++Ei) {
      RealType(T) absDi = ABS(*Di);
      if (absDi < dthresh) *Di = T(0);
      RealType(T) ethresh = eps * (absDi + absDim1);
      if (ethresh == RealType(T)(0)) ethresh = dthresh;
      if (ABS(*Ei) < ethresh) *Ei = T(0);
      absDim1 = absDi;
    }
  }

  template <class T, class TU> inline void BidiagonalZeroFirstRow(
      const MatrixView<TU>* U, const VectorView<T>& D, 
      const VectorView<T>& E)
  {
    // Input D,E form a bidiagonal matrix with the first element of D = 0:
    // (eg. for N = 5)
    //     [ 0 x 0 0 0 ]
    //     [ 0 x x 0 0 ]
    // B = [ 0 0 x x 0 ]
    //     [ 0 0 0 x x ]
    //     [ 0 0 0 0 x ]
    // Zero out the first row maintaining the constancy of U B
    // using Givens transformations.
    const size_t N = D.size();
    if (N <= 1) return; 
    TMVAssert(E.size() == N-1);
    if (U) TMVAssert(U->rowsize() == N);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);

    T* Di = D.ptr();
    T* Ei = E.ptr();
    TMVAssert(*Di == T(0));

    T x = *Ei;
    if (x != T(0)) {
      *Ei = T(0);
      ++Ei; ++Di;
      // Loop Invariant: x = B(0,i)
      for(size_t i=1; i<N; ++i,++Di,++Ei) {
	Givens<T> G = Givens_Rotate(*Di,x);
	// Make new B = G B
	if (i<N-1) G.Mult(*Ei,x);
	// Make new U = U Gt
	if (U) G.ConjMult(U->ColPair(i,0).Transpose());
      }
    }
  }

  template <class T, class TV> inline void BidiagonalZeroLastCol(
      const VectorView<T>& D, const VectorView<T>& E,
      const MatrixView<TV>* V)
  {
    // Input D,E form a bidiagonal matrix with the last element of D = 0:
    // (eg. for N = 5)
    //     [ x x 0 0 0 ]
    //     [ 0 x x 0 0 ]
    // B = [ 0 0 x x 0 ]
    //     [ 0 0 0 x x ]
    //     [ 0 0 0 0 0 ]
    // Zero out the last col maintaining the constancy of B V
    // using Givens transformations.
    const size_t N = D.size();
    if (N <= 1) return; 
    TMVAssert(E.size() == N-1);
    if (V) TMVAssert(V->colsize() == N);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(D(N-1) == T(0));

    T* Di = D.ptr()+N-2;
    T* Ei = E.ptr()+N-2;

    T x = *Ei;
    if (x != T(0)) {
      *Ei = T(0);
      // Loop Invariant: x = B(i,N-1)
      for(int i=N-2; i>=0; --i,--Di) {
	Givens<T> G = Givens_Rotate(*Di,x);
	// Make new B = B GT
	if (i>0) G.Mult(*(--Ei),x);
	// Make new V = G* V 
	if (V) G.ConjMult(V->RowPair(i,N-1));
      }
    }
  }

  template <class T> inline RealType(T) BidiagonalTrailingEigenValue(
      const VectorView<T>& D, const VectorView<T>& E)
  {
    // Return the Wilkinson choice for an eigenvalue of T = BtB, namely
    // the eigenvalue of the trailing 2x2 block of T which is closer
    // to the last diagonal element of T.
    //
    // Trailing 2x2 block =  (Use i = N-2, j = N-1)
    // [ a  b ] = [ |Di|^2 + |Ei-1|^2      Di Ei      ]
    // [ b* c ]   [     Di* Ei*       |Dj|^2 + |Ei|^2 ]
    // 
    // mu = c - d +- sqrt(d^2+|b|^2), where d = (c-a)/2
    // if d>0 we use +, if d<0 we use -.
    // 
    // For stability when |b| is small, we rearrange this to:
    // mu = c + |b|^2/(d +- sqrt(d^2+|b|^2))
    //    = c +- |b|^2/|d|(1 + sqrt(1+|b|^2/d^2))
    const size_t N = D.size();
    TMVAssert(E.size() == N-1);
    TMVAssert(N > 1);

    RealType(T) a = NORM(D(N-2)) + (N>2 ? NORM(E(N-3)) : RealType(T)(0));
    RealType(T) c = NORM(D(N-1)) + NORM(E(N-2));
    T b = D(N-2)*E(N-2);
    RealType(T) bsq = NORM(b);
    RealType(T) d = (c-a)/2;
    RealType(T) dsq = SQR(d);
    if (dsq < bsq) {
      RealType(T) x = SQRT(dsq+bsq);
      if (d > 0) return c - d + x;
      else return c - d - x;
    } else {
      RealType(T) x = bsq/ABS(d)/(1 + SQRT(1+bsq/dsq));
      if (d > 0) return c + x;
      else return c - x;
    }
  }

  template <class TUV, class T> inline void ReduceUnredBidiagonal(
      const MatrixView<TUV>* U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>* V)
  {
    // Reduce the superdiagonal elements of unreduced Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    // Note: the input B must be unreduced - ie. all entries are non-zero.
    const size_t N = D.size();
    TMVAssert(N>0);
    TMVAssert(E.size() == N-1);
    if (U) TMVAssert(U->rowsize() == N);
    if (V) TMVAssert(V->colsize() == N);
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);

    if (N == 1) return;

    // The reduction is based on the QR algorithm to diagonalize the
    // unreduced symmetric tridiagonal matrix T = BtB
    // The basic idea is as follows:
    // (see Golub and van Loan, chapter 8 for a derivation)
    //
    // if T is a symmetric tridiagonal matrix
    // and mu is (approximately) an eigenvalue of T
    // and the QR decomposition of T - mu I = V R
    // Then, T' = R V + mu I will be tridiagonal with the last 
    // subdiagonal element small.
    // (Note: T' = Vt (T-muI) V + muI = Vt T V.)
    //
    // Wilkinson (1968) suggested that a good choice for mu is
    // the eigenvalue of the trailing 2x2 block of T that is 
    // closer to the trailing diagonal element.
    //
    // Rather than explicitly forming T = BtB and doing this
    // procedure, Golub and van Load show that it can be done
    // in place.
    // If T' = Vt T V, 
    // then Bt'B' = Vt Bt B V
    // B' = U B V for some U
    //
    // So, start with the first Givens matrix in the QR algorithm for T:
    // G0 [ T00 - mu ] = [ x ]
    //    [   T10    ]   [ 0 ]
    // We apply this to the right of B which has the effect: (for N=5)
    //              [ x x 0 0 0 ]
    //              [ + x x 0 0 ]
    // B <- B G0T = [ 0 0 x x 0 ]
    //              [ 0 0 0 x x ]
    //              [ 0 0 0 0 x ]
    // The + is the element which screws up the bidiagonal structure.
    // The rest of the procedure simply involves chasing this + down
    // the diagonal using Givens rotations.
    // For each Givens rotation we use, we also multiply U or V by the
    // adjoint to maintain the constancy of U B V.
    //
    // At the end of this procedure, E(N-1) should be smaller than it was.
    // Note: This procedure works exactly if N=2.
    T* Di = D.ptr();
    T* Ei = E.ptr();

    T mu = BidiagonalTrailingEigenValue(D,E);
    T y = NORM(*Di) - mu;  // = T00 - mu
    T x = CONJ(*Di)*(*Ei);  // = T10
    Givens<T> G = Givens_Rotate(y,x);
    for(size_t i=1;i<N;++i) {
      G.Mult(*Di,*Ei);
      if (V) G.ConjMult(V->RowPair(i-1,i));
      TMVAssert(x==T(0));
      G.Mult(x,*(++Di)); // x = B(i,i-1)
      G = Givens_Rotate(*(Di-1),x);
      G.Mult(*Ei,*Di);
      if (U) G.ConjMult(U->ColPair(i-1,i).Transpose());
      if (i < N-1) {
	TMVAssert(x==T(0));
	G.Mult(x,*(++Ei)); // x = B(i-1,i+1)
	G = Givens_Rotate(*(Ei-1),x);
      } 
    }
  }

  template <class TUV, class T> inline void ReduceBidiagonal(
      const MatrixView<TUV>* U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>* V)
  {
    // Reduce the superdiagonal elements of Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    const size_t N = D.size();
    TMVAssert(E.size() == N-1);
    if (U) TMVAssert(U->rowsize() == N); 
    if (V) TMVAssert(V->colsize() == N);

    // The input E(i) are all assumed to be non-zero.
    // If there are any zeros in D, we can zero the corresponding
    // E's (above and right) directly, so look for these first.
    // Loop invariant: all D(i) with p<=i<q are non-zero.
    size_t p=0; 
    for(size_t q=0; q<N; ++q) {
      if (D(q) == T(0)) {
	if (p<q) {
	  if (V) {
	    MatrixView<TUV> V1 = V->Rows(0,q+1);
	    BidiagonalZeroLastCol(D.SubVector(p,q+1), E.SubVector(p,q), &V1);
	    MatrixView<TUV> V2 = V->Rows(p,q);
	    if (U) {
	      MatrixView<TUV> U1 = U->Cols(p,q);
	      ReduceUnredBidiagonal(&U1, D.SubVector(p,q), E.SubVector(p,q-1), 
		  &V2);
	    } else {
	      ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,q), 
		  E.SubVector(p,q-1), &V2);
	    }
	  } else {
	    BidiagonalZeroLastCol(D.SubVector(p,q+1), E.SubVector(p,q),
		ZMV<TUV>());
	    if (U) {
	      MatrixView<TUV> U1 = U->Cols(p,q);
	      ReduceUnredBidiagonal(&U1, D.SubVector(p,q), E.SubVector(p,q-1), 
		  ZMV<TUV>());
	    } else {
	      ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,q),
		  E.SubVector(p,q-1), ZMV<TUV>());
	    }
	  }
	}
	if (q<N-1) {
	  if (U) {
	    MatrixView<TUV> U1 = U->Cols(q,N);
	    BidiagonalZeroFirstRow(&U1, D.SubVector(q,N), E.SubVector(q,N-1));
	  } else {
	    BidiagonalZeroFirstRow(ZMV<TUV>(), D.SubVector(q,N), 
		E.SubVector(q,N-1));
	  }
	}
	p=q+1;
      }
    }
    if (p<N) {
      if (U) {
	MatrixView<TUV> U1 = U->Cols(p,N);
	if (V) {
	  MatrixView<TUV> V1 = V->Rows(p,N);
	  ReduceUnredBidiagonal(&U1, D.SubVector(p,N), E.SubVector(p,N-1), &V1);
	} else {
	  ReduceUnredBidiagonal(&U1, D.SubVector(p,N), E.SubVector(p,N-1),
	      ZMV<TUV>());
	}
      } else {
	if (V) {
	  MatrixView<TUV> V1 = V->Rows(p,N);
	  ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,N), 
	      E.SubVector(p,N-1), &V1);
	} else {
	  ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,N), 
	      E.SubVector(p,N-1), ZMV<TUV>());
	}
      }
    }
  }

  template <class T> inline void SV_Decompose_From_Bidiagonal_QR(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V)
  {
    TMVAssert(D.size() == E.size()+1);
    if (U) {
      TMVAssert(U->rowsize() == D.size());
    }
    if (V) {
      TMVAssert(V->colsize() == D.size());
    }

    // We successively reduce the superdiagonal of B (E) to 0
    // using a sequence of Givens rotations. 
    // The reduction procedure tends to push the values up and left, so it 
    // makes sense to start at the lower right and work back up the matrix.
    // We also set to zero any very small values based on machine precision.
    // Loop invariant: all E(i) with i>=q are 0.
    // Initially q = N-1. (ie. All E(i) are potentially non-zero.)
    // When q = 0, we are done.
    BidiagonalChopSmallElements(D,E);

    for(size_t q = E.size(); q>0; ) {
      if (E(q-1) == T(0)) --q;
      else {
	size_t p=q-1;
	while (p > 0 && (E(p-1) != T(0))) --p; 
	// Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.
	if (U)
	  if (V) {
	    MatrixView<T> U1 = U->Cols(p,q+1);
	    MatrixView<T> V1 = V->Rows(p,q+1);
	    ReduceBidiagonal(&U1,D.SubVector(p,q+1),E.SubVector(p,q),&V1);
	  } else {
	    MatrixView<T> U1 = U->Cols(p,q+1);
	    ReduceBidiagonal(&U1,D.SubVector(p,q+1),E.SubVector(p,q),ZMV<T>());
	  }
	else
	  if (V) {
	    MatrixView<T> V1 = V->Rows(p,q+1);
	    ReduceBidiagonal(ZMV<T>(),D.SubVector(p,q+1),E.SubVector(p,q),&V1);
	  } else
	    ReduceBidiagonal(ZMV<T>(),D.SubVector(p,q+1),E.SubVector(p,q),
		ZMV<T>());
	BidiagonalChopSmallElements(D,E);
      }
    }
  }

  // MJ: This is unfinished - write divide and conquer version.
  template <class T> inline void SV_Decompose_From_Bidiagonal_DC(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>* V)
  {
    // Solve the SVD of unreduced Bidiagonal Matrix B (given by D,E).
    // Note: the input B must be unreduced - ie. all entries are non-zero.
    // This routine implements the divide and conquer approach, calling
    // itself for the recursion, and ReduceUnredBidiagonal_QR when the
    // size gets too small for divide and conquer to be efficient.
    //
    // The basic idea of the divide and conquer algorithm is to split
    // the bidiagonal matrix B into two parts, B1 and B2 and a joining
    // element:
    //     [ D0 E0                    ]   [        |          ]
    //     [    D1 E1                 ]   [  B1    |          ] N1+1
    //     [       .. ..              ]   [     Dx | Ex       ]
    // B = [          Dx Ex           ] = [-------------------]
    //     [             D(x+1) ..    ]   [        |          ]
    //     [                    .. .. ]   [        |    B2    ] N2
    //     [                       DN ]   [        |          ]
    //                                       N1+1       N2
    //
    // The smaller bidiagonal matrices are first reduced recursively:
    // B1 = U1 S1 V1
    // B2 = U2 S2 V2
    //
    //
    //     [ B1    |     ]   [ U1 S1 V1    |            ]
    //     [    Dx | Ex  ]   [          Dx | Ex         ]
    // B = [-------------] = [--------------------------]
    //     [       |     ]   [             |            ]
    //     [       |  B2 ]   [             |   U2 S2 V2 ]
    //
    //     [ U1 0  0  ] [ S1 0  0  ] [  V1   0  ]
    //   = [ 0  1  0  ] [ w1 wx w2 ] [ 0  1  0  ]
    //     [ 0  0  U2 ] [ 0  0  S2 ] [ 0  0  V2 ]
    //
    // Note that U1, U2, V2 are all square, but V1 is not, since B1 is not.
    // The vector w = [ w1 wx w2 ] is such that w V = [ 0 Dx Ex 0 ]
    //
    //
    TMVAssert(D.size()>0);
    TMVAssert(E.size()+1 == D.size());
    if (U) TMVAssert(U->rowsize() == D.size()); 
    if (V) TMVAssert(V->colsize() == D.size()); 
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);

  }
   
  template <class T> inline void DoSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V)
  {
    SV_Decompose_From_Bidiagonal_QR(U,D,E,V);
  }

  template <class T> inline void NonLapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V, bool SetUV)
  {
#ifdef XDEBUG
    //cerr<<"Start Decompose from Bidiag:\n";
    //if (U) cerr<<"U = "<<Type(*U)<<"  "<<*U<<endl;
    //cerr<<"D = "<<Type(D)<<"  step "<<D.step()<<"  "<<D<<endl;
    //cerr<<"E = "<<Type(E)<<"  step "<<E.step()<<"  "<<E<<endl;
    //if (V) cerr<<"V = "<<Type(*V)<<"  "<<*V<<endl;
    //cerr<<"SetUV = "<<SetUV<<endl;
    Matrix<RealType(T)> B(D.size(),D.size(),RealType(T)(0));
    B.diag() = D;
    B.diag(1) = E;
    Matrix<T> A0(U&&V ? U->colsize() : D.size(),D.size());
    if (U && V && !SetUV) A0 = (*U) * B * (*V);
    else A0 = B;
    //cerr<<"A0 = "<<A0<<endl;
#endif

    const size_t N = D.size();

    if (SetUV) {
      TMVAssert(U && V);
      U->SetToIdentity();
      V->SetToIdentity();
    }

    // First chop any small elements in D,E
    BidiagonalChopSmallElements(D,E);

    // Find sub-problems to solve:
    for(size_t q = N-1; q>0; ) {
      if (E(q-1) == T(0)) --q;
      else {
	// Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.
	size_t p=q-1;
	while (p > 0 && (E(p-1) != T(0))) --p; 
	if (U)
	  if (V) {
	    MatrixView<T> U1 = U->Cols(p,q+1);
	    MatrixView<T> V1 = V->Rows(p,q+1);
	    DoSV_Decompose_From_Bidiagonal(&U1,D.SubVector(p,q+1),
		E.SubVector(p,q),&V1);
	  } else {
	    MatrixView<T> U1 = U->Cols(p,q+1);
	    DoSV_Decompose_From_Bidiagonal(&U1,D.SubVector(p,q+1),
		E.SubVector(p,q),(const MatrixView<T>*)(0));
	  }
	else
	  if (V) {
	    MatrixView<T> V1 = V->Rows(p,q+1);
	    DoSV_Decompose_From_Bidiagonal((const MatrixView<T>*)(0),
		D.SubVector(p,q+1),E.SubVector(p,q),&V1);
	  } else {
	    DoSV_Decompose_From_Bidiagonal((const MatrixView<T>*)(0),
		D.SubVector(p,q+1),E.SubVector(p,q),
		(const MatrixView<T>*)(0));
	  }
	q = p > 0 ? p-1 : 0;
      }
    }

    // Make all of the singular values positive
    RealType(T)* Di = D.ptr();
    for(size_t i=0;i<N;++i,++Di) if (*Di < 0) {
      *Di = -(*Di);
      if (V) V->row(i) = -V->row(i);
    }

    // Now A = U * S * V
    // Sort output singular values 
    auto_array<size_t> sortp(new size_t[N]);
    D.Sort(sortp.get(),DESCEND);
    if (U) U->PermuteCols(sortp.get());
    if (V) V->PermuteRows(sortp.get());

#ifdef XDEBUG
    if (U && V) {
      Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
      if (Norm(A0-AA) > 0.001*Norm(A0)) {
	cerr<<"SV_DecomposeFromBidiagonal: \n";
	cerr<<"input B = "<<B<<endl;
	cerr<<"UBV = "<<A0<<endl;
	cerr<<"USV = "<<AA<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"S = "<<D<<endl;
	cerr<<"V = "<<*V<<endl;
	abort();
      }
    }
#endif
  }

#ifdef LAP 
  template <class T> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV)
  { NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV); }
#ifdef INST_DOUBLE
  template <> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<double>* U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<double>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->ct()==NonConj);
      if (U && SetUV) TMVAssert(U->stor() == V->stor());
    }
    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      if (U->iscm()) {
	TMVAssert(V->iscm());
	int ldu = U->stepj();
	int ldv = V->stepj();
	LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	    LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(U->ptr()),LAPV(ldu),LAPP(V->ptr()),LAPV(ldv),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      } else {
	u = 'L';
	TMVAssert(U->isrm());
	TMVAssert(V->isrm());
	int ldu = U->stepi();
	int ldv = V->stepi();
	LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	    LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(V->ptr()),LAPV(ldv),LAPP(U->ptr()),LAPV(ldu),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      }
    } else if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
    }
    LAP_Results("dbdsdc");
  }
  template <> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<std::complex<double> >* U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<std::complex<double> >* V, 
      bool DEBUGPARAM(SetUV))
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->ct()==NonConj);
    }
    TMVAssert(!SetUV);

    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
    }
    LAP_Results("dbdsdc");
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<float>* U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<float>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->ct()==NonConj);
      if (U && SetUV) TMVAssert(U->stor() == V->stor());
    }
    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      if (U->iscm()) {
	TMVAssert(V->iscm());
	int ldu = U->stepj();
	int ldv = V->stepj();
	LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	    LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(U->ptr()),LAPV(ldu),LAPP(V->ptr()),LAPV(ldv),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      } else {
	u = 'L';
	TMVAssert(U->isrm());
	TMVAssert(V->isrm());
	int ldu = U->stepi();
	int ldv = V->stepi();
	LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(V->ptr()),LAPV(ldv),LAPP(U->ptr()),LAPV(ldu),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      }
    } else if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
    }
    LAP_Results("sbdsdc");
  }
  template <> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<std::complex<float> >* U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<std::complex<float> >* V, 
      bool DEBUGPARAM(SetUV))
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->ct()==NonConj);
    }
    TMVAssert(!SetUV);

    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
    }
    LAP_Results("sbdsdc");
  }
#endif // FLOAT
#endif // LAP

  template <class T> void SV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->ct()==NonConj);
    }
    TMVAssert((!U || U->iscm() || U->isrm()));
    TMVAssert((!V || V->iscm() || V->isrm()));
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(!U || !V || !SetUV || U->stor()==V->stor());

    if (D.size() > 0) {
#ifdef LAP
      LapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
#else 
      //NonLapSV_Decompose_From_Bidiagonal_DelayedU(U,D,E,V,SetUV);
      NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
#endif
    }
  }

  //
  // Main SVD Drivers
  //
  
  template <class T> inline void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det, bool StoreU)
  {
#ifdef XDEBUG
    Matrix<T> A0(U);
#endif
    // Decompose A (input as U) into U S V
    // where S is a diagonal real matrix, and U,V are unitary matrices.
    // A,U are M x N (M >= N)
    // S,V are N x N
    // The determinant is returned in det.
    // (Technically, det is multiplied by the determinant, so det should
    // be set to 1 on entry.)
    const size_t M = U.colsize();
    const size_t N = U.rowsize();
    if (N == 0) return;

    TMVAssert(N <= M);
    if (V) {
      TMVAssert(V->colsize() == N);
      TMVAssert(V->rowsize() == N);
    }
    TMVAssert(S.size() == N);
    TMVAssert(U.iscm() || U.isrm());

    // If M is much larger than N (technically M > 5/3 N), then it is quicker
    // to start by doing a QR decomposition and then do SVD on the square
    // R matrix.  Thus, the final U of the SVD is Q (from the QR decomp)
    // times U from R's SVD.
    if (M > 5*N/3) {
      if (StoreU) {
	Matrix<T,ColMajor> R(N,N);
	LowerTriMatrixViewOf(R).OffDiag().Zero();
	QR_Decompose(U,UpperTriMatrixViewOf(R),det);
	SV_Decompose(R.View(),S,V,det,StoreU);
	// Now R is a Unitary Matrix U'.  Need to multiply U by U'
	U = U*R;
      } else {
	Vector<T> Qbeta(N);
	QR_Decompose(U,Qbeta.View(),det);
	if (N > 1)
	  LowerTriMatrixViewOf(U.Rows(0,N)).OffDiag().Zero();
	SV_Decompose(U.Rows(0,N),S,V,det,StoreU);
      }
    } else {
      // First we reduce A to bidiagonal form: A = U * B * V
      // using a series of Householder transformations.
      // The diagonal of the Bidiagonal Matrix B is stored in D.
      // The superdiagonal is stored in E.
      Vector<RealType(T)> E(N-1);
      Vector<T> Ubeta(N);
      Vector<T> Vbeta(N-1);
      Bidiagonalize(U,Ubeta.View(),Vbeta.View(),S,E.View(),det);
      // The determinant of B is just the product of the diagonal elements:
      if (det!=T(0)) det *= DiagMatrixViewOf(S).Det();

      // Now UV stores Householder vectors for U in lower diagonal columns 
      // (HLi) and Householder vectors for V in upper diagonal rows (HRi)
      // The Householder matrices for U are actually the adjoints of the 
      // matrices that bidiagonalize A, and for V are the transposes:
      // U = HLn-1t ... HL1t HL0t A HR0T HR1T ... HRn-2T
      // Using the fact that H Ht = I, we get A = U B V with:
      // U = HL0 ... HLn-1 
      if (V) {
	V->row(0).MakeBasis(0);
	V->Rows(1,N) = U.Rows(0,N-1);
	V->col(0,1,N).Zero();
	GetQFromQR(V->SubMatrix(1,N,1,N).Transpose(),Vbeta);
      }
      if (StoreU) {
	GetQFromQR(U,Ubeta);
      }

      if (StoreU) SV_Decompose_From_Bidiagonal(&U,S,E.View(),V);
      else SV_Decompose_From_Bidiagonal((const MatrixView<T>*)(0),S,E.View(),V);

    }
#ifdef XDEBUG
    if (StoreU && V) {
      Matrix<T> A2 = U * DiagMatrixViewOf(S) * (*V);
      if (Norm(A0-A2) > 0.0001 * Norm(U) * Norm(S) * Norm(*V)) {
	cerr<<"SV_Decompose:\n";
	cerr<<"A = "<<A0<<endl;
	cerr<<"U = "<<U<<endl;
	cerr<<"S = "<<S<<endl;
	cerr<<"V = "<<*V<<endl;
	cerr<<"USV = "<<A2<<endl;
	abort();
      }
    }
#endif
  }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det, bool StoreU)
  { SV_Decompose(U,S,&V,det,StoreU); }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      T& det, bool StoreU)
  { SV_Decompose(U,S,(const MatrixView<T>*)(0),det,StoreU); }

  */

#define InstFile "TMV_Eigen_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


