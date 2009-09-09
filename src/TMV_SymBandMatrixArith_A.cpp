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
#include "TMV_SymBandMatrix.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_SymMatrixArith.h"
#include "TMV_BandMatrixArith.h"
#include "TMV_SymBandMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // MultMV
  //

  template <bool add, class T, class Ta, class Tx> inline void UnitAMultMV(
      const GenSymBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  {
    if (add) y += A.LowerBand() * x;
    else y = A.LowerBand() * x;

    const size_t N = A.size();
    if (N > 1 && A.nlo() > 0)
      y.SubVector(0,N-1) += A.UpperBandOff() * x.SubVector(1,N);
  }

  template <bool add, class T, class Ta, class Tx> inline void NonBlasMultMV(
      const T alpha, const GenSymBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
    // y (+)= alpha * A * x
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));

    if (A.uplo() == Upper) 
      if (A.isherm()) NonBlasMultMV<add>(alpha,A.Adjoint(),x,y);
      else NonBlasMultMV<add>(alpha,A.Transpose(),x,y);
    else if (y.isconj())
      NonBlasMultMV<add>(CONJ(alpha),A.Conjugate(),x.Conjugate(),
	  y.Conjugate());
    else {
      if (x.step() != 1) {
	if (IMAG(alpha) == RealType(T)(0)) {
	  Vector<Tx> xx = REAL(alpha)*x;
	  if (y.step() != 1) {
	    Vector<T> yy(y.size());
	    UnitAMultMV<false>(A,xx,yy.View());
	    if (!add) y = yy;
	    else y += yy;
	  }
	  else 
	    UnitAMultMV<add>(A,xx,y);
	} else {
	  Vector<T> xx = alpha*x;
	  if (y.step()!=1) {
	    Vector<T> yy(y.size());
	    UnitAMultMV<false>(A,xx,yy.View());
	    if (add) y += yy;
	    else y = yy;
	  }
	  else
	    UnitAMultMV<add>(A,xx,y);
	}
      } else if (y.step()!=1 || alpha!=T(1)) {
	Vector<T> yy(y.size());
	UnitAMultMV<false>(A,x,yy.View());
	if (add) y += alpha * yy;
	else y = alpha * yy;
      } else {
	TMVAssert(alpha == T(1));
	TMVAssert(y.step() == 1);
	TMVAssert(x.step() == 1);
	UnitAMultMV<add>(A,x,y);
      }
    }
  }

#ifdef BLAS
  template <class T, class Ta, class Tx> inline void BlasMultMV(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const GenVector<Tx>& x, const int beta, const VectorView<T>& y)
  {
    if (beta == 1) NonBlasMultMV<true>(alpha,A,x,y); 
    else NonBlasMultMV<false>(alpha,A,x,y); 
  }
#ifdef INST_DOUBLE
  template <> inline void BlasMultMV(const double alpha,
      const GenSymBandMatrix<double>& A, const GenVector<double>& x,
      const int beta, const VectorView<double>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    int n = A.size();
    int k = A.nlo();
    int lda = A.diagstep();
    int xs = x.step();
    int ys = y.step();
    double xbeta(beta);
    const double* Aptr = A.cptr();
    if (A.uplo() == Upper) Aptr -= A.nlo();
    BLASNAME(dsbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	BLASV(n),BLASV(k),BLASV(alpha),BLASP(Aptr),BLASV(lda),
	BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),
	BLASP(y.ptr()),BLASV(ys) BLAS1);
  }
  template <> inline void BlasMultMV(
      const std::complex<double> alpha,
      const GenSymBandMatrix<std::complex<double> >& A,
      const GenVector<std::complex<double> >& x,
      const int beta, const VectorView<std::complex<double> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    if (A.isherm()) {
      int n = A.size();
      int k = A.nlo();
      int lda = A.diagstep();
      int xs = x.step();
      int ys = y.step();
      std::complex<double> xbeta(beta);
      const std::complex<double>* Aptr = A.cptr();
      if (A.uplo() == Upper) Aptr -= A.nlo();
      BLASNAME(zhbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASP(&alpha),BLASP(Aptr),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),
	  BLASP(y.ptr()),BLASV(ys) BLAS1);
    } else {
      if (beta == 1) NonBlasMultMV<true>(alpha,A,x,y);
      else NonBlasMultMV<false>(alpha,A,x,y);
    }
  }
#endif
#ifdef INST_FLOAT
  template <> inline void BlasMultMV(const float alpha,
      const GenSymBandMatrix<float>& A, const GenVector<float>& x,
      const int beta, const VectorView<float>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    int n = A.size();
    int k = A.nlo();
    int lda = A.diagstep();
    int xs = x.step();
    int ys = y.step();
    float xbeta(beta);
    const float* Aptr = A.cptr();
    if (A.uplo() == Upper) Aptr -= A.nlo();
    BLASNAME(ssbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	BLASV(n),BLASV(k),BLASV(alpha),BLASP(Aptr),BLASV(lda),
	BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),
	BLASP(y.ptr()),BLASV(ys) BLAS1);
  }
  template <> inline void BlasMultMV(
      const std::complex<float> alpha,
      const GenSymBandMatrix<std::complex<float> >& A,
      const GenVector<std::complex<float> >& x,
      const int beta, const VectorView<std::complex<float> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    if (A.isherm()) {
      int n = A.size();
      int k = A.nlo();
      int lda = A.diagstep();
      int xs = x.step();
      int ys = y.step();
      std::complex<float> xbeta(beta);
      const std::complex<float>* Aptr = A.cptr();
      if (A.uplo() == Upper) Aptr -= A.nlo();
      BLASNAME(chbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASP(&alpha),BLASP(Aptr),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),
	  BLASP(y.ptr()),BLASV(ys) BLAS1);
    } else {
      if (beta == 1) NonBlasMultMV<true>(alpha,A,x,y);
      else NonBlasMultMV<false>(alpha,A,x,y);
    }
  }
#endif 
#endif // BLAS

  template <bool add, class T, class Ta, class Tx> inline void DoMultMV(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));
    //cerr<<"Start DoMultMV: alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y<<endl;

#ifdef BLAS
    if (IsComplex(T()) && (IsReal(Ta()) || IsReal(Tx())))
      BlasMultMV(alpha,A,x,add?1:0,y);
    else if (A.isrm()) {
      if (A.isherm()) DoMultMV<add>(alpha,A.Adjoint(),x,y);
      else DoMultMV<add>(alpha,A.Transpose(),x,y);
    }
    else if (A.isconj()) 
      DoMultMV<add>(CONJ(alpha),A.Conjugate(),x.Conjugate(),y.Conjugate());
    else {
      if (A.iscm()&&(A.nlo()==0 || A.stepj()>0)) {
	if (!y.isconj() && y.step() > 0) { 
	  if (!x.isconj() && x.step() > 0)
	    BlasMultMV(alpha,A,x,add?1:0,y);
	  else {
	    Vector<T> xx = alpha*x;
	    BlasMultMV(T(1),A,xx,add?1:0,y);
	  }
	} else {
	  Vector<T> yy(y.size());
	  if (!x.isconj() && x.step() > 0) {
	    BlasMultMV(T(1),A,x,0,yy.View());
	    if (add) y += alpha*yy;
	    else y = alpha*yy;
	  } else {
	    Vector<T> xx = alpha*x;
	    BlasMultMV(T(1),A,xx,0,yy.View());
	    if (add) y += yy;
	    else y = yy;
	  }
	}
      } else {
	if (IMAG(alpha) == RealType(T)(0)) {
	  if (A.isherm()) {
	    if (A.uplo() == Upper) {
	      HermBandMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    } else {
	      HermBandMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    }
	  } else {
	    if (A.uplo() == Upper) {
	      SymBandMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    } else {
	      SymBandMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    }
	  }
	} else {
	  if (A.isherm()) {
	    if (A.uplo() == Upper) {
	      HermBandMatrix<T,Upper,ColMajor> A2 = alpha*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    } else {
	      HermBandMatrix<T,Lower,ColMajor> A2 = alpha*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    }
	  } else {
	    if (A.uplo() == Upper) {
	      SymBandMatrix<T,Upper,ColMajor> A2 = alpha*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    } else {
	      SymBandMatrix<T,Lower,ColMajor> A2 = alpha*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    }
	  }
	}
      }
    }
#else
    NonBlasMultMV<add>(alpha,A,x,y);
#endif
  }

  template <bool add, class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenSymBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
    // y (+)= alpha * A * x
  { 
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
#ifdef XDEBUG
    //cerr<<"Start MultMV: alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y<<endl;
    Matrix<Ta> A0 = A;
    Vector<Tx> x0 = x;
    Vector<T> y0 = y;
    Vector<T> y2 = alpha*A0*x0;
    if (add) y2 += y0;
#endif

    if (y.size() > 0) {
      if (x.size()==0 || alpha==T(0)) {
	if (!add) y.Zero();
      } else if (SameStorage(x,y)) {
	Vector<T> yy(y.size());
	DoMultMV<false>(T(1),A,x,yy.View());
	if (add) y += alpha*yy;
	else y = alpha*yy;
      } else {
	DoMultMV<add>(alpha,A,x,y);
      } 
    }
#ifdef XDEBUG
    //cerr<<"--> y = "<<y<<endl;
    //cerr<<"y2 = "<<y2<<endl;
    if (Norm(y-y2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(x0)+
	  (add?Norm(y0):RealType(T)(0)))) {
      cerr<<"MultMV: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"--> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, const BandMatrixView<T>& B)
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif

#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> B2 = alpha*A0 + B0;
    //cerr<<"Start SymAddMM: alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());
    TMVAssert(B.nlo() >= A.nlo());
    TMVAssert(B.nhi() >= A.nlo());
    if (A.size() > 0) {
      if (B.isconj()) AddMM(CONJ(alpha),A.Conjugate(),B.Conjugate());
      else {
	if (SameStorage(A,B)) {
	  if (B.isrm()) {
	    BandMatrix<Ta,RowMajor> tempA(A);
	    AddMM(alpha,tempA,B);
	  } else {
	    BandMatrix<Ta,ColMajor> tempA(A);
	    AddMM(alpha,tempA,B);
	  }
	}
	else {
	  AddMM(alpha,A.UpperBand(),BandMatrixViewOf(B,0,A.nlo()));
	  if (A.nlo()>0)
	    AddMM(alpha,A.LowerBandOff(),
		BandMatrixViewOf(B,A.nlo(),0).Diags(-A.nlo(),0));
	}
      }
    }
#ifdef XDEBUG
    //cerr<<"Done\n";
    if (Norm(B-B2) > 0.001*(Norm(B0)+ABS(alpha)*Norm(A))) {
      cerr<<"SymAddMM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }


  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const T beta, const GenSymBandMatrix<Tb>& B, const BandMatrixView<T>& C)
  { 
#ifdef XTEST
    TMVAssert(A.HermOK());
    TMVAssert(B.HermOK());
#endif
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = alpha*A0 + beta*B0;
    //cerr<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.size());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());

    if (A.size() > 0) {
      if (SameStorage(A,C)) {
	if (SameStorage(B,C)) {
	  if (C.isrm()) {
	    BandMatrix<T,RowMajor> tempC = beta*B;
	    tempC += alpha*A;
	    C = tempC;
	  } else {
	    BandMatrix<T,ColMajor> tempC = beta*B;
	    tempC += alpha*A;
	    C = tempC;
	  }
	} else {
	  C = alpha*A;
	  AddMM(beta,B,C);
	}
      } else {
	C = beta*B;
	AddMM(alpha,A,C);
      }
    }

#ifdef XDEBUG
    //cerr<<"Done\n";
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)+ABS(beta)*Norm(B0))) {
      cerr<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  } 

  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const T beta, const GenSymBandMatrix<Tb>& B, const MatrixView<T>& C)
  { 
#ifdef XTEST
    TMVAssert(A.HermOK());
    TMVAssert(B.HermOK());
#endif
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = alpha*A0 + beta*B0;
    //cerr<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.size());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());

    const size_t N = A.size();
    const int k = std::max(A.nlo(),B.nlo());

    if (N > 0) {
      if (SameStorage(A,C) || SameStorage(B,C)) {
	AddMM(alpha,A,beta,B,BandMatrixViewOf(C,k,k));
	UpperTriMatrixViewOf(C.SubMatrix(0,N-k-1,k+1,N)).Zero();
	LowerTriMatrixViewOf(C.SubMatrix(k+1,N,0,N-k-1)).Zero();
      } else {
	C.Zero();
	AddMM(alpha,A,beta,B,BandMatrixViewOf(C,k,k));
      }
    }

#ifdef XDEBUG
    //cerr<<"Done\n";
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)+ABS(beta)*Norm(B0))) {
      cerr<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const T beta, const GenBandMatrix<Tb>& B, const BandMatrixView<T>& C)
  { 
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = alpha*A0 + beta*B0;
    //cerr<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());
    TMVAssert(C.nlo() >= A.nlo());
    TMVAssert(C.nlo() >= B.nlo());
    TMVAssert(C.nhi() >= A.nhi());
    TMVAssert(C.nhi() >= B.nhi());

    if (A.size() > 0) {
      if (SameStorage(A,C)) {
	if (SameStorage(B,C)) {
	  if (C.isrm()) {
	    BandMatrix<T,RowMajor> tempC(C.colsize(),C.rowsize(),
		C.nlo(),C.nhi());
	    tempC = beta*B;
	    tempC += alpha*A;
	    C = tempC;
	  } else {
	    BandMatrix<T,ColMajor> tempC(C.colsize(),C.rowsize(),
		C.nlo(),C.nhi());
	    tempC = beta*B;
	    tempC += alpha*A;
	    C = tempC;
	  }
	} else {
	  C = alpha*A;
	  C += beta*B;
	}
      } else {
	C = beta*B;
	AddMM(alpha,A,C);
      }
    }

#ifdef XDEBUG
    //cerr<<"Done\n";
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)+ABS(beta)*Norm(B0))) {
      cerr<<"SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  } 

  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const T beta, const GenMatrix<Tb>& B, const MatrixView<T>& C)
  { 
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = alpha*A0 + beta*B0;
    //cerr<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());

    const size_t N = A.size();
    const int k = A.nlo();
    if (N > 0) {
      AddMM(alpha,A,beta,BandMatrixViewOf(B,k,k),BandMatrixViewOf(C,k,k));
      UpperTriMatrixViewOf(C.SubMatrix(0,N-k-1,k+1,N)) =
	beta * UpperTriMatrixViewOf(B.SubMatrix(0,N-k-1,k+1,N)); 
      LowerTriMatrixViewOf(C.SubMatrix(k+1,N,0,N-k-1)) =
	beta * LowerTriMatrixViewOf(B.SubMatrix(k+1,N,0,N-k-1));
    }

#ifdef XDEBUG
    //cerr<<"Done\n";
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)+ABS(beta)*Norm(B0))) {
      cerr<<"SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymBandMatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


