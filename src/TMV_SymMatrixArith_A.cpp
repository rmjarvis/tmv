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
#include "TMV_SymMatrix.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_SymMatrixArith.h"

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

  template <bool add, class T, class Ta, class Tx> inline void DoUnitAMultMV(
      const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  {
    const size_t N = A.size();
    if (add) y += A.LowerTri() * x;
    else y = A.LowerTri() * x;

    if (N > 1)
      y.SubVector(0,N-1) += A.UpperTri().OffDiag() * x.SubVector(1,N);
  }

  template <bool add, class T, class Ta, class Tx> inline void UnitAMultMV(
      const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  {
    // Check for 0's in the beginning or end of x:
    //      [ A11 A12 A13 ] [ 0 ]          [ A12 ]
    // y += [ A21 A22 A23 ] [ x ] --> y += [ A22 ] x
    //      [ A31 A32 A33 ] [ 0 ]          [ A32 ]

    const size_t N = x.size(); // == A.size()
    size_t j2 = N;
    for(const Tx* x2=x.cptr()+N-1; j2>0 && *x2==Tx(0); --j2,--x2);
    if (j2 == 0) {
      if (!add) y.Zero();
      return;
    }
    size_t j1 = 0;
    for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1);
    if (j1 == 0 && j2 == N) DoUnitAMultMV<add>(A,x,y);
    else {
      if (j1 > 0)
	MultMV<add>(T(1),A.SubMatrix(0,j1,j1,j2),x.SubVector(j1,j2),
	    y.SubVector(0,j1));
      TMVAssert(j1 != j2);
      DoUnitAMultMV<add>(A.SubSymMatrix(j1,j2),x.SubVector(j1,j2),
	  y.SubVector(j1,j2));
      if (j2 < N)
	MultMV<add>(T(1),A.SubMatrix(j2,N,j1,j2),x.SubVector(j1,j2),
	    y.SubVector(j2,N));
    }
  }

  template <bool add, class T, class Ta, class Tx> inline void NonBlasMultMV(
      const T alpha, const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
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
      const T alpha, const GenSymMatrix<Ta>& A,
      const GenVector<Tx>& x, int beta, const VectorView<T>& y)
  { 
    if (beta==1) NonBlasMultMV<true>(alpha,A,x,y); 
    else NonBlasMultMV<false>(alpha,A,x,y); 
  }
#ifdef INST_DOUBLE
  template <> inline void BlasMultMV(const double alpha,
      const GenSymMatrix<double>& A, const GenVector<double>& x,
      int beta, const VectorView<double>& y)
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
    int lda = A.stepj();
    int xs = x.step();
    int ys = y.step();
    double xbeta(beta);
    BLASNAME(dsymv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),
	BLASP(y.ptr()),BLASV(ys) BLAS1);
  }
  template <> inline void BlasMultMV(
      const std::complex<double> alpha,
      const GenSymMatrix<std::complex<double> >& A,
      const GenVector<std::complex<double> >& x,
      int beta, const VectorView<std::complex<double> >& y)
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
      int lda = A.stepj();
      int xs = x.step();
      int ys = y.step();
      std::complex<double> xbeta(beta);
      BLASNAME(zhemv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),
	  BLASP(y.ptr()),BLASV(ys) BLAS1);
    } else {
#ifdef ELAP
      int n = A.size();
      int lda = A.stepj();
      int xs = x.step();
      int ys = y.step();
      std::complex<double> xbeta(beta);
      LAPNAMEX(zsymv) (LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
	  LAPV(n),LAPP(&alpha),LAPP(A.cptr()),LAPV(lda),
	  LAPP(x.cptr()),LAPV(xs),LAPP(&xbeta),
	  LAPP(y.ptr()),LAPV(ys) LAP1);
#else
      if (beta==1)
	NonBlasMultMV<true>(alpha,A,x,y);
      else
	NonBlasMultMV<false>(alpha,A,x,y);
#endif
    }
  }
#endif
#ifdef INST_FLOAT
  template <> inline void BlasMultMV(const float alpha,
      const GenSymMatrix<float>& A, const GenVector<float>& x,
      int beta, const VectorView<float>& y)
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
    TMVAssert(A.isherm());

    int n = A.size();
    int lda = A.stepj();
    int xs = x.step();
    int ys = y.step();
    float xbeta(beta);
    BLASNAME(ssymv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),
	BLASP(y.ptr()),BLASV(ys) BLAS1);
  }
  template <> inline void BlasMultMV(
      const std::complex<float> alpha,
      const GenSymMatrix<std::complex<float> >& A,
      const GenVector<std::complex<float> >& x,
      int beta, const VectorView<std::complex<float> >& y)
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
      int lda = A.stepj();
      int xs = x.step();
      int ys = y.step();
      std::complex<float> xbeta(beta);
      BLASNAME(chemv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),
	  BLASP(y.ptr()),BLASV(ys) BLAS1);
    } else {
#ifdef ELAP
      int n = A.size();
      int lda = A.stepj();
      int xs = x.step();
      int ys = y.step();
      std::complex<float> xbeta(beta);
      LAPNAMEX(csymv) (LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
	  LAPV(n),LAPP(&alpha),LAPP(A.cptr()),LAPV(lda),
	  LAPP(x.cptr()),LAPV(xs),LAPP(&xbeta),
	  LAPP(y.ptr()),LAPV(ys) LAP1);
#else
      if (beta==1)
	NonBlasMultMV<true>(alpha,A,x,y);
      else
	NonBlasMultMV<false>(alpha,A,x,y);
#endif
    }
  }
#endif 
#endif // BLAS

  template <bool add, class T, class Ta, class Tx> inline void DoMultMV(
      const T alpha, const GenSymMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));

#ifdef BLAS
    if (IsComplex(T()) && (IsReal(Ta()) || IsReal(Tx())))
      BlasMultMV(alpha,A,x,add?1:0,y);
    else if (A.isconj()) 
      DoMultMV<add>(CONJ(alpha),A.Conjugate(),x.Conjugate(),y.Conjugate());
    else if (A.isrm()) {
      if (A.isherm()) DoMultMV<add>(alpha,A.Adjoint(),x,y);
      else DoMultMV<add>(alpha,A.Transpose(),x,y);
    } else {
      if (A.iscm()&&A.stepj()>0) {
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
	      HermMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    } else {
	      HermMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    }
	  } else {
	    if (A.uplo() == Upper) {
	      SymMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    } else {
	      SymMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    }
	  }
	} else {
	  if (A.isherm()) {
	    if (A.uplo() == Upper) {
	      HermMatrix<T,Upper,ColMajor> A2 = alpha*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    } else {
	      HermMatrix<T,Lower,ColMajor> A2 = alpha*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    }
	  } else {
	    if (A.uplo() == Upper) {
	      SymMatrix<T,Upper,ColMajor> A2 = alpha*A;
	      DoMultMV<add>(T(1),A2,x,y);
	    } else {
	      SymMatrix<T,Lower,ColMajor> A2 = alpha*A;
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
      const T alpha, const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
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
      const T alpha, const GenSymMatrix<Ta>& A, const MatrixView<T>& B)
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
    if (A.size() > 0) {
      if (B.isconj()) AddMM(CONJ(alpha),A.Conjugate(),B.Conjugate());
      else {
	if (SameStorage(A,B)) {
	  if (B.isrm()) {
	    Matrix<Ta,RowMajor> tempA = alpha * A;
	    B += tempA;
	  } else {
	    Matrix<Ta,ColMajor> tempA = alpha * A;
	    B += tempA;
	  }
	}
	else {
	  UpperTriMatrixViewOf(B) += alpha * A.UpperTri();
	  if (A.size() > 1)
	    LowerTriMatrixViewOf(B).OffDiag() += alpha * A.LowerTri().OffDiag();
	}
      }
    }
#ifdef XDEBUG
    //cerr<<"Done\n";
    if (Norm(B-B2) > 0.001*Norm(B)) {
      cerr<<"SymAddMM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }


  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A,
      const T beta, const GenSymMatrix<Tb>& B, const MatrixView<T>& C)
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
	    Matrix<Ta,RowMajor> tempA = alpha * A;
	    C = beta*B;
	    C += tempA;
	  } else {
	    Matrix<Ta,ColMajor> tempA = alpha * A;
	    C = beta*B;
	    C += tempA;
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
      const T alpha, const GenSymMatrix<Ta>& A,
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

    if (A.size() > 0) {
      if (SameStorage(A,C)) {
	if (SameStorage(B,C)) {
	  if (C.isrm()) {
	    Matrix<Ta,RowMajor> tempA = alpha * A;
	    C = beta * B;
	    C += tempA;
	  } else {
	    Matrix<Ta,ColMajor> tempA = alpha * A;
	    C = beta * B;
	    C += tempA;
	  }
	} else {
	  C = alpha*A;
	  C += beta*B;
	}
      } else {
	C = beta * B;
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

#define InstFile "TMV_SymMatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


