///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_SymBandMatrixArithFunc.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"
#include "TMV_BandMatrix.h"
#include "TMV_BandMatrixArith.h"
#ifdef BLAS
#include "TMV_SymBandMatrixArith.h"
#endif
#include <iostream>

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <class T> const T* SymBandMatrixComposite<T>::cptr() const
  {                       
    if (!itsm1.get()) {         
      size_t s = this->size();        
      size_t lo = this->nlo();      
      size_t len = BandStorageLength(this->stor(),s,s,lo,0);
      itsm1.reset(new T[len]);                          
      itsm = this->isdm() ? itsm1.get() + lo*(s-1) : itsm1.get();
      this->AssignTosB(SymBandMatrixView<T>(itsm,s,lo,
	    stepi(),stepj(),diagstep(),
	    Sym,this->uplo(),this->stor(),NonConj                   
	    FIRSTLAST1(itsm.get(),itsm.get()+len) ));
    }                                     
    return itsm;
  }               

  template <class T> int SymBandMatrixComposite<T>::stepi() const
  {   
    return this->iscm() ? 1 : this->isrm() ? this->nlo() : 
      -int(this->size()) + 1;
  }

  template <class T> int SymBandMatrixComposite<T>::stepj() const
  { 
    return this->isrm() ? 1 : this->iscm() ? this->nlo() : 
      int(this->size());
  }

  template <class T> int SymBandMatrixComposite<T>::diagstep() const
  { 
    return this->isdm() ? 1 : this->nlo()+1;
  }


  //
  // MultMV
  //

  template <bool add, class T, class Ta, class Tx> static void UnitAMultMV(
      const GenSymBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  {
    if (add) y += A.LowerBand() * x;
    else y = A.LowerBand() * x;

    const int N = A.size();
    if (N > 1 && A.nlo() > 0)
      y.SubVector(0,N-1) += A.UpperBandOff() * x.SubVector(1,N);
  }

  template <bool add, class T, class Ta, class Tx> static void NonBlasMultMV(
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
  template <class T, class Ta, class Tx> static inline void BlasMultMV(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const GenVector<Tx>& x, const int beta, const VectorView<T>& y)
  {
    if (beta == 1) NonBlasMultMV<true>(alpha,A,x,y); 
    else NonBlasMultMV<false>(alpha,A,x,y); 
  }
#ifdef INST_DOUBLE
  template <> void BlasMultMV(const double alpha,
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
  template <> void BlasMultMV(
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
  template <> void BlasMultMV(
      const std::complex<double> alpha,
      const GenSymBandMatrix<std::complex<double> >& A,
      const GenVector<double>& x,
      const int beta, const VectorView<std::complex<double> >& y)
  { BlasMultMV(alpha,A,Vector<std::complex<double> >(x),beta,y); }
  template <> void BlasMultMV(
      const std::complex<double> alpha,
      const GenSymBandMatrix<double>& A, 
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

    if (beta == 0) {
      int n = A.size();
      int k = A.nlo();
      int lda = A.diagstep();
      int xs = 2*x.step();
      int ys = 2*y.step();
      double xalpha(1);
      double xbeta(beta);
      const double* Aptr = A.cptr();
      if (A.uplo() == Upper) Aptr -= A.nlo();
      BLASNAME(dsbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(xalpha),BLASP(Aptr),BLASV(lda),
	  BLASP((double*)x.cptr()),BLASV(xs),BLASV(xbeta),
	  BLASP((double*)y.ptr()),BLASV(ys) BLAS1);
      BLASNAME(dsbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(xalpha),BLASP(Aptr),BLASV(lda),
	  BLASP((double*)x.cptr()+1),BLASV(xs),BLASV(xbeta),
	  BLASP((double*)y.ptr()+1),BLASV(ys) BLAS1);
      if (x.isconj()) y.ConjugateSelf();
      y *= alpha;
    } else if (IMAG(alpha) == float(0) && !x.isconj()) {
      int n = A.size();
      int k = A.nlo();
      int lda = A.diagstep();
      int xs = 2*x.step();
      int ys = 2*y.step();
      double xalpha(REAL(alpha));
      double xbeta(beta);
      const double* Aptr = A.cptr();
      if (A.uplo() == Upper) Aptr -= A.nlo();
      BLASNAME(dsbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(xalpha),BLASP(Aptr),BLASV(lda),
	  BLASP((double*)x.cptr()),BLASV(xs),BLASV(xbeta),
	  BLASP((double*)y.ptr()),BLASV(ys) BLAS1);
      BLASNAME(dsbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(xalpha),BLASP(Aptr),BLASV(lda),
	  BLASP((double*)x.cptr()+1),BLASV(xs),BLASV(xbeta),
	  BLASP((double*)y.ptr()+1),BLASV(ys) BLAS1);
    } else {
      Vector<std::complex<double> > xx = alpha*x;
      BlasMultMV(std::complex<double>(1),A,xx,1,y);
    }
  }
  template <> void BlasMultMV(
      const std::complex<double> alpha,
      const GenSymBandMatrix<double>& A, 
      const GenVector<double>& x,
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

    int n = A.size();
    int k = A.nlo();
    int lda = A.diagstep();
    int xs = x.step();
    int ys = 2*y.step();
    double ar(REAL(alpha));
    double ai(IMAG(alpha));
    double xbeta(beta);
    const double* Aptr = A.cptr();
    if (A.uplo() == Upper) Aptr -= A.nlo();
    if (ar == double(0)) {
      if (beta == 0) y.Real().Zero();
    } else
      BLASNAME(dsbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(ar),BLASP(Aptr),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),
	  BLASP((double*)y.ptr()),BLASV(ys) BLAS1);
    if (ai == double(0)) {
      if (beta == 0) y.Imag().Zero();
    } else 
      BLASNAME(dsbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(ai),BLASP(Aptr),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),
	  BLASP((double*)y.ptr()+1),BLASV(ys) BLAS1);
  }
#endif
#ifdef INST_FLOAT
  template <> void BlasMultMV(const float alpha,
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
  template <> void BlasMultMV(
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
  template <> void BlasMultMV(
      const std::complex<float> alpha,
      const GenSymBandMatrix<std::complex<float> >& A,
      const GenVector<float>& x,
      const int beta, const VectorView<std::complex<float> >& y)
  { BlasMultMV(alpha,A,Vector<std::complex<float> >(x),beta,y); }
  template <> void BlasMultMV(
      const std::complex<float> alpha,
      const GenSymBandMatrix<float>& A, 
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

    if (beta == 0) {
      int n = A.size();
      int k = A.nlo();
      int lda = A.diagstep();
      int xs = 2*x.step();
      int ys = 2*y.step();
      float xalpha(1);
      float xbeta(beta);
      const float* Aptr = A.cptr();
      if (A.uplo() == Upper) Aptr -= A.nlo();
      BLASNAME(ssbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(xalpha),BLASP(Aptr),BLASV(lda),
	  BLASP((float*)x.cptr()),BLASV(xs),BLASV(xbeta),
	  BLASP((float*)y.ptr()),BLASV(ys) BLAS1);
      BLASNAME(ssbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(xalpha),BLASP(Aptr),BLASV(lda),
	  BLASP((float*)x.cptr()+1),BLASV(xs),BLASV(xbeta),
	  BLASP((float*)y.ptr()+1),BLASV(ys) BLAS1);
      if (x.isconj()) y.ConjugateSelf();
      y *= alpha;
    } else if (IMAG(alpha) == float(0) && !x.isconj()) {
      int n = A.size();
      int k = A.nlo();
      int lda = A.diagstep();
      int xs = 2*x.step();
      int ys = 2*y.step();
      float xalpha(REAL(alpha));
      float xbeta(beta);
      const float* Aptr = A.cptr();
      if (A.uplo() == Upper) Aptr -= A.nlo();
      BLASNAME(ssbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(xalpha),BLASP(Aptr),BLASV(lda),
	  BLASP((float*)x.cptr()),BLASV(xs),BLASV(xbeta),
	  BLASP((float*)y.ptr()),BLASV(ys) BLAS1);
      BLASNAME(ssbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(xalpha),BLASP(Aptr),BLASV(lda),
	  BLASP((float*)x.cptr()+1),BLASV(xs),BLASV(xbeta),
	  BLASP((float*)y.ptr()+1),BLASV(ys) BLAS1);
    } else {
      Vector<std::complex<float> > xx = alpha*x;
      BlasMultMV(std::complex<float>(1),A,xx,1,y);
    }
  }
  template <> void BlasMultMV(
      const std::complex<float> alpha,
      const GenSymBandMatrix<float>& A, 
      const GenVector<float>& x,
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

    int n = A.size();
    int k = A.nlo();
    int lda = A.diagstep();
    int xs = x.step();
    int ys = 2*y.step();
    float ar(REAL(alpha));
    float ai(IMAG(alpha));
    float xbeta(beta);
    const float* Aptr = A.cptr();
    if (A.uplo() == Upper) Aptr -= A.nlo();
    if (ar == float(0)) {
      if (beta == 0) y.Real().Zero();
    } else
      BLASNAME(ssbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(ar),BLASP(Aptr),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),
	  BLASP((float*)y.ptr()),BLASV(ys) BLAS1);
    if (ai == float(0)) {
      if (beta == 0) y.Imag().Zero();
    } else 
      BLASNAME(ssbmv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASV(k),BLASV(ai),BLASP(Aptr),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),
	  BLASP((float*)y.ptr()+1),BLASV(ys) BLAS1);
  }
#endif
#endif // BLAS

  template <bool add, class T, class Ta, class Tx> static void DoMultMV(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));
    //cout<<"Start DoMultMV: alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cout<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y<<endl;

#ifdef BLAS
    if (A.isrm()) {
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
    //cout<<"Start MultMV: alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cout<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y<<endl;
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
    //cout<<"--> y = "<<y<<endl;
    //cout<<"y2 = "<<y2<<endl;
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

#define InstFile "TMV_MultsBV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


