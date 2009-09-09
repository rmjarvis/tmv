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
#include "TMV_SymMatrixArithFunc.h"
#include "TMV_SymMatrixArith.h"
#include "TMV_SymMatrix.h"
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

  // 
  // Rank2Update
  //

  template <bool ha, bool add, class T, class Tx, class Ty> 
    static void UpperRank2Update(
	const GenVector<Tx>& x, const GenVector<Ty>& y,
	const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(A.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Upper);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(A.iscm());
      TMVAssert(x.ct() == NonConj);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(ha == !A.issym());
#ifdef XTEST
      RealType(T) NormA = Norm(A);
      RealType(T) Normx = Norm(x);
      RealType(T) Normy = Norm(y);
      RealType(T) eps = RealType(T)(2)*Normx*Normy;
      if (add) eps += NormA;
      eps *= A.size() * Epsilon<T>();
#endif

      const int sj = A.stepj();
      const int N = A.size();
      const Tx*const x0 = x.cptr();
      const Ty*const y0 = y.cptr();

      const Tx* xj = x0;
      const Ty* yj = y0;

      T A00;
      if (*xj == Tx(0) || *yj == Ty(0)) A00 = T(0);
      else if (ha) A00 = RealType(T)(2) * REAL(*xj * CONJ(*yj));
      else A00 = RealType(T)(2) * (*xj) * (*yj);
      ++xj; ++yj;
      T* Acolj = A.ptr()+sj;

      for (int j=1;j<N;++j,++xj,++yj,Acolj+=sj) {
	// A.col(j,0,j+1) += ax * y.SubVector(0,j+1) + ay * x.SubVector(0,j+1);
	if (*xj != Tx(0)) {
	  T* Aij = Acolj;
	  const Ty* yi = y0;
	  if (*yj != Ty(0)) {
	    const Tx* xi = x0;
	    for(int i=j+1;i>0;--i,++yi,++xi,++Aij) {
	      T temp = *xi * (ha ? CONJ(*yj) : *yj);
	      temp += *yi * (ha ? CONJ(*xj) : *xj);
#ifdef TMVFLDEBUG
	      TMVAssert(Aij >= A.first);
	      TMVAssert(Aij < A.last);
#endif
	      if (add) *Aij += temp;
	      else *Aij = temp;
	    }
	  } else {
	    for(int i=j+1;i>0;--j,++yi,++Aij) {
	      const T temp = *yi * (ha ? CONJ(*xj) : *xj);
#ifdef TMVFLDEBUG
	      TMVAssert(Aij >= A.first);
	      TMVAssert(Aij < A.last);
#endif
	      if (add) *Aij += temp;
	      else *Aij = temp;
	    }
	  }
	} else if (*yj != Ty(0)) {
	  T* Aij = Acolj;
	  const Tx* xi = x0;
	  for(int i=j+1;i>0;--i,++xi,++Aij) {
	    const T temp = *xi * (ha ? CONJ(*yj) : *yj);
#ifdef TMVFLDEBUG
	    TMVAssert(Aij >= A.first);
	    TMVAssert(Aij < A.last);
#endif
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	} else if (!add) {
	  T* Aij = Acolj;
	  memset(Aij,0,(j+1)*sizeof(T));
	}
      }
#ifdef TMVFLDEBUG
      TMVAssert(A.ptr() >= A.first);
      TMVAssert(A.ptr() < A.last);
#endif
      if (add) *A.ptr() += A00;
      else *A.ptr() = A00;
      if (ha && IsComplex(T())) {
#ifdef XTEST
        TMVAssert(NormInf(A.diag().Imag()) <= eps);
#endif
	A.diag().Imag().Zero();
      }
    }

  template <bool ha, bool add, class T, class Tx, class Ty> 
    static void LowerRank2Update(
	const GenVector<Tx>& x, const GenVector<Ty>& y,
	const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(A.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(A.iscm());
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(x.ct() == NonConj);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(ha == !A.issym());
#ifdef XTEST
      RealType(T) NormA = Norm(A);
      RealType(T) Normx = Norm(x);
      RealType(T) Normy = Norm(y);
      RealType(T) eps = RealType(T)(2)*Normx*Normy;
      if (add) eps += NormA;
      eps *= A.size() * Epsilon<T>();
#endif

      const int ds = A.stepj() + 1;
      const int N = A.size();
      const Tx* xj = x.cptr()+N-1;
      const Ty* yj = y.cptr()+N-1;
      T* Ajj = A.ptr()+(N-1)*ds;

      for (int jj=N,Nmj=1;jj>0;--jj,++Nmj,--xj,--yj,Ajj-=ds) {
	// Nmj = N-j
	// A.col(j,j,N) += (A.isherm() ? CONJ(*yj) : *yj) * x.SubVector(j,N);
	// A.col(j,j,N) += (A.isherm() ? CONJ(*xj) : *xj) * y.SubVector(j,N);
	if (*yj!=Ty(0)) {
	  T* Aij = Ajj;
	  const Tx* xi = xj;
	  if (*xj!=Tx(0)) {
	    const Ty* yi = yj;
	    for(int i=Nmj;i>0;--i,++xi,++yi,++Aij) {
	      T temp = *xi * (ha ? CONJ(*yj) : *yj);
	      temp += *yi * (ha ? CONJ(*xj) : *xj);
#ifdef TMVFLDEBUG
	      TMVAssert(Aij >= A.first);
	      TMVAssert(Aij < A.last);
#endif
	      if (add) *Aij += temp;
	      else *Aij = temp;
	    }
	  } else {
	    for(int i=Nmj;i>0;--i,++xi,++Aij) {
	      const T temp = *xi * (ha ? CONJ(*yj) : *yj);
#ifdef TMVFLDEBUG
	      TMVAssert(Aij >= A.first);
	      TMVAssert(Aij < A.last);
#endif
	      if (add) *Aij += temp;
	      else *Aij = temp;
	    }
	  }
	} else if (*xj!=Tx(0)) {
	  T* Aij = Ajj;
	  const Ty* yi = yj;
	  for(int i=Nmj;i>0;--i,++yi,++Aij) {
	    const T temp = *yi * (ha ? CONJ(*xj) : *xj);
#ifdef TMVFLDEBUG
	    TMVAssert(Aij >= A.first);
	    TMVAssert(Aij < A.last);
#endif
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	} else if (!add) {
	  T* Aij = Ajj;
	  memset(Aij,0,Nmj*sizeof(T));
	}
      }
      if (ha && IsComplex(T())) {
#ifdef XTEST
        TMVAssert(NormInf(A.diag().Imag()) <= eps);
#endif
	A.diag().Imag().Zero();
      }
    }

  template <bool add, class T, class Tx, class Ty> 
    struct UnitARank2Update
    {
      static void F(const GenVector<Tx>& x,
	const GenVector<Ty>& y, const SymMatrixView<T>& A)
      {
	TMVAssert(A.size() == x.size());
	TMVAssert(A.size() == y.size());
	TMVAssert(A.size() > 0);
	TMVAssert(A.ct() == NonConj);
	TMVAssert(A.iscm());
	TMVAssert(x.step() == 1);
	TMVAssert(y.step() == 1);
	TMVAssert(x.ct() == NonConj);
	TMVAssert(y.ct() == NonConj);

	if (A.isupper()) UpperRank2Update<false,add>(x,y,A);
	else LowerRank2Update<false,add>(x,y,A);
      }
    };

  template <bool add, class T, class Tx, class Ty> 
    struct UnitARank2Update<add,std::complex<T>,Tx,Ty>
    {
      static void F(const GenVector<Tx>& x,
	const GenVector<Ty>& y, const SymMatrixView<std::complex<T> >& A)
      {
	TMVAssert(A.size() == x.size());
	TMVAssert(A.size() == y.size());
	TMVAssert(A.size() > 0);
	TMVAssert(A.ct() == NonConj);
	TMVAssert(A.iscm());
	TMVAssert(x.step() == 1);
	TMVAssert(y.step() == 1);
	TMVAssert(x.ct() == NonConj);
	TMVAssert(y.ct() == NonConj);

	if (A.isherm())
	  if (A.isupper()) UpperRank2Update<true,add>(x,y,A);
	  else LowerRank2Update<true,add>(x,y,A);
	else
	  if (A.isupper()) UpperRank2Update<false,add>(x,y,A);
	  else LowerRank2Update<false,add>(x,y,A);
      }
    };
   

  template <bool add, class T, class Tx, class Ty> 
    static void NonBlasRank2Update(
	const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
	const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(A.iscm());
      TMVAssert(A.ct() == NonConj);
      TMVAssert(alpha != T(0));
      TMVAssert(A.size() > 0);

      if (x.step() != 1 || x.isconj()) {
	// Copy x to new storage
	if (y.step() != 1 || y.isconj()) {
	  // Copy x and y to new storage
	  if (x.size() <= y.size()) {
	    if (IMAG(alpha) == RealType(T)(0)) {
	      Vector<Tx> xx = REAL(alpha)*x;
	      Vector<Ty> yy = y;
	      UnitARank2Update<add,T,Tx,Ty>::F(xx,yy,A);
	    } else {
	      Vector<T> xx = alpha*x;
	      Vector<Ty> yy = y;
	      UnitARank2Update<add,T,T,Ty>::F(xx,yy,A);
	    }
	  } else {
	    if (IMAG(alpha) == RealType(T)(0)) {
	      Vector<Tx> xx = x;
	      Vector<Ty> yy = REAL(alpha)*y;
	      UnitARank2Update<add,T,Tx,Ty>::F(xx,yy,A);
	    } else {
	      Vector<Tx> xx = x;
	      Vector<T> yy = alpha*y;
	      UnitARank2Update<add,T,Tx,T>::F(xx,yy,A);
	    }
	  }
	} else {
	  // Copy only x to new storage
	  if (IMAG(alpha) == RealType(T)(0)) {
	    Vector<Tx> xx = REAL(alpha)*x;
	    UnitARank2Update<add,T,Tx,Ty>::F(xx,y,A);
	  } else {
	    Vector<T> xx = alpha*x;
	    UnitARank2Update<add,T,T,Ty>::F(xx,y,A);
	  }
	}
      } else if (y.step() != 1 || y.isconj()) {
	// Copy only y to new storage
	if (IMAG(alpha) == RealType(T)(0)) {
	  Vector<Ty> yy = REAL(alpha)*y;
	  UnitARank2Update<add,T,Tx,Ty>::F(x,yy,A);
	} else {
	  Vector<T> yy = alpha*y;
	  UnitARank2Update<add,T,Tx,T>::F(x,yy,A);
	}
      } else if (alpha != T(1)) {
	// Copy something to new storage to incorporate alpha
	if (x.size() <= y.size()) {
	  if (IMAG(alpha) == RealType(T)(0)) {
	    Vector<Tx> xx = REAL(alpha)*x;
	    UnitARank2Update<add,T,Tx,Ty>::F(xx,y,A);
	  } else {
	    Vector<T> xx = alpha*x;
	    UnitARank2Update<add,T,T,Ty>::F(xx,y,A);
	  }
	} else {
	  if (IMAG(alpha) == RealType(T)(0)) {
	    Vector<Ty> yy = REAL(alpha)*y;
	    UnitARank2Update<add,T,Tx,Ty>::F(x,yy,A);
	  } else {
	    Vector<T> yy = alpha*y;
	    UnitARank2Update<add,T,Tx,T>::F(x,yy,A);
	  }
	}
      } else {
	UnitARank2Update<add,T,Tx,Ty>::F(x,y,A);
      }
    }

#ifdef BLAS
  template <class T, class Tx, class Ty> static inline void BlasRank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  { NonBlasRank2Update<true>(alpha,x,y,A); }
#ifdef INST_DOUBLE
  template <> void BlasRank2Update(
      const double alpha, const GenVector<double>& x,
      const GenVector<double>& y, const SymMatrixView<double>& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());

    int n=A.size();
    int xs=x.step();
    int ys=y.step();
    int lda=A.stepj();
    BLASNAME(dsyr2) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	BLASV(n),BLASV(alpha),BLASP(x.cptr()),BLASV(xs),
	BLASP(y.cptr()),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
  }
  template <> void BlasRank2Update(
      const std::complex<double> alpha,
      const GenVector<std::complex<double> >& x, 
      const GenVector<std::complex<double> >& y, 
      const SymMatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    if (A.issym() && (x.step() != 1 || y.step() != 1)) {
      if (x.step() != 1) {
	Vector<std::complex<double> > xx = x;
	if (y.step() != 1) {
	  Vector<std::complex<double> > yy = y;
	  return BlasRank2Update(alpha,xx,yy,A);
	} 
	else return BlasRank2Update(alpha,xx,y,A);
      } else {
	TMVAssert(y.step() != 1);
	Vector<std::complex<double> > yy = y;
	return BlasRank2Update(alpha,x,yy,A);
      }
    } else {
      int n=A.size();
      int xs=x.step();
      int ys=y.step();
      int lda=A.stepj();
      if (A.isherm()) {
	BLASNAME(zher2) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	    BLASV(n),BLASP(&alpha),BLASP(x.cptr()),BLASV(xs),
	    BLASP(y.cptr()),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
      } else {
	int k=1;
	std::complex<double> beta(1);
	BLASNAME(zsyr2k) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	    BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
	    BLASP(x.cptr()),BLASV(n),BLASP(y.cptr()),BLASV(n),
	    BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
      }
    }
  }
  template <> void BlasRank2Update(
      const std::complex<double> alpha,
      const GenVector<std::complex<double> >& x, 
      const GenVector<double>& y, 
      const SymMatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    SymMatrix<double> A1(A.size(),double(0));
    BlasRank2Update(double(1),x.Real(),y,A1.View());
    A += alpha*A1;
    A1.Zero();
    BlasRank2Update(double(1),x.Imag(),y,A1.View());
    A += std::complex<double>(0,1)*alpha*A1;
  }
  template <> void BlasRank2Update(
      const std::complex<double> alpha,
      const GenVector<double>& x, 
      const GenVector<std::complex<double> >& y, 
      const SymMatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    SymMatrix<double> A1(A.size(),double(0));
    BlasRank2Update(double(1),x,y.Real(),A1.View());
    A += alpha*A1;
    A1.Zero();
    BlasRank2Update(double(1),x,y.Imag(),A1.View());
    A += std::complex<double>(0,1)*alpha*A1;
  }
  template <> void BlasRank2Update(
      const std::complex<double> alpha,
      const GenVector<double>& x, const GenVector<double>& y, 
      const SymMatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    SymMatrix<double> A1(A.size(),double(0));
    BlasRank2Update(double(1),x,y,A1.View());
    A += alpha*A1;
  }
#endif
#ifdef INST_FLOAT
  template <> void BlasRank2Update(
      const float alpha, const GenVector<float>& x,
      const GenVector<float>& y, const SymMatrixView<float>& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());

    int n=A.size();
    int xs=x.step();
    int ys=y.step();
    int lda=A.stepj();
    BLASNAME(ssyr2) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	BLASV(n),BLASV(alpha),BLASP(x.cptr()),BLASV(xs),
	BLASP(y.cptr()),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
  }
  template <> void BlasRank2Update(
      const std::complex<float> alpha,
      const GenVector<std::complex<float> >& x, 
      const GenVector<std::complex<float> >& y, 
      const SymMatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    if (A.issym() && (x.step() != 1 || y.step() != 1)) {
      if (x.step() != 1) {
	Vector<std::complex<float> > xx = x;
	if (y.step() != 1) {
	  Vector<std::complex<float> > yy = y;
	  return BlasRank2Update(alpha,xx,yy,A);
	} 
	else return BlasRank2Update(alpha,xx,y,A);
      } else {
	TMVAssert(y.step() != 1);
	Vector<std::complex<float> > yy = y;
	return BlasRank2Update(alpha,x,yy,A);
      }
    } else {
      int n=A.size();
      int xs=x.step();
      int ys=y.step();
      int lda=A.stepj();
      if (A.isherm()) {
	BLASNAME(cher2) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	    BLASV(n),BLASP(&alpha),BLASP(x.cptr()),BLASV(xs),
	    BLASP(y.cptr()),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
      } else {
	int k=1;
	std::complex<float> beta(1);
	BLASNAME(csyr2k) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	    BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
	    BLASP(x.cptr()),BLASV(n),BLASP(y.cptr()),BLASV(n),
	    BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
      }
    }
  }
  template <> void BlasRank2Update(
      const std::complex<float> alpha,
      const GenVector<std::complex<float> >& x, 
      const GenVector<float>& y, 
      const SymMatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    SymMatrix<float> A1(A.size(),float(0));
    BlasRank2Update(float(1),x.Real(),y,A1.View());
    A += alpha*A1;
    A1.Zero();
    BlasRank2Update(float(1),x.Imag(),y,A1.View());
    A += std::complex<float>(0,1)*alpha*A1;
  }
  template <> void BlasRank2Update(
      const std::complex<float> alpha,
      const GenVector<float>& x, 
      const GenVector<std::complex<float> >& y, 
      const SymMatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    SymMatrix<float> A1(A.size(),float(0));
    BlasRank2Update(float(1),x,y.Real(),A1.View());
    A += alpha*A1;
    A1.Zero();
    BlasRank2Update(float(1),x,y.Imag(),A1.View());
    A += std::complex<float>(0,1)*alpha*A1;
  }
  template <> void BlasRank2Update(
      const std::complex<float> alpha,
      const GenVector<float>& x, const GenVector<float>& y, 
      const SymMatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    SymMatrix<float> A1(A.size(),float(0));
    BlasRank2Update(float(1),x,y,A1.View());
    A += alpha*A1;
  }
#endif 
#endif // BLAS

  template <bool add, class T, class Tx, class Ty> void Rank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
    // if A is sym:  A (+)= alpha * (x ^ y + y ^ x)
    // if A is herm: A (+)= alpha * x ^ y* + conj(alpha) * y ^ x*
  {
#ifdef XTEST
    TMVAssert(!add || A.HermOK());
#endif
#ifdef XDEBUG
    Vector<Tx> x0 = x;
    Vector<Ty> y0 = y;
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    if (A.isherm()) {
      if (add) A2 += (alpha*x^y.Conjugate());
      else A2 = (alpha*x^y.Conjugate());
      A2 += (CONJ(alpha)*y^x.Conjugate());
    }
    else {
      if (add) A2 += alpha*(x^y);
      else A2 = alpha*(x^y);
      A2 += alpha*(y^x);
    }
    //cout<<"Start Rank2Update: alpha = "<<alpha<<endl;
    //cout<<"add = "<<add<<endl;
    //cout<<"x = "<<Type(x)<<"  step = "<<x.step()<<"  "<<x0<<endl;
    //cout<<"y = "<<Type(y)<<"  step = "<<y.step()<<"  "<<y0<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A0<<endl;
#endif

    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
	Rank2Update<add>(CONJ(alpha),x.Conjugate(),y.Conjugate(),
	    A.Conjugate());
      else if (A.isrm())
	if (A.isherm()) Rank2Update<add>(alpha,x,y,A.Adjoint());
	else Rank2Update<add>(alpha,x,y,A.Transpose());
      else if (!(A.iscm() 
#ifdef BLAS
	    && A.stepj()>0
#endif
	    )) {
	if (A.isherm()) {
	  HermMatrix<T,Lower,ColMajor> AA(A.size());
	  Rank2Update<false>(alpha,x,y,AA.View());
	  if (add) A += AA;
	  else A = AA;
	} else {
	  SymMatrix<T,Lower,ColMajor> AA(A.size());
	  Rank2Update<false>(alpha,x,y,AA.View());
	  if (add) A += AA;
	  else A = AA;
	}
      } else 
#ifdef BLAS
	if (x.isconj() || x.step()<0 || SameStorage(x,A)) {
	  if (y.isconj() || y.step()<0 || SameStorage(y,A)) {
	    if (IMAG(alpha) == RealType(T)(0)) {
	      Vector<Tx> xx = REAL(alpha)*x;
	      Vector<Ty> yy = y;
	      if (!add) A.Zero();
	      BlasRank2Update(T(1),xx,yy,A);
	    } else {
	      Vector<T> xx = alpha*x;
	      Vector<Ty> yy = y;
	      if (!add) A.Zero();
	      BlasRank2Update(T(1),xx,yy,A);
	    }
	  } else {
	    if (IMAG(alpha) == RealType(T)(0)) {
	      Vector<Tx> xx = REAL(alpha)*x;
	      if (!add) A.Zero();
	      BlasRank2Update(T(1),xx,y,A);
	    } else {
	      Vector<T> xx = alpha*x;
	      if (!add) A.Zero();
	      BlasRank2Update(T(1),xx,y,A);
	    }
	  }
	} else {
	  if (y.isconj() || y.step()<0 || SameStorage(y,A)) {
	    if (IMAG(alpha) == RealType(T)(0)) {
	      Vector<Ty> yy = REAL(alpha)*y;
	      if (!add) A.Zero();
	      BlasRank2Update(T(1),x,yy,A);
	    } else {
	      Vector<T> yy = CONJ(alpha)*y;
	      if (!add) A.Zero();
	      BlasRank2Update(T(1),x,yy,A);
	    }
	  } else {
	    if (!add) A.Zero();
	    BlasRank2Update(alpha,x,y,A);
	  }
	}
#else
      NonBlasRank2Update<add>(alpha,x,y,A);
#endif
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*(ABS(alpha)*Norm(x0)*Norm(y0)+Norm(A0))) {
      cerr<<"Rank2Update: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"x = "<<Type(x)<<"  step = "<<x.step()<<"  "<<x0<<endl;
      cerr<<"y = "<<Type(y)<<"  step = "<<y.step()<<"  "<<y0<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
  }

#define InstFile "TMV_Rank2_VVS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


