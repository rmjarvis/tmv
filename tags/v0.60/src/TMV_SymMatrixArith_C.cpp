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
#include "TMV_SymMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  // 
  // Rank2Update
  //

  template <bool cx, bool cy, bool ha, bool rm, bool add, class T, class Tx, class Ty> 
    inline void RowRank2Update(
	const GenVector<Tx>& x, const GenVector<Ty>& y,
	const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(A.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(rm == A.isrm());
      TMVAssert(cx == x.isconj());
      TMVAssert(cy == y.isconj());
      TMVAssert(ha == A.isherm());
#ifdef XTEST
      RealType(T) NormA = Norm(A);
      RealType(T) Normx = Norm(x);
      RealType(T) Normy = Norm(y);
      RealType(T) eps = RealType(T)(2)*Normx*Normy;
      if (add) eps += NormA;
      eps *= A.size() * Epsilon<T>();
#endif

      const size_t si = A.stepi();
      const size_t sj = A.stepj();
      const size_t N = A.size();
      const Tx*const x0 = x.cptr();
      const Ty*const y0 = y.cptr();

      const Tx* xi = x0;
      const Ty* yi = y0;

      T A00;
      if (*xi == Tx(0) || *yi == Ty(0)) A00 = T(0);
      else if (ha) A00 = RealType(T)(2) * REAL(*xi * CONJ(*yi));
      else A00 = RealType(T)(2) * (cx?CONJ(*xi):*xi) * (cy?CONJ(*yi):*yi);
      ++xi; ++yi;
      T* Arowi = A.ptr()+si;

      for (size_t i=1;i<N;++i,++xi,++yi,Arowi+=si) {
	// A.row(i,0,i+1) += ax * y.SubVector(0,i+1) + ay * x.SubVector(0,i+1);
	if (*xi != Tx(0)) {
	  const Tx xival = cx ? CONJ(*xi) : *xi;
	  T* Aij = Arowi;
	  const Ty* yj = y0;
	  if (*yi != Ty(0)) {
	    const Ty yival = cy ? CONJ(*yi) : *yi;
	    const Tx* xj = x0;
	    for(size_t j=i+1;j>0;--j,++yj,++xj,(rm?++Aij:Aij+=sj)) {
	      T temp = xival * (ha==cy ? *yj : CONJ(*yj));
	      temp += yival * (ha==cx ? *xj : CONJ(*xj));
	      if (add) *Aij += temp;
	      else *Aij = temp;
	    }
	  } else {
	    for(size_t j=i+1;j>0;--j,++yj,(rm?++Aij:Aij+=sj)) {
	      const T temp = xival * (ha==cy ? *yj : CONJ(*yj));
	      if (add) *Aij += temp;
	      else *Aij = temp;
	    }
	  }
	} else if (*yi != Ty(0)) {
	  const Ty yival = cy ? CONJ(*yi) : *yi;
	  T* Aij = Arowi;
	  const Tx* xj = x0;
	  for(size_t j=i+1;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
	    const T temp = yival * (ha==cx ? *xj : CONJ(*xj));
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	}
      }
      if (add) *A.ptr() += A00;
      else *A.ptr() = A00;
      if (ha && IsComplex(T())) {
#ifdef XTEST
        TMVAssert(NormInf(A.diag().Imag()) <= eps);
#endif
	A.diag().Imag().Zero();
      }
    }

  template <bool cx, bool cy, bool ha, bool cm, bool add, class T, class Tx, class Ty> 
    inline void ColRank2Update(
	const GenVector<Tx>& x, const GenVector<Ty>& y,
	const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(A.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(cm == A.iscm());
      TMVAssert(cx == x.isconj());
      TMVAssert(cy == y.isconj());
      TMVAssert(ha == A.isherm());
#ifdef XTEST
      RealType(T) NormA = Norm(A);
      RealType(T) Normx = Norm(x);
      RealType(T) Normy = Norm(y);
      RealType(T) eps = RealType(T)(2)*Normx*Normy;
      if (add) eps += NormA;
      eps *= A.size() * Epsilon<T>();
#endif

      const size_t si = cm ? 1 : A.stepi();
      const size_t ds = A.stepj() + si;
      const size_t N = A.size();
      const Tx* xj = x.cptr()+N-1;
      const Ty* yj = y.cptr()+N-1;
      T* Ajj = A.ptr()+(N-1)*ds;

      for (size_t jj=N,Nmj=1;jj>0;--jj,++Nmj,--xj,--yj,Ajj-=ds) {
	// Nmj = N-j
	// A.col(j,j,N) += (A.isherm() ? CONJ(*yj) : *yj) * x.SubVector(j,N);
	// A.col(j,j,N) += (A.isherm() ? CONJ(*xj) : *xj) * y.SubVector(j,N);
	if (*yj!=Tx(0)) {
	  const Ty yjval = (ha==cy) ? *yj : CONJ(*yj);
	  T* Aij = Ajj;
	  const Tx* xi = xj;
	  if (*xj!=Tx(0)) {
	    const Tx xjval = (ha==cx) ? *xj : CONJ(*xj);
	    const Ty* yi = yj;
	    for(size_t i=Nmj;i>0;--i,++xi,++yi,(cm?++Aij:Aij+=si)) {
	      T temp = (cx ? CONJ(*xi) : *xi) * yjval;
	      temp += (cy ? CONJ(*yi) : *yi) * xjval;
	      if (add) *Aij += temp;
	      else *Aij = temp;
	    }
	  } else {
	    for(size_t i=Nmj;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
	      const T temp = (cx ? CONJ(*xi) : *xi) * yjval;
	      if (add) *Aij += temp;
	      else *Aij = temp;
	    }
	  }
	} else if (*xj!=Tx(0)) {
	  const Tx xjval = (ha==cx) ? *xj : CONJ(*xj);
	  T* Aij = Ajj;
	  const Ty* yi = yj;
	  for(size_t i=Nmj;i>0;--i,++yi,(cm?++Aij:Aij+=si)) {
	    const T temp = (cy ? CONJ(*yi) : *yi) * xjval;
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	}
      }
      if (ha && IsComplex(T())) {
#ifdef XTEST
        TMVAssert(NormInf(A.diag().Imag()) <= eps);
#endif
	A.diag().Imag().Zero();
      }
    }

  template <bool cx, bool cy, bool add, class T, class Tx, class Ty> 
    inline void UnitARank2Update(const GenVector<Tx>& x,
	const GenVector<Ty>& y, const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(A.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(x.step() == 1);
      TMVAssert(y.step() == 1);
      TMVAssert(cx == x.isconj());
      TMVAssert(cy == y.isconj());

      if (A.isherm())
	if (A.isrm()) RowRank2Update<cx,cy,true,true,add>(x,y,A);
	else if (A.iscm()) ColRank2Update<cx,cy,true,true,add>(x,y,A);
	else RowRank2Update<cx,cy,true,false,add>(x,y,A);
      else
	if (A.isrm()) RowRank2Update<cx,cy,false,true,add>(x,y,A);
	else if (A.iscm()) ColRank2Update<cx,cy,false,true,add>(x,y,A);
	else RowRank2Update<cx,cy,false,false,add>(x,y,A);
    }

  template <bool add, class T, class Tx, class Ty> inline void NonBlasRank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(A.size() > 0);

    if (A.uplo() == Upper) 
      return NonBlasRank2Update<add>(alpha,x,y,A.issym()?A.Transpose():A.Adjoint());
    else if (A.isconj())
      return NonBlasRank2Update<add>(CONJ(alpha),x.Conjugate(),y.Conjugate(),
	  A.Conjugate());
    else {
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);

      if (x.step() != 1 || alpha != T(1)) {
	if (x.step() == 1 && y.step() != 1) {
	  if (IMAG(alpha) == RealType(T)(0)) {
	    Vector<Ty> yy = REAL(alpha)*y;
	    if (x.isconj())
	      UnitARank2Update<true,false,add>(x,yy,A);
	    else
	      UnitARank2Update<false,false,add>(x,yy,A);
	  } else {
	    Vector<T> yy = alpha*y;
	    if (x.isconj())
	      UnitARank2Update<true,false,add>(x,yy,A);
	    else
	      UnitARank2Update<false,false,add>(x,yy,A);
	  }
	} else {
	  if (IMAG(alpha) == RealType(T)(0)) {
	    Vector<Tx> xx = REAL(alpha)*x;
	    if (y.step() == 1)
	      if (y.isconj())
		UnitARank2Update<false,true,add>(xx,y,A);
	      else
		UnitARank2Update<false,false,add>(xx,y,A);
	    else {
	      Vector<Ty> yy = y;
	      UnitARank2Update<false,false,add>(xx,yy,A);
	    }
	  } else {
	    Vector<T> xx = alpha*x;
	    if (y.step() == 1)
	      if (y.isconj())
		UnitARank2Update<false,true,add>(xx,y,A);
	      else
		UnitARank2Update<false,false,add>(xx,y,A);
	    else {
	      Vector<Ty> yy = y;
	      UnitARank2Update<false,false,add>(xx,yy,A);
	    }
	  }
	}
      } else {
	if (x.isconj())
	  if (y.isconj())
	    UnitARank2Update<true,true,add>(x,y,A);
	  else
	    UnitARank2Update<true,false,add>(x,y,A);
	else
	  if (y.isconj())
	    UnitARank2Update<false,true,add>(x,y,A);
	  else
	    UnitARank2Update<false,false,add>(x,y,A);
      }
    }
  }

#ifdef BLAS
  template <class T, class Tx, class Ty> inline void BlasRank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  { NonBlasRank2Update<true>(alpha,x,y,A); }
#ifdef INST_DOUBLE
  template <> inline void BlasRank2Update(
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
    TMVAssert(A.isherm());

    int n=A.size();
    int xs=x.step();
    int ys=y.step();
    int lda=A.stepj();
    BLASNAME(dsyr2) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	BLASV(n),BLASV(alpha),BLASP(x.cptr()),BLASV(xs),
	BLASP(y.cptr()),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
  }
  template <> inline void BlasRank2Update(
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
#endif
#ifdef INST_FLOAT
  template <> inline void BlasRank2Update(
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
    TMVAssert(A.isherm());

    int n=A.size();
    int xs=x.step();
    int ys=y.step();
    int lda=A.stepj();
    BLASNAME(ssyr2) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	BLASV(n),BLASV(alpha),BLASP(x.cptr()),BLASV(xs),
	BLASP(y.cptr()),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
  }
  template <> inline void BlasRank2Update(
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
    //cerr<<"Start Rank2Update: alpha = "<<alpha<<endl;
    //cerr<<"add = "<<add<<endl;
    //cerr<<"x = "<<Type(x)<<"  step = "<<x.step()<<"  "<<x0<<endl;
    //cerr<<"y = "<<Type(y)<<"  step = "<<y.step()<<"  "<<y0<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
#endif

    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    if (alpha != T(0) && A.size() > 0) {
#ifdef BLAS
      if (IsComplex(T()) && (IsReal(Tx()) || IsReal(Ty())) && add)
	BlasRank2Update(alpha,x,y,A);
      else if (A.isconj()) 
	Rank2Update<add>(CONJ(alpha),x.Conjugate(),y.Conjugate(),
	    A.Conjugate());
      else if (A.isrm())
	if (A.isherm()) Rank2Update<add>(alpha,x,y,A.Adjoint());
	else Rank2Update<add>(alpha,x,y,A.Transpose());
      else if (A.iscm() && A.stepj()>0) {
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
      } else {
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

#define InstFile "TMV_SymMatrixArith_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


