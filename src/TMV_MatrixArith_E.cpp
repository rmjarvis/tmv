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
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  // 
  // Rank1Update
  //

  template <bool cx, bool cy, bool y1, bool cm, bool add, class T, class Tx, class Ty> 
    inline void ColRank1Update(
	const GenVector<Tx>& x, const GenVector<Ty>& y,
	const MatrixView<T>& A)
    {
      TMVAssert(A.colsize() == x.size());
      TMVAssert(A.rowsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(!A.isrm());
      TMVAssert(x.step() == 1);
      TMVAssert(cm == A.iscm());
      TMVAssert(y1 == (y.step() == 1));
      TMVAssert(cx == x.isconj());
      TMVAssert(cy == y.isconj());

      const Ty* yj = y.cptr();
      const Tx*const xptr = x.cptr();
      T* Acolj = A.ptr();
      const int sj = A.stepj();
      const int si = (cm ? 1 : A.stepi());
      const int ys = y.step();
      const size_t M = A.colsize();
      const size_t N = A.rowsize();

      for (size_t j=N; j>0; --j,(y1?++yj:yj+=ys),Acolj+=sj) if (*yj!=Ty(0)) {
	T* Aij = Acolj;
	const Tx* xi = xptr;
	for (size_t i=M; i>0; --i,++xi,(cm?++Aij:Aij+=si)) {
	  const T temp = (cx ? CONJ(*xi) : *xi) * (cy ? CONJ(*yj) : *yj);
	  if (add) *Aij += temp;
	  else *Aij = temp;
	}
      }
    }

  template <bool cx, bool y1, bool add, class T, class Tx, class Ty> 
    inline void UnitARank1Update(
	const GenVector<Tx>& x,
	const GenVector<Ty>& y, const MatrixView<T>& A)
    {
      TMVAssert(A.colsize() == x.size());
      TMVAssert(A.rowsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(!A.isrm());
      TMVAssert(x.step() == 1);
      TMVAssert(cx == x.isconj());
      TMVAssert(y1 == (y.step() == 1));

      if (A.iscm()) 
	if (y.isconj())
	  ColRank1Update<cx,true,y1,true,add>(x,y,A);
	else
	  ColRank1Update<cx,false,y1,true,add>(x,y,A);
      else
	if (y.isconj())
	  ColRank1Update<cx,true,y1,false,add>(x,y,A);
	else
	  ColRank1Update<cx,false,y1,false,add>(x,y,A);
    }

  template <bool add, class T, class Tx, class Ty> 
    inline void NonBlasRank1Update(
	const T alpha, const GenVector<Tx>& x,
	const GenVector<Ty>& y, const MatrixView<T>& A)
    {
      TMVAssert(A.colsize() == x.size());
      TMVAssert(A.rowsize() == y.size());
      TMVAssert(alpha != T(0));
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(!A.isrm());

      if (x.step() != 1 || alpha != T(1)) {
	if (x.step() == 1 && y.size() < x.size()) {
	  if (IMAG(alpha) == RealType(T)(0)) {
	    Vector<Ty> yy = REAL(alpha)*y;
	    if (x.isconj())
	      UnitARank1Update<true,true,add>(x,yy,A);
	    else
	      UnitARank1Update<false,true,add>(x,yy,A);
	  } else {
	    Vector<T> yy = alpha*y;
	    if (x.isconj())
	      UnitARank1Update<true,true,add>(x,yy,A);
	    else
	      UnitARank1Update<false,true,add>(x,yy,A);
	  }
	} else {
	  if (IMAG(alpha) == RealType(T)(0)) {
	    Vector<Tx> xx = REAL(alpha)*x;
	    if (y.step() == 1)
	      UnitARank1Update<false,true,add>(xx,y,A);
	    else
	      UnitARank1Update<false,false,add>(xx,y,A);
	  } else {
	    Vector<T> xx = alpha*x;
	    if (y.step() == 1)
	      UnitARank1Update<false,true,add>(xx,y,A);
	    else
	      UnitARank1Update<false,false,add>(xx,y,A);
	  }
	}
      } else {
	if (x.isconj())
	  if (y.step() == 1)
	    UnitARank1Update<true,true,add>(x,y,A);
	  else
	    UnitARank1Update<true,false,add>(x,y,A);
	else
	  if (y.step() == 1)
	    UnitARank1Update<false,true,add>(x,y,A);
	  else
	    UnitARank1Update<false,false,add>(x,y,A);
      }
    }

#ifdef BLAS
  template <class T, class Tx, class Ty> inline void BlasRank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  { NonBlasRank1Update<true>(alpha,x,y,A); }
#ifdef INST_DOUBLE
  template <> inline void BlasRank1Update(
      const double alpha, const GenVector<double>& x,
      const GenVector<double>& y, const MatrixView<double>& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0); 
    TMVAssert(y.step() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = A.colsize();
    int n = A.rowsize();
    int xs = x.step();
    int ys = y.step();
    int lda = A.stepj();
    BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(alpha),
	BLASP(x.cptr()),BLASV(xs),BLASP(y.cptr()),BLASV(ys),
	BLASP(A.ptr()),BLASV(lda));
  }
  template <> inline void BlasRank1Update(
      const std::complex<double> alpha,
      const GenVector<std::complex<double> >& x, 
      const GenVector<std::complex<double> >& y,
      const MatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj || y.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = A.colsize();
    int n = A.rowsize();
    int xs = x.step();
    int ys = y.step();
    int lda = A.stepj();
    if (x.isconj()) {
#ifdef CBLAS
      BLASNAME(zgerc) (BLASRM BLASV(n),BLASV(m),BLASP(&alpha),
	  BLASP(y.cptr()),BLASV(ys),BLASP(x.cptr()),BLASV(xs),
	  BLASP(A.ptr()),BLASV(lda));
#else
      Vector<std::complex<double> > xx = alpha*x;
      xs = 1;
      std::complex<double> alpha2(1);
      BLASNAME(zgeru) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha2),
	  BLASP(xx.cptr()),BLASV(xs),BLASP(y.cptr()),BLASV(ys),
	  BLASP(A.ptr()),BLASV(lda));
#endif
    }
    else if (y.isconj())
      BLASNAME(zgerc) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(xs),BLASP(y.cptr()),BLASV(ys),
	  BLASP(A.ptr()),BLASV(lda));
    else
      BLASNAME(zgeru) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(xs),BLASP(y.cptr()),BLASV(ys),
	  BLASP(A.ptr()),BLASV(lda));
  }
#endif
#ifdef INST_FLOAT
  template <> inline void BlasRank1Update(
      const float alpha, const GenVector<float>& x,
      const GenVector<float>& y, const MatrixView<float>& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0); 
    TMVAssert(y.step() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = A.colsize();
    int n = A.rowsize();
    int xs = x.step();
    int ys = y.step();
    int lda = A.stepj();
    BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(alpha),
	BLASP(x.cptr()),BLASV(xs),BLASP(y.cptr()),BLASV(ys),
	BLASP(A.ptr()),BLASV(lda));
  }
  template <> inline void BlasRank1Update(
      const std::complex<float> alpha,
      const GenVector<std::complex<float> >& x, 
      const GenVector<std::complex<float> >& y,
      const MatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj || y.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = A.colsize();
    int n = A.rowsize();
    int xs = x.step();
    int ys = y.step();
    int lda = A.stepj();
    if (x.isconj()) {
#ifdef CBLAS
      BLASNAME(cgerc) (BLASRM BLASV(n),BLASV(m),BLASP(&alpha),
	  BLASP(y.cptr()),BLASV(ys),BLASP(x.cptr()),BLASV(xs),
	  BLASP(A.ptr()),BLASV(lda));
#else
      Vector<std::complex<float> > xx = alpha*x;
      xs = 1;
      std::complex<float> alpha2(1);
      BLASNAME(cgeru) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha2),
	  BLASP(xx.cptr()),BLASV(xs),BLASP(y.cptr()),BLASV(ys),
	  BLASP(A.ptr()),BLASV(lda));
#endif
    }
    else if (y.isconj())
      BLASNAME(cgerc) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(xs),BLASP(y.cptr()),BLASV(ys),
	  BLASP(A.ptr()),BLASV(lda));
    else
      BLASNAME(cgeru) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(xs),BLASP(y.cptr()),BLASV(ys),
	  BLASP(A.ptr()),BLASV(lda));
  }
#endif
#endif // BLAS

  template <bool add, class T, class Tx, class Ty> void Rank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
    // A (+)= beta + alpha * x * yT
  {
#ifdef XDEBUG
    //cerr<<"Rank1Update: alpha = "<<alpha<<endl;
    //cerr<<"add = "<<add<<endl;
    //cerr<<"x = "<<Type(x)<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<"  "<<y<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    Vector<Tx> x0 = x;
    Vector<Ty> y0 = y;
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    for(size_t i=0;i<x.size();i++) for(size_t j=0;j<y.size();j++) 
      if (add)
	A2(i,j) += alpha*x0(i)*y0(j);
      else
	A2(i,j) = alpha*x0(i)*y0(j);
#endif

    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());

    if (A.colsize() > 0 && A.rowsize() > 0) {
      if (alpha == T(0)) {
	if (!add) A.Zero();
      } else {
	if (A.rowsize() == 1) {
	  Ty y0 = y.isconj() ? CONJ(*y.cptr()) : *y.cptr();
	  if (add)
	    // A.col(0) = alpha * y(0) * x;
	    if (alpha == T(1)) A.col(0) += y0 * x;
	    else A.col(0) += alpha * y0 * x;
	  else
	    // A.col(0) += alpha * y(0) * x;
	    if (alpha == T(1)) A.col(0) = y0 * x;
	    else A.col(0) = alpha * y0 * x;
	} else if (A.colsize() == 1)
	  Rank1Update<add>(alpha,y,x,A.Transpose());
	else if (A.isconj()) 
	  Rank1Update<add>(CONJ(alpha),x.Conjugate(),y.Conjugate(),
	      A.Conjugate());
	else if (A.isrm())
	  Rank1Update<add>(alpha,y,x,A.Transpose());
#ifdef BLAS
	else if (IsComplex(T()) && (IsReal(Tx()) || IsReal(Ty())) && add)
	  BlasRank1Update(alpha,x,y,A);
	else if (!((A.iscm() && A.stepj()>0))) {
	  Matrix<T,ColMajor> A2(A);
	  Rank1Update<add>(alpha,x,y,A2.View());
	  A = A2;
	} else {
	  if (x.step() < 0 || SameStorage(x,A)) {
	    if (y.step() < 0 || SameStorage(y,A)) {
	      if (IMAG(alpha) == RealType(T)(0)) {
		if (x.size() <= y.size()) {
		  Vector<Tx> xx = REAL(alpha)*x;
		  Vector<Ty> yy = y;
		  if (!add) A.Zero();
		  BlasRank1Update(T(1),xx,yy,A);
		} else {
		  Vector<Tx> xx = x;
		  Vector<Ty> yy = REAL(alpha)*y;
		  if (!add) A.Zero();
		  BlasRank1Update(T(1),xx,yy,A);
		}
	      } else {
		if (x.size() <= y.size()) {
		  Vector<T> xx = alpha*x;
		  Vector<Ty> yy = y;
		  if (!add) A.Zero();
		  BlasRank1Update(T(1),xx,yy,A);
		} else {
		  Vector<Tx> xx = x;
		  Vector<T> yy = alpha*y;
		  if (!add) A.Zero();
		  BlasRank1Update(T(1),xx,yy,A);
		}
	      }
	    } else {
	      if (IMAG(alpha) == RealType(T)(0)) {
		Vector<Tx> xx = REAL(alpha)*x;
		if (!add) A.Zero();
		BlasRank1Update(T(1),xx,y,A);
	      } else {
		Vector<T> xx = alpha*x;
		if (!add) A.Zero();
		BlasRank1Update(T(1),xx,y,A);
	      }
	    }
	  } else {
	    if (y.step() < 0 || SameStorage(A,y)) {
	      if (IMAG(alpha) == RealType(T)(0)) {
		Vector<Ty> yy = REAL(alpha)*y;
		if (!add) A.Zero();
		BlasRank1Update(T(1),x,yy,A);
	      } else {
		Vector<T> yy = alpha*y;
		if (!add) A.Zero();
		BlasRank1Update(T(1),x,yy,A);
	      }
	    } else {
	      if (!add) A.Zero();
	      if (x.isconj() && y.isconj()) {
		if (IMAG(alpha) == RealType(T)(0)) {
		  if (x.size() <= y.size()) {
		    Vector<Tx> xx = REAL(alpha)*x;
		    BlasRank1Update(T(1),xx,y,A);
		  } else {
		    Vector<Ty> yy = REAL(alpha)*y;
		    BlasRank1Update(T(1),x,yy,A);
		  }
		} else {
		  if (x.size() <= y.size()) {
		    Vector<T> xx = alpha*x;
		    BlasRank1Update(T(1),xx,y,A);
		  } else {
		    Vector<T> yy = alpha*y;
		    BlasRank1Update(T(1),x,yy,A);
		  }
		}
	      } else {
		BlasRank1Update(alpha,x,y,A);
	      }
	    }
	  }
	}
#else
	else NonBlasRank1Update<add>(alpha,x,y,A);
#endif
      }
    }
  
#ifdef XDEBUG
    //cerr<<"Done Rank1Update: A->"<<A<<endl;
    if (Norm(A-A2) > 0.001*(ABS(alpha)*Norm(x0)*Norm(y0)+
	  (add?Norm(A0):RealType(T)(0)))) {
      cerr<<"Rank1Update: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"x = "<<Type(x)<<"  step = "<<x.step()<<"  "<<x0<<endl;
      cerr<<"y = "<<Type(y)<<"  step = "<<y.step()<<"  "<<y0<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MatrixArith_E.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


