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
#include "TMV_SymMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  // 
  // Rank1Update
  //

  template <bool a1, bool cx, bool ha, bool rm, bool add, class T, class Tx, class Ta> 
    inline void RowRank1Update(const Ta alpha,
	const GenVector<Tx>& x, const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(x.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(rm == A.isrm());
      TMVAssert(a1 == (alpha == Ta(1)));
      TMVAssert(x.step()==1);
      TMVAssert(rm == A.isrm());
      TMVAssert(cx == x.isconj());
      TMVAssert(ha == A.isherm());
#ifdef XTEST
      RealType(T) NormA = Norm(A);
      RealType(T) Normx = Norm(x);
      RealType(T) eps = std::abs(alpha)*Normx*Normx;
      if (add) eps += NormA;
      eps *= A.size() * Epsilon<T>();
#endif

      const size_t si = A.stepi();
      const size_t sj = A.stepj();
      const size_t N = A.size();
      const Tx* xi = x.cptr();
      const Tx*const x0 = x.cptr();

      T A00;
      if (*xi == Tx(0)) A00 = T(0);
      else if (ha) A00 = alpha * NORM(*xi);
      else if (cx) A00 = alpha * CONJ(*xi * *xi);
      else A00 = alpha * (*xi * *xi);
      ++xi;
      T* Arowi = A.ptr()+si;

      for (size_t i=1;i<N;++i,++xi,Arowi+=si) if (*xi != Tx(0)) {
	// A.row(i,0,i+1) += ax * x.SubVector(0,i+1);
	T* Aij = Arowi;
	const Tx* xj = x0;

	if (a1 || IsReal(alpha)) {
	  Tx axi = cx ? CONJ(*xi) : *xi;
	  if (!a1) axi *= REAL(alpha);
	  for(size_t j=i+1;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
	    const T temp = axi * (ha==cx ? *xj : CONJ(*xj));
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	} else {
	  T axi = alpha * (cx ? CONJ(*xi) : *xi);
	  for(size_t j=i+1;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
	    const T temp = axi * (ha==cx ? *xj : CONJ(*xj));
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	}
      }

      // Do this at end in case A.ptr = x.ptr, so don't mess up x
      if (add) *A.ptr() += A00; 
      else *A.ptr() = A00;
      
      if (ha && IsComplex(T())) {
#ifdef XTEST
        TMVAssert(NormInf(A.diag().Imag()) <= eps);
#endif
	A.diag().Imag().Zero();
      }
    }

  template <bool a1, bool cx, bool ha, bool cm, bool add, class T, class Tx, class Ta> 
    inline void ColRank1Update(const Ta alpha,
	const GenVector<Tx>& x, const SymMatrixView<T>& A)
    {
      TMVAssert(x.step()==1);
      TMVAssert(x.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(a1 == (alpha == Ta(1)));
      TMVAssert(x.step() == 1);
      TMVAssert(cx == x.isconj());
      TMVAssert(cm == A.iscm());
      TMVAssert(ha == A.isherm());
#ifdef XTEST
      RealType(T) NormA = Norm(A);
      RealType(T) Normx = Norm(x);
      RealType(T) eps = std::abs(alpha)*Normx*Normx;
      if (add) eps += NormA;
      eps *= A.size() * Epsilon<T>();
#endif

      const size_t si = cm ? 1 : A.stepi();
      const size_t ds = A.stepj()+si;
      const size_t N = A.size();
      const Tx* xj = x.cptr()+N-1;
      T* Ajj = A.ptr()+(N-1)*ds;

      for (size_t jj=N,Nmj=1;jj>0;--jj,++Nmj,--xj,Ajj-=ds) if (*xj!=Tx(0)) {
	// Nmj = N-j
	// A.col(j,j,N) += *xj * x.SubVector(j,N);
	T* Aij = Ajj;
	const Tx* xi = xj;
	if (a1 || IsReal(alpha)) {
	  Tx axj = (ha!=cx) ? CONJ(*xj) : *xj;
	  if (!a1) axj *= REAL(alpha);
	  for(size_t i=Nmj;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
	    const T temp = axj * (cx ? CONJ(*xi) : *xi);
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	} else {
	  T axj = alpha * ((ha!=cx) ? CONJ(*xj) : *xj);
	  for(size_t i=Nmj;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
	    const T temp = axj * (cx ? CONJ(*xi) : *xi);
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

  template <bool a1, bool cx, bool add, class T, class Ta, class Tx>
    inline void DoRank1Update(const Ta alpha,
	const GenVector<Tx>& x, const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(x.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(a1 == (alpha == Ta(1)));
      TMVAssert(x.step() == 1);
      TMVAssert(cx == x.isconj());

      if (A.isherm()) 
	if (A.isrm()) RowRank1Update<a1,cx,true,true,add>(alpha,x,A);
	else if (A.iscm()) ColRank1Update<a1,cx,true,true,add>(alpha,x,A);
	else RowRank1Update<a1,cx,true,false,add>(alpha,x,A);
      else
	if (A.isrm()) RowRank1Update<a1,cx,false,true,add>(alpha,x,A);
	else if (A.iscm()) ColRank1Update<a1,cx,false,true,add>(alpha,x,A);
	else RowRank1Update<a1,cx,false,false,add>(alpha,x,A);
    }

  template <bool add, class T, class Tx> inline void NonBlasRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(x.size() > 0);

    if (A.uplo() == Upper) 
      return NonBlasRank1Update<add>(alpha,x,A.issym()?A.Transpose():A.Adjoint());
    else if (A.isconj()) 
      return NonBlasRank1Update<add>(CONJ(alpha),x.Conjugate(),A.Conjugate());
    else if (x.step() != 1) {
      Vector<T> xx = x;
      if (IMAG(alpha) == RealType(T)(0)) 
	if (REAL(alpha) == RealType(T)(1)) 
	  DoRank1Update<true,false,add>(REAL(alpha),xx,A);
	else
	  DoRank1Update<false,false,add>(REAL(alpha),xx,A);
      else 
	DoRank1Update<false,false,add>(alpha,xx,A);
    } else {
      if (x.isconj())
	if (IMAG(alpha) == RealType(T)(0)) 
	  if (REAL(alpha) == RealType(T)(1)) 
	    DoRank1Update<true,true,add>(REAL(alpha),x,A);
	  else
	    DoRank1Update<false,true,add>(REAL(alpha),x,A);
	else 
	  DoRank1Update<false,true,add>(alpha,x,A);
      else
	if (IMAG(alpha) == RealType(T)(0)) 
	  if (REAL(alpha) == RealType(T)(1)) 
	    DoRank1Update<true,false,add>(REAL(alpha),x,A);
	  else
	    DoRank1Update<false,false,add>(REAL(alpha),x,A);
	else 
	  DoRank1Update<false,false,add>(alpha,x,A);
    }
  }

#ifdef BLAS
  template <class T, class Tx> inline void BlasRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  { NonBlasRank1Update<true>(alpha,x,A); }
#ifdef INST_DOUBLE
  template <> inline void BlasRank1Update(
      const double alpha, const GenVector<double>& x,
      const SymMatrixView<double>& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(x.step() > 0);
    TMVAssert(A.iscm());

    int n=A.size();
    int xs=x.step();
    int lda=A.stepj();
    BLASNAME(dsyr) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	BLASV(n),BLASV(alpha),BLASP(x.cptr()),BLASV(xs),
	BLASP(A.ptr()),BLASV(lda) BLAS1);
  }
  template <> inline void BlasRank1Update(
      const std::complex<double> alpha,
      const GenVector<std::complex<double> >& x, 
      const SymMatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != double(0));
    TMVAssert(IMAG(alpha)==double(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(x.step() > 0);
    TMVAssert(A.iscm());

#ifndef ELAP
    if (A.issym() && x.step() != 1) {
      Vector<std::complex<double> > xx = x;
      return BlasRank1Update(alpha,xx,A);
    } else
#endif
    if (A.isherm()) {
      TMVAssert(IMAG(alpha)==double(0));
      int n=A.size();
      int xs=x.step();
      int lda=A.stepj();
      double ralpha = REAL(alpha);
      BLASNAME(zher) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	  BLASV(n),BLASV(ralpha),BLASP(x.cptr()),BLASV(xs),
	  BLASP(A.ptr()),BLASV(lda) BLAS1);
    } else {
      int n = A.size();
      int lda = A.stepj();
#ifdef ELAP
      int xs = x.step();
      LAPNAME(zsyr) (LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
	  LAPV(n),LAPP(&alpha),LAPP(x.cptr()),LAPV(xs),
	  LAPP(A.ptr()),LAPV(lda) LAP1);
#else
      int k=1;
      std::complex<double> beta(1);
      BLASNAME(zsyrk) (BLASCM 
	  A.uplo()==Upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
	  BLASV(n),BLASV(k),BLASP(&alpha),BLASP(x.cptr()),BLASV(n),
	  BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
#endif
    }
  }
#endif
#ifdef INST_FLOAT
  template <> inline void BlasRank1Update(
      const float alpha, const GenVector<float>& x,
      const SymMatrixView<float>& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(x.step() > 0);
    TMVAssert(A.iscm());

    int n=A.size();
    int xs=x.step();
    int lda=A.stepj();
    BLASNAME(ssyr) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	BLASV(n),BLASV(alpha),BLASP(x.cptr()),BLASV(xs),
	BLASP(A.ptr()),BLASV(lda) BLAS1);
  }
  template <> inline void BlasRank1Update(
      const std::complex<float> alpha,
      const GenVector<std::complex<float> >& x, 
      const SymMatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != float(0));
    TMVAssert(IMAG(alpha)==float(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(x.step() > 0);
    TMVAssert(A.iscm());

#ifndef ELAP
    if (A.issym() && x.step() != 1) {
      Vector<std::complex<float> > xx = x;
      return BlasRank1Update(alpha,xx,A);
    } else
#endif
    if (A.isherm()) {
      TMVAssert(IMAG(alpha)==float(0));
      int n=A.size();
      int xs=x.step();
      int lda=A.stepj();
      float ralpha = REAL(alpha);
      BLASNAME(cher) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	  BLASV(n),BLASV(ralpha),BLASP(x.cptr()),BLASV(xs),
	  BLASP(A.ptr()),BLASV(lda) BLAS1);
    } else {
      int n = A.size();
      int lda = A.stepj();
#ifdef ELAP
      int xs = x.step();
      LAPNAME(csyr) (LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
	  LAPV(n),LAPP(&alpha),LAPP(x.cptr()),LAPV(xs),
	  LAPP(A.ptr()),LAPV(lda) LAP1);
#else
      int k=1;
      std::complex<float> beta(1);
      BLASNAME(csyrk) (BLASCM 
	  A.uplo()==Upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
	  BLASV(n),BLASV(k),BLASP(&alpha),BLASP(x.cptr()),BLASV(n),
	  BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
#endif
    }
  }
#endif 
#endif // BLAS

  template <bool add, class T, class Tx> void Rank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
    // A (+)= alpha * x * xT
  {
#ifdef XDEBUG
    Vector<Tx> x0 = x;
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    if (!add) A2.Zero();
    if (A.isherm())
      A2 += (alpha*x0^x0.Conjugate());
    else 
      A2 += (alpha*x0^x0);
    //cerr<<"Start Rank1Update: alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<"  "<<x<<endl;
#endif
#ifdef XTEST
    TMVAssert(!add || A.HermOK());
#endif

    TMVAssert(A.size() == x.size());
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    if (alpha != T(0) && A.size() > 0) {
#ifdef BLAS
      if (IsComplex(T()) && IsReal(Tx()) && add)
	BlasRank1Update(alpha,x,A);
      else if (A.isrm())
	return Rank1Update<add>(alpha,x,A.issym()?A.Transpose():A.Adjoint());
      else if (A.isconj()) 
	return Rank1Update<add>(CONJ(alpha),x.Conjugate(),A.Conjugate());
      else if (A.iscm() && A.stepj()>0) {
	if (x.isconj() || x.step() < 0 || SameStorage(x,A)) {
	  Vector<Tx> xx = x;
	  if (!add) A.Zero();
	  BlasRank1Update(alpha,xx,A);
	} else {
	  if (!add) A.Zero();
	  BlasRank1Update(alpha,x,A);
	}
      } else {
	if (A.isherm()) {
	  HermMatrix<T,Lower,ColMajor> AA(A.size());
	  Rank1Update<false>(alpha,x,AA.View());
	  if (add) A += AA;
	  else A = AA;
	} else {
	  SymMatrix<T,Lower,ColMajor> AA(A.size(),T(0));
	  Rank1Update<false>(alpha,x,AA.View());
	  if (add) A += AA;
	  else A = AA;
	}
      }
#else
      NonBlasRank1Update<add>(alpha,x,A);
#endif
    }
  
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
#ifdef XDEBUG
    //cerr<<"Done Rank1\n";
    if (Norm(A-A2) > 0.001*(ABS(alpha)*SQR(Norm(x0))+Norm(A0))) {
      cerr<<"Rank1Update: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"x = "<<Type(x)<<"  step = "<<x.step()<<"  "<<x<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymMatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv

