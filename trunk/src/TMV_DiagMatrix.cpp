///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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


//#define XDEBUG


#include <iostream>
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixIO.h"

namespace tmv {

  const int XX = UNKNOWN;

#if 0
  template <class T> 
  T Det(const ConstDiagMatrixView<T>& m) const
  {
    T signdet(1);
    RT logdet = LogDet(m,&signdet);
    if (signdet == T(0)) return T(0);
    else return signdet * EXP(logdet);
  }

  template <class T> 
  RT LogDet(const ConstDiagMatrixView<T>& m, T* sign) const
  {
    T s(1);
    RT logdet(0);
    const int n=m.size();
    for(int i=0;i<n;++i) {
      const T mi = m(i);
      if (mi == T(0)) { 
        logdet = LOG(REAL(mi)); // i.e. -inf
        if (sign) s = T(0); 
      }
      else {
        RT a = ABS(mi);
        logdet += LOG(a);
        if (sign) {
          if (IsReal(T())) {
            if (REAL(mi) < RT(0)) s = -s;
          } else {
            s *= (mi/a);
          }
        }
      }
    }
    if (sign) *sign = s;
    return logdet;
  }

  template <class T>
  void InvertSelf(DiagMatrixView<T> m) 
  {
    const int n = m.size();
    for(int i=0;i<n;++i) {
      const T mi = m(i);
      if (mi == T(0)) {
#ifdef NOTHROW
        { std::cerr<<"Singular DiagMatrix found\n"; exit(1); }
#else
        throw SingularDiagMatrix<T>(*this);
#endif
      }
      if (IMAG(mi) == RT(0))
        mi = RT(1) / REAL(mi);
      else
        mi = RT(1) / mi;
    }
    return *this;
  }

  template <class T> void GenDiagMatrix<T>::DoInverseATA(
      const DiagMatrixView<T>& ata) const
  {
    Inverse(ata);
    T* mi = ata.diag().ptr();
    const int ds = ata.diag().step();
    if (ds==1)
      for(int i=size();i>0;--i,++mi) {
#ifdef TMVFLDEBUG
        TMVAssert(mi >= ata.diag().first);
        TMVAssert(mi < ata.diag().last);
#endif
        *mi = NORM(*mi);
      }
    else
      for(int i=size();i>0;--i,mi+=ds) {
#ifdef TMVFLDEBUG
        TMVAssert(mi >= ata.diag().first);
        TMVAssert(mi < ata.diag().last);
#endif
        *mi = NORM(*mi);
      }
  }

#define CT std::complex<T>

  template <bool cd, class T, class Td> static void DoDiagLDivEq1(
      const GenDiagMatrix<Td>& d, const VectorView<T>& v)
  {
    TMVAssert(v.size() == d.size());
    TMVAssert(v.ct()==NonConj);
    TMVAssert(v.size() > 0);

    const Td* di = d.diag().cptr();
    T* vi = v.ptr();
    const int dstep = d.diag().step();
    const int vstep = v.step();

    if (dstep == 1 && vstep == 1)
      for(int i=v.size();i>0;--i,++di,++vi) {
#ifdef TMVFLDEBUG
        TMVAssert(vi >= v.first);
        TMVAssert(vi < v.last);
#endif
        if (*di == Td(0))
#ifdef NOTHROW
        { std::cerr<<"Singular DiagMatrix found\n"; exit(1); }
#else
        throw SingularDiagMatrix<Td>(d);
#endif
        if (IMAG(*di) == RT(0)) {
          if (REAL(*di) != RT(1))
            *vi /= REAL(*di);
        } else
          *vi /= (cd?CONJ(*di):*di);
      }
    else
      for(int i=v.size();i>0;--i,di+=dstep,vi+=vstep) {
#ifdef TMVFLDEBUG
        TMVAssert(vi >= v.first);
        TMVAssert(vi < v.last);
#endif
        if (*di == Td(0))
#ifdef NOTHROW
        { std::cerr<<"Singular DiagMatrix found\n"; exit(1); }
#else
        throw SingularDiagMatrix<Td>(d);
#endif
        if (IMAG(*di) == RT(0)) {
          if (REAL(*di) != RT(1))
            *vi /= REAL(*di);
        } else
          *vi /= (cd?CONJ(*di):*di);
      }
  }

  template <class T, class Td> static inline void DoDiagLDivEq(
      const GenDiagMatrix<Td>& d, const VectorView<T>& v)
  { 
    if (d.diag().isconj()) DoDiagLDivEq1<true>(d,v);
    else DoDiagLDivEq1<false>(d,v);
  }

  template <class T> template <class T1> void GenDiagMatrix<T>::DoLDivEq(
      const VectorView<T1>& v) const
  {
#ifdef XDEBUG
    DiagMatrix<T> d0(*this);
    Vector<T1> v0(v);
#endif

    TMVAssert(v.size() == size());

    if (v.size() > 0) {
      if (v.isconj()) DoDiagLDivEq(Conjugate(),v.Conjugate());
      else DoDiagLDivEq(*this,v);
    }

#ifdef XDEBUG
    Vector<T1> v1 = d0*v;
    if (Norm(v1-v0) > 0.001*Norm(v0)) {
      cerr<<"DiagLDivEq v: \n";
      cerr<<"d = "<<TypeText(*this)<<"  "<<d0<<endl;
      cerr<<"v = "<<TypeText(v)<<"  "<<v0<<endl;
      cerr<<"-> v/d = "<<v<<endl;
      cerr<<"d*(v/d) = "<<v1<<endl;
      abort();
    }
#endif
  }

  template <class T> template <class T1, class T0> 
  void GenDiagMatrix<T>::DoLDiv(
      const ConstVectorView<T1>& v1, const VectorView<T0>& v0) const
  {
    TMVAssert(v1.size() == size());
    TMVAssert(v0.size() == size());
    if (SameStorage(diag(),v0)) {
      DiagMatrix<T> temp = *this;
      temp.DoLDivEq(v0=v1);
    } else {
      DoLDivEq(v0=v1);
    }
  }

  template <bool rm, bool cd, class T, class Td> static void RowDiagLDivEq(
      const GenDiagMatrix<Td>& d, const MatrixView<T>& m)
  {
    TMVAssert(d.size() == m.colsize());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(rm == m.isrm());
    TMVAssert(cd == d.diag().isconj());

    const Td* di = d.diag().cptr();
    T* mrowi = m.ptr();
    const int dstep = d.diag().step();
    const int stepj = m.stepj();
    const int stepi = m.stepi();
    const int M = m.colsize();
    const int N = m.rowsize();

    for(int i=M;i>0;--i,di+=dstep,mrowi+=stepi) {
      T* mij = mrowi;
      if (*di == Td(0))
#ifdef NOTHROW
      { std::cerr<<"Singular DiagMatrix found\n"; exit(1); }
#else
      throw SingularDiagMatrix<Td>(d);
#endif
      else if (IMAG(*di) == RT(0)) {
        RT invdi = RT(1)/REAL(*di);
        for(int j=N;j>0;--j,(rm?++mij:mij+=stepj)) {
#ifdef TMVFLDEBUG
          TMVAssert(mij >= m.first);
          TMVAssert(mij < m.last);
#endif
          *mij *= invdi;
        }
      }
      else {
        Td invdi = RT(1)/(cd?CONJ(*di):*di);
        for(int j=N;j>0;--j,(rm?++mij:mij+=stepj)) {
#ifdef TMVFLDEBUG
          TMVAssert(mij >= m.first);
          TMVAssert(mij < m.last);
#endif
          *mij *= invdi;
        }
      }
    }
  }

  template <bool cm, bool cd, class T, class Td> static void ColDiagLDivEq(
      const GenDiagMatrix<Td>& d, const MatrixView<T>& m)
  {
    TMVAssert(d.size() == m.colsize());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(m.ct()==NonConj);
    TMVAssert(cm == m.iscm());
    TMVAssert(cd == d.diag().isconj());

    DiagMatrix<Td> invd(d.size());
    const Td* di = d.diag().cptr();
    const int step = d.diag().step();
    Td* invdi = invd.diag().ptr();

    if (step == 1)
      for(int i=d.size();i>0;--i,++di,++invdi) {
        if (*di == Td(0))
#ifdef NOTHROW
        { std::cerr<<"Singular DiagMatrix found\n"; exit(1); }
#else
        throw SingularDiagMatrix<Td>(d);
#endif
#ifdef TMVFLDEBUG
        TMVAssert(invdi >= invd.diag().first);
        TMVAssert(invdi < invd.diag().last);
#endif
        *invdi = RT(1)/(cd?CONJ(*di):*di);
      }
    else
      for(int i=d.size();i>0;--i,di+=step,++invdi) {
        if (*di == Td(0))
#ifdef NOTHROW
        { std::cerr<<"Singular DiagMatrix found\n"; exit(1); }
#else
        throw SingularDiagMatrix<Td>(d);
#endif
#ifdef TMVFLDEBUG
        TMVAssert(invdi >= invd.diag().first);
        TMVAssert(invdi < invd.diag().last);
#endif
        *invdi = RT(1)/(cd?CONJ(*di):*di);
      }
    m = invd*m;
  }

  template <class T, class Td> static void DoDiagLDivEq(
      const GenDiagMatrix<Td>& d, const MatrixView<T>& m)
  {
    if (d.diag().isconj())
      if (m.isrm()) RowDiagLDivEq<true,true>(d,m);
      else if (m.iscm()) ColDiagLDivEq<true,true>(d,m);
      else if (m.colsize() > m.rowsize()) ColDiagLDivEq<false,true>(d,m);
      else RowDiagLDivEq<false,true>(d,m);
    else
      if (m.isrm()) RowDiagLDivEq<true,false>(d,m);
      else if (m.iscm()) ColDiagLDivEq<true,false>(d,m);
      else if (m.colsize() > m.rowsize()) ColDiagLDivEq<false,false>(d,m);
      else RowDiagLDivEq<false,false>(d,m);
  }

  template <class T> template <class T1> void GenDiagMatrix<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(size() == m.colsize());

#ifdef XDEBUG
    DiagMatrix<T> d0(*this);
    Matrix<T1> m0(m);
#endif

    if (m.colsize() > 0 && m.rowsize() > 0) {
      if (m.isconj()) Conjugate().DoLDivEq(m.Conjugate());
      else if (m.rowsize() == 1) DoLDivEq(m.col(0));
      else DoDiagLDivEq(*this,m);
    }

#ifdef XDEBUG
    Matrix<T1> m1 = d0*m;
    if (Norm(m1-m0) > 0.001*Norm(m0)) {
      cerr<<"DiagLDivEq m: \n";
      cerr<<"d = "<<TypeText(*this)<<"  "<<d0<<endl;
      cerr<<"m = "<<TypeText(m)<<"  "<<m0<<endl;
      cerr<<"-> m/d = "<<m<<endl;
      cerr<<"d*(m/d) = "<<m1<<endl;
      abort();
    }
#endif
  }

  template <class T> template <class T1, class T0> 
  void GenDiagMatrix<T>::DoLDiv(
      const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
  {
    TMVAssert(m1.rowsize() == m0.rowsize());
    TMVAssert(m1.colsize() == size());
    TMVAssert(m0.colsize() == size());
    if (SameStorage(diag(),m0)) {
      DiagMatrix<T> temp = *this;
      temp.DoLDivEq(m0=m1);
    } else {
      DoLDivEq(m0=m1);
    }
  }

#undef CT
#endif

  template <class T, bool C>
  void InstWrite(std::ostream& os, const ConstDiagMatrixView<T,XX,C>& m)
  {
    if (m.step() == 1)
      InlineWrite(os,m.CMView()); 
    else
      InlineWrite(os,m); 
  }

  template <class T, bool C>
  void InstWrite(std::ostream& os, const ConstDiagMatrixView<T,XX,C>& m,
      RealType(T) thresh)
  {
    if (m.step() == 1)
      InlineWrite(os,m.CMView(),thresh); 
    else
      InlineWrite(os,m,thresh); 
  }

#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


