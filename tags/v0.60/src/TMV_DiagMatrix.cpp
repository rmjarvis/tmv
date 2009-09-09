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



#include "TMV_Vector.h"
#include "TMV_VIt.h"
#include "TMV_VectorRE.h"
#include "TMV_Matrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_DiagMatrixArith.h"
#include "TMV_TriMatrix.h"
#include <ostream>

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_VectorArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <class T> void GenDiagMatrix<T>::DoAssignToM(
      const MatrixView<RealType(T)>& m2) const
  {
    if (size() > 1) {
      if (SameStorage(diag(),m2)) {
	UpperTriMatrixViewOf(m2).OffDiag().Zero();
	LowerTriMatrixViewOf(m2).OffDiag().Zero();
      } else {
	m2.Zero();
      }
    }
    m2.diag() = diag();
  }

  template <class T> void GenDiagMatrix<T>::DoAssignToM(
      const MatrixView<ComplexType(T)>& m2) const
  {
    if (size() > 1) {
      if (SameStorage(diag(),m2)) {
	UpperTriMatrixViewOf(m2).OffDiag().Zero();
	LowerTriMatrixViewOf(m2).OffDiag().Zero();
      } else {
	m2.Zero();
      }
    }
    m2.diag() = diag();
  }

  template <class T> void GenDiagMatrix<T>::DoAssignToU(
      const UpperTriMatrixView<RealType(T)>& m2) const
  {
    if (size() > 1) {
      if (SameStorage(diag(),m2)) {
	m2.OffDiag().Zero();
      } else {
	m2.Zero();
      }
    }
    m2.diag() = diag();
  }

  template <class T> void GenDiagMatrix<T>::DoAssignToU(
      const UpperTriMatrixView<ComplexType(T)>& m2) const
  {
    if (size() > 1) {
      if (SameStorage(diag(),m2)) {
	m2.OffDiag().Zero();
      } else {
	m2.Zero();
      }
    }
    m2.diag() = diag();
  }

  template <class T> void GenDiagMatrix<T>::DoAssignToL(
      const LowerTriMatrixView<RealType(T)>& m2) const
  {
    if (size() > 1) {
      if (SameStorage(diag(),m2)) {
	m2.OffDiag().Zero();
      } else {
	m2.Zero();
      }
    }
    m2.diag() = diag();
  }

  template <class T> void GenDiagMatrix<T>::DoAssignToL(
      const LowerTriMatrixView<ComplexType(T)>& m2) const
  {
    if (size() > 1) {
      if (SameStorage(diag(),m2)) {
	m2.OffDiag().Zero();
      } else {
	m2.Zero();
      }
    }
    m2.diag() = diag();
  }

  template <class T> T GenDiagMatrix<T>::Det() const
  {
    const T* di = diag().cptr();
    const int ds = diag().step();
    T det(1);
    if (ds == 1)
      for(size_t i=size();i>0;--i,++di) det *= *di;
    else
      for(size_t i=size();i>0;--i,di+=ds) det *= *di;
    if (diag().isconj()) return CONJ(det);
    else return det;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::NewCopy() const
  {
    auto_ptr<BaseMatrix<T> > a(new DiagMatrix<T>(*this));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::NewView() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(View()));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::NewTranspose() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(Transpose()));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::NewConjugate() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(Conjugate()));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::NewAdjoint() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(Adjoint()));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::NewInverse() const
  {
    auto_ptr<DiagMatrix<T> > minv(new DiagMatrix<T>(*this));
    minv->InvertSelf();
    BaseMatrix<T>* ret1 = minv.release();
    auto_ptr<BaseMatrix<T> > ret(ret1);
    return ret;
  }


  template <class T> class SingularDiagMatrix :
    public Singular
  {
    public:
      DiagMatrix<T> A;

      SingularDiagMatrix(const GenDiagMatrix<T>& _A) :
	Singular("DiagMatrix"), A(_A) {}
      ~SingularDiagMatrix() throw() {}
      void Write(std::ostream& os) const throw()
      {
	Singular::Write(os);
	os<<A<<std::endl;
      }
  };

  template <class T, IndexStyle I> 
    const DiagMatrixView<T,I>& DiagMatrixView<T,I>::InvertSelf() const
    {
      T* di = diag().ptr();
      const int dstep = diag().step();

      if (dstep == 1)
	for(size_t i=size();i>0;--i,++di) {
	  if (*di == T(0))
	    throw SingularDiagMatrix<T>(*this);
	  if (IMAG(*di) == RealType(T)(0))
	    *di = RealType(T)(1) / REAL(*di);
	  else
	    *di = RealType(T)(1) / *di;
	}
      else {
	for(size_t i=size();i>0;--i,di+=dstep) {
	  if (*di == T(0))
	    throw SingularDiagMatrix<T>(*this);
	  if (IMAG(*di) == RealType(T)(0))
	    *di = RealType(T)(1) / REAL(*di);
	  else
	    *di = RealType(T)(1) / *di;
	}
      }
      return *this;
    }


  template <class T> template <class T1> void GenDiagMatrix<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    bool ss = SameStorage(diag(),minv);
    if (!ss) minv.Zero();
    (DiagMatrixViewOf(minv.diag()) = *this).InvertSelf();
    if (ss && size() > 1) {
      UpperTriMatrixViewOf(minv).OffDiag().Zero();
      LowerTriMatrixViewOf(minv).OffDiag().Zero();
    }
  }

  template <class T> template <class T1> void GenDiagMatrix<T>::DoInverse(
      const DiagMatrixView<T1>& minv) const
  { (minv = *this).InvertSelf(); }

  template <class T> void GenDiagMatrix<T>::DoInverseATA(
      const DiagMatrixView<T>& ata) const
  {
    Inverse(ata);
    T* mi = ata.diag().ptr();
    const int ds = ata.diag().step();
    if (ds==1)
      for(size_t i=size();i>0;--i,++mi) *mi = NORM(*mi);
    else
      for(size_t i=size();i>0;--i,mi+=ds) *mi = NORM(*mi);
  }

  template <class T> void GenDiagMatrix<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    ata.Zero();
    InverseATA(DiagMatrixViewOf(ata.diag()));
  }

  template <class T> QuotXD<T,T> GenDiagMatrix<T>::QInverse() const
  { return QuotXD<T,T>(T(1),*this); }
  
#define CT std::complex<T>

  template <bool cd, class T, class Td> inline void DoDiagLDivEq1(
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
      for(size_t i=v.size();i>0;--i,++di,++vi) {
	if (*di == Td(0))
	  throw SingularDiagMatrix<Td>(d);
	if (IMAG(*di) == RealType(Td)(0)) {
	  if (REAL(*di) != RealType(Td)(1))
	    *vi /= REAL(*di);
	} else
	  *vi /= (cd?CONJ(*di):*di);
      }
    else
      for(size_t i=v.size();i>0;--i,di+=dstep,vi+=vstep) {
	if (*di == Td(0))
	  throw SingularDiagMatrix<Td>(d);
	if (IMAG(*di) == RealType(Td)(0)) {
	  if (REAL(*di) != RealType(Td)(1))
	    *vi /= REAL(*di);
	} else
	  *vi /= (cd?CONJ(*di):*di);
      }
  }

  template <class T, class Td> inline void DoDiagLDivEq(
      const GenDiagMatrix<Td>& d, const VectorView<T>& v)
  { 
    if (d.diag().isconj()) DoDiagLDivEq1<true>(d,v);
    else DoDiagLDivEq1<false>(d,v);
  }
  template <class T> inline void DoDiagLDivEq(
      const GenDiagMatrix<CT>& , const VectorView<T>& )
  { TMVAssert(FALSE); }

  template <class T> template <class T1> void GenDiagMatrix<T>::DoLDivEq(
      const VectorView<T1>& v) const
  {
#ifdef XDEBUG
    DiagMatrix<T> d0(*this);
    Vector<T1> v0(v);
#endif

    TMVAssert(v.size() == size());
    TMVAssert(v.ct() == NonConj);

    if (v.size() > 0)
      if (v.isconj()) Conjugate().DoLDivEq(v.Conjugate());
      else DoDiagLDivEq(*this,v);

#ifdef XDEBUG
    Vector<T1> v1 = d0*v;
    if (tmv::Norm(v1-v0) > 0.001*tmv::Norm(v0)) {
      cerr<<"DiagLDivEq v: \n";
      cerr<<"d = "<<Type(*this)<<"  "<<d0<<endl;
      cerr<<"v = "<<Type(v)<<"  "<<v0<<endl;
      cerr<<"-> v/d = "<<v<<endl;
      cerr<<"d*(v/d) = "<<v1<<endl;
      abort();
    }
#endif
  }

  template <class T> template <class T1, class T0> 
    void GenDiagMatrix<T>::DoLDiv(
	const GenVector<T1>& v1, const VectorView<T0>& v0) const
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

  template <bool rm, bool cd, class T, class Td> inline void RowDiagLDivEq(
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
    const size_t M = m.colsize();
    const size_t N = m.rowsize();
    const RealType(Td) one(1);
    const RealType(Td) zero(0);

    for(size_t i=M;i>0;--i,di+=dstep,mrowi+=stepi) {
      T* mij = mrowi;
      if (*di == Td(0))
	throw SingularDiagMatrix<Td>(d);
      else if (IMAG(*di) == zero) {
	RealType(Td) invdi = one/REAL(*di);
	for(size_t j=N;j>0;--j,(rm?++mij:mij+=stepj))
	  *mij *= invdi;
      }
      else {
	Td invdi = one/(cd?CONJ(*di):*di);
	for(size_t j=N;j>0;--j,(rm?++mij:mij+=stepj))
	  *mij *= invdi;
      }
    }
  }

  template <bool cm, bool cd, class T, class Td> inline void ColDiagLDivEq(
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
    const RealType(Td) one(1);

    if (step == 1)
      for(size_t i=d.size();i>0;--i,++di,++invdi) {
	if (*di == Td(0))
	  throw SingularDiagMatrix<Td>(d);
	*invdi = one/(cd?CONJ(*di):*di);
      }
    else
      for(size_t i=d.size();i>0;--i,di+=step,++invdi) {
	if (*di == Td(0))
	  throw SingularDiagMatrix<Td>(d);
	*invdi = one/(cd?CONJ(*di):*di);
      }
    m = invd*m;
  }

  template <class T, class Td> inline void DoDiagLDivEq(
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
  template <class T> inline void DoDiagLDivEq(
      const GenDiagMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T> template <class T1> void GenDiagMatrix<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(size() == m.colsize());

#ifdef XDEBUG
    DiagMatrix<T> d0(*this);
    Matrix<T1> m0(m);
#endif

    if (m.colsize() > 0 && m.rowsize() > 0)
      if (m.isconj()) Conjugate().DoLDivEq(m.Conjugate());
      else if (m.rowsize() == 1) DoLDivEq(m.col(0));
      else DoDiagLDivEq(*this,m);

#ifdef XDEBUG
    Matrix<T1> m1 = d0*m;
    if (tmv::Norm(m1-m0) > 0.001*tmv::Norm(m0)) {
      cerr<<"DiagLDivEq m: \n";
      cerr<<"d = "<<Type(*this)<<"  "<<d0<<endl;
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
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

  template <class T> void GenDiagMatrix<T>::Write(std::ostream& os) const
  {
    const int sd = diag().step();
    const T* di = diag().cptr();
    os << size() <<' '<< size() << std::endl;
    for (size_t i=0,nmi=size();nmi>0;++i,--nmi,di+=sd) {
      os << "( ";
      for(size_t k=i;k>0;--k) os <<' '<<T(0)<<' ';
      if (diag().isconj())
	os <<' '<<CONJ(*di)<<' ';
      else
	os <<' '<<*di<<' ';
      for(size_t k=nmi-1;k>0;--k) os <<' '<<T(0)<<' ';
      os << " )\n";
    }
  }

  template <class T> void GenDiagMatrix<T>::Write(std::ostream& os,
      RealType(T) thresh) const
  {
    const int sd = diag().step();
    const T* di = diag().cptr();
    os << size() <<' '<< size() << std::endl;
    for (size_t i=0,nmi=size();nmi>0;++i,--nmi,di+=sd) {
      os << "( ";
      for(size_t k=i;k>0;--k) os <<' '<<T(0)<<' ';
      if (diag().isconj())
	os <<' '<<(ABS(*di)<thresh ? T(0) : CONJ(*di))<<' ';
      else
	os <<' '<<(ABS(*di)<thresh ? T(0) : *di)<<' ';
      for(size_t k=nmi-1;k>0;--k) os <<' '<<T(0)<<' ';
      os << " )\n";
    }
  }

  template <class T> class DiagMatrixReadError :
    public ReadError
  {
    public :
      size_t i;
      mutable auto_ptr<DiagMatrix<T> > m;
      char exp,got;
      size_t s;
      bool is, iseof, isbad;

      inline DiagMatrixReadError(std::istream& _is) throw() :
	ReadError("DiagMatrix"),
	i(0), m(0), exp(0), got(0), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline DiagMatrixReadError(size_t _i, const GenDiagMatrix<T>& _m,
	  char _e, char _g, size_t _s,
	  bool _is, bool _iseof, bool _isbad) throw() :
	ReadError("DiagMatrix"),
	i(_i), m(new DiagMatrix<T>(_m)), exp(_e), got(_g), s(_s),
	is(_is), iseof(_iseof), isbad(_isbad) {}
      inline DiagMatrixReadError(const GenDiagMatrix<T>& _m,
	  std::istream& _is, size_t _s) throw() :
	ReadError("DiagMatrix"),
	i(0), m(new DiagMatrix<T>(_m)), s(_s),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline DiagMatrixReadError(std::istream& _is, char _e, char _g) throw() :
	ReadError("DiagMatrix"),
	i(0), m(0), exp(_e), got(_g), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline DiagMatrixReadError(const DiagMatrixReadError<T>& rhs) :
	i(rhs.i), m(rhs.m), exp(rhs.exp), got(rhs.got), s(rhs.s),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
      virtual inline ~DiagMatrixReadError() throw() {}

      virtual void Write(std::ostream& os) const throw();
  };

  template <class T> void DiagMatrixReadError<T>::Write(
      std::ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for DiagMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (m.get() && s != m->size()) {
      os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
    }
    if (!is) {
      if (iseof) {
	os<<"Input stream reached end-of-file prematurely.\n";
      } else if (isbad) {
	os<<"Input stream is corrupted.\n";
      } else {
	os<<"Input stream cannot read next character.\n";
      }
    }
    if (m.get()) {
      os<<"The portion of the DiagMatrix which was successfully read is: \n";
      ConstDiagMatrixView<T> mm = m->View();
      os<<"( ";
      for(size_t ii=0;ii<i;++ii)
	os<<' '<<mm(ii,ii)<<' ';
      os<<" )\n";
    }
  }

  template <class T, IndexStyle I> std::istream& operator>>(std::istream& is,
      auto_ptr<DiagMatrix<T,I> >& m)
  {
    char d;
    is >> d;
    if (!is || d != 'D') 
      throw DiagMatrixReadError<T>(is,'D',d);
    size_t size;
    is >> size;
    if (!is) 
      throw DiagMatrixReadError<T>(is);
    m.reset(new DiagMatrix<T,I>(size));
    try {
      m->diag().Read(is);
    }
    catch (VectorReadError<T>& ve) {
      throw DiagMatrixReadError<T>(ve.i,*m,ve.exp,ve.got,ve.s,
	  ve.is,ve.iseof,ve.isbad);
    }
    return is;
  }

  template <class T> std::istream& operator>>(std::istream& is,
      const DiagMatrixView<T>& m)
  {
    char d;
    is >> d;
    if (!is || d != 'D') 
      throw DiagMatrixReadError<T>(is,'D',d);
    size_t s;
    is >> s;
    if (!is) 
      throw DiagMatrixReadError<T>(is);
    if (m.size() != s)
      throw DiagMatrixReadError<T>(m,is,s);
    TMVAssert(m.size() == s);
    try {
      m.diag().Read(is);
    }
    catch (VectorReadError<T>& ve) {
      throw DiagMatrixReadError<T>(ve.i,m,ve.exp,ve.got,ve.s,
	  ve.is,ve.iseof,ve.isbad);
    }
    return is;
  }

#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


