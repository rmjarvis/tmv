
#include "TMV.h"
#include "TMV_Diag.h"
#include "TMV_DiagDiv.h"

//#define XDEBUG

namespace tmv {

  template <class T> inline bool DiagDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenDiagMatrix<T>* dm = dynamic_cast<const GenDiagMatrix<T>*>(&m);
    TMVAssert(dm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*dm)<<"  "<<*dm<<endl;
      *fout << "D = "<<*itsm<<endl;
    }
    RealType(T) nm = Norm(*itsm-*dm);
    nm /= Norm(*itsm);
    if (fout) {
      *fout << "Norm(M-D)/Norm(D) = "<<nm<<endl;
    }
    return nm < Epsilon<T>();
  }

  template <bool cd, class T, class Td> inline void DoDiagLDivEq(
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

  template <class T, class Td> void DiagLDivEq(
      const GenDiagMatrix<Td>& d, const VectorView<T>& v)
  { 
#ifdef XDEBUG
    DiagMatrix<Td> d0 = d;
    Vector<T> v0 = v;
#endif

    TMVAssert(v.size() == d.size());
    TMVAssert(v.ct()==NonConj);
    if (v.size() > 0)  
      if (d.diag().isconj())
	DoDiagLDivEq<true>(d,v);
      else
	DoDiagLDivEq<false>(d,v);
#ifdef XDEBUG
    Vector<T> v1 = d0*v;
    if (Norm(v1-v0) > 0.001*Norm(v0)) {
      cerr<<"DiagLDivEq v: \n";
      cerr<<"d = "<<Type(d)<<"  "<<d<<endl;
      cerr<<"v = "<<Type(v)<<"  "<<v0<<endl;
      cerr<<"-> v/d = "<<v<<endl;
      cerr<<"d*(v/d) = "<<v1<<endl;
      abort();
    }
#endif
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
    TMVAssert(d.diag().step() == 1);
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

  template <class T, class Td> void DiagLDivEq(
      const GenDiagMatrix<Td>& d, const MatrixView<T>& m)
  {
    TMVAssert(d.size() == m.colsize());

#ifdef XDEBUG
    DiagMatrix<Td> d0 = d;
    Matrix<T> m0 = m;
#endif

    if (m.isconj()) DiagLDivEq(d.Conjugate(),m.Conjugate());
    else if (m.colsize() > 0 && m.rowsize() > 0) 
      if (m.rowsize() == 1) DiagLDivEq(d,m.col(0));
      else if (d.diag().isconj()) 
        if (m.isrm()) RowDiagLDivEq<true,true>(d,m);
	else if (m.iscm()) ColDiagLDivEq<true,true>(d,m);
	else if (m.colsize() > m.rowsize()) ColDiagLDivEq<false,true>(d,m);
	else RowDiagLDivEq<false,true>(d,m);
      else 
        if (m.isrm()) RowDiagLDivEq<true,false>(d,m);
	else if (m.iscm()) ColDiagLDivEq<true,false>(d,m);
	else if (m.colsize() > m.rowsize()) ColDiagLDivEq<false,false>(d,m);
	else RowDiagLDivEq<false,false>(d,m);
#ifdef XDEBUG
    Matrix<T> m1 = d0*m;
    if (Norm(m1-m0) > 0.001*Norm(m0)) {
      cerr<<"DiagLDivEq m: \n";
      cerr<<"d = "<<Type(d)<<"  "<<d0<<endl;
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m/d = "<<m<<endl;
      cerr<<"d*(m/d) = "<<m1<<endl;
      abort();
    }
#endif
  }

  template <class T> T DiagDiv<T>::Det() const 
  {
    if (!calcdet) {
      const T* di = itsm->diag().cptr();
      const int ds = itsm->diag().step();
      if (ds == 1)
	for(size_t i=itsm->diag().size();i>0;--i,++di) det *= *di;
      else 
	for(size_t i=itsm->diag().size();i>0;--i,di+=ds) det *= *di;
      if (itsm->diag().isconj()) det = CONJ(det);
      calcdet = true;
    }
    return det;
  }

  template <class T> void DiagDiv<T>::Inverse(
      const DiagMatrixView<T>& minv) const 
  {
    TMVAssert(minv.size() == itsm->size());
    minv.SetToIdentity();
    if (minv.diag().isconj()) 
      DiagLDivEq(itsm->Conjugate(),minv.diag().Conjugate());
    else
      DiagLDivEq(*itsm,minv.diag());
  }

  template <class T> void DiagDiv<T>::InverseATA(
      const DiagMatrixView<T>& minv) const
  {
    TMVAssert(minv.size() == itsm->size());
    Inverse(minv);
    T* mi = minv.diag().ptr();
    const int ds = minv.diag().step();
    if (ds==1) 
      for(size_t i=minv.size();i>0;--i,++mi) *mi = NORM(*mi);
    else
      for(size_t i=minv.size();i>0;--i,mi+=ds) *mi = NORM(*mi);
  }

#define InstFile "TMV_DiagDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


