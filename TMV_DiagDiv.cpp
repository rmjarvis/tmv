
#include "TMV.h"
#include "TMV_Diag.h"
#include "TMV_VectorArith_Inline.h"

namespace tmv {

  template <class T> bool DiagDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenDiagMatrix<T>* dm = dynamic_cast<const GenDiagMatrix<T>*>(&m);
    TMVAssert(dm);
    if (fout) {
      *fout << "M = "<<*dm<<endl;
      *fout << "D = "<<*itsm<<endl;
    }
    RealType(T) nm = Norm(*itsm-*dm);
    nm /= Norm(*itsm);
    if (fout) {
      *fout << "Norm(M-D)/Norm(D) = "<<nm<<endl;
    }
    return nm < Epsilon<T>();
  }

  template <class T1, StepItType S1, ConjItType C1, class T2, StepItType S2> 
    void DoDiagLDivEq(CVIt<T1,S1,C1> dit, VIt<T2,S2,NonConj> vit,
	const size_t size)
  {
    for(size_t i=size;i>0;--i,++dit,++vit) {
      if (*dit == T1(0)) tmv_error("Division by 0 in DiagMatrix /=");
      if (*dit != T1(1)) *vit /= *dit;
    }
  }

  template <class T1, class T2> void DiagLDivEq(
      const GenDiagMatrix<T1>& d, const VectorView<T2>& v)
  { 
    TMVAssert(v.size() == d.size());
    TMVAssert(v.ct()==NonConj);
    if (v.size() > 0)  
      if (d.diag().step()==1)
	if (d.diag().isconj())
	  if (v.step()==1) 
	    DoDiagLDivEq(CVIt<T1,Unit,Conj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(v.begin()),v.size());
	  else
	    DoDiagLDivEq(CVIt<T1,Unit,Conj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(v.begin()),v.size());
	else
	  if (v.step()==1) 
	    DoDiagLDivEq(CVIt<T1,Unit,NonConj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(v.begin()),v.size());
	  else
	    DoDiagLDivEq(CVIt<T1,Unit,NonConj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(v.begin()),v.size());
      else
	if (d.diag().isconj())
	  if (v.step()==1) 
	    DoDiagLDivEq(CVIt<T1,Step,Conj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(v.begin()),v.size());
	  else
	    DoDiagLDivEq(CVIt<T1,Step,Conj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(v.begin()),v.size());
	else
	  if (v.step()==1) 
	    DoDiagLDivEq(CVIt<T1,Step,NonConj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(v.begin()),v.size());
	  else
	    DoDiagLDivEq(CVIt<T1,Step,NonConj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(v.begin()),v.size());
  }

  template <class T1, StepItType S1, ConjItType C1, class T2> 
    void RowMajorDiagLDivEq(CVIt<T1,S1,C1> dit, const MatrixView<T2>& m)
    {
      TMVAssert(m.isrm());
      TMVAssert(m.ct()==NonConj);
      T2* mptr = m.ptr();
      for(size_t i=0;i<m.colsize();++i,++dit,mptr+=m.stepi()) {
	if (*dit == T1(0)) tmv_error("Division by 0 in DiagMatrix /=");
	// m.row(i) /= *dit;
	if (*dit != T1(1)) 
	  DoMultXV(RealType(T1)(1)/(*dit),VIt<T2,Unit,NonConj>(mptr,1
		FIRSTLAST1(m.first,m.last) ), m.rowsize());
      }
    }

  template <class T1, class T2> void DoRowMajorDiagLDivEq(
      const GenDiagMatrix<T1>& d, const MatrixView<T2>& m)
  {
    TMVAssert(d.size()==m.colsize());
    TMVAssert(m.isrm());
    TMVAssert(m.ct()==NonConj);

    if (d.diag().step()==1) 
      if (d.diag().isconj())
	RowMajorDiagLDivEq(CVIt<T1,Unit,Conj>(d.diag().begin()),m);
      else
	RowMajorDiagLDivEq(CVIt<T1,Unit,NonConj>(d.diag().begin()),m);
    else
      if (d.diag().isconj())
	RowMajorDiagLDivEq(CVIt<T1,Step,Conj>(d.diag().begin()),m);
      else
	RowMajorDiagLDivEq(CVIt<T1,Step,NonConj>(d.diag().begin()),m);
  }

  template <StepItType S2, class T1, class T2> void DoColDiagLDivEq(
      const GenDiagMatrix<T1>& d, const MatrixView<T2>& m)
  {
    TMVAssert(d.size()==m.colsize());
    TMVAssert(m.iscm() || S2 == Step);
    TMVAssert(m.ct() == NonConj);

    T2* mptr = m.ptr();
    if (d.diag().step()==1)
      if (d.diag().isconj())
	for(size_t j=0;j<m.rowsize();++j,mptr+=m.stepj()) 
	  DoDiagLDivEq(CVIt<T1,Unit,Conj>(d.diag().begin()),
	      VIt<T2,S2,NonConj>(mptr,m.stepi() FIRSTLAST1(m.first,m.last) ),
	      d.size());
      else
	for(size_t j=0;j<m.rowsize();++j,mptr+=m.stepj()) 
	  DoDiagLDivEq(CVIt<T1,Unit,NonConj>(d.diag().begin()),
	      VIt<T2,S2,NonConj>(mptr,m.stepi() FIRSTLAST1(m.first,m.last) ),
	      d.size());
    else
      if (d.diag().isconj())
	for(size_t j=0;j<m.rowsize();++j,mptr+=m.stepj()) 
	  DoDiagLDivEq(CVIt<T1,Step,Conj>(d.diag().begin()),
	      VIt<T2,S2,NonConj>(mptr,m.stepi() FIRSTLAST1(m.first,m.last) ),
	      d.size());
      else
	for(size_t j=0;j<m.rowsize();++j,mptr+=m.stepj()) 
	  DoDiagLDivEq(CVIt<T1,Step,NonConj>(d.diag().begin()),
	      VIt<T2,S2,NonConj>(mptr,m.stepi() FIRSTLAST1(m.first,m.last) ),
	      d.size());
  }

  template <class T1, class T2> void DiagLDivEq(
      const GenDiagMatrix<T1>& d, const MatrixView<T2>& m)
  { 
    TMVAssert(m.colsize() == d.size());
    if (m.isconj()) DiagLDivEq(d.Conjugate(),m.Conjugate());
    else if (m.colsize() > 0 && m.rowsize() > 0)
      if (m.isconj()) DiagLDivEq(d.Conjugate(),m.Conjugate());
      else if (m.rowsize() == 1) DiagLDivEq(d,m.col(0));
      else if (m.isrm()) DoRowMajorDiagLDivEq(d,m);
      else if (m.iscm()) DoColDiagLDivEq<Unit>(d,m);
      else DoColDiagLDivEq<Step>(d,m);
  }

  template <class T> T DiagDiv<T>::Det() const 
  {
    if (!calcdet) {
      if (itsm->diag().step() == 1) 
	if (itsm->diag().isconj()) {
	  CVIt<T,Unit,Conj> it=itsm->diag().begin();
	  for(size_t i=itsm->diag().size();i>0;--i,++it) det *= *it;
	} else {
	  CVIt<T,Unit,NonConj> it=itsm->diag().begin();
	  for(size_t i=itsm->diag().size();i>0;--i,++it) det *= *it;
	}
      else 
	if (itsm->diag().isconj()) {
	  CVIt<T,Step,Conj> it=itsm->diag().begin();
	  for(size_t i=itsm->diag().size();i>0;--i,++it) det *= *it;
	} else {
	  CVIt<T,Step,NonConj> it=itsm->diag().begin();
	  for(size_t i=itsm->diag().size();i>0;--i,++it) det *= *it;
	}
      calcdet = true;
    }
    return det;
  }

  template <class T> void DiagDiv<T>::DInverse(
      const DiagMatrixView<T>& minv) const 
  {
    TMVAssert(minv.size() == itsm->size());
    minv.SetToIdentity();
    if (minv.diag().isconj()) 
      DiagLDivEq(itsm->Conjugate(),minv.diag().Conjugate());
    else
      DiagLDivEq(*itsm,minv.diag());
  }

  template <class T> void DiagDiv<T>::DInverseATA(
      const DiagMatrixView<T>& minv) const
  {
    TMVAssert(minv.size() == itsm->size());
    DInverse(minv);
    if (minv.diag().step()==1) {
      VIt<T,Unit,NonConj> mit = minv.diag().isconj() ? 
	minv.diag().Conjugate().begin() : minv.diag().begin();
      for(size_t i=minv.size();i>0;--i,++mit) *mit = NORM(*mit);
    } else {
      VIt<T,Step,NonConj> mit = minv.diag().isconj() ? 
	minv.diag().Conjugate().begin() : minv.diag().begin();
      for(size_t i=minv.size();i>0;--i,++mit) *mit = NORM(*mit);
    }
  }

#define InstFile "TMV_DiagDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


