
#include "TMV.h"
#include "TMV_Diag.h"
#include "TMV_VectorArith_Inline.h"

namespace tmv {

  template <class T1, StepItType S1, ConjItType C1, class T2, StepItType S2> 
    void DoDiagLDivEq(CVIt<T1,S1,C1> dit, VIt<T2,S2,NonConj> vit,
	const size_t size)
  {
    const VIt<T2,S2,NonConj> _end = vit + size;
    for(; vit!=_end; ++dit,++vit) {
      if (*dit == T1(0)) tmv_error("Division by 0 in DiagMatrix /=");
      *vit /= *dit;
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
	  const CVIt<T,Unit,Conj> _end = itsm->diag().end();
	  for(CVIt<T,Unit,Conj> it=itsm->diag().begin(); it!=_end; ++it) 
	    det *= *it;
	} else {
	  const CVIt<T,Unit,NonConj> _end = itsm->diag().end();
	  for(CVIt<T,Unit,NonConj> it=itsm->diag().begin(); it!=_end; ++it) 
	    det *= *it;
	}
      else 
	if (itsm->diag().isconj()) {
	  const CVIt<T,Step,Conj> _end = itsm->diag().end();
	  for(CVIt<T,Step,Conj> it=itsm->diag().begin(); it!=_end; ++it) 
	    det *= *it;
	} else {
	  const CVIt<T,Step,NonConj> _end = itsm->diag().end();
	  for(CVIt<T,Step,NonConj> it=itsm->diag().begin(); it!=_end; ++it) 
	    det *= *it;
	}
      calcdet = true;
    }
    return det;
  }

  template <class T> DiagMatrix<T> DiagDiv<T>::DInverse() const 
  {
    DiagMatrix<T> temp(itsm->size());
    if (itsm->size() == 0) return temp;
    const VIt<T,Unit,NonConj> _end = temp.diag().end();
    VIt<T,Unit,NonConj> mit = temp.diag().begin();
    if (itsm->diag().step()==1) 
      if (itsm->diag().isconj())
	for(CVIt<T,Unit,Conj> dit=itsm->diag().begin(); mit!=_end;
	    ++dit,++mit) {
	  if (*dit == T(0)) tmv_error("Division by 0 in DiagMatrix DInverse");
	  *mit = RealType(T)(1)/(*dit);
	}
      else
	for(CVIt<T,Unit,NonConj> dit=itsm->diag().begin(); mit!=_end;
	    ++dit,++mit) {
	  if (*dit == T(0)) tmv_error("Division by 0 in DiagMatrix DInverse");
	  *mit = RealType(T)(1)/(*dit);
	}
    else
      if (itsm->diag().isconj())
	for(CVIt<T,Step,Conj> dit=itsm->diag().begin(); mit!=_end;
	    ++dit,++mit) {
	  if (*dit == T(0)) tmv_error("Division by 0 in DiagMatrix DInverse");
	  *mit = RealType(T)(1)/(*dit);
	}
      else
	for(CVIt<T,Step,NonConj> dit=itsm->diag().begin(); mit!=_end;
	    ++dit,++mit) {
	  if (*dit == T(0)) tmv_error("Division by 0 in DiagMatrix DInverse");
	  *mit = RealType(T)(1)/(*dit);
	}
    return temp;
  }

  template <class T> DiagMatrix<T> DiagDiv<T>::DInverseATA() const
  {
    DiagMatrix<T> temp(itsm->size());
    if (itsm->size() == 0) return temp;
    const VIt<T,Unit,NonConj> _end = temp.diag().end();
    VIt<T,Unit,NonConj> mit = temp.diag().begin();
    if (itsm->diag().step()==1) 
      if (itsm->diag().isconj())
	for(CVIt<T,Unit,Conj> dit=itsm->diag().begin(); mit!=_end;
	    ++dit,++mit) {
	  if (*dit == T(0)) tmv_error("Division by 0 in DiagMatrix DInvATA");
	  *mit = RealType(T)(1)/((*dit)*(*dit));
	}
      else
	for(CVIt<T,Unit,NonConj> dit=itsm->diag().begin(); mit!=_end;
	    ++dit,++mit) {
	  if (*dit == T(0)) tmv_error("Division by 0 in DiagMatrix DInvATA");
	  *mit = RealType(T)(1)/((*dit)*(*dit));
	}
    else 
      if (itsm->diag().isconj())
	for(CVIt<T,Step,Conj> dit=itsm->diag().begin(); mit!=_end;
	    ++dit,++mit) {
	  if (*dit == T(0)) tmv_error("Division by 0 in DiagMatrix DInvATA");
	  *mit = RealType(T)(1)/((*dit)*(*dit));
	}
      else
	for(CVIt<T,Step,NonConj> dit=itsm->diag().begin(); mit!=_end;
	    ++dit,++mit) {
	  if (*dit == T(0)) tmv_error("Division by 0 in DiagMatrix DInvATA");
	  *mit = RealType(T)(1)/((*dit)*(*dit));
	}
    return temp;
  }

#define InstFile "TMV_DiagDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


