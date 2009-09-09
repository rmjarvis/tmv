
#include "TMV.h"
#include "TMV_Diag.h"

namespace tmv {

  template <class T1, StepItType S1, ConjItType C1, class T2, StepItType S2> 
    void DoLU_LDivEq(CVIt<T1,S1,C1> dit, VIt<T2,S2,NonConj> vit,
	const size_t size)
  {
    const VIt<T2,S2,NonConj> _end = vit + size;
    for(; vit!=_end; ++dit,++vit) {
      if (*dit == T1(0)) tmv_error("Division by 0 in DiagMatrix /=");
      *vit /= *dit;
    }
  }

  template <class T1, class T2> void LU_LDivEq(
      const GenDiagMatrix<T1>& d, const VectorView<T2>& v)
  { 
    TMVAssert(v.size() == d.size());
    TMVAssert(v.ct()==NonConj);
    if (v.size() > 0)  
      if (d.diag().step()==1)
	if (d.diag().isconj())
	  if (v.step()==1) 
	    DoLU_LDivEq(CVIt<T1,Unit,Conj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(v.begin()),v.size());
	  else
	    DoLU_LDivEq(CVIt<T1,Unit,Conj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(v.begin()),v.size());
	else
	  if (v.step()==1) 
	    DoLU_LDivEq(CVIt<T1,Unit,NonConj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(v.begin()),v.size());
	  else
	    DoLU_LDivEq(CVIt<T1,Unit,NonConj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(v.begin()),v.size());
      else
	if (d.diag().isconj())
	  if (v.step()==1) 
	    DoLU_LDivEq(CVIt<T1,Step,Conj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(v.begin()),v.size());
	  else
	    DoLU_LDivEq(CVIt<T1,Step,Conj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(v.begin()),v.size());
	else
	  if (v.step()==1) 
	    DoLU_LDivEq(CVIt<T1,Step,NonConj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(v.begin()),v.size());
	  else
	    DoLU_LDivEq(CVIt<T1,Step,NonConj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(v.begin()),v.size());
  }

  template <class T1, StepItType S1, ConjItType C1, class T2> void RowMajorLU_LDivEq(
      CVIt<T1,S1,C1> dit, const MatrixView<T2>& m)
  {
    for(size_t i=0;i<m.colsize();++i,++dit) {
      if (*dit == T1(0)) tmv_error("Division by 0 in DiagMatrix /=");
      m.row(i) /= *dit;
    }
  }

  template <class T1, class T2> void LU_LDivEq(
      const GenDiagMatrix<T1>& d, const MatrixView<T2>& m)
  { 
    TMVAssert(m.colsize() == d.size());
    if (m.colsize() > 0 && m.rowsize() > 0)
      if (m.isconj()) LU_LDivEq(d.Conjugate(),m.Conjugate());
      else if (m.rowsize() == 1) LU_LDivEq(d,m.col(0));
      else if (m.isrm()) 
	if (d.diag().step()==1) 
	  if (d.diag().isconj())
	    RowMajorLU_LDivEq(CVIt<T1,Unit,Conj>(d.diag().begin()),m);
	  else
	    RowMajorLU_LDivEq(CVIt<T1,Unit,NonConj>(d.diag().begin()),m);
	else
	  if (d.diag().isconj())
	    RowMajorLU_LDivEq(CVIt<T1,Step,Conj>(d.diag().begin()),m);
	  else
	    RowMajorLU_LDivEq(CVIt<T1,Step,NonConj>(d.diag().begin()),m);
      else if (m.iscm())
	// MJ: Add Bptr stuff
	if (d.diag().step()==1)
	  if (d.diag().isconj())
	    for(size_t j=0;j<m.rowsize();++j) 
	      DoLU_LDivEq(CVIt<T1,Unit,Conj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(m.col(j).begin()),d.size());
	  else
	    for(size_t j=0;j<m.rowsize();++j) 
	      DoLU_LDivEq(CVIt<T1,Unit,NonConj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(m.col(j).begin()),d.size());
	else
	  if (d.diag().isconj())
	    for(size_t j=0;j<m.rowsize();++j) 
	      DoLU_LDivEq(CVIt<T1,Step,Conj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(m.col(j).begin()),d.size());
	  else
	    for(size_t j=0;j<m.rowsize();++j) 
	      DoLU_LDivEq(CVIt<T1,Step,NonConj>(d.diag().begin()),
		  VIt<T2,Unit,NonConj>(m.col(j).begin()),d.size());
      else
	if (d.diag().step()==1)
	  if (d.diag().isconj())
	    for(size_t j=0;j<m.rowsize();++j) 
	      DoLU_LDivEq(CVIt<T1,Unit,Conj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(m.col(j).begin()),d.size());
	  else
	    for(size_t j=0;j<m.rowsize();++j) 
	      DoLU_LDivEq(CVIt<T1,Unit,NonConj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(m.col(j).begin()),d.size());
	else
	  if (d.diag().isconj())
	    for(size_t j=0;j<m.rowsize();++j) 
	      DoLU_LDivEq(CVIt<T1,Step,Conj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(m.col(j).begin()),d.size());
	  else
	    for(size_t j=0;j<m.rowsize();++j) 
	      DoLU_LDivEq(CVIt<T1,Step,NonConj>(d.diag().begin()),
		  VIt<T2,Step,NonConj>(m.col(j).begin()),d.size());
  }

  template <class T> T DiagLUDiv<T>::Det() const 
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

  template <class T> DiagMatrix<T> DiagLUDiv<T>::DInverse() const 
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

  template <class T> DiagMatrix<T> DiagLUDiv<T>::DInverseATA() const
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

#define InstFile "TMV_DiagLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


