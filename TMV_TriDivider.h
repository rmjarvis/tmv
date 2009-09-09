//---------------------------------------------------------------------------
//


#ifndef TMV_TriDivider_H
#define TMV_TriDivider_H

#include "TMV_Divider.h"
#include "TMV_TriMatrix.h"

namespace tmv {

#define RT RealType(T)
#define CT ComplexType(T)

  template <class T> class UpperTriDivider : virtual public Divider<T>
  {

    public :

      UpperTriDivider() {}
      virtual ~UpperTriDivider() {}

      virtual UpperTriMatrix<T,NonUnitDiag,ColMajor> TInverse() const =0;

#define DefDivEq(T1) \
      virtual void LDivEq(const UpperTriMatrixView<T1>&) const =0; \
      virtual void RDivEq(const UpperTriMatrixView<T1>&) const =0; \

      DefDivEq(RT);
      DefDivEq(CT);

#undef DefDivEq

#define DefDiv(T1,T2) \
      virtual void LDiv(const GenUpperTriMatrix<T1>& b, \
	  const UpperTriMatrixView<T2>& x) const =0; \
      virtual void RDiv(const GenUpperTriMatrix<T1>& b, \
	  const UpperTriMatrixView<T2>& x) const =0; \

      DefDiv(RT,T);
      DefDiv(CT,CT);

#undef DefDiv

      using Divider<T>::LDivEq;
      using Divider<T>::RDivEq;
      using Divider<T>::LDiv;
      using Divider<T>::RDiv;

  };

  template <class T> class LowerTriDivider : virtual public Divider<T>
  {

    public :

      LowerTriDivider() {}
      virtual ~LowerTriDivider() {}

      virtual LowerTriMatrix<T,NonUnitDiag,ColMajor> TInverse() const =0;

#define DefDivEq(T1) \
      virtual void LDivEq(const LowerTriMatrixView<T1>&) const =0; \
      virtual void RDivEq(const LowerTriMatrixView<T1>&) const =0; \

      DefDivEq(RT);
      DefDivEq(CT);

#undef DefDivEq

#define DefDiv(T1,T2) \
      virtual void LDiv(const GenLowerTriMatrix<T1>& b, \
	  const LowerTriMatrixView<T2>& x) const =0; \
      virtual void RDiv(const GenLowerTriMatrix<T1>& b, \
	  const LowerTriMatrixView<T2>& x) const =0; \

      DefDiv(RT,T);
      DefDiv(CT,CT);

#undef DefDiv

      using Divider<T>::LDivEq;
      using Divider<T>::RDivEq;
      using Divider<T>::LDiv;
      using Divider<T>::RDiv;

  };

#undef RT
#undef CT

  template <class T> class MetaUpperTriDivider : 
    public UpperTriDivider<T>, 
    private MetaDivider<T>
  {
    public :
      MetaUpperTriDivider(bool _trans, bool _conj, const BaseMatrix<T>* _m) :
	MetaDivider<T>(_trans,_conj,_m) {}

      MetaUpperTriDivider(const MetaUpperTriDivider<T>& rhs) :
	MetaDivider<T>(rhs) {}

      ~MetaUpperTriDivider() {}

      using MetaDivider<T>::IsSV;
      using MetaDivider<T>::MetaM;

      inline Divider<T>* MetaCopy() const 
      { return new MetaUpperTriDivider<T>(*this); }

      using MetaDivider<T>::istrans;
      using MetaDivider<T>::isconj;
      using MetaDivider<T>::GetM;

      inline UpperTriMatrix<T,NonUnitDiag,ColMajor> TInverse() const
      {
	if (istrans()) {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) 
	    return UpperTriMatrix<T,NonUnitDiag,ColMajor>(
		tm->TInverse().QuickAdjoint());
	  else
	    return UpperTriMatrix<T,NonUnitDiag,ColMajor>(
		tm->TInverse().QuickTranspose());
	} else {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj())
	    return UpperTriMatrix<T,NonUnitDiag,ColMajor>(
		tm->TInverse().QuickConjugate());
	  else 
	    return tm->TInverse();
	}
      }

      inline Matrix<T,ColMajor> Inverse() const 
      { return Matrix<T,ColMajor>(TInverse()); }

      template <class T1> inline void DoLDivEq(
	  const UpperTriMatrixView<T1>& m) const
      {
	if (istrans()) {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  // m = MT^-1 m --> mT = mT M^-1
	  if (isconj()) tm->RDivEq(m.QuickAdjoint());
	  else tm->RDivEq(m.QuickTranspose());
	} else {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->LDivEq(m.QuickConjugate());
	  else tm->LDivEq(m);
	}
      }

      template <class T1> inline void DoRDivEq(
	  const UpperTriMatrixView<T1>& m) const
      {
	if (istrans()) {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  // m = m MT^-1 --> mT = M^-1 mT
	  if (isconj()) tm->LDivEq(m.QuickAdjoint());
	  else tm->LDivEq(m.QuickTranspose());
	} else {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->RDivEq(m.QuickConjugate());
	  else tm->RDivEq(m);
	}
      }

      template <class T1, class T2> inline void DoLDiv(
	  const GenUpperTriMatrix<T1>& m,
	  const UpperTriMatrixView<T2>& x) const
      {
	if (istrans()) {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->RDiv(m.QuickAdjoint(),x.QuickAdjoint());
	  else tm->RDiv(m.QuickTranspose(),x.QuickTranspose());
	} else {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->LDiv(m.QuickConjugate(),x.QuickConjugate());
	  else tm->LDiv(m,x);
	}
      }

      template <class T1, class T2> inline void DoRDiv(
	  const GenUpperTriMatrix<T1>& m,
	  const UpperTriMatrixView<T2>& x) const
      {
	if (istrans()) {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->LDiv(m.QuickAdjoint(),x.QuickAdjoint());
	  else tm->LDiv(m.QuickTranspose(),x.QuickTranspose());
	} else {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->RDiv(m.QuickConjugate(),x.QuickConjugate());
	  else tm->RDiv(m,x);
	}
      }

      using MetaDivider<T>::DoLDivEq;
      using MetaDivider<T>::DoRDivEq;
      using MetaDivider<T>::DoLDiv;
      using MetaDivider<T>::DoRDiv;

#include "TMV_AuxAllDiv.h"
#include "TMV_AuxUpperTriAllDiv.h"

      using MetaDivider<T>::Det;
      using MetaDivider<T>::InverseATA;
      using MetaDivider<T>::Norm2;
      using MetaDivider<T>::Singular;

      inline std::string Type() const
      {
	std::string s = std::string("MetaUpperTriDivider<") + tmv::Type(T());
	if (istrans()) s += ",trans";
	if (isconj()) s += ",conj";
	s += ">";
	return s;
      }
  }; // MetaUpperTriDivider

  template <class T> class MetaLowerTriDivider : 
    public LowerTriDivider<T>, 
    private MetaDivider<T>
  {
    public :
      MetaLowerTriDivider(bool _trans, bool _conj, const BaseMatrix<T>* _m) :
	MetaDivider<T>(_trans,_conj,_m) {}

      MetaLowerTriDivider(const MetaLowerTriDivider<T>& rhs) :
	MetaDivider<T>(rhs) {}

      ~MetaLowerTriDivider() {}

      using MetaDivider<T>::IsSV;
      using MetaDivider<T>::MetaM;

      inline Divider<T>* MetaCopy() const 
      { return new MetaLowerTriDivider<T>(*this); }

      using MetaDivider<T>::istrans;
      using MetaDivider<T>::isconj;
      using MetaDivider<T>::GetM;

      inline LowerTriMatrix<T,NonUnitDiag,ColMajor> TInverse() const
      {
	if (istrans()) {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) 
	    return LowerTriMatrix<T,NonUnitDiag,ColMajor>(
		tm->TInverse().QuickAdjoint());
	  else
	    return LowerTriMatrix<T,NonUnitDiag,ColMajor>(
		tm->TInverse().QuickTranspose());
	} else {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj())
	    return LowerTriMatrix<T,NonUnitDiag,ColMajor>(
		tm->TInverse().QuickConjugate());
	  else 
	    return tm->TInverse();
	}
      }

      inline Matrix<T,ColMajor> Inverse() const 
      { return Matrix<T,ColMajor>(TInverse()); }

      template <class T1> inline void DoLDivEq(
	  const LowerTriMatrixView<T1>& m) const
      {
	if (istrans()) {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  // m = MT^-1 m --> mT = mT M^-1
	  if (isconj()) tm->RDivEq(m.QuickAdjoint());
	  else tm->RDivEq(m.QuickTranspose());
	} else {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->LDivEq(m.QuickConjugate());
	  else tm->LDivEq(m);
	}
      }

      template <class T1> inline void DoRDivEq(
	  const LowerTriMatrixView<T1>& m) const
      {
	if (istrans()) {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  // m = m MT^-1 --> mT = M^-1 mT
	  if (isconj()) tm->LDivEq(m.QuickAdjoint());
	  else tm->LDivEq(m.QuickTranspose());
	} else {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->RDivEq(m.QuickConjugate());
	  else tm->RDivEq(m);
	}
      }

      template <class T1, class T2> inline void DoLDiv(
	  const GenLowerTriMatrix<T1>& m,
	  const LowerTriMatrixView<T2>& x) const
      {
	if (istrans()) {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->RDiv(m.QuickAdjoint(),x.QuickAdjoint());
	  else tm->RDiv(m.QuickTranspose(),x.QuickTranspose());
	} else {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->LDiv(m.QuickConjugate(),x.QuickConjugate());
	  else tm->LDiv(m,x);
	}
      }

      template <class T1, class T2> inline void DoRDiv(
	  const GenLowerTriMatrix<T1>& m,
	  const LowerTriMatrixView<T2>& x) const
      {
	if (istrans()) {
	  const GenUpperTriMatrix<T>* tm =
	    dynamic_cast<const GenUpperTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->LDiv(m.QuickAdjoint(),x.QuickAdjoint());
	  else tm->LDiv(m.QuickTranspose(),x.QuickTranspose());
	} else {
	  const GenLowerTriMatrix<T>* tm =
	    dynamic_cast<const GenLowerTriMatrix<T>*>(GetM());
	  TMVAssert(tm);
	  if (isconj()) tm->RDiv(m.QuickConjugate(),x.QuickConjugate());
	  else tm->RDiv(m,x);
	}
      }

      using MetaDivider<T>::DoLDivEq;
      using MetaDivider<T>::DoRDivEq;
      using MetaDivider<T>::DoLDiv;
      using MetaDivider<T>::DoRDiv;

#include "TMV_AuxAllDiv.h"
#include "TMV_AuxLowerTriAllDiv.h"

      using MetaDivider<T>::Det;
      using MetaDivider<T>::InverseATA;
      using MetaDivider<T>::Norm2;
      using MetaDivider<T>::Singular;

      inline std::string Type() const
      {
	std::string s = std::string("MetaLowerTriDivider<") + tmv::Type(T());
	if (istrans()) s += ",trans";
	if (isconj()) s += ",conj";
	s += ">";
	return s;
      }
  }; // MetaLowerTriDivider

  template <class T> MetaUpperTriDivider<T>* MakeMetaUpperTriDivider(
      bool trans, bool conj, const BaseMatrix<T>* m)
  {
    if (trans) {
      const MetaLowerTriDivider<T>* md1 = 
	dynamic_cast<const MetaLowerTriDivider<T>*>(m->GetDiv());
      if (md1) {
	return new MetaUpperTriDivider<T>(!md1->istrans(),
	  conj ? !md1->isconj() : md1->isconj(), md1->GetM());
      } else {
	TMVAssert(!dynamic_cast<const MetaUpperTriDivider<T>*>(m->GetDiv()));
	return 0;
      }
    } else {
      const MetaUpperTriDivider<T>* md1 = 
	dynamic_cast<const MetaUpperTriDivider<T>*>(m->GetDiv());
      if (md1) {
	return new MetaUpperTriDivider<T>(md1->istrans(),
	  conj ? !md1->isconj() : md1->isconj(), md1->GetM());
      } else {
	TMVAssert(!dynamic_cast<const MetaLowerTriDivider<T>*>(m->GetDiv()));
	return 0;
      }
    }
  }

  template <class T> MetaLowerTriDivider<T>* MakeMetaLowerTriDivider(
      bool trans, bool conj, const BaseMatrix<T>* m)
  {
    if (trans) {
      const MetaUpperTriDivider<T>* md1 = 
	dynamic_cast<const MetaUpperTriDivider<T>*>(m->GetDiv());
      if (md1) {
	return new MetaLowerTriDivider<T>(!md1->istrans(),
	  conj ? !md1->isconj() : md1->isconj(), md1->GetM());
      } else {
	TMVAssert(!dynamic_cast<const MetaLowerTriDivider<T>*>(m->GetDiv()));
	return 0;
      }
    } else {
      const MetaLowerTriDivider<T>* md1 = 
	dynamic_cast<const MetaLowerTriDivider<T>*>(m->GetDiv());
      if (md1) {
	return new MetaLowerTriDivider<T>(md1->istrans(),
	  conj ? !md1->isconj() : md1->isconj(), md1->GetM());
      } else {
	TMVAssert(!dynamic_cast<const MetaUpperTriDivider<T>*>(m->GetDiv()));
	return 0;
      }
    }
  }

} // namespace tmv

#endif
