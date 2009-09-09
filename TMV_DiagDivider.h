//---------------------------------------------------------------------------
//


#ifndef TMV_DiagDivider_H
#define TMV_DiagDivider_H

#include "TMV_Divider.h"
#include "TMV_DiagMatrix.h"

namespace tmv {

#define RT RealType(T)
#define CT ComplexType(T)

  template <class T> class DiagDivider : virtual public Divider<T>
  {

    public :

      DiagDivider() {}
      virtual ~DiagDivider() {}

      virtual DiagMatrix<T> DInverse() const =0;
      virtual DiagMatrix<T> DInverseATA() const =0;

      using Divider<T>::LDivEq;
      using Divider<T>::RDivEq;
      using Divider<T>::LDiv;
      using Divider<T>::RDiv;

#define DefDivEq(T1) \
      inline void DivEq(const DiagMatrixView<T1>& m) const\
      { this->LDivEq(ColVectorViewOf(m.diag())); }

      DefDivEq(RT);
      DefDivEq(CT);

#undef DefDivEq

#define DefDiv(T1,T2) \
      inline void Div(const GenDiagMatrix<T1>& b, \
	  const DiagMatrixView<T2>& x) const \
      { this->LDiv(ColVectorViewOf(b.diag()),ColVectorViewOf(x.diag())); }

      DefDiv(RT,T);
      DefDiv(CT,CT);

#undef DefDiv

  };

#undef RT
#undef CT

  template <class T> class MetaDiagDivider : 
    public DiagDivider<T>,
    private MetaDivider<T>
  {
    public :
      MetaDiagDivider(bool _trans, bool _conj, const BaseMatrix<T>* _m) :
	MetaDivider<T>(_trans,_conj,_m) {}

      MetaDiagDivider(const MetaDiagDivider<T>& rhs) :
	MetaDivider<T>(rhs) {}

      ~MetaDiagDivider() {}

      inline Divider<T>* MetaCopy() const 
      { return new MetaDiagDivider<T>(*this); }
      inline bool IsMeta() const { return true; }
      inline void MetaDivideUsing(DivType dt) const 
      { MetaDivider<T>::MetaDivideUsing(dt); }

      using MetaDivider<T>::istrans;
      using MetaDivider<T>::isconj;
      using MetaDivider<T>::GetM;

      inline DiagMatrix<T> DInverse() const
      {
	const DiagMatrix<T>* dm =
	  dynamic_cast<const DiagMatrix<T>*>(GetM());
	TMVAssert(dm);
	if (istrans()) {
	  if (isconj()) 
	    return DiagMatrix<T>(Adjoint(dm->DInverse()));
	  else
	    return DiagMatrix<T>(Transpose(dm->DInverse()));
	} else {
	  if (isconj())
	    return DiagMatrix<T>(Conjugate(dm->DInverse()));
	  else 
	    return dm->DInverse();
	}
      }

      inline DiagMatrix<T> DInverseATA() const
      {
	const DiagMatrix<T>* dm =
	  dynamic_cast<const DiagMatrix<T>*>(GetM());
	TMVAssert(dm);
	if (istrans()) {
	  if (isconj()) 
	    return DiagMatrix<T>(Adjoint(dm->DInverseATA()));
	  else
	    return DiagMatrix<T>(Transpose(dm->DInverseATA()));
	} else {
	  if (isconj())
	    return DiagMatrix<T>(Conjugate(dm->DInverseATA()));
	  else 
	    return dm->DInverseATA();
	}
      }

      inline Matrix<T,ColMajor> Inverse() const 
      { return Matrix<T,ColMajor>(DInverse()); }

      using MetaDivider<T>::DoLDivEq;
      using MetaDivider<T>::DoRDivEq;
      using MetaDivider<T>::DoLDiv;
      using MetaDivider<T>::DoRDiv;

#include "TMV_AuxAllDiv.h"

      inline T Det() const
      { return MetaDivider<T>::Det(); }

      inline Matrix<T,ColMajor> InverseATA() const
      { return MetaDivider<T>::InverseATA(); }

      inline RealType(T) Norm2() const
      { return GetM()->Norm2(); }

      inline bool Singular() const 
      { return GetM()->Singular(); }

      inline std::string Type() const
      {
	std::string s = std::string("MetaDiagDivider<") + tmv::Type(T());
	if (istrans()) s += ",trans";
	if (isconj()) s += ",conj";
	s += ">";
	return s;
      }

  };

  template <class T> MetaDiagDivider<T>* MakeMetaDiagDivider(
      bool trans, bool conj, const BaseMatrix<T>* m)
  {
    const MetaDiagDivider<T>* md1 = 
      dynamic_cast<const MetaDiagDivider<T>*>(m->GetDiv());
    if (md1) return new MetaDiagDivider<T>(
	trans ? !md1->istrans() : md1->istrans(),
	conj ? !md1->isconj() : md1->isconj(),
	md1->GetM());
    else return 0;
  }

} // namespace tmv

#endif
