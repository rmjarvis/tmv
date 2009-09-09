//---------------------------------------------------------------------------
//
// This file defines the TMV Divider class.
//
// There are currently 4 algorithms for doing division (and Inverse and Det)
//
// LU Decomposition
// Singular Value Decomposition
// QR Decomposition with or without Permutation
//
// To tell a Matrix to use a particular algorithm, use the command:
// m.DivideUsing(ALG)
// where ALG is 'LU', 'SV', 'QR' or 'QRP' for the algorithms above.
//
// The default algorithm is LU for square matrices or SV for non-square.
//
// By default, the appropriate Divider class is created the first
// time it is needed (eg. when the statement v = b/m is called).
// However, you can also setup the Divider class beforehand manually
// by calling m.SetupDiv().
//
// You can also query whether the Divider class is already set up.
// This will only be true, if it was previously set up, _and_ the 
// Matrix hasn't been modified since then.
//
// If you want access to the various Divider functions directly,
// They can be accessed by:
//
// m.LUD()
// m.SVD()
// m.QRD()
// m.QRPD()
//
// The one of these that is probably most useful to access is SVD(),
// since it is generally a good idea to look for small singular values
// and zero them out before using SVD for division.
//
// To set to zero all singular value which are less than thresh * 
// the largest singular value use:
//
// m.SVD()->SetThresh(thresh);
//
// To use only the largest nsv singular values use:
//
// m.SVD()->SetTop(nsv);
//
// Also, the singular value decomposition can be used for principal
// component analysis of a Matrix.  The principal component vectors
// are the rows of V.  You can access the decomposition using:
//
// m.SVD()->GetU();
// m.SVD()->GetS(); // The diagonal vector of the (diagonal) matrix W
// m.SVD()->GetV();
//
//


#ifndef TMV_Divider_H
#define TMV_Divider_H

#include "TMV_Matrix.h"

namespace tmv {

#define RT RealType(T)
#define CT ComplexType(T)

  // This is defined in TMV_BaseMatrix.h
  //enum DivType { XXX, LU, QR, QRP, SV, CH };

  template <class T> class Divider 
  {

    public :

      Divider() {}
      virtual ~Divider() {}

      virtual inline bool IsSV() const 
      { return false; }
      virtual inline Divider<T>* MetaCopy() const { return 0; }
      virtual inline const BaseMatrix<T>* MetaM() const { return 0; }

      virtual T Det() const =0;
      virtual Matrix<T,ColMajor> Inverse() const =0;
      virtual Matrix<T,ColMajor> InverseATA() const =0;
      virtual bool Singular() const =0;
      virtual RealType(T) Norm2() const 
      { TMVAssert(false); return RealType(T)(0); }

#define DefDivEq(T) \
      virtual void LDivEq(const MatrixView<T>&) const =0; \
      virtual void RDivEq(const MatrixView<T>&) const =0; 
      DefDivEq(RT);
      DefDivEq(CT);
#undef DefDivEq

#define DefDiv(T1,T2) \
      virtual void LDiv(const GenMatrix<T1>& b, \
	  const MatrixView<T2>& x) const =0; \
      virtual void RDiv(const GenMatrix<T1>& b, \
	  const MatrixView<T2>& x) const =0; 
      DefDiv(RT,T);
      DefDiv(CT,CT);
#undef DefDiv

      virtual std::string Type() const = 0;
  };

#undef RT
#undef CT

  template <class T> class BaseLUDiv : virtual public Divider<T> 
  {
    public :
      BaseLUDiv() {}
      virtual ~BaseLUDiv() {}

      virtual Matrix<T,ColMajor> LU_GetL() const = 0;
      virtual Matrix<T,ColMajor> LU_GetU() const = 0;
      virtual Permutation LU_GetP() const = 0;
      virtual Permutation LU_GetQ() const = 0;
      // Normally, decomp is A = PLU, but allow for permutation
      // on other side, or both sides (eg. for Hermitian matrices)
      // So A = PLUQ
  };

  template <class T> class BaseQRDiv : virtual public Divider<T> 
  {
    public :
      BaseQRDiv() {}
      virtual ~BaseQRDiv() {}

      virtual Matrix<T,ColMajor> QR_GetQ() const = 0;
      virtual Matrix<T,ColMajor> QR_GetR() const = 0;
      virtual Permutation QR_GetP() const = 0;
      virtual bool QR_IsTrans() const = 0;
      // Allow for R to be permuted, so general decomposition is
      // A = QRP
  };

  template <class T> class BaseSVDiv : virtual public Divider<T> 
  {
    public :
      BaseSVDiv() {}
      virtual ~BaseSVDiv() {}

      inline bool IsSV() const { return true; }
      virtual void Thresh(RealType(T) toler, ostream* debugout=0) const =0;
      virtual void Top(size_t neigen, ostream* debugout=0) const =0;
      virtual size_t GetKMax() const = 0;

      virtual Matrix<T,ColMajor> SV_GetU() const = 0;
      virtual Vector<RealType(T)> SV_GetS() const = 0;
      virtual Matrix<T,ColMajor> SV_GetV() const = 0;
  };

  template <class T> class BaseCHDiv : virtual public Divider<T> 
  {
    public :
      BaseCHDiv() {}
      virtual ~BaseCHDiv() {}

      virtual Matrix<T,ColMajor> CH_GetL() const = 0;
  };

  template <class T> class MetaDivider : virtual public Divider<T>
  {
    public :
      MetaDivider(bool _trans, bool _conj, const BaseMatrix<T>* _m) :
	trans(_trans), conj(_conj), itsm(_m) {}

      MetaDivider(const MetaDivider<T>& rhs) :
	trans(rhs.trans), conj(rhs.conj), itsm(rhs.itsm) { }

      ~MetaDivider() { }

      inline bool IsSV() const 
      {
	return (itsm->GetDivType() == tmv::SV || 
	    (itsm->GetDiv() && itsm->GetDiv()->IsSV())); 
      }
      inline Divider<T>* MetaCopy() const 
      { return new MetaDivider<T>(*this); }
      inline const BaseMatrix<T>* MetaM() const { return GetM(); }

      bool istrans() const { return trans; }
      bool isconj() const { return conj; }
      const BaseMatrix<T>* GetM() const { return itsm; }

      inline T Det() const
      { 
	if (conj) return tmv::CONJ(itsm->Det());
	else return itsm->Det();
      }

      inline Matrix<T,ColMajor> Inverse() const
      {
	if (conj) 
	  if (trans)
	    return Matrix<T,ColMajor>(itsm->Inverse().QuickAdjoint());
	  else
	    return Matrix<T,ColMajor>(itsm->Inverse().QuickConjugate());
	else
	  if (trans)
	    return Matrix<T,ColMajor>(itsm->Inverse().QuickTranspose());
	  else 
	    return itsm->Inverse();
      }

      inline Matrix<T,ColMajor> InverseATA() const
      {
	// AtA is Hermitian, so conj, trans are same thing
	// != is xor
	if (conj != trans) 
	  return Matrix<T,ColMajor>(itsm->InverseATA().QuickTranspose());
	else
	  return itsm->InverseATA();
      }

      inline RealType(T) Norm2() const
      { return itsm->Norm2(); }

      inline bool Singular() const 
      { return itsm->Singular(); }

      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const
      {
	if (trans) {
	  // m = MT^-1 m --> mT = mT M^-1
	  if (conj) itsm->RDivEq(m.QuickAdjoint());
	  else itsm->RDivEq(m.QuickTranspose());
	} else {
	  if (conj) itsm->LDivEq(m.QuickConjugate());
	  else itsm->LDivEq(m);
	}
      }

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const
      {
	if (trans) {
	  // m = m MT^-1 --> mT = M^-1 mT
	  if (conj) itsm->LDivEq(m.QuickAdjoint());
	  else itsm->LDivEq(m.QuickTranspose());
	} else {
	  if (conj) itsm->RDivEq(m.QuickConjugate());
	  else itsm->RDivEq(m);
	}
      }

      template <class T1, class T2> inline void DoLDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const
      {
	if (trans) {
	  if (conj) itsm->RDiv(m.QuickAdjoint(),x.QuickAdjoint());
	  else itsm->RDiv(m.QuickTranspose(),x.QuickTranspose());
	} else {
	  if (conj) itsm->LDiv(m.QuickConjugate(),x.QuickConjugate());
	  else itsm->LDiv(m,x);
	}
      }

      template <class T1, class T2> inline void DoRDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const
      {
	if (trans) {
	  if (conj) itsm->LDiv(m.QuickAdjoint(),x.QuickAdjoint());
	  else itsm->LDiv(m.QuickTranspose(),x.QuickTranspose());
	} else {
	  if (conj) itsm->RDiv(m.QuickConjugate(),x.QuickConjugate());
	  else itsm->RDiv(m,x);
	}
      }

#include "TMV_AuxAllDiv.h"

      inline std::string Type() const 
      { 
	std::string s = std::string("MetaDivider<") + tmv::Type(T());
	if (trans) s += ",trans";
	if (conj) s += ",conj";
	s += ">";
	return s;
      }

    private :

      bool trans;
      bool conj;
      const BaseMatrix<T>* itsm;
  };

  template <class T> MetaDivider<T>* MakeMetaDivider(bool trans, bool conj,
      const BaseMatrix<T>* m)
  {
    const MetaDivider<T>* md1 = 
      dynamic_cast<const MetaDivider<T>*>(m->GetDiv());
    if (md1) return new MetaDivider<T>(
	trans ? !md1->istrans() : md1->istrans(),
	conj ? !md1->isconj() : md1->isconj(), md1->GetM());
    else return 0;
  }


  inline std::string Text(DivType dt)
  {
    switch (dt) {
      case tmv::XXX : return "DT: XXX";
      case tmv::LU : return "DT: LU";
      case tmv::SV : return "DT: SV";
      case tmv::QR : return "DT: QR";
      case tmv::QRP : return "DT: QRP";
      case tmv::CH : return "DT: CH";
      default : return "DT: unknown";
    }
  }

  template <class T> inline std::string Type(const Divider<T>& d)
  { return std::string("Divider<")+tmv::Type(T())+">{"+d.Type()+"}"; }

} // namespace tmv

#endif
