//---------------------------------------------------------------------------
//
// This file defines the TMV SymDivider class.
//
// It adds Inverse(SymMatrix) capability to the regular Divier class.
//

#ifndef TMV_SymDivider_H
#define TMV_SymDivider_H

#include "TMV_Divider.h"

namespace tmv {

  template <class T> class SymDivider : public Divider<T>
  {

    public :

      SymDivider() {}
      virtual ~SymDivider() {}

      using Divider<T>::Inverse;
      virtual void Inverse(const SymMatrixView<T>& sinv) const =0;
  };

  template <class T> inline string Type(const SymDivider<T>& d)
  { return string("SymDivider<")+tmv::Type(T())+">{"+d.Type()+"}"; }

} // namespace tmv

#endif
