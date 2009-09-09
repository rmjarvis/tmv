#ifndef TMV_SymHouseholder_H
#define TMV_SymHouseholder_H

#include "TMV_Sym.h"
#include "TMV_Householder.h"

namespace tmv {

  template <class T1, class T2> void Householder_LRMult(
      const GenVector<T1>& v, T1 beta, const SymMatrixView<T2>& M);
  // The input vector, v, is taken to be the vector for a  
  // Householder matrix, H.  This routine takes M <- H M Ht
  // if M is Hermitian or H M HT if M is symmetric.

} // namespace tmv

#endif
