#ifndef TMV_SYMLUDiv_A_H
#define TMV_SYMLUDiv_A_H

namespace tmv {

  template <bool herm, class T> void SymInvert_2x2(
      T& a, T& b, T& c, T* dd=0);
  template <bool herm, class T, class T1> void PseudoDiag_LDivEq(
      const GenVector<T1>& D, const GenVector<T1>& xD,
      const MatrixView<T>& m);
  template <bool herm, class T, class T1> void PseudoDiag_LMultEq(
      const GenVector<T1>& D, const GenVector<T1>& xD,
      const MatrixView<T>& m);

}

#endif
