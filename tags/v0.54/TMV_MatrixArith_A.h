#ifndef TMV_MATRIXARITH_A_H
#define TMV_MATRIXARITH_A_H

namespace tmv {

  template <bool b0, bool cx, class T, class Ta, class Tx> 
    void UnitAMultMV1(
	const GenMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y);

} // namespace mv

#endif
