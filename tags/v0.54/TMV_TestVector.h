#ifndef TMV_TESTVECTOR_H
#define TMV_TESTVECTOR_H

#include "TMV_Test.h"
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"
using tmv::Vector;

template <class T> void TestVectorReal();
template <class T> void TestVectorComplex();
template <class V1, class T> void DoTestVa(
    const V1& a, const Vector<T>& v, string label);
template <class V1, class T> void DoTestV(
    const V1& a, const Vector<T>& v, string label);
template <class V1, class T, class T2> void DoTestVXa(
    const V1& a, const Vector<T>& v, T2 x, string label);
template <class V1, class T, class T2> void DoTestVX(
    const V1& a, const Vector<T>& v, T2 x, string label);
template <class V1, class V2, class T, class T2> void DoTestVVa(
    const V1& a, const V2& b, const Vector<T>& v1, const Vector<T2>& v2,
    string label);
template <class V1, class V2, class T, class T2> void DoTestVV(
    const V1& a, const V2& b, const Vector<T>& v1, const Vector<T2>& v2,
    string label);
template <class T> void TestVectorArith();

#endif
