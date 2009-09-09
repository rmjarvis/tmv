#ifndef TMV_TESTVECTOR_H
#define TMV_TESTVECTOR_H

#include "TMV_Test.h"
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"
using tmv::Vector;

template <class T, class V> inline void TestV(
    const V& a, std::string label)
{
  tmv::Vector<T> v = a;
  Assert(Norm(a-v) <= EPS*Norm(v),label+" a != v");

  Assert(std::abs(Norm1(a)-Norm1(v)) <= EPS*std::abs(Norm1(v)),label+" Norm1");
  Assert(std::abs(Norm2(a)-Norm2(v)) <= EPS*std::abs(Norm2(v)),label+" Norm2");
  Assert(std::abs(NormInf(a)-NormInf(v)) <= EPS*std::abs(NormInf(v)),
      label+" NormInf");
  Assert(std::abs(NormSq(a)-NormSq(v)) <= EPS*std::abs(NormSq(v)),
      label+" NormSq");
}

template <class T, class V, class T2> inline void TestVX(
    const V& a, T2 x, std::string label)
{
  tmv::Vector<T> v = a;
  Assert(Norm(a-v) <= EPS*Norm(v),label+" a != v");

  Assert(Norm((x*a)-(x*v)) <= EPS*Norm(v)*std::abs(x),label+" x*a");
  Assert(Norm((a*x)-(x*v)) <= EPS*Norm(v)*std::abs(x),label+" a*x");
  Assert(Norm((a/x)-(v/x)) <= EPS*Norm(v)*std::abs(x),label+" v/x");
}

template <class T, class V, class T2> inline void TestVX2(
    const V& a, T2 x, std::string label)
{
  tmv::Vector<T> v = a;
  Assert(Norm(a-v) <= EPS*Norm(v),label+" a != v");
  double normv = std::abs(x)*Norm(v);
  tmv::Vector<T> v0 = v;

  a *= x;
  v = tmv::Vector<T>(v*x);
  Assert(Norm(a-v) <= EPS*normv,label+" a *= x");
  Assert(Norm((a*=x)-(v*=x)) <= EPS*normv,label+" a *= x (2)");
  a = v = v0;

  a /= x;
  v = tmv::Vector<T>(v/x);
  Assert(Norm(a-v) <= EPS*normv,label+" a /= x");
  a = v = v0;

#ifdef ALIASOK
  a = a*x;
  v = tmv::Vector<T>(v*x);
  Assert(Norm(a-v) <= EPS*normv,label+" a = a*x");
  a = v = v0;

  a = x*a;
  v = tmv::Vector<T>(v*x);
  Assert(Norm(a-v) <= EPS*normv,label+" a = x*a");
  a = v = v0;

  a = a/x;
  v = tmv::Vector<T>(v/x);
  Assert(Norm(a-v) <= EPS*normv,label+" a = a/x");
  a = v0;
#endif
}

template <class T, class T2, class V1, class V2> inline void TestVV(
    const V1& a, const V2& b, std::string label)
{
  tmv::Vector<T> v1 = a;
  tmv::Vector<T2> v2 = b;
  Assert(Norm(a-v1) <= EPS*Norm(v1),label+" a != v1");
  Assert(Norm(b-v2) <= EPS*Norm(v2),label+" b != v2");

  Assert(Norm((a+v2)-(v1+v2)) <= EPS*Norm(v1+v2),label+" a+v");
  Assert(Norm((v1+b)-(v1+v2)) <= EPS*Norm(v1+v2),label+" v+b");
  Assert(Norm((a-v2)-(v1-v2)) <= EPS*Norm(v1+v2),label+" a-v");
  Assert(Norm((v1-b)-(v1-v2)) <= EPS*Norm(v1+v2),label+" v-b");
  Assert(std::abs((a*v2)-(v1*v2)) <= EPS*std::abs(v1*v2),label+" a*v");
  Assert(std::abs((v1*b)-(v1*v2)) <= EPS*std::abs(v1*v2),label+" v*b");
  Assert(Norm((a+b)-(v1+v2)) <= EPS*Norm(v1+v2),label+" a+b");
  Assert(Norm((a-b)-(v1-v2)) <= EPS*Norm(v1+v2),label+" a-b");
  Assert(std::abs((a*b)-(v1*v2)) <= EPS*std::abs(v1*v2),label+" a*b");
}

template <class T, class T2, class V1, class V2> inline void TestVV2(
    const V1& a, const V2& b, std::string label)
{
  tmv::Vector<T> v1 = a;
  tmv::Vector<T2> v2 = b;
  Assert(Norm(a-v1) <= EPS*Norm(v1),label+" a != v1");
  Assert(Norm(b-v2) <= EPS*Norm(v2),label+" b != v2");

  double normv = Norm(v1)+Norm(v2);

  tmv::Vector<T> v3 = v1;
  tmv::Vector<T> v4 = v1;
  v3 += b;
  v4 = v1+v2;
  Assert(Norm(v3-v4) <= EPS*normv,label+" v += b");
  v3 = v4 = v1;
  v3 -= b;
  v4 = v1-v2;
  Assert(Norm(v3-v4) <= EPS*normv,label+" v -= b");
  v3 = v4 = v1;
  v3 = v3+b;
  v4 = v1+v2;
  Assert(Norm(v3-v4) <= EPS*normv,label+" v = v+b");
  v3 = v4 = v1;
  v3 = v3-b;
  v4 = v1-v2;
  Assert(Norm(v3-v4) <= EPS*normv,label+" v = v-b");
  v3 = v4 = v1;
  v3 = b+v3;
  v4 = v1+v2;
  Assert(Norm(v3-v4) <= EPS*normv,label+" v = b+v");
  v3 = v4 = v1;
  v3 = b-v3;
  v4 = v2-v1;
  Assert(Norm(v3-v4) <= EPS*normv,label+" v = b-v");
  v3 = v4 = v1;

  a += v2;
  v4 = v1+v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a += v");
  a = v4 = v1;
  a -= v2;
  v4 = v1-v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a -= v");
  a = v4 = v1;
#ifdef ALIASOK
  a = a+v2;
  v4 = v1+v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a = a+v");
  a = v4 = v1;
  a = a-v2;
  v4 = v1-v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a = a-v");
  a = v4 = v1;
  a = v2+a;
  v4 = v1+v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a = v+a");
  a = v4 = v1;
  a = v2-a;
  v4 = v2-v1;
  Assert(Norm(a-v4) <= EPS*normv,label+" a = v-a");
  a = v4 = v1;
#endif

  a += b;
  v4 = v1+v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a += b");
  a = v4 = v1;
  a -= b;
  v4 = v1-v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a -= b");
  a = v4 = v1;
#ifdef ALIASOK
  a = a+b;
  v4 = v1+v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a = a+b");
  a = v4 = v1;
  a = a-b;
  v4 = v1-v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a = a-b");
  a = v4 = v1;
  a = b+a;
  v4 = v1+v2;
  Assert(Norm(a-v4) <= EPS*normv,label+" a = b+a");
  a = v4 = v1;
  a = b-a;
  v4 = v2-v1;
  Assert(Norm(a-v4) <= EPS*normv,label+" a = b-a");
  a = v4 = v1;
#endif
}

template <class T, class V, class CV> inline void TestVectorArith1(
    const V& a, const CV& ca, std::string label)
{
  TestV<T>(a,label+" R");
  TestV<std::complex<T> >(ca,label+" C");

  T x = 12;
  std::complex<T> z(9.-2);

  TestVX<T>(a,x,label+" R,R");
  TestVX2<T>(a,x,label+" R,R");
  TestVX<T>(a,z,label+" R,C");
  TestVX<std::complex<T> >(ca,x,label+" C,R");
  TestVX2<std::complex<T> >(ca,x,label+" C,R");
  TestVX<std::complex<T> >(ca,z,label+" C,C");
  TestVX2<std::complex<T> >(ca,z,label+" C,C");
#ifdef ALIASOK
  TestVV<T,T>(a,a,label+" R,R");
  TestVV2<T,T>(a,a,label+" R,R");
  TestVV<std::complex<T>,std::complex<T> >(ca,ca,label+" C,C");
  TestVV2<std::complex<T>,std::complex<T> >(ca,ca,label+" C,C");
#endif

}

template <class T, class V1, class CV1, class V2, class CV2> 
inline void TestVectorArith2(
    const V1& a, const CV1& ca, const V2& b, const CV2& cb, std::string label)
{
  TestVV<T,T>(a,b,label+" R,R");
  TestVV2<T,T>(a,b,label+" R,R");
  TestVV<T,std::complex<T> >(a,cb,label+" R,C");
  TestVV<std::complex<T>,T>(ca,b,label+" C,R");
  TestVV2<std::complex<T>,T>(ca,b,label+" C,R");
  TestVV<std::complex<T>,std::complex<T> >(ca,cb,label+" C,C");
  TestVV2<std::complex<T>,std::complex<T> >(ca,cb,label+" C,C");
}

#endif
