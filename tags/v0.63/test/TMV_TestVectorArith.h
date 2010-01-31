#ifndef TMV_TESTVECTOR_H
#define TMV_TESTVECTOR_H

#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "TMV_Test.h"
using tmv::Vector;

#ifdef NOVIEWS
#define CONST
#endif

#ifndef CONST
#define CONST const
#endif

#ifdef USETEMP
#define IFTEMP(x) x
#define IFTEMP1(x) x,
#else
#define IFTEMP(x)
#define IFTEMP1(x)
#endif

template <class T1, class T2> struct ProdType
{ typedef T1 Tprod; };

template <class T> struct ProdType<T,std::complex<T> >
{ typedef std::complex<T> Tprod; };

#define ProductType(T1,T2) typename ProdType<T1,T2>::Tprod

#ifdef NOMIX
#define VEC(T,v) (tmv::Vector<T>(v))
#define VEC2(T1,T2,v) (tmv::Vector<ProductType(T1,T2)>(v))
#else
#define VEC(T,v) (v)
#define VEC2(T1,T2,v) (v)
#endif

template <class T, class V> inline void TestV(
    CONST V& a, std::string label)
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
    CONST V& a, T2 x, std::string label)
{
  tmv::Vector<T> v = a;
  Assert(Norm(a-v) <= EPS*Norm(v),label+" a != v");

  Assert(Norm(VEC2(T,T2,x*a)-(x*v)) <= EPS*Norm(v)*std::abs(x),label+" x*a");
  Assert(Norm(VEC2(T,T2,a*x)-(x*v)) <= EPS*Norm(v)*std::abs(x),label+" a*x");
  Assert(Norm(VEC2(T,T2,a/x)-(v/x)) <= EPS*Norm(v)*std::abs(x),label+" a/x");
}

template <class T, class V, class T2> inline void TestVX2(
    CONST V& a, T2 x, std::string label)
{
  tmv::Vector<T> v = a;
  Assert(Norm(VEC(T,a)-v) <= EPS*Norm(v),label+" a != v");
  double normv = std::abs(x)*Norm(v);
  tmv::Vector<T> v0 = v;

  a *= x;
  v = tmv::Vector<T>(v*x);
  Assert(Norm(VEC(T,a)-v) <= EPS*normv,label+" a *= x");
  Assert(Norm(VEC(T,a*=x)-(v*=x)) <= EPS*normv,label+" a *= x (2)");
  a = v = v0;

  a /= x;
  v = tmv::Vector<T>(v/x);
  Assert(Norm(VEC(T,a)-v) <= EPS*normv,label+" a /= x");
  a = v = v0;

#ifdef ALIASOK
  a = a*x;
  v = tmv::Vector<T>(v*x);
  Assert(Norm(VEC(T,a)-v) <= EPS*normv,label+" a = a*x");
  a = v = v0;

  a = x*a;
  v = tmv::Vector<T>(v*x);
  Assert(Norm(VEC(T,a)-v) <= EPS*normv,label+" a = x*a");
  a = v = v0;

  a = a/x;
  v = tmv::Vector<T>(v/x);
  Assert(Norm(VEC(T,a)-v) <= EPS*normv,label+" a = a/x");
  a = v0;
#endif
}

template <class T, class T2, IFTEMP1(class V0) IFTEMP1(class CV0) class V1, class V2> inline void TestVV(
    IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) const V1& a, const V2& b, std::string label)
{
  tmv::Vector<T> v1 = a;
  tmv::Vector<T2> v2 = b;
  TMV_RealType(T) eps = EPS * (Norm(v1) + Norm(v2));
  TMV_RealType(T) eps2 = EPS * Norm(v1) * Norm(v2);
  Assert(Norm(VEC(T,a)-v1) <= EPS*Norm(v1),label+" a != v1");
  Assert(Norm(VEC(T2,b)-v2) <= EPS*Norm(v2),label+" b != v2");

  Assert(Norm(VEC2(T,T2,IFTEMP(temp=)a+b)-(v1+v2)) <= eps,label+" a+b");
  Assert(Norm(VEC2(T,T2,IFTEMP(temp=)a-b)-(v1-v2)) <= eps,label+" a-b");
  Assert(std::abs((a*b)-(v1*v2)) <= eps2,label+" a*b");
#ifndef NOMIX
  Assert(Norm((a+v2)-(v1+v2)) <= eps,label+" a+v");
  Assert(Norm((v1+b)-(v1+v2)) <= eps,label+" v+b");
  Assert(Norm((a-v2)-(v1-v2)) <= eps,label+" a-v");
  Assert(Norm((v1-b)-(v1-v2)) <= eps,label+" v-b");
  Assert(std::abs((a*v2)-(v1*v2)) <= eps2,label+" a*v");
  Assert(std::abs((v1*b)-(v1*v2)) <= eps2,label+" v*b");
#endif

  TMV_RealType(T) x(5);
  TMV_ComplexType(T) z(3,4);
  Assert(Norm(VEC2(T,T2,IFTEMP(temp=)a-x*b)-(v1-x*v2)) <= x*eps,label+" a-x*b");
  Assert(Norm(VEC2(T,T2,IFTEMP(temp=)a+x*b)-(v1+x*v2)) <= x*eps,label+" a+x*b");
  Assert(Norm(VEC2(T,T2,IFTEMP(temp=)x*a-b)-(x*v1-v2)) <= x*eps,label+" x*a-b");
  Assert(Norm(VEC2(T,T2,IFTEMP(temp=)x*a+b)-(x*v1+v2)) <= x*eps,label+" x*a+b");
  Assert(Norm(VEC2(T,T2,IFTEMP(temp=)x*a-x*b)-(x*v1-x*v2)) <= x*eps,label+" x*a-x*b");
  Assert(Norm(VEC2(T,T2,IFTEMP(temp=)x*a+x*b)-(x*v1+x*v2)) <= x*eps,label+" x*a+x*b");

  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)a-z*b)-(v1-z*v2)) <= x*eps,label+" a-z*b");
  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)a+z*b)-(v1+z*v2)) <= x*eps,label+" a+z*b");
  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*a-b)-(z*v1-v2)) <= x*eps,label+" z*a-b");
  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*a+b)-(z*v1+v2)) <= x*eps,label+" z*a+b");
  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*a-z*b)-(z*v1-z*v2)) <= x*eps,label+" z*a-z*b");
  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*a+z*b)-(z*v1+z*v2)) <= x*eps,label+" z*a+z*b");

  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*a-x*b)-(z*v1-x*v2)) <= x*eps,label+" z*a-x*b");
  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*a+x*b)-(z*v1+x*v2)) <= x*eps,label+" z*a+x*b");
  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)x*a-z*b)-(x*v1-z*v2)) <= x*eps,label+" x*a-z*b");
  Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)x*a+z*b)-(x*v1+z*v2)) <= x*eps,label+" x*a+z*b");

  Assert(std::abs(((x*a)*b)-(x*v1*v2)) <= x*eps2,label+" (x*a)*b");
  Assert(std::abs((a*(x*b))-(x*v1*v2)) <= x*eps2,label+" a*(x*b)");
  Assert(std::abs((x*(a*b))-(x*v1*v2)) <= x*eps2,label+" x*(a*b)");

  Assert(std::abs(((z*a)*b)-(z*v1*v2)) <= x*eps2,label+" (z*a)*b");
  Assert(std::abs((a*(z*b))-(z*v1*v2)) <= x*eps2,label+" a*(z*b)");
  Assert(std::abs((z*(a*b))-(z*v1*v2)) <= x*eps2,label+" z*(a*b)");

  Assert(std::abs(((x*a)*(x*b))-(x*x*v1*v2)) <= x*x*eps2,label+" (x*a)*(x*b)");
  Assert(std::abs(((z*a)*(x*b))-(z*x*v1*v2)) <= x*x*eps2,label+" (z*a)*(x*b)");
  Assert(std::abs(((x*a)*(z*b))-(x*z*v1*v2)) <= x*x*eps2,label+" (x*a)*(z*b)");
  Assert(std::abs(((z*a)*(z*b))-(z*z*v1*v2)) <= x*x*eps2,label+" (z*a)*(z*b)");
}

template <class T, class T2, class V1, class V2> inline void TestVV2(
    CONST V1& a, const V2& b, std::string label)
{

  tmv::Vector<T> v1 = a;
  tmv::Vector<T2> v2 = b;
  Assert(Norm(a-v1) <= EPS*Norm(v1),label+" a != v1");
  Assert(Norm(b-v2) <= EPS*Norm(v2),label+" b != v2");

  double normv = Norm(v1)+Norm(v2);

  tmv::Vector<T> v3 = v1;
  tmv::Vector<T> v4 = v1;
#ifndef NOMIX
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
#endif

  a += b;
  v4 = v1+v2;
  Assert(Norm(VEC(T,a)-v4) <= EPS*normv,label+" a += b");
  a = v4 = v1;
  a -= b;
  v4 = v1-v2;
  Assert(Norm(VEC(T,a)-v4) <= EPS*normv,label+" a -= b");
  a = v4 = v1;
  a += -b;
  v4 = v1-v2;
  Assert(Norm(VEC(T,a)-v4) <= EPS*normv,label+" a += -b");
  a = v4 = v1;
  a -= -b;
  v4 = v1+v2;
  Assert(Norm(VEC(T,a)-v4) <= EPS*normv,label+" a -= -b");
  a = v4 = v1;
#ifdef ALIASOK
  a = a+b;
  v4 = v1+v2;
  Assert(Norm(VEC(T,a)-v4) <= EPS*normv,label+" a = a+b");
  a = v4 = v1;
  a = a-b;
  v4 = v1-v2;
  Assert(Norm(VEC(T,a)-v4) <= EPS*normv,label+" a = a-b");
  a = v4 = v1;
  a = b+a;
  v4 = v1+v2;
  Assert(Norm(VEC(T,a)-v4) <= EPS*normv,label+" a = b+a");
  a = v4 = v1;
  a = b-a;
  v4 = v2-v1;
  Assert(Norm(VEC(T,a)-v4) <= EPS*normv,label+" a = b-a");
  a = v4 = v1;
#endif
}

template <class T, IFTEMP1(class V0) IFTEMP1(class CV0) class V, class CV> inline void TestVectorArith1(
#ifdef ALIASOK
    IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) 
#else
    IFTEMP1(V0& ) IFTEMP1(CV0& ) 
#endif
    CONST V& a, CONST CV& ca, std::string label)
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
  TestVV<T,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,a,label+" R,R");
  TestVV2<T,T>(a,a,label+" R,R");
  TestVV<std::complex<T>,std::complex<T> >(IFTEMP1(ctemp) IFTEMP1(ctemp) 
      ca,ca,label+" C,C");
  TestVV2<std::complex<T>,std::complex<T> >(ca,ca,label+" C,C");
#endif

}
template <class T, IFTEMP1(class V0) IFTEMP1(class CV0) class V1, class CV1, class V2, class CV2> 
inline void TestVectorArith2(
    IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) 
    CONST V1& a, CONST CV1& ca, const V2& b, const CV2& cb,
    std::string label)
{
  TestVV<T,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label+" R,R");
  TestVV2<T,T>(a,b,label+" R,R");
  TestVV<T,std::complex<T> >(IFTEMP1(ctemp) IFTEMP1(ctemp) a,cb,label+" R,C");
  TestVV<std::complex<T>,T>(IFTEMP1(ctemp) IFTEMP1(ctemp) ca,b,label+" C,R");
  TestVV2<std::complex<T>,T>(ca,b,label+" C,R");
  TestVV<std::complex<T>,std::complex<T> >(IFTEMP1(ctemp) IFTEMP1(ctemp) 
      ca,cb,label+" C,C");
  TestVV2<std::complex<T>,std::complex<T> >(ca,cb,label+" C,C");
}

#endif
