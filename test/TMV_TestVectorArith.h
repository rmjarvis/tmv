#ifndef TMV_TESTVECTOR_H
#define TMV_TESTVECTOR_H

#include "TMV_Test.h"

template <class T, class V> 
inline void TestV(const V& a, std::string label)
{
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start V "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
    }

    tmv::Vector<T> v = a;
    RT eps = EPS*v.size();

    if (XXDEBUG1) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"v = "<<tmv::TMV_Text(v)<<" = "<<v<<std::endl;
        std::cout<<"a-v = "<<a-v<<std::endl;
        std::cout<<"Norm(a-v) = "<<Norm(a-v)<<std::endl;
        std::cout<<"Norm1(a) = "<<Norm1(a)<<"  "<<Norm1(v)<<std::endl;
        std::cout<<"Norm2(a) = "<<Norm2(a)<<"  "<<Norm2(v)<<std::endl;
        std::cout<<"NormInf(a) = "<<NormInf(a)<<"  "<<NormInf(v)<<std::endl;
        std::cout<<"NormSq(a) = "<<NormSq(a)<<"  "<<NormSq(v)<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS(NormSq(a)-NormSq(v))<<std::endl;
        std::cout<<"eps*normsq = "<<eps*NormSq(v)<<std::endl;
    }

    Assert(Equal(a,v,eps),label+" a != v");

    if (!(std::numeric_limits<RT>::is_integer)) {
        Assert(Equal2(Norm2(a),Norm2(v),eps*Norm2(v)),label+" Norm2");
    }
    if (!(std::numeric_limits<RT>::is_integer && tmv::Traits<T>::iscomplex)) {
        Assert(Equal2(Norm1(a),Norm1(v),eps*Norm1(v)),label+" Norm1");
        Assert(Equal2(NormInf(a),NormInf(v),eps*NormInf(v)),label+" NormInf");
    }
    Assert(Equal2(NormSq(a),NormSq(v),eps*NormSq(v)),label+" NormSq");

    if (showstartdone) {
        std::cout<<"Done V "<<std::endl;
    }
}

template <class T, class V, class T2> 
inline void TestVX(const V& a, T2 x, std::string label)
{
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VX "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"x "<<x<<std::endl;
    }

    tmv::Vector<T> v = a;
    RT eps = EPS*v.size();
    if (!(std::numeric_limits<RT>::is_integer)) 
        eps *= Norm(v) * tmv::TMV_ABS2(x);

    if (XXDEBUG2) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"v = "<<tmv::TMV_Text(v)<<" = "<<v<<std::endl;
        std::cout<<"a-v = "<<a-v<<std::endl;
        std::cout<<"Norm(a-v) = "<<Norm(a-v)<<std::endl;
        std::cout<<"x*a = "<<(x*a)<<std::endl;
        std::cout<<"x*v = "<<(x*v)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm((x*a)-(x*v))<<std::endl;
        std::cout<<"eps*norm = "<<eps<<std::endl;
    }
    Assert(Equal(a,v,eps),label+" a != v");

    Assert(Equal((x*a),(x*v),eps),label+" x*a");
    Assert(Equal((a*x),(x*v),eps),label+" a*x");
    if (!(std::numeric_limits<RT>::is_integer)) {
        Assert(Equal((a/x),(v/x),eps),label+" a/x");
    }

    if (showstartdone) {
        std::cout<<"Done VX "<<std::endl;
    }
}

template <class T, class V, class T2> 
inline void TestVX2(V& a, T2 x, std::string label)
{
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VX2 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"x "<<x<<std::endl;
    }

    tmv::Vector<T> v = a;
    RT eps = EPS*v.size();
    if (!(std::numeric_limits<RT>::is_integer)) 
        eps *= Norm(v) * tmv::TMV_ABS2(x);

    Assert(Equal(a,v,eps),label+" a != v");

    typename V::copy_type a0 = a;

    if (XXDEBUG3) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"v = "<<tmv::TMV_Text(v)<<" = "<<v<<std::endl;
        std::cout<<"x = "<<x<<std::endl;
        std::cout<<"a*=x = "<<(a*=x)<<std::endl;
        std::cout<<"x*v = "<<(x*v)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(a-(x*v))<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
        a = a0;
    }

    a *= x;
    v = tmv::Vector<T>(v*x);
    Assert(Equal(a,v,eps),label+" a *= x");
    Assert(Equal((a*=x),(v*=x),eps),label+" a *= x (2)");
    a = v = a0;

    if (!(std::numeric_limits<RT>::is_integer)) {
        a /= x;
        v = tmv::Vector<T>(v/x);
        Assert(Equal(a,v,eps),label+" a /= x");
        a = v = a0;
    }

#ifdef ALIASOK
    a = a*x;
    v = tmv::Vector<T>(v*x);
    Assert(Equal(a,v,eps),label+" a = a*x");
    a = v = a0;

    a = x*a;
    v = tmv::Vector<T>(v*x);
    Assert(Equal(a,v,eps),label+" a = x*a");
    a = v = a0;

    if (!(std::numeric_limits<RT>::is_integer)) {
        a = a/x;
        v = tmv::Vector<T>(v/x);
        Assert(Equal(a,v,eps),label+" a = a/x");
        a = a0;
    }
#endif

    if (showstartdone) {
        std::cout<<"Done VX2 "<<std::endl;
    }
}

template <class T, class T2, class V1, class V2> 
inline void TestVV(const V1& a, const V2& b, std::string label)
{
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VV "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }

    tmv::Vector<T> v1 = a;
    tmv::Vector<T2> v2 = b;
    RT eps = EPS;
    RT eps2 = EPS;
    Assert(Equal(a,v1,eps),label+" a != v1");
    Assert(Equal(b,v2,eps),label+" b != v2");
    if (!(std::numeric_limits<RT>::is_integer)) {
        eps *= Norm(v1) + Norm(v2);
        eps2 *= Norm(v1) * Norm(v2);
    }

    Assert(Equal(a,v1,eps),label+" a != v1");
    Assert(Equal(b,v2,eps),label+" b != v2");

    if (XXDEBUG4) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<" = "<<b<<std::endl;
        std::cout<<"v1 = "<<tmv::TMV_Text(v1)<<" = "<<v1<<std::endl;
        std::cout<<"v2 = "<<tmv::TMV_Text(v2)<<" = "<<v2<<std::endl;
        std::cout<<"a+b = "<<(a+b)<<std::endl;
        std::cout<<"v1+v2 = "<<(v1+v2)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm((a+b)-(v1+v2))<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
    }

    Assert(Equal((a+b),(v1+v2),eps),label+" a+b");
    Assert(Equal((a-b),(v1-v2),eps),label+" a-b");
    Assert(Equal2((a*b),(v1*v2),eps2),label+" a*b");
    Assert(Equal((a+v2),(v1+v2),eps),label+" a+v");
    Assert(Equal((v1+b),(v1+v2),eps),label+" v+b");
    Assert(Equal((a-v2),(v1-v2),eps),label+" a-v");
    Assert(Equal((v1-b),(v1-v2),eps),label+" v-b");
    Assert(Equal2((a*v2),(v1*v2),eps2),label+" a*v");
    Assert(Equal2((v1*b),(v1*v2),eps2),label+" v*b");

    RT x(5);
    std::complex<RT> z(3,4);
    if (XXDEBUG4) {
        std::cout<<"a-x*b = "<<(a-x*b)<<std::endl;
        std::cout<<"v1-x*v2 = "<<(v1-x*v2)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm((a-x*b)-(v1-x*v2))<<std::endl;
        std::cout<<"eps*x = "<<x*eps<<std::endl;
    }
    Assert(Equal((a-x*b),(v1-x*v2),x*eps),label+" a-x*b");
    Assert(Equal((a+x*b),(v1+x*v2),x*eps),label+" a+x*b");
    Assert(Equal((x*a-b),(x*v1-v2),x*eps),label+" x*a-b");
    Assert(Equal((x*a+b),(x*v1+v2),x*eps),label+" x*a+b");
    Assert(Equal((x*a-x*b),(x*v1-x*v2),x*eps),label+" x*a-x*b");
    Assert(Equal((x*a+x*b),(x*v1+x*v2),x*eps),label+" x*a+x*b");

    Assert(Equal((a-z*b),(v1-z*v2),x*eps),label+" a-z*b");
    Assert(Equal((a+z*b),(v1+z*v2),x*eps),label+" a+z*b");
    if (XXDEBUG4) {
        std::cout<<"z*a-b = "<<(z*a-b)<<std::endl;
        std::cout<<"z*v1-v2 = "<<(z*v1-v2)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm((z*a-b)-(z*v1-v2))<<std::endl;
        std::cout<<"eps*x = "<<x*eps<<std::endl;
    }
    Assert(Equal((z*a-b),(z*v1-v2),x*eps),label+" z*a-b");
    Assert(Equal((z*a+b),(z*v1+v2),x*eps),label+" z*a+b");
    Assert(Equal((z*a-z*b),(z*v1-z*v2),x*eps),label+" z*a-z*b");
    Assert(Equal((z*a+z*b),(z*v1+z*v2),x*eps),label+" z*a+z*b");

    Assert(Equal((z*a-x*b),(z*v1-x*v2),x*eps),label+" z*a-x*b");
    Assert(Equal((z*a+x*b),(z*v1+x*v2),x*eps),label+" z*a+x*b");
    Assert(Equal((x*a-z*b),(x*v1-z*v2),x*eps),label+" x*a-z*b");
    Assert(Equal((x*a+z*b),(x*v1+z*v2),x*eps),label+" x*a+z*b");

    Assert(Equal2(((x*a)*b),(x*v1*v2),x*eps2),label+" (x*a)*b");
    Assert(Equal2((a*(x*b)),(x*v1*v2),x*eps2),label+" a*(x*b)");
    Assert(Equal2((x*(a*b)),(x*v1*v2),x*eps2),label+" x*(a*b)");

    Assert(Equal2(((z*a)*b),(z*v1*v2),x*eps2),label+" (z*a)*b");
    Assert(Equal2((a*(z*b)),(z*v1*v2),x*eps2),label+" a*(z*b)");
    Assert(Equal2((z*(a*b)),(z*v1*v2),x*eps2),label+" z*(a*b)");

    Assert(Equal2(((x*a)*(x*b)),(x*x*v1*v2),x*x*eps2),label+" (x*a)*(x*b)");
    Assert(Equal2(((z*a)*(x*b)),(z*x*v1*v2),x*x*eps2),label+" (z*a)*(x*b)");
    Assert(Equal2(((x*a)*(z*b)),(x*z*v1*v2),x*x*eps2),label+" (x*a)*(z*b)");
    Assert(Equal2(((z*a)*(z*b)),(z*z*v1*v2),x*x*eps2),label+" (z*a)*(z*b)");

    if (showstartdone) {
        std::cout<<"Done VV "<<std::endl;
    }
}

template <class T, class T2, class V1, class V2> 
inline void TestVV2(V1& a, const V2& b, std::string label)
{
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VV2 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }

    tmv::Vector<T> v1 = a;
    tmv::Vector<T2> v2 = b;
    RT eps = EPS;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(v1) + Norm(v2);
    }
    Assert(Equal(a,v1,eps),label+" a != v1");
    Assert(Equal(b,v2,eps),label+" b != v2");

    tmv::Vector<T> v4 = v1;
    tmv::Vector<T> v3 = v1;
    if (XXDEBUG5) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<" = "<<b<<std::endl;
        std::cout<<"v3+=b = "<<(v3+=b)<<std::endl;
        std::cout<<"v1+v2 = "<<(v1+v2)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(v3-(v1+v2))<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
        v3 = v1;
    }

    v3 += b;
    v4 = v1+v2;
    Assert(Equal(v3,v4,eps),label+" v += b");
    v3 = v4 = v1;
    v3 -= b;
    v4 = v1-v2;
    Assert(Equal(v3,v4,eps),label+" v -= b");
    v3 = v4 = v1;
    v3 = v3+b;
    v4 = v1+v2;
    Assert(Equal(v3,v4,eps),label+" v = v+b");
    v3 = v4 = v1;
    v3 = v3-b;
    v4 = v1-v2;
    Assert(Equal(v3,v4,eps),label+" v = v-b");
    v3 = v4 = v1;
    v3 = b+v3;
    v4 = v1+v2;
    Assert(Equal(v3,v4,eps),label+" v = b+v");
    v3 = v4 = v1;
    v3 = b-v3;
    v4 = v2-v1;
    Assert(Equal(v3,v4,eps),label+" v = b-v");
    v3 = v4 = v1;
    a += v2;
    v4 = v1+v2;
    Assert(Equal(a,v4,eps),label+" a += v");
    a = v4 = v1;
    a -= v2;
    v4 = v1-v2;
    Assert(Equal(a,v4,eps),label+" a -= v");
    a = v4 = v1;
#ifdef ALIASOK
    a = a+v2;
    v4 = v1+v2;
    Assert(Equal(a,v4,eps),label+" a = a+v");
    a = v4 = v1;
    a = a-v2;
    v4 = v1-v2;
    Assert(Equal(a,v4,eps),label+" a = a-v");
    a = v4 = v1;
    a = v2+a;
    v4 = v1+v2;
    Assert(Equal(a,v4,eps),label+" a = v+a");
    a = v4 = v1;
    a = v2-a;
    v4 = v2-v1;
    Assert(Equal(a,v4,eps),label+" a = v-a");
    a = v4 = v1;
#endif

    if (XXDEBUG5) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<" = "<<b<<std::endl;
        std::cout<<"a+=b = "<<(a+=b)<<std::endl;
        std::cout<<"a-=b = "<<(a-=b)<<std::endl;
        std::cout<<"v1+v2 = "<<(v1+v2)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(a-(v1+v2))<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
        a = v1;
    }

    a += b;
    v4 = v1+v2;
    Assert(Equal(a,v4,eps),label+" a += b");
    a = v4 = v1;
    a -= b;
    v4 = v1-v2;
    Assert(Equal(a,v4,eps),label+" a -= b");
    a = v4 = v1;
    a += -b;
    v4 = v1-v2;
    Assert(Equal(a,v4,eps),label+" a += -b");
    a = v4 = v1;
    a -= -b;
    v4 = v1+v2;
    Assert(Equal(a,v4,eps),label+" a -= -b");
    a = v4 = v1;
#ifdef ALIASOK
    a = a+b;
    v4 = v1+v2;
    Assert(Equal(a,v4,eps),label+" a = a+b");
    a = v4 = v1;
    a = a-b;
    v4 = v1-v2;
    Assert(Equal(a,v4,eps),label+" a = a-b");
    a = v4 = v1;
    a = b+a;
    v4 = v1+v2;
    Assert(Equal(a,v4,eps),label+" a = b+a");
    a = v4 = v1;
    a = b-a;
    v4 = v2-v1;
    Assert(Equal(a,v4,eps),label+" a = b-a");
    a = v4 = v1;
#endif

    if (showstartdone) {
        std::cout<<"Done VV2 "<<std::endl;
    }
}

template <class T, class V, class CV> 
inline void TestVectorArith1(V& a, CV& ca, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start VectorArith1 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
    }

    TestV<T>(a,label+" R");
    TestV<std::complex<T> >(ca,label+" C");

    T x = 12;
    std::complex<T> z(9,-2);

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
    V1& a, CV1& ca, const V2& b, const CV2& cb, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start VectorArith2 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
    }

    TestVV<T,T>(a,b,label+" R,R");
    TestVV2<T,T>(a,b,label+" R,R");
    TestVV<T,std::complex<T> >(a,cb,label+" R,C");
    TestVV<std::complex<T>,T>(ca,b,label+" C,R");
    TestVV2<std::complex<T>,T>(ca,b,label+" C,R");
    TestVV<std::complex<T>,std::complex<T> >(ca,cb,label+" C,C");
    TestVV2<std::complex<T>,std::complex<T> >(ca,cb,label+" C,C");
}

#endif
