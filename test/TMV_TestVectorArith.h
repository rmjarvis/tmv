#ifndef TMV_TESTVECTOR_H
#define TMV_TESTVECTOR_H

#include "TMV_Test.h"

template <class V> 
inline void TestV(const V& a, std::string label)
{
    typedef typename V::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start V "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
    }

    const tmv::Vector<T> v = a;
    const int N = a.size();
    RT eps = EPS*N;

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

    RT normsq(0);
    RT norm1(0);
    RT norminf(0);
    for(int i=0;i<N;++i) {
        normsq += tmv::TMV_NORM(v(i));
        norm1 += tmv::TMV_ABS(v(i));
        if (tmv::TMV_ABS(v(i)) > norminf) norminf = tmv::TMV_ABS(v(i));
    }

    if (!(std::numeric_limits<RT>::is_integer)) {
        Assert(Equal2(Norm2(a),tmv::TMV_SQRT(normsq),eps*tmv::TMV_SQRT(normsq)),
               label+" Norm2");
    }
    if (!(std::numeric_limits<RT>::is_integer && tmv::Traits<T>::iscomplex)) {
        Assert(Equal2(Norm1(a),norm1,eps*norm1),label+" Norm1");
        Assert(Equal2(NormInf(a),norminf,eps*norminf),label+" NormInf");
    }
    Assert(Equal2(NormSq(a),normsq,eps*normsq),label+" NormSq");

    if (showstartdone) {
        std::cout<<"Done V "<<std::endl;
    }
}

template <class V, class T2> 
inline void TestVX(const V& a, T2 x, std::string label)
{
    typedef typename V::value_type T;
    typedef typename tmv::Traits2<T,T2>::type PT;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VX "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"x "<<x<<std::endl;
    }

    const tmv::Vector<T> v = a;
    const int N = a.size();
    RT eps = EPS*N;
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
        std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Equal(a,v,eps),label+" a != v");

    tmv::Vector<PT> v2 = v;

    // v2 = x*v
    for(int i=0;i<N;++i) v2(i) = x*v(i);
    Assert(Equal((x*a),v2,eps),label+" x*a");
    Assert(Equal((a*x),v2,eps),label+" a*x");

    if (!(std::numeric_limits<RT>::is_integer)) {
        // v2 = v / x
        for(int i=0;i<N;++i) v2(i) = v(i)/x;
        Assert(Equal((a/x),v2,eps),label+" a/x");
    }

    if (showstartdone) {
        std::cout<<"Done VX "<<std::endl;
    }
}

template <class V, class T2> 
inline void TestVX2(V& a, T2 x, std::string label)
{
    typedef typename V::value_type T;
    typedef typename tmv::Traits2<T,T2>::type PT;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VX2 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"x "<<x<<std::endl;
    }

    const tmv::Vector<T> v = a;
    const int N = a.size();
    RT eps = EPS*N;
    if (!(std::numeric_limits<RT>::is_integer)) 
        eps *= Norm(v) * tmv::TMV_ABS2(x);

    Assert(Equal(a,v,eps),label+" a != v");

    tmv::Vector<PT> v2 = v;

    if (XXDEBUG3) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"v = "<<tmv::TMV_Text(v)<<" = "<<v<<std::endl;
        std::cout<<"x = "<<x<<std::endl;
        std::cout<<"a*=x = "<<(a*=x)<<std::endl;
        std::cout<<"x*v = "<<(x*v)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(a-(x*v))<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
        a = v;
    }

    // v2 = x*v
    for(int i=0;i<N;++i) v2(i) = x*v(i);
    (a=v) *= x;
    Assert(Equal(a,v2,eps),label+" a *= x");
    a = v;
    Assert(Equal((a*=x),v2,eps),label+" a *= x (2)");
#ifdef ALIASOK
    (a=v) = a*x;
    Assert(Equal(a,v2,eps),label+" a = a*x");
    (a=v) = x*a;
    Assert(Equal(a,v2,eps),label+" a = x*a");
#endif

    if (!(std::numeric_limits<RT>::is_integer)) {
        // v2 = v/x
        for(int i=0;i<N;++i) v2(i) = v(i)/x;
        (a=v) /= x;
        Assert(Equal(a,v2,eps),label+" a /= x");
#ifdef ALIASOK
        (a=v) = a/x;
        Assert(Equal(a,v2,eps),label+" a = a/x");
#endif
    }
    a = v;

    if (showstartdone) {
        std::cout<<"Done VX2 "<<std::endl;
    }
}

template <class V1, class V2> 
inline void TestVV(const V1& a, const V2& b, std::string label)
{
    typedef typename V1::value_type T;
    typedef typename V2::value_type T2;
    typedef typename tmv::Traits2<T,T2>::type PT;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VV "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }

    const tmv::Vector<T> v1 = a;
    const tmv::Vector<T2> v2 = b;
    const int N=a.size();
    RT eps = EPS;
    RT eps2 = EPS;
    Assert(Equal(a,v1,eps),label+" a != v1");
    Assert(Equal(b,v2,eps),label+" b != v2");
    if (!(std::numeric_limits<RT>::is_integer)) {
        eps *= Norm(v1) + Norm(v2);
        eps2 *= Norm(v1) * Norm(v2);
    }
    tmv::Vector<std::complex<RT> > v3 = v1;

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

    for(int i=0;i<N;++i) v3(i) = v1(i) + v2(i);
    Assert(Equal((a+v2),v3,eps),label+" a+v");
    Assert(Equal((v1+b),v3,eps),label+" v+b");
    Assert(Equal((a+b),v3,eps),label+" a+b");

    for(int i=0;i<N;++i) v3(i) = v1(i) - v2(i);
    Assert(Equal((a-v2),v3,eps),label+" a-v");
    Assert(Equal((v1-b),v3,eps),label+" v-b");
    Assert(Equal((a-b),v3,eps),label+" a-b");

    for(int i=0;i<N;++i) v3(i) = v1(i) * v2(i);
    Assert(Equal(ElemProd(a,v2),v3,eps),label+" ElemProd(a,v)");
    Assert(Equal(ElemProd(v1,b),v3,eps),label+" ElemProd(v,b)");
    Assert(Equal(ElemProd(a,b),v3,eps),label+" ElemProd(a,b)");

    PT prod = T(0);
    for(int i=0;i<N;++i) prod += v1(i) * v2(i);
    Assert(Equal2((a*v2),prod,eps2),label+" a*v");
    Assert(Equal2((v1*b),prod,eps2),label+" v*b");
    Assert(Equal2((a*b),prod,eps2),label+" a*b");

    RT x(5);
    std::complex<RT> z(3,4);
    if (XXDEBUG4) {
        std::cout<<"a-x*b = "<<(a-x*b)<<std::endl;
        std::cout<<"v1-x*v2 = "<<(v1-x*v2)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm((a-x*b)-(v1-x*v2))<<std::endl;
        std::cout<<"eps*x = "<<x*eps<<std::endl;
    }

    for(int i=0;i<N;++i) v3(i) = v1(i) + x*v2(i);
    Assert(Equal((a+x*b),v3,x*eps),label+" a+x*b");
    for(int i=0;i<N;++i) v3(i) = v1(i) - x*v2(i);
    Assert(Equal((a-x*b),v3,x*eps),label+" a-x*b");
    for(int i=0;i<N;++i) v3(i) = x*v1(i) + v2(i);
    Assert(Equal((x*a+b),v3,x*eps),label+" x*a+b");
    for(int i=0;i<N;++i) v3(i) = x*v1(i) - v2(i);
    Assert(Equal((x*a-b),v3,x*eps),label+" x*a-b");
    for(int i=0;i<N;++i) v3(i) = x*v1(i) + x*v2(i);
    Assert(Equal((x*a+x*b),v3,x*eps),label+" x*a+x*b");
    for(int i=0;i<N;++i) v3(i) = x*v1(i) - x*v2(i);
    Assert(Equal((x*a-x*b),v3,x*eps),label+" x*a-x*b");

    for(int i=0;i<N;++i) v3(i) = v1(i) + z*v2(i);
    Assert(Equal((a+z*b),v3,x*eps),label+" a+z*b");
    for(int i=0;i<N;++i) v3(i) = v1(i) - z*v2(i);
    Assert(Equal((a-z*b),v3,x*eps),label+" a-z*b");
    if (XXDEBUG4) {
        std::cout<<"z*a-b = "<<(z*a-b)<<std::endl;
        std::cout<<"z*v1-v2 = "<<(z*v1-v2)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm((z*a-b)-(z*v1-v2))<<std::endl;
        std::cout<<"eps*x = "<<x*eps<<std::endl;
    }
    for(int i=0;i<N;++i) v3(i) = z*v1(i) + v2(i);
    Assert(Equal((z*a+b),v3,x*eps),label+" z*a+b");
    for(int i=0;i<N;++i) v3(i) = z*v1(i) - v2(i);
    Assert(Equal((z*a-b),v3,x*eps),label+" z*a-b");
    for(int i=0;i<N;++i) v3(i) = z*v1(i) + z*v2(i);
    Assert(Equal((z*a+z*b),v3,x*eps),label+" z*a+z*b");
    for(int i=0;i<N;++i) v3(i) = z*v1(i) - z*v2(i);
    Assert(Equal((z*a-z*b),v3,x*eps),label+" z*a-z*b");

    for(int i=0;i<N;++i) v3(i) = z*v1(i) + x*v2(i);
    Assert(Equal((z*a+x*b),v3,x*eps),label+" z*a+x*b");
    for(int i=0;i<N;++i) v3(i) = z*v1(i) - x*v2(i);
    Assert(Equal((z*a-x*b),v3,x*eps),label+" z*a-x*b");
    for(int i=0;i<N;++i) v3(i) = x*v1(i) + z*v2(i);
    Assert(Equal((x*a+z*b),v3,x*eps),label+" x*a+z*b");
    for(int i=0;i<N;++i) v3(i) = x*v1(i) - z*v2(i);
    Assert(Equal((x*a-z*b),v3,x*eps),label+" x*a-z*b");

    Assert(Equal2(((x*a)*b),x*prod,x*eps2),label+" (x*a)*b");
    Assert(Equal2((a*(x*b)),x*prod,x*eps2),label+" a*(x*b)");
    Assert(Equal2((x*(a*b)),x*prod,x*eps2),label+" x*(a*b)");

    Assert(Equal2(((z*a)*b),z*prod,x*eps2),label+" (z*a)*b");
    Assert(Equal2((a*(z*b)),z*prod,x*eps2),label+" a*(z*b)");
    Assert(Equal2((z*(a*b)),z*prod,x*eps2),label+" z*(a*b)");

    Assert(Equal2(((x*a)*(x*b)),x*x*prod,x*x*eps2),label+" (x*a)*(x*b)");
    Assert(Equal2(((z*a)*(x*b)),x*z*prod,x*x*eps2),label+" (z*a)*(x*b)");
    Assert(Equal2(((x*a)*(z*b)),x*z*prod,x*x*eps2),label+" (x*a)*(z*b)");
    Assert(Equal2(((z*a)*(z*b)),z*z*prod,x*x*eps2),label+" (z*a)*(z*b)");

    if (showstartdone) {
        std::cout<<"Done VV "<<std::endl;
    }
}

template <class V1, class V2> 
inline void TestVV2(V1& a, const V2& b, std::string label)
{
    typedef typename V1::value_type T;
    typedef typename V2::value_type T2;
    typedef typename tmv::Traits2<T,T2>::type PT;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VV2 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }

    const tmv::Vector<T> v1 = a;
    const tmv::Vector<T2> v2 = b;
    const int N = a.size();
    RT eps = EPS;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(v1) + Norm(v2);
    }
    Assert(Equal(a,v1,eps),label+" a != v1");
    Assert(Equal(b,v2,eps),label+" b != v2");

    tmv::Vector<PT> v3 = v1;
    tmv::Vector<PT> v4 = v1;
    if (XXDEBUG5) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<" = "<<b<<std::endl;
        std::cout<<"v3+=b = "<<(v3+=b)<<std::endl;
        std::cout<<"v1+v2 = "<<(v1+v2)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(v3-(v1+v2))<<std::endl;
        std::cout<<"a+=b = "<<(a+=b)<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(a-(v1+v2))<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
        a = v1;
    }

    // v4 = v1 + v2
    for(int i=0;i<N;++i) v4(i) = v1(i) + v2(i);
    (v3=v1) += b;
    Assert(Equal(v3,v4,eps),label+" v += b");
    (v3=v1) = v3+b;
    Assert(Equal(v3,v4,eps),label+" v = v+b");
    (v3=v1) = b+v3;
    Assert(Equal(v3,v4,eps),label+" v = b+v");
    (a=v1) += v2;
    Assert(Equal(a,v4,eps),label+" a += v");
    (a=v1) += b;
    Assert(Equal(a,v4,eps),label+" a += b");
    (a=v1) -= -b;
    Assert(Equal(a,v4,eps),label+" a -= -b");
#ifdef ALIASOK
    (a=v1) = a+v2;
    Assert(Equal(a,v4,eps),label+" a = a+v");
    (a=v1) = v2+a;
    Assert(Equal(a,v4,eps),label+" a = v+a");
    (a=v1) = a+b;
    Assert(Equal(a,v4,eps),label+" a = a+b");
    (a=v1) = b+a;
    Assert(Equal(a,v4,eps),label+" a = b+a");
#endif
    a = v1;

    // v4 = v1 - v2
    for(int i=0;i<N;++i) v4(i) = v1(i) - v2(i);
    (v3=v1) -= b;
    Assert(Equal(v3,v4,eps),label+" v -= b");
    (v3=v1) = v3-b;
    Assert(Equal(v3,v4,eps),label+" v = v-b");
    (a=v1) -= v2;
    Assert(Equal(a,v4,eps),label+" a -= v");
    (a=v1) -= b;
    Assert(Equal(a,v4,eps),label+" a -= b");
    (a=v1) += -b;
    Assert(Equal(a,v4,eps),label+" a += -b");
#ifdef ALIASOK
    (a=v1) = a-v2;
    Assert(Equal(a,v4,eps),label+" a = a-v");
    (a=v1) = a-b;
    Assert(Equal(a,v4,eps),label+" a = a-b");
#endif
    a = v1;

    // v4 = v2 - v1
    for(int i=0;i<N;++i) v4(i) = v2(i) - v1(i);
    (v3=v1) = b-v3;
    Assert(Equal(v3,v4,eps),label+" v = b-v");
#ifdef ALIASOK
    (a=v1) = v2-a;
    Assert(Equal(a,v4,eps),label+" a = v-a");
    (a=v1) = b-a;
    Assert(Equal(a,v4,eps),label+" a = b-a");
#endif
    a = v1;

    // v4 = ElemProd(v1,v2)
    for(int i=0;i<N;++i) v4(i) = v1(i) * v2(i);
    (v3=v1) = ElemProd(v3,b);
    Assert(Equal(v3,v4,eps),label+" v = ElemProd(v,b)");
    (v3=v1) = ElemProd(b,v3);
    Assert(Equal(v3,v4,eps),label+" v = ElemProd(b,v)");
#ifdef ALIASOK
    (a=v1) = ElemProd(a,v2);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(a,v)");
    (a=v1) = ElemProd(v2,a);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(v,a)");
    (a=v1) = ElemProd(a,b);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(a,b)");
    (a=v1) = ElemProd(b,a);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(b,a)");
#endif
    a = v1;

    // v4 = v1 + ElemProd(v1,v2)
    for(int i=0;i<N;++i) v4(i) = v1(i) + v1(i) * v2(i);
    (v3=v1) += ElemProd(v3,b);
    Assert(Equal(v3,v4,eps),label+" v += ElemProd(v,b)");
    (v3=v1) += ElemProd(b,v3);
    Assert(Equal(v3,v4,eps),label+" v += ElemProd(b,v)");
#ifdef ALIASOK
    (a=v1) += ElemProd(a,v2);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(a,v)");
    (a=v1) += ElemProd(v2,a);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(v,a)");
    (a=v1) += ElemProd(a,b);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(a,b)");
    (a=v1) += ElemProd(b,a);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(b,a)");
#endif
    a = v1;

    // v4 = v1 - ElemProd(v1,v2)
    for(int i=0;i<N;++i) v4(i) = v1(i) - v1(i) * v2(i);
    (v3=v1) -= ElemProd(v3,b);
    Assert(Equal(v3,v4,eps),label+" v -= ElemProd(v,b)");
    (v3=v1) -= ElemProd(b,v3);
    Assert(Equal(v3,v4,eps),label+" v -= ElemProd(b,v)");
#ifdef ALIASOK
    (a=v1) -= ElemProd(a,v2);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(a,v)");
    (a=v1) -= ElemProd(v2,a);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(v,a)");
    (a=v1) -= ElemProd(a,b);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(a,b)");
    (a=v1) -= ElemProd(b,a);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(b,a)");
#endif
    a = v1;

    // v4 = x*ElemProd(v1,v2)
    RT x = 5;
    for(int i=0;i<N;++i) v4(i) = x * v1(i) * v2(i);
    (v3=v1) = x*ElemProd(v3,b);
    Assert(Equal(v3,v4,eps),label+" v = x*ElemProd(v,b)");
    (v3=v1) = x*ElemProd(b,v3);
    Assert(Equal(v3,v4,eps),label+" v = x*ElemProd(b,v)");
    (v3=v1) = ElemProd(x*v3,b);
    Assert(Equal(v3,v4,eps),label+" v = ElemProd(x*v,b)");
    (v3=v1) = ElemProd(v3,x*b);
    Assert(Equal(v3,v4,eps),label+" v = ElemProd(v,x*b)");
    (v3=v1) = ElemProd(x*b,v3);
    Assert(Equal(v3,v4,eps),label+" v = ElemProd(x*b,v)");
    (v3=v1) = ElemProd(b,x*v3);
    Assert(Equal(v3,v4,eps),label+" v = ElemProd(b,x*v)");
#ifdef ALIASOK
    (a=v1) = x*ElemProd(a,v2);
    Assert(Equal(a,v4,eps),label+" a = x*ElemProd(a,v)");
    (a=v1) = x*ElemProd(v2,a);
    Assert(Equal(a,v4,eps),label+" a = x*ElemProd(v,a)");
    (a=v1) = ElemProd(x*a,v2);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(x*a,v)");
    (a=v1) = ElemProd(a,x*v2);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(a,x*v)");
    (a=v1) = ElemProd(x*v2,a);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(x*v,a)");
    (a=v1) = ElemProd(v2,x*a);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(v,x*a)");
    (a=v1) = x*ElemProd(a,b);
    Assert(Equal(a,v4,eps),label+" a = x*ElemProd(a,b)");
    (a=v1) = x*ElemProd(b,a);
    Assert(Equal(a,v4,eps),label+" a = x*ElemProd(b,a)");
    (a=v1) = ElemProd(x*a,b);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(x*a,b)");
    (a=v1) = ElemProd(a,x*b);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(a,x*b)");
    (a=v1) = ElemProd(x*b,a);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(x*b,a)");
    (a=v1) = ElemProd(b,x*a);
    Assert(Equal(a,v4,eps),label+" a = ElemProd(b,x*a)");
#endif
    a = v1;

    // v4 = v1 + x*ElemProd(v1,v2)
    for(int i=0;i<N;++i) v4(i) = v1(i) + x * v1(i) * v2(i);
    (v3=v1) += x*ElemProd(v3,b);
    Assert(Equal(v3,v4,eps),label+" v += x*ElemProd(v,b)");
    (v3=v1) += x*ElemProd(b,v3);
    Assert(Equal(v3,v4,eps),label+" v += x*ElemProd(b,v)");
    (v3=v1) += ElemProd(x*v3,b);
    Assert(Equal(v3,v4,eps),label+" v += ElemProd(x*v,b)");
    (v3=v1) += ElemProd(v3,x*b);
    Assert(Equal(v3,v4,eps),label+" v += ElemProd(v,x*b)");
    (v3=v1) += ElemProd(x*b,v3);
    Assert(Equal(v3,v4,eps),label+" v += ElemProd(x*b,v)");
    (v3=v1) += ElemProd(b,x*v3);
    Assert(Equal(v3,v4,eps),label+" v += ElemProd(b,x*v)");
#ifdef ALIASOK
    (a=v1) += x*ElemProd(a,v2);
    Assert(Equal(a,v4,eps),label+" a += x*ElemProd(a,v)");
    (a=v1) += x*ElemProd(v2,a);
    Assert(Equal(a,v4,eps),label+" a += x*ElemProd(v,a)");
    (a=v1) += ElemProd(x*a,v2);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(x*a,v)");
    (a=v1) += ElemProd(a,x*v2);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(a,x*v)");
    (a=v1) += ElemProd(x*v2,a);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(x*v,a)");
    (a=v1) += ElemProd(v2,x*a);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(v,x*a)");
    (a=v1) += x*ElemProd(a,b);
    Assert(Equal(a,v4,eps),label+" a += x*ElemProd(a,b)");
    (a=v1) += x*ElemProd(b,a);
    Assert(Equal(a,v4,eps),label+" a += x*ElemProd(b,a)");
    (a=v1) += ElemProd(x*a,b);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(x*a,b)");
    (a=v1) += ElemProd(a,x*b);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(a,x*b)");
    (a=v1) += ElemProd(x*b,a);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(x*b,a)");
    (a=v1) += ElemProd(b,x*a);
    Assert(Equal(a,v4,eps),label+" a += ElemProd(b,x*a)");
#endif
    a = v1;

    // v4 = v1 - x*ElemProd(v1,v2)
    for(int i=0;i<N;++i) v4(i) = v1(i) - x * v1(i) * v2(i);
    (v3=v1) -= x*ElemProd(v3,b);
    Assert(Equal(v3,v4,eps),label+" v -= x*ElemProd(v,b)");
    (v3=v1) -= x*ElemProd(b,v3);
    Assert(Equal(v3,v4,eps),label+" v -= x*ElemProd(b,v)");
    (v3=v1) -= ElemProd(x*v3,b);
    Assert(Equal(v3,v4,eps),label+" v -= ElemProd(x*v,b)");
    (v3=v1) -= ElemProd(v3,x*b);
    Assert(Equal(v3,v4,eps),label+" v -= ElemProd(v,x*b)");
    (v3=v1) -= ElemProd(x*b,v3);
    Assert(Equal(v3,v4,eps),label+" v -= ElemProd(x*b,v)");
    (v3=v1) -= ElemProd(b,x*v3);
    Assert(Equal(v3,v4,eps),label+" v -= ElemProd(b,x*v)");
#ifdef ALIASOK
    (a=v1) -= x*ElemProd(a,v2);
    Assert(Equal(a,v4,eps),label+" a -= x*ElemProd(a,v)");
    (a=v1) -= x*ElemProd(v2,a);
    Assert(Equal(a,v4,eps),label+" a -= x*ElemProd(v,a)");
    (a=v1) -= ElemProd(x*a,v2);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(x*a,v)");
    (a=v1) -= ElemProd(a,x*v2);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(a,x*v)");
    (a=v1) -= ElemProd(x*v2,a);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(x*v,a)");
    (a=v1) -= ElemProd(v2,x*a);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(v,x*a)");
    (a=v1) -= x*ElemProd(a,b);
    Assert(Equal(a,v4,eps),label+" a -= x*ElemProd(a,b)");
    (a=v1) -= x*ElemProd(b,a);
    Assert(Equal(a,v4,eps),label+" a -= x*ElemProd(b,a)");
    (a=v1) -= ElemProd(x*a,b);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(x*a,b)");
    (a=v1) -= ElemProd(a,x*b);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(a,x*b)");
    (a=v1) -= ElemProd(x*b,a);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(x*b,a)");
    (a=v1) -= ElemProd(b,x*a);
    Assert(Equal(a,v4,eps),label+" a -= ElemProd(b,x*a)");
#endif
    a = v1;

    if (showstartdone) {
        std::cout<<"Done VV2 "<<std::endl;
    }
}

template <class V, class CV> 
inline void TestVectorArith1(V& a, CV& ca, std::string label)
{
    typedef typename V::value_type T;
    if (showstartdone) {
        std::cout<<"Start VectorArith1 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
    }

    TestV(a,label+" R");
    TestV(ca,label+" C");

    T x = 12;
    std::complex<T> z(9,-2);

    TestVX(a,x,label+" R,R");
    TestVX2(a,x,label+" R,R");
    TestVX(a,z,label+" R,C");
    TestVX(ca,x,label+" C,R");
    TestVX2(ca,x,label+" C,R");
    TestVX(ca,z,label+" C,C");
    TestVX2(ca,z,label+" C,C");
    TestVV(a,a,label+" R,R");
    TestVV(ca,ca,label+" C,C");
#ifdef ALIASOK
    TestVV2(a,a,label+" R,R self_arith");
    TestVV2(ca,ca,label+" C,C self_arith");
#endif
}

template <class V1, class CV1, class V2, class CV2> 
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

    TestVV(a,b,label+" R,R");
    TestVV2(a,b,label+" R,R");
    TestVV(a,cb,label+" R,C");
    TestVV(ca,b,label+" C,R");
    TestVV2(ca,b,label+" C,R");
    TestVV(ca,cb,label+" C,C");
    TestVV2(ca,cb,label+" C,C");
}

#endif
