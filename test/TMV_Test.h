#ifndef TMV_TEST_H
#define TMV_TEST_H

#ifndef XTEST
#define XTEST 0
#endif

#define TMV_TEXT

#include <iostream>
#include <typeinfo>
#include <cmath>
#include "tmv/TMV_Base.h"

#define EPS (10*tmv::TMV_Epsilon<FT>())

#ifndef NO_TEST_DOUBLE
#define TEST_DOUBLE
#endif

#ifndef NO_TEST_FLOAT
#define TEST_FLOAT
#endif

#ifndef NO_TEST_COMPLEX
#define TEST_COMPLEX
#endif

extern bool showtests;
extern bool showacc;
extern bool showdiv;
extern bool donorm2;
extern bool showstartdone;
extern bool aliasok;
extern bool symoprod;
extern bool dontthrow;
extern std::string lastsuccess;

void PreAssert(std::string s);
void DoAssert(bool x, std::string s);

#define Assert(x,s) \
    do {  \
        PreAssert(s);  \
        DoAssert(x,s); \
    } while (false)

template <class M1, class M2, class T>
inline bool Equal(const M1& a, const M2& b, T eps)
{
    T normdiff = Norm(a-b);
#ifdef XXD
    if (showacc && !(normdiff <= eps)) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"a-b = "<<a-b<<std::endl;
        std::cout<<"Norm(a-b) = "<<normdiff<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
    }
    if (showtests) std::cout<<"  "<<normdiff<<" <=? "<<eps<<"  ";
#endif
    return normdiff <= eps; 
}
template <class X1, class X2, class T>
inline bool Equal2(const X1& a, const X2& b, T eps)
{
    T absdiff = tmv::TMV_ABS2(a-b);
#ifdef XXD
    if (showacc && !(absdiff <= eps)) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"abs2(a-b) = "<<absdiff<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
    }
    if (showtests) std::cout<<"  "<<absdiff<<" <=? "<<eps<<"  ";
#endif
    return absdiff <= eps;
}

template <class M1, class M2>
inline bool Equal(const M1& a, const M2& b, int )
{
    bool eq = (a == b);
#ifdef XXD
    if (showacc && !(eq)) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"Norm(a-b) = "<<Norm(a-b)<<std::endl;
        std::cout<<"a == b = "<<eq<<std::endl;
    }
#endif
    return eq;
}
template <class X1, class X2>
inline bool Equal2(const X1& a, const X2& b, int )
{
    bool eq = (a == b);
#ifdef XXD
    if (showacc && !(eq)) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"abs(a-b) = "<<tmv::TMV_ABS2(a-b)<<std::endl;
        std::cout<<"a == b = "<<eq<<std::endl;
    }
#endif
    return eq;
}

template <class T>
static inline std::string Text(const T&)
{ return std::string("Unknown (") + typeid(T).name() + ")"; }

static inline std::string Text(const double&)
{ return "double"; }

static inline std::string Text(const float&)
{ return "float"; }

static inline std::string Text(const int&)
{ return "int"; }

static inline std::string Text(const long double&)
{ return "long double"; }

template <class T>
static inline std::string Text(std::complex<T>)
{ return std::string("complex<") + Text(T()) + ">"; }

static inline std::string Text(tmv::DivType d)
{
    return
        d==tmv::LU ? "LU" :
        d==tmv::CH ? "CH" :
        d==tmv::QR ? "QR" :
        d==tmv::QRP ? "QRP" :
        d==tmv::SV ? "SV" : "XX";
}

static inline std::string Text(tmv::StorageType s)
{ 
    return
        s==tmv::ColMajor ? "ColMajor" :
        s==tmv::RowMajor ? "RowMajor" :
        s==tmv::DiagMajor ? "DiagMajor" :
        "NonMajor";
}

extern bool XXDEBUG1;
extern bool XXDEBUG2;
extern bool XXDEBUG3;
extern bool XXDEBUG4;
extern bool XXDEBUG5;
extern bool XXDEBUG6;
extern bool XXDEBUG7;
extern bool XXDEBUG8;
extern bool XXDEBUG9;

#endif
