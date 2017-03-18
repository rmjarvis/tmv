#ifndef TMV_TEST_H
#define TMV_TEST_H

#ifndef XTEST
#define XTEST 0
#endif

#include <iostream>
#include <cmath>
#include "tmv/TMV_Base.h"

#ifndef NO_TEST_DOUBLE
#define TEST_DOUBLE
#endif

#ifndef NO_TEST_FLOAT
#define TEST_FLOAT
#endif

#ifndef NO_TEST_COMPLEX
#define TEST_COMPLEX
#endif

#define EPS (10*tmv::TMV_Epsilon<T>())

extern bool showtests;
extern bool showacc;
extern bool showdiv;
extern bool donorm2;
extern bool showstartdone;
extern bool aliasok;
extern bool symoprod;
extern bool dontthrow;
extern std::string lastsuccess;

#ifdef TMV_DEBUG
inline void PreAssert(const std::string& s)
{
    if (showtests) { 
        std::cout<<"Trying: "<<s;  
        std::cout.flush(); 
    } 
}
#else
inline void PreAssert(const std::string& ) {}
#endif

inline void DoAssert(bool x, std::string s)
{
    if (x) { 
#ifdef TMV_DEBUG
        if (showtests) std::cout<<"  Passed"<<std::endl;
        lastsuccess = s; 
#endif
    } else { 
#ifdef TMV_DEBUG
        if (showtests) std::cout<<"  Failed"<<std::endl;
        if (dontthrow) std::cout<<"Failed test: "<<s<<std::endl;  
        else {
#endif
#ifdef NOTHROW
            std::cerr<<"Error in test: "<<s<<std::endl; exit(1); 
#else
            throw tmv::Error("Error in test: "+s);  
#endif
#ifdef TMV_DEBUG
        }
#endif
    } 
}

#define Assert(x,s) \
    do {  \
        PreAssert(s);  \
        DoAssert(x,s); \
    } while (false)

template <class M1, class M2, class T>
inline bool Equal(const M1& a, const M2& b, T eps)
{ 
    T normdiff = Norm(a-b);
    if (showtests) std::cout<<"  "<<normdiff<<" <=? "<<eps<<"  ";
    return normdiff <= eps; 
}
template <class X1, class X2, class T>
inline bool Equal2(const X1& a, const X2& b, T eps)
{
    T absdiff = tmv::TMV_ABS2(a-b);
    if (showtests) std::cout<<"  "<<absdiff<<" <=? "<<eps<<"  ";
    return absdiff <= eps;
}

template <class M1, class M2>
inline bool Equal(const M1& a, const M2& b, int )
{ return a == b; }
template <class X1, class X2>
inline bool Equal2(const X1& a, const X2& b, int )
{ return a == b; }

// C++ I/O doesn't seem capable of reading in at higher accuracy than double precision.
// If you do fin >> x; where x is a long double variable and the text in the file is 
// 1.234, say, then (x-1.234) will not be zero after this.  Instead it will be something
// of order 1.e-17.  So we only check the IO at double precision for long doubles.
template <class M1, class M2, class T>
inline bool EqualIO(const M1& a, const M2& b, T eps)
{ return Equal(a,b,eps); }
template <class M1, class M2>
inline bool EqualIO(const M1& a, const M2& b, long double eps )
{ return Equal(a,b,10*tmv::TMV_Epsilon<double>()); }


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
