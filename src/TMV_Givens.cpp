///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "tmv/TMV_Givens.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VIt.h"

//#include <iostream>

namespace tmv {

#define RT TMV_RealType(T)

    // 
    // GivensRotate
    //
    template <class T> Givens<T> GivensRotate(T& x, T& y)
    {
        //std::cout<<"Start GivensRotate: x,y = "<<x<<','<<y<<std::endl;
        // Use a Givens matrix G to rotate the vector so y = 0:
        // G [ x ] = [ r ]
        //   [ y ]   [ 0 ]
        // Also, return the Givens matrix used to do so.
        const RT sqrteps = TMV_SQRT(TMV_Epsilon<T>());

        RT maxabs_x = TMV_MAXABS(x);
        RT maxabs_y = TMV_MAXABS(y);
        if (maxabs_y == RT(0)) {
            y = RT(0);
            //std::cout<<"Simple1: return 1,0\n";
            return Givens<T>(RT(1),T(0));
        } else if (maxabs_x == RT(0)) {
            x = RT(0);
            RT absy = TMV_ABS(y);
            T s = TMV_SIGN(TMV_CONJ(y),absy);
            x = absy; y = RT(0);
            //std::cout<<"Simple2: return 0,"<<s<<std::endl;
            return Givens<T>(RT(0),s);
        } else {
            // This test is better than abs(x) > abs(y) for complex,
            // since it saves 2 sqrt's and has roughly the same effect
            // of preventing overflow and reducing rounding errors.
            if (maxabs_x > maxabs_y) {
                if (maxabs_y <= sqrteps*maxabs_x) {
                    // Then everything simplifies:
                    // c = 1
                    // s = (y/x)* 
                    // r = f
                    T s = TMV_CONJ(y)/TMV_CONJ(x);
                    y = RT(0);
                    //std::cout<<"Simple3: return 1,"<<s<<std::endl;
                    return Givens<T>(RT(1),s);
                } else {
                    // c = 1/sqrt(1+|y/x|^2)
                    // s = (y/x)*/sqrt(1+|y/x|^2)
                    // r = x sqrt(1+|y/x|^2)
                    // We get a slightly more accurate calculation of r
                    // if we calculate r-x and add this to x:
                    // r-x = x (sqrt(1+|y/x|^2)-1) = x |y/x|^2 / (1 + sqrt(1+|y/x|^2))
                    T yoverx = y/x;
                    RT n = TMV_NORM(yoverx);
                    RT sqrtfactor = TMV_SQRT(RT(1)+n);
                    RT c = RT(1)/sqrtfactor;
                    T s = TMV_CONJ(yoverx)*c;
                    x += x*(n/(RT(1)+sqrtfactor));
                    y = RT(0);
                    //std::cout<<"x>y: return "<<c<<','<<s<<std::endl;
                    return Givens<T>(c,s);
                }
            } else {
                // As above, we store c-1 rather than c
                // even though it isn't small here.
                T xovery = x/y;
                RT n = TMV_NORM(xovery);
                RT absxovery = TMV_SQRT(n);
                if (n <= TMV_Epsilon<T>()) {
                    // c = |x/y|
                    // s = (x/y)/|x/y|
                    // r = x/|x/y|
                    T s = TMV_SIGN(xovery,absxovery);
                    x = s*y;
                    y = RT(0);
                    //std::cout<<"Simple4: return "<<absxovery<<','<<s<<std::endl;
                    return Givens<T>(absxovery,s);
                } else {
                    // c = |x/y|/sqrt(1+|x/y|^2)
                    // s = (x/y)/|x/y|/sqrt(1+|x/y|^2)
                    // r = x/|x/y| sqrt(|x|^2+|y|^2)
                    RT sqrtfactor = TMV_SQRT(RT(1)+n);
                    RT invsqrtfactor = RT(1)/sqrtfactor;
                    T signxovery = TMV_SIGN(xovery,absxovery); // (x/y)/|x/y|
                    T s = signxovery*invsqrtfactor;
                    x = y * signxovery * sqrtfactor; 
                    y = T(0);
                    RT c = absxovery*invsqrtfactor;
                    //std::cout<<"x<y: return "<<c<<','<<s<<std::endl;
                    return Givens<T>(c,s);
                }
            }
        }
    }

    //
    // GivensMult Scalars
    //
    template <class T, class Tx> 
    void GivensMult(RT c, T s, Tx& x, Tx& y) 
    {
        // [ x' ] = [  c  s  ] [ x ] = [  cx+sy  ]
        // [ y' ]   [ -s* c* ] [ y ]   [ c*y-s*x ]
        Tx xx = c*x+s*y;
        y = c*y-TMV_CONJ(s)*x;
        x = xx;
    }
#if 0
    template <class T, class Tx> 
    void GivensMult(RT c, T s, ConjRef<Tx> x, ConjRef<Tx> y)
    {
        Tx xx = c*x+s*y;
        y = c*y-TMV_CONJ(s)*x;
        x = xx;
    }
    template <class T, class Tx> 
    void GivensMult(RT c, T s, VarConjRef<Tx> x, VarConjRef<Tx> y)
    {
        Tx xx = c*x+s*y;
        y = c*y-TMV_CONJ(s)*x;
        x = xx;
    }
#endif

    //
    // GivensMult vectors
    //

    template <bool c0, class T, class Tx> 
    static void NonBlasGivensMult(
        RT c, T s, VectorView<Tx> v1, VectorView<Tx> v2)
    {
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() > 0);
        TMVAssert(v1.step() != 0);
        TMVAssert(v2.step() != 0);
        TMVAssert(!v1.isSameAs(v2));
        TMVAssert(v2.ct() == NonConj);
        TMVAssert(v2.step() != -1);
        TMVAssert(v1.step() != -1 || v2.step() == 1);
        TMVAssert(v2.step() > 0 || v1.step() == 1);
        TMVAssert(c0 == v1.isconj());

        Tx* v1ptr = v1.ptr();
        Tx* v2ptr = v2.ptr();
        const ptrdiff_t step0 = v1.step();
        const ptrdiff_t step1 = v2.step();

        if (step0 == 1 && step1 == 1)
            for(ptrdiff_t i=v1.size();i>0;--i,++v1ptr,++v2ptr) {
#ifdef TMVFLDEBUG
                TMVAssert(v1ptr >= v1._first);
                TMVAssert(v1ptr < v1._last);
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr < v2._last);
#endif
                GivensMult(c,s,*v1ptr,*v2ptr); 
            }
        else
            for(ptrdiff_t i=v1.size();i>0;--i,v1ptr+=step0,v2ptr+=step1) {
#ifdef TMVFLDEBUG
                TMVAssert(v1ptr >= v1._first);
                TMVAssert(v1ptr < v1._last);
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr < v2._last);
#endif
                GivensMult(c,s,*v1ptr,*v2ptr); 
            }
    }

    template <class T, class Tx> 
    static void DoGivensMult(
        RT c, T s, VectorView<Tx> v1, VectorView<Tx> v2)
    {
        if (v1.isconj()) NonBlasGivensMult<true>(c,s,v1,v2); 
        else NonBlasGivensMult<false>(c,s,v1,v2); 
    }

#ifdef BLAS
#ifdef INST_DOUBLE
    static void DoGivensMult(
        double c, double s,
        VectorView<double> v1, VectorView<double> v2)
    { 
        TMVAssert(v1.size()==v2.size());
        TMVAssert(v1.ct()==NonConj);
        TMVAssert(v2.ct()==NonConj);
        int n=v1.size();
        int s1=v1.step();
        int s2=v2.step();
        double* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        double* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(drot) (
            BLASV(n),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
    }
#ifdef BLASZDROT
    static void DoGivensMult(
        double c, double s,
        VectorView<std::complex<double> > v1,
        VectorView<std::complex<double> > v2)
    { 
        TMVAssert(v1.size()==v2.size());
        TMVAssert(v1.ct()==NonConj);
        TMVAssert(v2.ct()==NonConj);
        int n=v1.size();
        int s1=v1.step();
        int s2=v2.step();
        std::complex<double>* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        std::complex<double>* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(zdrot) (
            BLASV(n),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
    }
#endif
#ifdef ELAP
    static void DoGivensMult(
        double c, std::complex<double> s,
        VectorView<std::complex<double> > v1,
        VectorView<std::complex<double> > v2)
    {
        TMVAssert(v1.size()==v2.size());
        TMVAssert(v1.ct()==NonConj);
        TMVAssert(v2.ct()==NonConj);
        int n = v1.size();
        int s1 = v1.step();
        int s2 = v2.step();
        std::complex<double>* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        std::complex<double>* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        LAPNAME(zrot) (LAPV(n),LAPP(v1p),LAPV(s1),LAPP(v2p),LAPV(s2),
                       LAPV(c),LAPP(&s)); 
    }
#endif // ELAP
#endif // DOUBLE
#ifdef INST_FLOAT
    static void DoGivensMult(
        float c, float s,
        VectorView<float> v1, VectorView<float> v2)
    { 
        TMVAssert(v1.size()==v2.size());
        TMVAssert(v1.ct()==NonConj);
        TMVAssert(v2.ct()==NonConj);
        int n=v1.size();
        int s1=v1.step();
        int s2=v2.step();
        float* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        float* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(srot) (
            BLASV(n),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
    }
#ifdef BLASZDROT
    static void DoGivensMult(
        float c, float s,
        VectorView<std::complex<float> > v1,
        VectorView<std::complex<float> > v2)
    { 
        TMVAssert(v1.size()==v2.size());
        TMVAssert(v1.ct()==NonConj);
        TMVAssert(v2.ct()==NonConj);
        int n=v1.size();
        int s1=v1.step();
        int s2=v2.step();
        std::complex<float>* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        std::complex<float>* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(csrot) (
            BLASV(n),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
    }
#endif
#ifdef ELAP
    static void DoGivensMult(
        float c, std::complex<float> s,
        VectorView<std::complex<float> > v1,
        VectorView<std::complex<float> > v2)
    {
        TMVAssert(v1.size()==v2.size());
        TMVAssert(v1.ct()==NonConj);
        TMVAssert(v2.ct()==NonConj);
        int n = v1.size();
        int s1 = v1.step();
        int s2 = v2.step();
        std::complex<float>* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        std::complex<float>* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        LAPNAME(crot) (LAPV(n),LAPP(v1p),LAPV(s1),LAPP(v2p),LAPV(s2),
                       LAPV(c),LAPP(&s)); 
    }
#endif // ELAP
#endif // DOUBLE
#endif // BLAS

    template <class T, class Tx> 
    void GivensMult(
        RT c, T s, VectorView<Tx> v1, VectorView<Tx> v2)
    { 
        TMVAssert(v1.size() == v2.size());
        TMVAssert(!v1.isSameAs(v2));

        if (v1.size() > 0 && s != T(0)) {
            if (shouldReverse(v1.step(),v2.step()))
                GivensMult(c,s,v1.reverse(),v2.reverse());
            else if (v2.isconj()) 
                GivensMult(c,TMV_CONJ(s),v1.conjugate(),v2.conjugate());
#ifdef BLAS
            else if (v1.isconj()) {
                Vector<Tx> v1x = v1;
                DoGivensMult(c,s,v1x.view(),v2);
                v1 = v1x;
            }
#endif
            else DoGivensMult(c,s,v1,v2); 
        }
    }

    // 
    // Symmetric Givens Mult
    //
    template <class T, class Tx> 
    void GivensHermMult(RT c, T s, Tx& d0, Tx& d1, Tx& e0) 
    {
        // [ d0 e0* ] = [  c  s ] [ d0 e0* ] [ c  -s ]
        // [ e0 d1  ]   [ -s* c ] [ e0 d1  ] [ s*  c ]
        // = [ c^2 d0 + 2c Re(s e0) + |s|^2 d1   cs(d1-d0) + c^2 e0* - s^2 e0    ]
        //   [ cs*(d1-d0) + c^2 e0 - s*^2 e0*    c^2 d1 - 2c Re(s e0) + |s|^2 d0 ]
        // (using c^2 = 1-|s|^2:)
        // d0' = d0 + 2c Re(s e0) + |s|^2 (d1-d0)
        // e0' = e0 - 2s* Re(s e0) + c s* (d1-d0)
        // d1' = d1 - 2c Re(s e0) - |s|^2 (d1-d0)

        TMV_RealType(Tx) Rese0 = TMV_REAL(s*e0);
        Tx d1md0 = d1-d0;
        // Note: s might be very small, so don't explicitly form |s|^2.
        Tx dd = RT(2)*c*Rese0 + (s*d1md0)*TMV_CONJ(s);
        d0 += dd;
        d1 -= dd;
        e0 += TMV_CONJ(s)*(c*d1md0 - RT(2)*Rese0);
    }
    template <class T, class Tx> 
    void GivensSymMult(RT c, T s, Tx& d0, Tx& d1, Tx& e0) 
    {
        // [ d0 e0 ] = [  c  s ] [ d0 e0 ] [ c -s* ]
        // [ e0 d1 ]   [ -s* c ] [ e0 d1 ] [ s  c  ]
        // d0' = d0 + 2c (s e0) + s^2 (d1-d0)
        // e0' = e0 - 2|s|^2 e0 + c (s d1 - s* d0)
        // d1' = d1 - 2c (s e0) - s^2 (d1-d0)

        Tx se0 = s*e0;
        Tx d1md0 = d1-d0;
        Tx dd = RT(2)*c*se0 + (s*d1md0)*s;
        d0 += dd;
        d1 -= dd;
        if (isReal(T()))
            e0 += s*(c*d1md0 - RT(2)*se0);
        else
            e0 += c*(s*d1-TMV_CONJ(s)*d0) - RT(2)*TMV_CONJ(s)*se0;
    }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_Givens.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


