///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
//
// This file defines some basic helper functions used by the TMV Matrix
// and Vector classes.  
//
// Probably the only thing someone might care about here is the
// class tmv::tmv_exception which gets thrown on run-time errors,
// so you could catch it if you don't want the program to crash.
//
// Also, it is a good idea to compile any code that uses TMV with 
// TMV_DEBUG defined.  (defined below in this file)
//
// This will turn on some basic debugging in the code and will help 
// you catch programming mistakes.
// For example, if you have two vectors: 
// v1 of length 5 and v2  of lengths 10,
// then v1 = v2 will only produce a runtime error if TMV_DEBUG is 
// turned on.  Otherwise some random data will be overwritten and 
// who knows what will happen then.
//
// Once you know your code is working properly, you can recompile without
// the TMV_DEBUG flag to speed up the code.
// (And actually, the slow down isn't much, so it is probably
// a good idea to just leave it on to help diagnose any bug that 
// might turn up.)
//
// A comment here about the letters used to identify each type of Matrix.
// The writeCompact routines for sparse Matrix types have an identifying
// letter(s) before the output.  The same letter (or combination) is 
// used in the names of the composite arithmetic types like SumDD, etc.
// Here is the list of letters used for reference:
//
// (* indicates that the type is not yet implemented.)
//
//  X   Scalar 
//  V   Vector
//  M   Matrix
//  D   DiagMatrix
//  U   UpperTriMatrix
//  L   LowerTriMatrix
//  B   BandMatrix
//  S   SymMatrix
//  H   HermMatrix
// sB   SymBandMatrix 
// hB   HermBandMatrix  
// *kD  BlockDiagMatrix
// *skD SymBlockDiagMatrix 
// *hkD HermBlockDiagMatrix 
// *W   SparseVector
// *Q   SparseMatrix
// *sQ  SymSparseMatrix
// *hQ  HermSparseMatrix
// *kQ  BlockSparseMatrix 
// *skQ SymBlockSparseMatrix 
// *hkQ HermBlockSparseMatrix 
//
// TODO:  BlockDiagMatrix and related varieties
// TODO:  SparseMatrix and related varieties


#ifndef TMV_Base_H
#define TMV_Base_H

#if (!defined(NDEBUG) && !defined(TMVNDEBUG))
#define TMV_DEBUG
#endif

#ifdef TMV_INLINE
#define TMV_NO_WARN
#define TMV_NO_INST_DOUBLE
#define TMV_NO_INST_FLOAT
#endif

#ifndef TMV_NO_WARN
#define TMV_WARN
#endif

// Default optimization level is 2 if not specified with a -D flag.
#ifndef TMV_OPT
#define TMV_OPT 2
#endif

#include <iosfwd>
#include <limits>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <algorithm>

#if defined(__SSE2__) || defined(__SSE__)
#include "xmmintrin.h"
#endif

#ifdef TMV_DEBUG
#include <typeinfo>
#include <iostream>
#endif

#include <typeinfo>

#ifdef TMV_MEM_TEST
#include "util/mmgr.h"
#endif

#ifndef TMV_NO_INST_DOUBLE
#define TMV_INST_DOUBLE
#endif

#ifndef TMV_NO_INST_FLOAT
#define TMV_INST_FLOAT
#endif

#ifndef TMV_NO_INST_COMPLEX
#define TMV_INST_COMPLEX
#endif

// These next three are not instantiated by default:

//#define TMV_INST_MIX

//#define TMV_INST_LONGDOUBLE

//#define TMV_INST_INT


namespace tmv {

    static inline std::string TMV_Version() { return "0.70"; }
#define TMV_MAJOR_VERSION 0
#define TMV_MINOR_VERSION 70
#define TMV_VERSION_AT_LEAST(major,minor) \
    ( (major > TMV_MAJOR_VERSION) || \
      (major == TMV_MAJOR_VERSION && minor >= TMV_MINOR_VERSION) )

    // StorageType defines the order to store the elements of a matrix
    // TODO: I haven't implemented RowPacked and ColPakced yet.
    enum StorageType { 
        RowMajor, ColMajor, DiagMajor, NoMajor, RowPacked, ColPacked };

    // IndexStyle defines which kind of indexing to use for a vector
    enum IndexStyle { CStyle, FortranStyle=789234 };

    // UNKNOWN is the value of _size, _step, etc. whenever it
    // is not known at compile time.
    // We use for this value the maximally negative int.
    // In binary, this is a 1 followed by all zeros.
    const int UNKNOWN = (1<<(sizeof(int)*8-1));

    enum DivType { 
        XXX=1, LU=2, CH=4, QR=8, QRP=16, SV=32,
        // We store the divtype in a binary field integer.
        // In addition to the above, we also use the same object to 
        // store the following other flags related to division.
        // So these values must not clash with the above DivType values.
        // These aren't technically DivType's but since they are 
        // stored together, I think this adds to type-safety.
        DivInPlaceFlag = 64,
        SaveDivFlag = 128,
        // And finally shorthand for "one of the real DivType's":
        DivTypeFlags = 62
    };
    // I use things like &, |, |= to manipulate DivType's.  These are legal
    // in C, but not C++, so write overrides for these functions:
    static inline DivType operator|(DivType a, DivType b) 
    { return static_cast<DivType>(static_cast<int>(a) | static_cast<int>(b)); }
    static inline DivType operator&(DivType a, DivType b) 
    { return static_cast<DivType>(static_cast<int>(a) & static_cast<int>(b)); }
    static inline DivType& operator|=(DivType& a, DivType b) 
    { a = (a|b); return a; }
    static inline DivType& operator&=(DivType& a, DivType b) 
    { a = (a&b); return a; }
    static inline DivType operator~(DivType a) 
    { return static_cast<DivType>(~static_cast<int>(a)); }

    // TODO: I haven't implemented ZeroDiag yet.
    enum DiagType { UnitDiag, NonUnitDiag, ZeroDiag, UnknownDiag };

    template <class T>
    struct InstType
    { enum { inst = false }; };
#ifdef TMV_INST_INT
    template <>
    struct InstType<int>
    { enum { inst = true }; };
#endif
#ifdef TMV_INST_FLOAT
    template <>
    struct InstType<float>
    { enum { inst = true }; };
#endif
#ifdef TMV_INST_DOUBLE
    template <>
    struct InstType<double>
    { enum { inst = true }; };
#endif
#ifdef TMV_INST_LONGDOUBLE
    template <>
    struct InstType<double>
    { enum { inst = true }; };
#endif

    template <class T> class ConjRef;
    template <class T, bool C> class TriRef;

    // This helper class acts as a ? : operator for a typedef
    template <bool first, class T1, class T2>
    struct TypeSelect // first = true
    { typedef T1 type; };
    template <class T1, class T2>
    struct TypeSelect<false,T1,T2> // first = false
    { typedef T2 type; };
 
    template <class T>
    struct Traits
    {
        enum { isreal = true };
        enum { iscomplex = false };
        enum { isinst = InstType<T>::inst };
        enum { isinteger = std::numeric_limits<T>::is_integer };
        enum { twoifcomplex = 1 };

        typedef T real_type;
        typedef std::complex<T> complex_type;
        typedef typename TypeSelect<isinteger,double,T>::type float_type;

        typedef T& conj_reference;

        static const T& convert(const T& x) { return x; }
        template <class T2>
        static T convert(const T2& x) { return T(x); }
    };

    template <class T>
    struct Traits<std::complex<T> >
    {
        enum { isreal = false };
        enum { iscomplex = true };
        enum { isinst = InstType<T>::inst };
        enum { isinteger = Traits<T>::isinteger };
        enum { twoifcomplex = 2 };

        typedef T real_type;
        typedef std::complex<T> complex_type;
        typedef std::complex<typename Traits<T>::float_type> float_type;

        typedef ConjRef<std::complex<T> >& conj_reference;

        static const std::complex<T>& convert(const std::complex<T>& x) 
        { return x; }
        template <class T2>
        static std::complex<T> convert(const T2& x) 
        { return std::complex<T>(x); }
        // This last one is the real reason to have a convert function.
        // Because std::complex doesn't have a constructor from a 
        // std::complex of a different base type.
        // So this makes it easier to do the conversion.
        template <class T2>
        static std::complex<T> convert(const std::complex<T2>& x) 
        { return std::complex<T>(x.real(),x.imag()); }
    };

    template <class T>
    struct Traits<T&> : public Traits<T> {};
    template <class T>
    struct Traits<const T> : public Traits<T> {};
    template <class T>
    struct Traits<const T&> : public Traits<T> {};

    template <class T1, class T2>
    struct Traits2
    {
        enum { sametype = false };  // Same type
        enum { samebase = false };  // Same type ignoring complex-ness
        typedef T1 type;            // type of product, sum, etc.
    };
    template <class T1, class T2>
    struct Traits2<T1,std::complex<T2> >
    {
        enum { sametype = false };
        enum { samebase = false };
        typedef std::complex<typename Traits2<T1,T2>::type> type; 
    };
    template <class T1, class T2>
    struct Traits2<std::complex<T1>,T2>
    {
        enum { sametype = false };
        enum { samebase = false };
        typedef std::complex<typename Traits2<T1,T2>::type> type; 
    };
    template <class T1, class T2>
    struct Traits2<std::complex<T1>,std::complex<T2> >
    {
        enum { sametype = false };
        enum { samebase = false };
        typedef std::complex<typename Traits2<T1,T2>::type> type; 
    };
    // Specialize when both are the same base type
    template <class T>
    struct Traits2<T,T>
    {
        enum { sametype = true }; 
        enum { samebase = true }; 
        typedef T type;
    };
    template <class T>
    struct Traits2<std::complex<T>,T>
    {
        enum { sametype = false }; 
        enum { samebase = true }; 
        typedef std::complex<T> type;
    };
    template <class T>
    struct Traits2<T,std::complex<T> >
    {
        enum { sametype = false }; 
        enum { samebase = true }; 
        typedef std::complex<T> type;
    };
    template <class T>
    struct Traits2<std::complex<T>,std::complex<T> >
    {
        enum { sametype = true }; 
        enum { samebase = true }; 
        typedef std::complex<T> type;
    };
    // Specialize the real pairs for which the second value is the type to
    // use for sums and products rather than the first.
    template <>
    struct Traits2<int,float>
    {
        enum { sametype = false }; 
        enum { samebase = false }; 
        typedef float type;
    };
    template <>
    struct Traits2<int,double>
    {
        enum { sametype = false }; 
        enum { samebase = false }; 
        typedef double type;
    };
    template <>
    struct Traits2<int,long double>
    {
        enum { sametype = false }; 
        enum { samebase = false }; 
        typedef long double type;
    };
    template <>
    struct Traits2<float,double>
    {
        enum { sametype = false }; 
        enum { samebase = false }; 
        typedef double type;
    };
    template <>
    struct Traits2<float,long double>
    {
        enum { sametype = false }; 
        enum { samebase = false }; 
        typedef long double type;
    };
    template <>
    struct Traits2<double,long double>
    {
        enum { sametype = false }; 
        enum { samebase = false }; 
        typedef long double type;
    };

    template <class M1, class M2>
    struct ProdType
    { 
        typedef typename Traits<M1>::value_type T1;
        typedef typename Traits<M2>::value_type T2;
        typedef typename Traits2<T1,T2>::type type;
    };

    // A helper structure that conjugates a value if necessary:
    template <bool C, class T>
    struct DoConj_Helper // C = false
    { static T apply(const T& x) { return x; } };
    template <class T>
    struct DoConj_Helper<true,std::complex<T> >
    {
        static std::complex<T> apply(std::complex<T> x) 
        { return std::conj(x); } 
    };
    template <bool C, class T>
    static T DoConj(const T& x) { return DoConj_Helper<C,T>::apply(x); }

    template <int S>
    struct IntTraits
    {
        enum { negS = -S };
        enum { twoS = S<<1 };
        enum { halfS = S>>1 };
        enum { Sm1 = S-1 };
        enum { Sp1 = S+1 };
        enum { Sp2 = S+2 };
        enum { ispowerof2 = S > 0 && !(S & (S-1)) };
        enum { log = (
                S == 1 ? 0 : 
                S == 2 ? 1 : 
                S == 4 ? 2 : 
                S == 8 ? 3 : 
                S == 16 ? 4 : 
                S == 32 ? 5 : 
                S == 64 ? 6 : 
                S == 128 ? 7 : 
                S == 256 ? 8 : 
                S == 512 ? 9 : 
                S == 1024 ? 10 : 
                S == 2048 ? 11 : 
                S == 4096 ? 12 :  
                // Probably this is high enough for any conceivable use.
                // I only use up to 64 right now, but I provided a few
                // extra values to be safe.
                UNKNOWN ) 
        };
        enum { half_roundup = (
                // For very large S, just call it UNKNOWN to keep from 
                // having big complicated recursive structures.
                S > 128 ? UNKNOWN :
                S > 16 ? ((((S-1)>>5)+1)<<4) :
                (S>>1)  )
        };
        static int text() { return S; }
    };
    template <>
    struct IntTraits<UNKNOWN>
    {
        enum { negS = UNKNOWN };
        enum { twoS = UNKNOWN };
        enum { halfS = UNKNOWN };
        enum { Sm1 = UNKNOWN };
        enum { Sp1 = UNKNOWN };
        enum { Sp2 = UNKNOWN };
        enum { ispowerof2 = false };
        enum { log = UNKNOWN };
        enum { half_roundup = UNKNOWN };
        static const char* text() { return "UNKNOWN"; }
    };

    template <int S1, int S2, bool safe>
    struct SafeIntTraits2;

    template <int S1, int S2>
    struct IntTraits2
    {
        enum { sum = S1 + S2 };
        enum { diff = S1 - S2 };
        enum { prod = S1 * S2 };
        enum { safeprod = SafeIntTraits2<S1,S2,(S1<300 && S2<300)>::prod };
        enum { min = S1 < S2 ? S1 : S2 };
        enum { max = S1 > S2 ? S1 : S2 };
    };
    template <int S1>
    struct IntTraits2<S1,UNKNOWN>
    {
        enum { sum = UNKNOWN };
        enum { diff = UNKNOWN };
        enum { prod = UNKNOWN };
        enum { safeprod = UNKNOWN };
        enum { min = UNKNOWN };
        enum { max = UNKNOWN };
    };
    template <int S2>
    struct IntTraits2<UNKNOWN,S2>
    {
        enum { sum = UNKNOWN };
        enum { diff = UNKNOWN };
        enum { prod = UNKNOWN };
        enum { safeprod = UNKNOWN };
        enum { min = UNKNOWN };
        enum { max = UNKNOWN };
    };
    template <>
    struct IntTraits2<UNKNOWN,UNKNOWN>
    {
        enum { sum = UNKNOWN };
        enum { diff = UNKNOWN };
        enum { prod = UNKNOWN };
        enum { safeprod = UNKNOWN };
        enum { min = UNKNOWN };
        enum { max = UNKNOWN };
    };

    template <int S1, int S2>
    struct SafeIntTraits2<S1,S2,true>
    {
        enum { prod = IntTraits2<S1,S2>::prod };
    };
    template <int S1, int S2>
    struct SafeIntTraits2<S1,S2,false>
    {
        enum { prod = UNKNOWN };
    };

#if 1
    // This is based on the BOOST implementation BOOST_STATIC_ASSERT:
    template<bool> 
    struct tmv_static_assert;
    template<> 
    struct tmv_static_assert<true>
    { static void call() {} };
#define TMVStaticAssert(e) \
    do { tmv_static_assert<(e) != 0>::call(); } while (false)
#else
    // This is copied from the CERT secure coding practices version:
#define TMV_JOIN(x, y) TMV_JOIN_AGAIN(x, y)
#define TMV_JOIN_AGAIN(x, y) x ## y

#define TMVStaticAssert(e) \
    typedef char TMV_JOIN(assertion_failed_at_line_, __LINE__) [(e) ? 1 : -1]
#endif

    template <class T> 
    static T TMV_SQR(const T& x) 
    { return x*x; }

    template <class T> 
    static typename Traits<T>::float_type TMV_SQRT(const T& x) 
    { return std::sqrt(typename Traits<T>::float_type(x)); }

    template <class T> 
    static typename Traits<T>::float_type TMV_EXP(const T& x) 
    { return std::exp(typename Traits<T>::float_type(x)); }

    template <class T> 
    static typename Traits<T>::float_type TMV_LOG(const T& x) 
    { return std::log(typename Traits<T>::float_type(x)); }

    template <class T> 
    static T TMV_NORM(const T& x) 
    { return x*x; }

    template <class T> 
    static T TMV_NORM(const std::complex<T>& x) 
    { return std::norm(x); }

    template <class T> 
    static T TMV_CONJ(const T& x)
    { return x; }

    template <class T> 
    static std::complex<T> TMV_CONJ(const std::complex<T>& x)
    { return std::conj(x); }

    template <class T> 
    static T TMV_REAL(const T& x)
    { return x; }

    template <class T> 
    static T TMV_REAL(const std::complex<T>& x)
    { return x.real(); }

    template <class T> 
    static T TMV_IMAG(const T& )
    { return T(0); }

    template <class T> 
    static T TMV_IMAG(const std::complex<T>& x)
    { return x.imag(); }

    // Many implemenations of complex allow x.real() and x.imag() to return 
    // references to the actual locations of these elements.
    // However, the standard technically says that these return by 
    // value, not reference.  So we need to use a reinterpret_cast to 
    // really make sure we get the right thing.
    template <class T> 
    static T& TMV_REAL_PART(std::complex<T>& x)
    { return reinterpret_cast<T&>(x); }

    template <class T> 
    static T& TMV_REAL_PART(T& x)
    { return x; }

    template <class T> 
    static T& TMV_IMAG_PART(std::complex<T>& x)
    { return *(reinterpret_cast<T*>(&x)+1); }

    template <class T> 
    static T TMV_ARG(const T& x)
    { return x >= T(0) ? T(1) : T(-1); }

    template <class T> 
    static typename Traits<T>::float_type TMV_ARG(const std::complex<T>& x)
    {
        typedef typename Traits<T>::float_type FT;
        return arg(Traits<std::complex<FT> >::convert(x)); 
    }

    template <class T> 
    static T TMV_ABS(const T& x)
    { return std::abs(x); }

    template <class T> 
    static typename Traits<T>::float_type TMV_ABS(const std::complex<T>& x)
    { 
        typedef typename Traits<T>::float_type FT;
        // This is the same as the usual std::abs algorithm.
        // However, I have come across implementations that don't 
        // protext against overflow, so I just duplicate it here.

        FT xr = x.real();
        FT xi = x.imag();
        const FT s = std::max(std::abs(xr),std::abs(xi));
        if (s == T(0)) return s; // Check for s == 0
        xr /= s;
        xi /= s;
        return s * std::sqrt(xr*xr + xi*xi);
    }

    template <class T> 
    static T TMV_ABS2(const T& x)
    { return std::abs(x); }

    template <class T> 
    static T TMV_ABS2(const std::complex<T>& x)
    { return std::abs(x.real()) + std::abs(x.imag()); }

    template <class T> 
    static T TMV_SIGN(const T& x, const T& )
    { return x > 0 ? T(1) : T(-1); }

    template <class T> 
    static std::complex<T> TMV_SIGN(const std::complex<T>& x, const T& absx)
    { return absx > 0 ? x/absx : std::complex<T>(1); }

    template <class T> 
    static T TMV_MIN(const T& x, const T& y)
    { return x > y ? y : x; }

    template <class T> 
    static T TMV_MAX(const T& x, const T& y)
    { return x > y ? x : y; }

    template <class T> 
    static void TMV_SWAP(T& x, T& y)
    { T z = x; x = y; y = z; }

  // Useful for debugging SSE routines:
#ifdef __SSE__
    static inline std::ostream& operator<<(std::ostream& os, const __m128& xf)
    {
        os << &xf<<": ";
        float ff[4];
        _mm_storeu_ps(ff,xf);
        os << ff[0]<<","<<ff[1]<<","<<ff[2]<<","<<ff[3]<<" ";
        return os;
    }
#endif

#ifdef __SSE2__
    static inline std::ostream& operator<<(std::ostream& os, const __m128d& xd)
    {
        os << &xd<<": ";
        double dd[2];
        _mm_storeu_pd(dd,xd);
        os << dd[0]<<","<<dd[1]<<" ";
        return os;
    }
#endif

    template <int ix, class T> 
    struct Scaling;

    template <bool conj1, bool conj2> 
    struct ZProd;

    struct ZSum;

   // A helper struct to pick one of two possibile behaviors
    // according to a template bool argument.  
    // It is defined below ZProd.
    template <bool yn> 
    struct Maybe;

    // The way the STL (typically) implements std::complex means
    // that the exporession z = x * y when x,y,z are complex
    // is not very amenable to optimization by most compilers.
    // It turns out that doing the real and the imaginary parts 
    // separately usually presents the calculation in a form that 
    // the compiler can optimize, so we do this using a helper
    // structure ZProd.

    // First another helper class to do the correct complex multiplication
    // depending on the conj status of x and y:
    template <bool conj1, bool conj2> 
    struct ZProd2;

    template <> 
    struct ZProd2<false,false>
    {
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type rprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*ry - ix*iy); }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type iprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*iy + ix*ry); }
    };
    template <> 
    struct ZProd2<true,false>
    {
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type rprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*ry + ix*iy); }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type iprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*iy - ix*ry); }
    };
    template <> 
    struct ZProd2<false,true>
    {
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type rprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*ry + ix*iy); }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type iprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (ry*ix - iy*rx); }
    };
    template <> 
    struct ZProd2<true,true>
    {
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type rprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*ry - ix*iy); }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type iprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return -(rx*iy + ix*ry); }
    };

    // Now the real ZProd struct that does the right thing for each real
    // or complex possibility, including the option of x being a Scaling.
    template <bool conj1, bool conj2> 
    struct ZProd
    {
        // complex * complex
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type rprod(
            const std::complex<T1>& x, const std::complex<T2>& y)
        { 
            return ZProd2<conj1,conj2>::rprod(
                x.real(),x.imag(),y.real(),y.imag()); 
        }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type iprod(
            const std::complex<T1>& x, const std::complex<T2>& y)
        {
            return ZProd2<conj1,conj2>::iprod(
                x.real(),x.imag(),y.real(),y.imag()); 
        }

        // real * complex
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type rprod(
            const T1& x, const std::complex<T2>& y)
        { return x*y.real(); }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type iprod(
            const T1& x, const std::complex<T2>& y)
        { return Maybe<conj2>::neg(x*y.imag()); }

        // complex * real
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type rprod(
            const std::complex<T1>& x, const T2& y)
        { return y*x.real(); }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type iprod(
            const std::complex<T1>& x, const T2& y)
        { return Maybe<conj1>::neg(y*x.imag()); }

        // real * real
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type rprod(
            const T1& x, const T2& y)
        { return x*y; }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type iprod(
            const T1& x, const T2& y)
        { return typename Traits2<T1,T2>::type(0); }

        // scaling * real
        template <class T1, class T2> 
        static typename Traits2<typename Traits<T1>::real_type,T2>::type rprod(
            const Scaling<0,T1>& x, const T2& y)
        { return rprod(x.x,y); }
        template <class T1, class T2> 
        static typename Traits2<typename Traits<T1>::real_type,T2>::type iprod(
            const Scaling<0,T1>& x, const T2& y)
        { return iprod(x.x,y); }
        template <class T1, class T2> 
        static T2 rprod(const Scaling<1,T1>& x, const T2& y)
        { return y; }
        template <class T1, class T2> 
        static T2 iprod(const Scaling<1,T1>& x, const T2& y)
        { return T2(0); }
        template <class T1, class T2> 
        static T2 rprod(const Scaling<-1,T1>& x, const T2& y)
        { return -y; }
        template <class T1, class T2> 
        static T2 iprod(const Scaling<-1,T1>& x, const T2& y)
        { return T2(0); }

        // scaling * complex
        template <class T1, class T2> 
        static typename Traits2<typename Traits<T1>::real_type,T2>::type rprod(
            const Scaling<0,T1>& x, const std::complex<T2>& y)
        { return rprod(x.x,y); }
        template <class T1, class T2> 
        static typename Traits2<typename Traits<T1>::real_type,T2>::type iprod(
            const Scaling<0,T1>& x, const std::complex<T2>& y)
        { return iprod(x.x,y); }
        template <class T1, class T2> 
        static T2 rprod(
            const Scaling<1,T1>& x, const std::complex<T2>& y)
        { return y.real(); }
        template <class T1, class T2> 
        static T2 iprod(
            const Scaling<1,T1>& x, const std::complex<T2>& y)
        { return Maybe<conj2>::neg(y.imag()); }
        template <class T1, class T2> 
        static T2 rprod(
            const Scaling<-1,T1>& x, const std::complex<T2>& y)
        { return -y.real(); }
        template <class T1, class T2> 
        static T2 iprod(
            const Scaling<-1,T1>& x, const std::complex<T2>& y)
        { return Maybe<!conj2>::neg(y.imag()); }

        // full prod
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type prod(
            const T1& x, const T2& y)
        {
            typedef typename Traits2<T1,T2>::type PT;
            const bool iscomplex = Traits<PT>::iscomplex;
            return makeprod<iscomplex,T1,T2>::call(x,y);
        }
        template <class T1, class T2> 
        static T2 prod(const Scaling<1,T1>& x, const T2& y)
        { return Maybe<conj2>::conj(y); }
        template <int ix, class T1, class T2> 
        static T1 prod(const T1& x, const Scaling<1,T2>& y)
        { return Maybe<conj1>::conj(x); }
        template <class T1, class T2> 
        static T2 prod(const Scaling<-1,T1>& x, const T2& y)
        { return -Maybe<conj2>::conj(y); }
        template <int ix, class T1, class T2> 
        static T1 prod(const T1& x, const Scaling<-1,T2>& y)
        { return -Maybe<conj1>::conj(x); }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type prod(
            const Scaling<0,T1>& x, const T2& y)
        { return prod(x.x , y); }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type prod(
            const T1& x, const Scaling<0,T2>& y)
        { return prod(x , y.x); }

        template <bool iscomplex, class T1, class T2> 
        struct makeprod;
        template <class T1, class T2>
        struct makeprod<true,T1,T2>
        {
            typedef typename Traits2<T1,T2>::type PT;
            static PT call(const T1& x, const T2& y)
            { return PT(rprod(x,y),iprod(x,y)); }
        };
        template <class T1, class T2>
        struct makeprod<false,T1,T2>
        {
            typedef typename Traits2<T1,T2>::type PT;
            static PT call(const T1& x, const T2& y)
            { return x*y; }
        };

        // quot
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type quot(
            const T1& x, const T2& y)
        {
            const bool c2 = Traits<T2>::iscomplex;
            return makequot<c2,T1,T2>::call(x,y);
        }
        template <int ix, class T1, class T2> 
        static typename Traits2<T1,T2>::type quot(
            const Scaling<ix,T1>& x, const T2& y)
        { return quot(T1(x),y); }
        template <int ix, class T1, class T2> 
        static typename Traits2<T1,T2>::type quot(
            const T1& x, const Scaling<ix,T2>& y)
        { return prod(x,y); }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type quot(
            const T1& x, const Scaling<0,T2>& y)
        { return quot(x , y.x); }

        template <bool c2, class T1, class T2> 
        struct makequot;
        template <class T1, class T2>
        struct makequot<true,T1,T2>
        {
            typedef typename Traits2<T1,T2>::type PT;
            typedef typename Traits<PT>::real_type RT;
            static PT call(const T1& x, const T2& y)
            { return ZProd<conj1,!conj2>::prod(x,y) / std::norm(y); }
        };
        template <class T1, class T2>
        struct makequot<false,T1,T2>
        {
            typedef typename Traits2<T1,T2>::type PT;
            static PT call(const T1& x, const T2& y)
            { return Maybe<conj1>::conj(x/y); }
        };

    };

    template <> 
    struct Maybe<true>
    {
        // real_type or T
        template <class T>
        struct RealType { typedef typename Traits<T>::real_type type; };

        // Type of T1*T2 or T2
        template <class T1, class T2>
        struct ProdType { typedef typename Traits2<T1,T2>::type type; };

        // x or y
        template <class T>
        static T select(const T& x, const T& /*y*/) { return x; }

        // abs(x) or x
        template <class T>
        static typename Traits<T>::real_type abs(const T& x) 
        { return TMV_ABS(x); }

        // real(x) or x
        template <class T>
        static typename Traits<T>::real_type real(const T& x) 
        { return TMV_REAL(x); }

        // conj(x) or x
        template <class T>
        static T conj(const T& x) { return TMV_CONJ(x); }

        // x = conj(x) or nothing
        template <class T>
        static void conjval(T& x) { x = TMV_CONJ(x); }

        // -x or x
        template <class T>
        static T neg(const T& x) { return -x; }

        // x<y or x>y
        template <class T>
        static bool less(const T& x, const T& y) { return x<y; }

        // 2*x or x
        template <class T>
        static T twox(const T& x)
        { return 2*x; }
        static int twox(const int& x)
        { return x>>1; }

        // ++x or nothing
        template <class T>
        static void increment(T& x) { ++x; }

        // --x or nothing
        template <class T>
        static void decrement(T& x) { --x; }

        // 2*x or x  (intended for integer types)
        template <class T>
        static T multby2(const T& x) { return x<<1; }
        // x/2 or x  (intended for integer types)
        template <class T>
        static T divby2(const T& x) { return x>>1; }

        // x += y or x = y
        template <class T1, class T2>
        static void add(T1& x, const T2& y) { x += y; }
        template <class T1, class T2>
        static void add(ConjRef<T1> x, const T2& y) { x += y; }
        template <class T1, bool C, class T2>
        static void add(TriRef<T1,C> x, const T2& y) { x += y; }

        // x *= y or nothing
        template <class T1, class T2>
        static void scale(T1& x, const T2& y) { x *= y; }
        template <class T1, class T2>
        static void scale(ConjRef<T1> x, const T2& y) { x *= y; }
        template <class T1, bool C, class T2>
        static void scale(TriRef<T1,C> x, const T2& y) { x *= y; }

        // x /= y or nothing
        template <class T1, class T2>
        static void invscale(T1& x, const T2& y) { x /= y; }
        template <class T1, class T2>
        static void invscale(ConjRef<T1> x, const T2& y) { x /= y; }
        template <class T1, bool C, class T2>
        static void invscale(TriRef<T1,C> x, const T2& y) { x /= y; }

        // x + y or y
        template <class T1, class T2>
        static typename Traits2<T1,T2>::type sum(
            const T1& x, const T2& y) 
        { return x + y; }

        template <class T1, int ix2, class T2>
        static typename Traits2<T1,T2>::type sum(
            const T1& x, const Scaling<ix2,T2>& y) 
        { return x + T2(y); }

        // x * y or y
        template <class T1, class T2>
        static typename Traits2<T1,T2>::type prod(
            const T1& x, const T2& y) 
        { return ZProd<false,false>::prod(x , y); }

        template <class T1, int ix2, class T2>
        static typename Traits2<T1,T2>::type prod(
            const T1& x, const Scaling<ix2,T2>& y) 
        { return ZProd<false,false>::prod(x , y); }

        // x^-1 * y or y
        template <class T1, class T2>
        static typename Traits2<T1,T2>::type invprod(
            const T1& x, const T2& y) 
        { return ZProd<false,false>::quot(y , x); }

        template <class T1, int ix2, class T2>
        static typename Traits2<T1,T2>::type invprod(
            const T1& x, const Scaling<ix2,T2>& y) 
        { return ZProd<false,false>::quot(y , x); }

        // x * y or y
        template <bool c1, bool c2, class T1, class T2>
        static typename Traits2<T1,T2>::type zprod(
            const T1& x, const T2& y)
        { return ZProd<c1,c2>::prod(x , y); }

        // x - y or x + y
        template <class T>
        static T diff(const T& x, const T& y) { return x - y; }

        // x = y or nothing
        template <class T1, class T2>
        static void set(T1& x, const T2& y) { x = y; }
        template <class T1, class T2>
        static void set(ConjRef<T1> x, const T2& y) { x = y; }
        template <class T1, bool C, class T2>
        static void set(TriRef<T1,C> x, const T2& y) { x = y; }

        // m.ref(i,i) = x or nothing
        template <class M, class T>
        static void setdiag(M& m, int i, const T& x) { m.ref(i,i) = x; }
        template <class M, class T>
        static void setdiag2(M m, int i, const T& x) { m.ref(i,i) = x; }

        // m.diag() = x or nothing
        template <class M, class T>
        static void setdiag(M& m, const T& x) { m.diag().setAllTo(x); }
        template <class M, class T>
        static void setdiag2(M m, const T& x) { m.diag().setAllTo(x); }

        // m.setZero() or nothing
        template <class M>
        static void zero(M& m) { m.setZero(); }
        template <class M>
        static void zero2(M m) { m.setZero(); }

        // m.conjugateSelf() or nothing
        template <class M>
        static void conjself(M& m) { m.conjugateSelf(); }
        template <class M>
        static void conjself2(M m) { m.conjugateSelf(); }

        // m.conjugate() or m
        template <class M>
        static typename M::const_conjugate_type conjugate(const M& m) 
        { return m.conjugate(); }
        template <class M>
        static typename M::conjugate_type conjugate(M& m) 
        { return m.conjugate(); }

        // m.transpose() or m
        template <class M>
        static typename M::const_transpose_type transpose(const M& m) 
        { return m.transpose(); }
        template <class M>
        static typename M::transpose_type transpose(M& m) 
        { return m.transpose(); }

        // v.addToAll(x) or v.setAllTo(x)
        template <class V, class T>
        static void addtoall(V& v, const T& x) { v.addToAll(x); }
        template <class V, class T>
        static void addtoall2(V v, const T& x) { v.addToAll(x); }

        // m.diag().setAllTo(1) or nothing
        template <class M>
        static void makeunitdiag(M& m) 
        { m.diag().setAllTo(typename M::value_type(1)); }
        template <class M>
        static void makeunitdiag2(M m) 
        { m.diag().setAllTo(typename M::value_type(1)); }

        // m.offDiag().setZero() or nothing
        template <class M>
        static void zero_offdiag(M& m) { m.offDiag().setZero(); }
        template <class M>
        static void zero_offdiag2(M m) { m.offDiag().setZero(); }

        // m.offDiag() or m.view()
        template <class M>
        static typename M::offdiag_type offdiag(M& m) 
        { return m.offDiag(); }
        template <class M>
        static typename M::offdiag_type offdiag2(M m) 
        { return m.offDiag(); }

        // m2.viewAsUnitDiag() or m2
        template <class M>
        static typename M::const_unitdiag_type unitview(const M& m) 
        { return m.viewAsUnitDiag(); }
        template <class M>
        static typename M::unitdiag_type unitview(M& m) 
        { return m.viewAsUnitDiag(); }

        // m2.UnitUpperTri() or m2.UpperTri()
        template <class M>
        static typename M::const_unit_uppertri_type unit_uppertri(
            const M& m) 
        { return m.unitUpperTri(); }
        template <class M>
        static typename M::unit_uppertri_type unit_uppertri(M& m) 
        { return m.unitUpperTri(); }
        template <class M>
        static typename M::unit_uppertri_type unit_uppertri2(M m) 
        { return m.unitUpperTri(); }

        // m2.UnitLowerTri() or m2.LowerTri()
        template <class M>
        static typename M::const_unit_lowertri_type unit_lowertri(
            const M& m) 
        { return m.unitLowerTri(); }
        template <class M>
        static typename M::unit_lowertri_type unit_lowertri(M& m) 
        { return m.unitLowerTri(); }
        template <class M>
        static typename M::unit_lowertri_type unit_lowertri2(M m) 
        { return m.unitLowerTri(); }

        // m2.UpperTri() or m2.LowerTri()
        template <class M>
        static typename M::const_uppertri_type uppertri(const M& m) 
        { return m.upperTri(); }
        template <class M>
        static typename M::uppertri_type uppertri(M& m) 
        { return m.upperTri(); }
        template <class M>
        static typename M::uppertri_type uppertri2(M m) 
        { return m.upperTri(); }

#ifdef __SSE__
        // _mm_load_ps or _mm_set_ps
        static void sse_load(
            __m128& m,
            const float* x, const float*, const float*, const float*)
        { m = _mm_load_ps(x); }

        // _mm_loadu_ps or _mm_set_ps
        static void sse_loadu(
            __m128& m,
            const float* x, const float*, const float*, const float*)
        { m = _mm_loadu_ps(x); }

        // _mm_store_ps or cast and assign
        static void sse_store(
            float* x, float*, float*, float*, const __m128& m)
        { _mm_store_ps(x,m); }

        // _mm_storeu_ps or cast and assign
        static void sse_storeu(
            float* x, float*, float*, float*, const __m128& m)
        { _mm_storeu_ps(x,m); }

        // _mm_add_ps or cast and add
        static void sse_add(
            float* x, float*, float*, float*, const __m128& m)
        { _mm_store_ps(x,_mm_add_ps(_mm_load_ps(x),m)); }
        static void sse_addu(
            float* x, float*, float*, float*, const __m128& m)
        { _mm_storeu_ps(x,_mm_add_ps(_mm_loadu_ps(x),m)); }

        // _mm_load_ps or two _mm_load_pi's
        static void sse_load(
            __m128& m,
            const std::complex<float>* x, const std::complex<float>* )
        { m = _mm_load_ps( (const float*) x ); }

        // _mm_loadu_ps or two _mm_load?_pi's
        static void sse_loadu(
            __m128& m,
            const std::complex<float>* x, const std::complex<float>* )
        { m = _mm_loadu_ps( (const float*) x ); }

        // _mm_store_ps or two _mm_store?_pi's
        static void sse_store(
            std::complex<float>* x, std::complex<float>* , const __m128& m)
        { _mm_store_ps((float*)(x) , m); }

        // _mm_storeu_ps or two _mm_store?_pi's
        static void sse_storeu(
            std::complex<float>* x, std::complex<float>* , const __m128& m)
        { _mm_storeu_ps((float*)(x) , m); }

        // _mm_add_ps or cast and add
        static void sse_add(
            std::complex<float>* x, std::complex<float>* , const __m128& m)
        { _mm_store_ps((float*)x,_mm_add_ps(_mm_load_ps((float*)x),m)); }
        static void sse_addu(
            std::complex<float>* x, std::complex<float>* , const __m128& m)
        { _mm_storeu_ps((float*)x,_mm_add_ps(_mm_loadu_ps((float*)x),m)); }

        // A = _mm_mul_ps(x,A) or nothing
        static void sse_mult(const __m128& x, __m128& A)
        { A = _mm_mul_ps(x,A); }
#endif
#ifdef __SSE2__
        static void sse_load(
            __m128d& m, const double* x, const double* )
        { m = _mm_load_pd(x); }
        static void sse_loadu(
            __m128d& m, const double* x, const double* )
        { m = _mm_loadu_pd(x); }

        static void sse_store(double* x, double*, const __m128d& m)
        { _mm_store_pd(x,m); }
        static void sse_storeu(double* x, double*, const __m128d& m)
        { _mm_storeu_pd(x,m); }

        static void sse_add(double* x, double*, const __m128d& m)
        { _mm_store_pd(x,_mm_add_pd(_mm_load_pd(x),m)); }
        static void sse_addu(double* x, double*, const __m128d& m)
        { _mm_storeu_pd(x,_mm_add_pd(_mm_loadu_pd(x),m)); }

        static void sse_load(
            __m128d& m, const std::complex<double>* x)
        { m = _mm_load_pd((const double*) x); }
        static void sse_loadu(
            __m128d& m, const std::complex<double>* x)
        { m = _mm_loadu_pd((const double*) x); }

        static void sse_store(std::complex<double>* x, const __m128d& m)
        { _mm_store_pd((double*) x,m); }
        static void sse_storeu(std::complex<double>* x, const __m128d& m)
        { _mm_storeu_pd((double*) x,m); }

        static void sse_add(std::complex<double>* x, const __m128d& m)
        { _mm_store_pd((double*)x,_mm_add_pd(_mm_load_pd((double*)x),m)); }
        static void sse_addu(std::complex<double>* x, const __m128d& m)
        { _mm_storeu_pd((double*)x,_mm_add_pd(_mm_loadu_pd((double*)x),m)); }

        static void sse_mult(const __m128d& x, __m128d& A)
        { A = _mm_mul_pd(x,A); }
#endif
    };
    template <> 
    struct Maybe<false>
    {
        template <class T>
        struct RealType { typedef T type; };

        template <class T1, class T2>
        struct ProdType { typedef T2 type; };

        template <class T>
        static T select(const T& /*x*/, const T& y) { return y; }

        template <class T>
        static T abs(const T& x) { return x; }

        template <class T>
        static T real(const T& x) { return x; }

        template <class T>
        static T conj(const T& x) { return x; }

        template <class T>
        static void conjval(T& ) { }

        template <class T>
        static T neg(const T& x) { return x; }

        template <class T>
        static bool less(const T& x, const T& y) { return x>y; }

        template <class T>
        static T twox(const T& x) { return x; }

        template <class T>
        static void increment(const T& /*x*/) { }

        template <class T>
        static void decrement(const T& /*x*/) { }

        template <class T>
        static const T& multby2(const T& x) { return x; }
        template <class T>
        static const T& divby2(const T& x) { return x; }

        template <class T1, class T2>
        static void add(T1& x, const T2& y) { x = y; }
        template <class T1, class T2>
        static void add(ConjRef<T1> x, const T2& y) { x = y; }
        template <class T1, bool C, class T2>
        static void add(TriRef<T1,C> x, const T2& y) { x = y; }

        template <class T1, class T2>
        static void scale(T1& /*x*/, const T2& /*y*/) { }
        template <class T1, class T2>
        static void scale(ConjRef<T1> /*x*/, const T2& /*y*/) { }
        template <class T1, bool C, class T2>
        static void scale(TriRef<T1,C> /*x*/, const T2& /*y*/) { }

        template <class T1, class T2>
        static void invscale(T1& /*x*/, const T2& /*y*/) { }
        template <class T1, class T2>
        static void invscale(ConjRef<T1> /*x*/, const T2& /*y*/) { }
        template <class T1, bool C, class T2>
        static void invscale(TriRef<T1,C> /*x*/, const T2& /*y*/) { }

        template <class T1, class T2>
        static T2 sum(const T1& /*x*/, const T2& y) { return y; }

        template <class T1, class T2>
        static T2 prod(const T1& /*x*/, const T2& y) { return y; }

        template <class T1, int ix2, class T2>
        static Scaling<ix2,T2> prod(
            const T1& /*x*/, const Scaling<ix2,T2>& y) 
        { return y; }

        template <class T1, class T2>
        static T2 invprod(const T1& /*x*/, const T2& y)  { return y; }

        template <class T1, int ix2, class T2>
        static Scaling<ix2,T2> invprod(
            const T1& /*x*/, const Scaling<ix2,T2>& y) { return y; }

        template <bool c1, bool c2, class T1, class T2>
        static T2 zprod(const T1& /*x*/, const T2& y)
        { return Maybe<c2>::conj(y); }

        template <class T>
        static T diff(const T& x, const T& y) { return x + y; }

        template <class T1, class T2>
        static void set(T1& /*x*/, const T2& /*y*/) { }
        template <class T1, class T2>
        static void set(ConjRef<T1> /*x*/, const T2& /*y*/) { }
        template <class T1, bool C, class T2>
        static void set(TriRef<T1,C> /*x*/, const T2& /*y*/) { }

        template <class M, class T>
        static void setdiag(M& /*m*/, int /*i*/, const T& /*x*/) { }
        template <class M, class T>
        static void setdiag2(M /*m*/, int /*i*/, const T& /*x*/) { }

        template <class M, class T>
        static void setdiag(M& /*m*/, const T& /*x*/) { }
        template <class M, class T>
        static void setdiag2(M /*m*/, const T& /*x*/) { }

        template <class M>
        static void zero(M& /*m*/) { }
        template <class M>
        static void zero2(M /*m*/) { }

        template <class M>
        static void conjself(M& /*m*/) { }
        template <class M>
        static void conjself2(M /*m*/) { }

        template <class M>
        static const M& conjugate(const M& m) { return m; }
        template <class M>
        static M& conjugate(M& m) { return m; }

        template <class M>
        static const M& transpose(const M& m) { return m; }
        template <class M>
        static M& transpose(M& m) { return m; }

        template <class V, class T>
        static void addtoall(V& v, const T& x) { v.setAllTo(x); }
        template <class V, class T>
        static void addtoall2(V v, const T& x) { v.setAllTo(x); }

        template <class M>
        static void makeunitdiag(M& /*m*/) { }
        template <class M>
        static void makeunitdiag2(M /*m*/) { }

        template <class M1, class M2>
        static void offdiag_copy(const M1& m1, M2& m2) { m2 = m1; }
        template <class M1, class M2>
        static void offdiag_copy2(const M1 m1, M2& m2) { m2 = m1; }

        template <class M>
        static void zero_offdiag(M& /*m*/) { }
        template <class M>
        static void zero_offdiag2(M /*m*/) { }

        template <class M>
        static typename M::view_type offdiag(M& m) { return m.view(); }
        template <class M>
        static typename M::view_type offdiag2(M m) { return m.view(); }

        template <class M>
        static const M& unitview(const M& m) 
        { return m; }
        template <class M>
        static M& unitview(M& m) 
        { return m; }

        template <class M>
        static typename M::const_uppertri_type unit_uppertri(const M& m) 
        { return m.upperTri(); }
        template <class M>
        static typename M::uppertri_type unit_uppertri(M& m) 
        { return m.upperTri(); }
        template <class M>
        static typename M::uppertri_type unit_uppertri2(M m) 
        { return m.upperTri(); }

        template <class M>
        static typename M::const_lowertri_type unit_lowertri(const M& m) 
        { return m.lowerTri(); }
        template <class M>
        static typename M::lowertri_type unit_lowertri(M& m) 
        { return m.lowerTri(); }
        template <class M>
        static typename M::lowertri_type unit_lowertri2(M m) 
        { return m.lowerTri(); }

        template <class M>
        static typename M::const_lowertri_type uppertri(const M& m) 
        { return m.lowerTri(); }
        template <class M>
        static typename M::lowertri_type uppertri(M& m) 
        { return m.lowerTri(); }
        template <class M>
        static typename M::lowertri_type uppertri2(M m) 
        { return m.lowerTri(); }

#ifdef __SSE__
        static void sse_load(
            __m128& m,
            const float* x1, const float* x2, const float* x3, const float* x4)
        { m = _mm_set_ps(*x4,*x3,*x2,*x1); }
        static void sse_loadu(
            __m128& m,
            const float* x1, const float* x2, const float* x3, const float* x4)
        { m = _mm_set_ps(*x4,*x3,*x2,*x1); }

        static void sse_store(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 = mf[0]; *x2 = mf[1]; *x3 = mf[2]; *x4 = mf[3];
        }
        static void sse_storeu(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 = mf[0]; *x2 = mf[1]; *x3 = mf[2]; *x4 = mf[3];
        }

        static void sse_add(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 += mf[0]; *x2 += mf[1]; *x3 += mf[2]; *x4 += mf[3];
        }
        static void sse_addu(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 += mf[0]; *x2 += mf[1]; *x3 += mf[2]; *x4 += mf[3];
        }

        static void sse_load(
            __m128& m,
            const std::complex<float>* x1, const std::complex<float>* x2)
        {
            // The junk variable is to avoid a warning about m being
            // used uninitialized.
            const __m128 junk = _mm_set1_ps(0.F);
            m = _mm_loadl_pi(junk,(const __m64*)x1);
            m = _mm_loadh_pi(m,(const __m64*)x2);
        }
        static void sse_loadu(
            __m128& m,
            const std::complex<float>* x1, const std::complex<float>* x2)
        {
            const __m128 junk = _mm_set1_ps(0.F);
            m = _mm_loadl_pi(junk,(const __m64*)x1);
            m = _mm_loadh_pi(m,(const __m64*)x2);
        }

        static void sse_store(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        {
            _mm_storel_pi((__m64*)x1,m);
            _mm_storeh_pi((__m64*)x2,m);
        }
        static void sse_storeu(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        {
            _mm_storel_pi((__m64*)x1,m);
            _mm_storeh_pi((__m64*)x2,m);
        }

        static void sse_add(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 += std::complex<float>(mf[0],mf[1]);
            *x2 += std::complex<float>(mf[2],mf[3]);
        }
        static void sse_addu(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 += std::complex<float>(mf[0],mf[1]);
            *x2 += std::complex<float>(mf[2],mf[3]);
        }

        static void sse_mult(const __m128& , __m128& ) {}
#endif
#ifdef __SSE2__
        static void sse_load(
            __m128d& m, const double* x1, const double* x2)
        { m = _mm_set_pd(*x2,*x1); }
        static void sse_loadu(
            __m128d& m, const double* x1, const double* x2)
        { m = _mm_set_pd(*x2,*x1); }

        static void sse_store(double* x1, double* x2, const __m128d& m)
        {
            const double* md= (const double*)(&m);
            *x1 = md[0]; *x2 = md[1];
        }
        static void sse_storeu(double* x1, double* x2, const __m128d& m)
        {
            const double* mf= (const double*)(&m);
            *x1 = mf[0]; *x2 = mf[1];
        }

        static void sse_add(double* x1, double* x2, const __m128d& m)
        {
            const double* md= (const double*)(&m);
            *x1 += md[0]; *x2 += md[1];
        }
        static void sse_addu(double* x1, double* x2, const __m128d& m)
        {
            const double* mf= (const double*)(&m);
            *x1 += mf[0]; *x2 += mf[1];
        }

        // These are the same for Maybe<true> or Maybe<false>.
        // I include them to make it easier to write the code that uses
        // these loads and stores, since it makes the syntax more similar
        // between all the different varieties (float/double, real/complex).
        static void sse_load(
            __m128d& m, const std::complex<double>* x)
        { m = _mm_load_pd((const double*)x); }
        static void sse_loadu(
            __m128d& m, const std::complex<double>* x)
        { m = _mm_loadu_pd((const double*)x); }

        static void sse_store(std::complex<double>* x, const __m128d& m)
        { _mm_store_pd((double*)x,m); }
        static void sse_storeu(std::complex<double>* x, const __m128d& m)
        { _mm_storeu_pd((double*)x,m); }

        static void sse_add(std::complex<double>* x, const __m128d& m)
        { _mm_store_pd((double*)x,_mm_add_pd(_mm_load_pd((double*)x),m)); }
        static void sse_addu(std::complex<double>* x, const __m128d& m)
        { _mm_storeu_pd((double*)x,_mm_add_pd(_mm_loadu_pd((double*)x),m)); }

        static void sse_mult(const __m128d& , __m128d& ) {}
#endif
    };

    // A couple things benefit from having two compile time booleans
    template <bool yn1, bool yn2> 
    struct Maybe2;

    template <bool yn2> 
    struct Maybe2<true,yn2>
    {
        // Maybe<add>::add or nothing
        template <class T1, class T2>
        static void add(T1& x, const T2& y) 
        { Maybe<yn2>::add(x,y); }
        template <class T1, class T2>
        static void add(ConjRef<T1> x, const T2& y) 
        { Maybe<yn2>::add(x,y); }
        template <class T1, bool C, class T2>
        static void add(TriRef<T1,C> x, const T2& y) 
        { Maybe<yn2>::add(x,y); }

#ifdef __SSE__
        // Maybe<unit>::sse_load or Maybe<unit>::sse_loadu
        static void sse_load(
            __m128& m,
            const float* x1, const float* x2, const float* x3, const float* x4)
        { Maybe<yn2>::sse_load(m,x1,x2,x3,x4); }
        static void sse_load(
            __m128& m,
            const std::complex<float>* x1, const std::complex<float>* x2)
        { Maybe<yn2>::sse_load(m,x1,x2); }

        // Maybe<unit>::sse_add or Maybe<unit>::sse_store
        static void sse_add(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        { Maybe<yn2>::sse_add(x1,x2,x3,x4,m); }
        static void sse_add(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        { Maybe<yn2>::sse_add(x1,x2,m); }

        // Maybe<unit>::sse_addu or Maybe<unit>::sse_storeu
        static void sse_addu(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        { Maybe<yn2>::sse_addu(x1,x2,x3,x4,m); }
        static void sse_addu(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        { Maybe<yn2>::sse_addu(x1,x2,m); }
#endif

#ifdef __SSE2__
        static void sse_load(
            __m128d& m, const double* x1, const double* x2)
        { Maybe<yn2>::sse_load(m,x1,x2); }

        static void sse_add(
            double* x1, double* x2, const __m128d& m)
        { Maybe<yn2>::sse_add(x1,x2,m); }
        static void sse_add(
            std::complex<double>* x, const __m128d& m)
        { Maybe<yn2>::sse_add(x,m); }

        static void sse_addu(
            double* x1, double* x2, const __m128d& m)
        { Maybe<yn2>::sse_addu(x1,x2,m); }
        static void sse_addu(
            std::complex<double>* x, const __m128d& m)
        { Maybe<yn2>::sse_addu(x,m); }
#endif
    };
        
    template <bool yn2> 
    struct Maybe2<false,yn2>
    {
        template <class T1, class T2>
        static void add(T1& , const T2& ) {}
        template <class T1, class T2>
        static void add(ConjRef<T1> , const T2& ) {}
        template <class T1, bool C, class T2>
        static void add(TriRef<T1,C> , const T2& ) {}

#ifdef __SSE2__
        static void sse_load(
            __m128& m,
            const float* x1, const float* x2, const float* x3, const float* x4)
        { Maybe<yn2>::sse_loadu(m,x1,x2,x3,x4); }
        static void sse_load(
            __m128& m,
            const std::complex<float>* x1, const std::complex<float>* x2)
        { Maybe<yn2>::sse_loadu(m,x1,x2); }

        static void sse_add(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        { Maybe<yn2>::sse_store(x1,x2,x3,x4,m); }
        static void sse_add(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        { Maybe<yn2>::sse_store(x1,x2,m); }

        static void sse_addu(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        { Maybe<yn2>::sse_storeu(x1,x2,x3,x4,m); }
        static void sse_addu(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        { Maybe<yn2>::sse_storeu(x1,x2,m); }
#endif

#ifdef __SSE2__
        static void sse_load(
            __m128d& m, const double* x1, const double* x2)
        { Maybe<yn2>::sse_loadu(m,x1,x2); }

        static void sse_add(
            double* x1, double* x2, const __m128d& m)
        { Maybe<yn2>::sse_store(x1,x2,m); }
        static void sse_add(
            std::complex<double>* x, const __m128d& m)
        { Maybe<yn2>::sse_store(x,m); }

        static void sse_addu(
            double* x1, double* x2, const __m128d& m)
        { Maybe<yn2>::sse_storeu(x1,x2,m); }
        static void sse_addu(
            std::complex<double>* x, const __m128d& m)
        { Maybe<yn2>::sse_storeu(x,m); }
#endif
    };
       
    // A simpler structure that is mostly just to handing things like
    // int + complex<double>
    // which don't have a definition normally.
    struct ZSum
    {
        // complex + complex
        template <class T>
        static std::complex<T> sum(
            const std::complex<T>& x, const std::complex<T>& y)
        { return x + y; }
        template <class T1, class T2> 
        static std::complex<typename Traits2<T1,T2>::type> sum(
            const std::complex<T1>& x, const std::complex<T2>& y)
        {
            typedef typename Traits2<T1,T2>::type T12;
            return 
                Traits<std::complex<T12> >::convert(x) +
                Traits<std::complex<T12> >::convert(y);
        }

        // real + complex
        template <class T>
        static std::complex<T> sum(const T& x, const std::complex<T>& y)
        { return x + y; }
        template <class T1, class T2> 
        static std::complex<typename Traits2<T1,T2>::type> sum(
            const T1& x, const std::complex<T2>& y)
        {
            typedef typename Traits2<T1,T2>::type T12;
            return T12(x) + Traits<std::complex<T12> >::convert(y);
        }

        // complex + real
        template <class T>
        static std::complex<T> sum(const std::complex<T>& x, const T& y)
        { return x + y; }
        template <class T1, class T2> 
        static std::complex<typename Traits2<T1,T2>::type> sum(
            const std::complex<T1>& x, const T2& y)
        {
            typedef typename Traits2<T1,T2>::type T12;
            return Traits<std::complex<T12> >::convert(x) + T12(y);
        }

        // real + real
        template <class T>
        static T sum(const T& x, const T& y)
        { return x + y; }
        template <class T1, class T2> 
        static typename Traits2<T1,T2>::type sum(
            const T1& x, const T2& y)
        {
            typedef typename Traits2<T1,T2>::type T12;
            return T12(x) + T12(y);
        }
    };

#ifndef TMV_NO_THROW
    class Error : public std::runtime_error
    {
    public :
        std::string s1;
        std::string s2;
        Error(std::string _s1, std::string _s2="") throw() :
            std::runtime_error("TMV Error"), s1(_s1), s2(_s2) {}
        virtual ~Error() throw() {}
        virtual void write(std::ostream& os) const throw()
        { os << "TMV Error: " << s1 << s2 << std::endl; }
        virtual const char* what() const throw()
        { return s1.c_str(); }
    };

    static inline std::ostream& operator<<(
        std::ostream& os, const Error& e) throw()
    { e.write(os); return os; }

    class FailedAssert : public Error
    {
    public :
        std::string failed_assert;
        unsigned long line;
        std::string file;

        FailedAssert(std::string s, unsigned long l, 
                            std::string f) throw() :
            Error("Failed Assert statement ",s),
            failed_assert(s), line(l), file(f) {}
        ~FailedAssert() throw() {}
        void write(std::ostream& os) const throw()
        {
            os<<"TMV Failed Assert: "<<failed_assert<<std::endl<< 
                "on line "<<line<<" in file "<<file<<std::endl;
        }
    };
#endif

#ifndef TMV_DEBUG
#define TMVAssert(x) do {} while (0)
#else
#ifdef TMV_NO_THROW
#define TMVAssert(x) do { if(!(x)) { \
    std::cerr<<"Failed Assert: "<<#x<<" on line " \
    <<__LINE__<<" in file "<<__FILE__<<std::endl; exit(1); } } while (false)
#else
#define TMVAssert(x) do { if(!(x)) { \
    throw tmv::FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif
#endif
#ifdef TMV_NO_THROW
#define TMVAssert2(x) do { if(!(x)) { \
    std::cerr<<"Failed Assert: "<<#x<<" on line " \
    <<__LINE__<<" in file "<<__FILE__<<std::endl; exit(1); } } while (false)
#else
#define TMVAssert2(x) do { if(!(x)) { \
    throw tmv::FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif

#ifndef TMV_NO_THROW
    class ReadError : public Error
    {
    public :
        ReadError() throw() :
            Error("Invalid istream input encountered.") {}
        ReadError(std::string s) throw() :
            Error("Invalid istream input encountered while reading ",s) {}
        ~ReadError() throw() {}
        void write(std::ostream& os) const throw()
        { 
            os << "TMV Read Error: " << Error::s1 << ' ' << Error::s2 << 
                std::endl; 
        }
    };

    // Each Vector, Matrix type has it's own derived ReadError class
    // which outputs whatever has been read so far, and possibly
    // what was expected instead of what was actually read in.

    class Singular : public Error
    {
    public :
        Singular() throw() :
            Error("Encountered singular matrix.") {}
        Singular(std::string s) throw() :
            Error("Encountered singular ",s) {}
        ~Singular() throw() {}
        void write(std::ostream& os) const throw()
        {
            os << "TMV Singular: " << Error::s1 << ' ' << Error::s2 <<
                std::endl; 
        }
    };

    template <class M>
    class SingularMatrix : public Singular
    {
    public:
        typename M::copy_type m;

        SingularMatrix(const M& _m) : 
            Singular(TMV_Text(_m)), m(_m) {}
        ~SingularMatrix() throw() {}
        void write(std::ostream& os) const throw()
        {
            Singular::write(os);
            m.write(os);
            os << std::endl;
        }
    };

    class NonPosDef : public Error
    {
    public:
        NonPosDef() throw() :
            Error("Invalid non-positive-definite matrix found.") {}
        NonPosDef(std::string s) throw() :
            Error("Non-positive-definite matrix found in ",s) {}
        ~NonPosDef() throw() {}
        void write(std::ostream& os) const throw()
        {
            os << "TMV NonPosDef: " << Error::s1 << ' ' << Error::s2 << 
                std::endl; 
        }
    };

    template <class M>
    class NonPosDefMatrix : public NonPosDef
    {
    public:
        typename M::copy_type m;

        NonPosDefMatrix(const M& _m) : 
            NonPosDef(TMV_Text(_m)), m(_m) {}
        ~NonPosDefMatrix() throw() {}
        void write(std::ostream& os) const throw()
        {
            NonPosDef::write(os);
            m.write(os);
            os << std::endl;
        }
    };
#endif

#ifdef TMV_WARN
    // defined in TMV_Vector.cpp
    // initialized with &std::cout;
    extern std::ostream* warn_out;
    static inline void TMV_Warning(std::string s)
    {
        if (warn_out) { 
            *warn_out << "Warning:\n" << s << std::endl;
        }
    }

    static inline std::ostream* WriteWarningsTo(std::ostream* newos)
    { std::ostream* temp = warn_out; warn_out = newos; return temp; }

    static inline void NoWarnings()
    { warn_out = 0; }
#else
    static inline void TMV_Warning(std::string ) {}
    static inline std::ostream* WriteWarningsTo(std::ostream* os) { return os; }
    static inline void NoWarnings() {}
#endif

    // A helper structure that acts like an int,
    // but only bothers to make the integer if S == UNKNOWN.
    // It also checks the constructor if S != UNKNOWN
    template <int S>
    struct CheckedInt
    {
        CheckedInt(int s) 
        { TMVAssert(s == S); }
        operator int() const { return S; }
    };
    template <>
    struct CheckedInt<FortranStyle> 
    // This helps catch programming errors from the change in the 
    // template parameters of VectorView, MatrixView, etc. from
    // version 0.63 to 0.70.
    // It used to be that to make a fortran-style view, you would write:
    // VectorView<T,FortranStyle> v = m.row(i);
    // This is now wrong, since the first template parameter is now the step
    // size, so this would treat FortranStyle as the step size of the view.
    // One thing I did to help catch this was to make FortranStyle a 
    // crazy number that will probably never be used as an actual step
    // (789234).  This specialization of CheckedInt catches when FortranStyle
    // is used as a step and gives an error.
    {
        CheckedInt(int s) 
        { 
#if 0 
            TMVStaticAssert(false);
#else
            std::string str = 
                "FortranStyle used as a step size is probably an error...\n"
                "Now you should use F at the end of the name "
                "(e.g. VectorViewF)\n"
                "instead of using I=FotranStyle as a template parameter.\n";
            TMV_Warning(str);
            TMVAssert(s == FortranStyle); 
#endif
        }
        operator int() const { return FortranStyle; }
    };
    template <>
    struct CheckedInt<UNKNOWN>
    {
        int step;
        CheckedInt(int s) : step(s) {}
        operator int() const 
        { 
#ifdef TMV_DEBUG
            TMVAssert(step !=  -987234);
#endif
            return step; 
        }
        ~CheckedInt() {
#ifdef TMV_DEBUG
            step = -987234;
#endif
        }
    };

    // A similar helper to account for possibly unknown DiagType
    template <DiagType D>
    struct DiagInt
    {
        DiagInt(DiagType d) { TMVAssert(d == D); }
        operator DiagType() const { return D; }
    };
    template <>
    struct DiagInt<UnknownDiag>
    {
        DiagType dt;
        DiagInt(DiagType d) : dt(d) {}
        operator DiagType() const { return dt; }
    };

    template <class T> 
    static typename Traits<T>::real_type TMV_Epsilon() 
    { return std::numeric_limits<typename Traits<T>::real_type>::epsilon(); }

    template <class T>
    static bool TMV_Underflow(T x)
    {
        typedef typename Traits<T>::real_type RT;
        return TMV_ABS2(x) < 
            std::numeric_limits<RT>::min() * RT(tmv::Traits<T>::twoifcomplex); 
    }

    static inline bool TMV_Underflow(int )
    { return false; }

    template <class T> 
    static std::string TMV_Text(const T&)
    { return std::string("Unknown (") + typeid(T).name() + ")"; }

    static inline std::string TMV_Text(const double&)
    { return "double"; }

    static inline std::string TMV_Text(const float&)
    { return "float"; }

    static inline std::string TMV_Text(const int&)
    { return "int"; }

    static inline std::string TMV_Text(const long double&)
    { return "long double"; }

    template <class T> 
    static std::string TMV_Text(std::complex<T>)
    { return std::string("complex<") + TMV_Text(T()) + ">"; }

    static inline std::string TMV_Text(StorageType s)
    {
        return 
            s == RowMajor ? "RowMajor" :
            s == ColMajor ? "ColMajor" :
            s == DiagMajor ? "DiagMajor" :
            s == NoMajor ? "NoMajor" :
            s == RowPacked ? "RowPacked" :
            s == ColPacked ? "ColPacked" :
            "unkown StorageType";
    }

    static inline std::string TMV_Text(IndexStyle i)
    { return i == CStyle ? "CStyle" : "FortranStyle"; }

    static inline std::string TMV_Text(DivType d)
    { 
        return 
            d==XXX ? "XXX" :
            d==LU ? "LU" :
            d==CH ? "CH" :
            d==QR ? "QR" :
            d==QRP ? "QRP" :
            d==SV ? "SV" :
            ( std::string("unkown DivType: ") + 
              // This is a quick and dirty int->string for int < 1000
              char('0' + (d/100)) + 
              char('0' + (d%100)/10) + 
              char('0' + (d%10)) );
    }

    static inline std::string TMV_Text(DiagType d)
    { 
        return 
            d == UnitDiag ? "UnitDiag" :
            d == NonUnitDiag ? "NonUnitDiag" :
            "UnknownDiag"; 
    }

#if 0
    static inline std::string TMV_Text(UpLoType u)
    { return u == Upper ? "Upper" : "Lower"; }

    static inline std::string TMV_Text(SymType s)
    { return s == Sym ? "Sym" : "Herm"; }
#endif


#ifdef TMV_NO_STL_AUTO_PTR
    // This is copied more or less verbatim from gcc's auto_ptr implementation.
    template <class X> 
    class auto_ptr {

    private:
        X* ptr;
        template <class Y> 
        struct auto_ptr_ref {
            Y* ptr;
            explicit auto_ptr_ref(Y* p) : ptr(p) {}
        };

    public:
        typedef X element_type;
        explicit auto_ptr(X* p = 0) throw() : ptr(p) {}
        auto_ptr(auto_ptr& a) throw() : ptr(a.release()) {}
        auto_ptr(auto_ptr_ref<X> ref) throw() : ptr(ref.ptr) {}
        template <class Y> 
        auto_ptr(auto_ptr<Y>& a) throw() :
            ptr(a.release()) {}

        auto_ptr& operator=(auto_ptr& a) throw() 
        {
            reset(a.release());
            return *this;
        }
        template <class Y> 
        auto_ptr& operator=(auto_ptr<Y>& a) throw() 
        {
            reset(a.release());
            return *this;
        }
        auto_ptr& operator=(auto_ptr_ref<X> ref) throw() 
        {
            if (ref.ptr != this->get())
            { 
                delete ptr;
                ptr = ref.ptr;
            }
            return *this;
        }

        ~auto_ptr() throw() { delete ptr; }

        template <class Y> 
        operator auto_ptr_ref<Y>() throw()
        { return auto_ptr_ref<Y>(this->release()); }
        template <class Y> 
        operator auto_ptr<Y>() throw()
        { return auto_ptr<Y>(this->release()); }

        X& operator*() const throw() { return *ptr; }
        X* operator->() const throw() { return ptr; }
        X* get() const throw() { return ptr; }
        X* release() throw() 
        {
            X* tmp = ptr;
            ptr = 0;
            return tmp;
        }
        void reset(X* p = 0) throw() 
        {
            if (p != ptr) {
                delete ptr;
                ptr = p;
            }
        }
    }
#else
    using std::auto_ptr;
#endif

    // Use DEBUGPARAM(x) for parameters that are only used in TMVAssert
    // statements.  So then they don't give warnings when compiled with 
    // -DNDEBUG
#ifdef TMV_DEBUG
#define TMV_DEBUG_PARAM(x) x
#else
#define TMV_DEBUG_PARAM(x)
#endif

} // namespace tmv

#endif
