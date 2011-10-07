

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

#ifdef TMV_NO_LIB
#define TMV_NO_WARN
#define TMV_NO_INST_DOUBLE
#define TMV_NO_INST_FLOAT
#endif

#ifndef TMV_NO_WARN
#define TMV_WARN
#endif

// The meanings of the different optimization levels are:
// 0 = No optimizations.  Use the simplest code possible.  
//     However, we still respect the majority of the matrices, so we
//     at least pick an algorithm that strides unit steps when possible.
//     This could be considered similar to the reference BLAS.
// 1 = Use optimizations that significantly increase the performance of
//     a large number of matrix sizes without a large increase in code size.
// 2 = Add in optimizations that may increase the code size somewhat,
//     but are deemed worthwhile because they significantly improve
//     performance for many matrix sizes.
// 3 = Add in more optimizations that are only useful for particular 
//     sizes, or that increase code size by a lot.
//
// The default optimization level is 2 if not specified with a -D flag.

#ifndef TMV_OPT
#define TMV_OPT 2
#endif

// The L1_CACHE and L2_CACHE parameters can be set to the cache values 
// for your machine if you know them (in KB).  
// If don't know these or don't want to bother defining them, it's ok.
// The below values will produce reasonable code for all machines,
// even if it may not be exactly optimal.
#ifndef TMV_L1_CACHE
#define TMV_L1_CACHE 32
#endif
#ifndef TMV_L2_CACHE
#define TMV_L2_CACHE 256
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
#include <sstream>
#include <iostream>
#endif

#include "../util/portable_platform.h"

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

#ifndef TMV_NO_INST_MIX
#define TMV_INST_MIX
#endif

// These next two are not instantiated by default:

//#define TMV_INST_LONGDOUBLE

//#define TMV_INST_INT

#ifdef __GNUC__
#define TMV_DEPRECATED(func) func __attribute__ ((deprecated))
#define TMV_INLINE inline __attribute__ ((__always_inline__))
#elif defined(_MSC_VER)
#define TMV_DEPRECATED(func) __declspec(deprecated) func
#define TMV_INLINE inline __forceinline
#elif defined(__PGI)
#define TMV_DEPRECATED(func) func __attribute__ ((deprecated))
#define TMV_INLINE inline __attribute__ ((__always_inline__))
#else
#define TMV_DEPRECATED(func) func
#define TMV_INLINE inline
#endif
#ifdef TMV_DEBUG
#define TMV_INLINE_ND inline
#else
#define TMV_INLINE_ND TMV_INLINE
#endif

namespace tmv {

    TMV_INLINE std::string TMV_Version() { return "0.70"; }
#define TMV_MAJOR_VERSION 0
#define TMV_MINOR_VERSION 70
#define TMV_VERSION_AT_LEAST(major,minor) \
    ( (major > TMV_MAJOR_VERSION) || \
      (major == TMV_MAJOR_VERSION && minor >= TMV_MINOR_VERSION) )

    // This is based on the BOOST implementation BOOST_STATIC_ASSERT:
    template<bool>
    struct tmv_static_assert;
    template<>
    struct tmv_static_assert<true>
    { static TMV_INLINE void call() {} };
#define TMVStaticAssert(e) \
    do { tmv_static_assert<(e) != 0>::call(); } while (false)

    // Attributes are stored as a binary field, so they can be |'ed
    // together to calculate the full set of attributes.
    //
    // For each object we have a default set of attributes that is 
    // implied by A=0.  This way we can write, for example:
    // UpperTriMatrix<T,UnitDiag> U;
    // and this will imply (UnitDiag | ColMajor | NonPacked | CStyle).
    //
    // Attributes that are always the default (whenever they are allowed)
    // can be defined as 0, since the absence of the converse is sufficient.
    // e.g. CStyle, NonConj, NonUnit, NonMajor, NonPacked, UnknownDiag.
    //
    // Some pairs of attributes are logical opposites, but we can't
    // set either to 0, since some objects have one as the default, and
    // others have the other.  e.g. WithDivider and NoDivider.
    // So these each need to have a non-zero bit.
    //
    // (Not implemented yet: Packed, ZeroDiag)
    //
    enum ConjType { NonConj = 0, Conj = 0x1 };
    enum IndexStyle { CStyle = 0, FortranStyle = 0x2 };
    enum StepType { NonUnit = 0, Unit = 0x4 };
    enum StorageType {
        NonMajor = 0, ColMajor = 0x4, RowMajor = 0x8, DiagMajor = 0x10,
        AllStorageType = 0x1c };
    enum DiagType {
        UnknownDiag = 0, NonUnitDiag = 0x20, UnitDiag = 0x40,
        ZeroDiag = 0x80, AllDiagType = 0xe0 };
    enum PackType { NonPacked = 0, Packed = 0x100 };
    enum DivStatus { 
        NoDivider = 0x200, WithDivider = 0x400, AllDivStatus = 0x600 };
    enum AliasStatus { 
        NoAlias = 0x800, CheckAlias = 0x1000, AllAliasStatus = 0x1800 };
    enum UpLoType { Upper = 0x2000, Lower = 0x4000 };

    template <bool iszero>
    struct DoNonZeroMessage;
    template <>
    struct DoNonZeroMessage<true>
    { 
        template <class Dummy>
        static TMV_INLINE void call() {} 
    };
#if 1
    template <>
    struct DoNonZeroMessage<false>
    { 
        template <class Dummy>
        static TMV_INLINE void call() 
        {
            // I think this might be portable.  We'll see...
            // At the very least, if it doesn't emit the #pragma message,
            // I think it will just ignore it, so at least I don't
            // think this will break anything...
//#pragma message ("Use of multiple Attributes is deprecated.\ne.g. Matrix<T,ColMajor,FortranStyle> should now be \nMatrix<T,ColMajor | FortranStyle>.")
            TMVStaticAssert(Dummy() == 100);
        }
    };
#endif
    template <bool iszero>
    struct NonZeroMessage
    { 
        static TMV_INLINE void call() 
        { DoNonZeroMessage<iszero>::template call<int>(); } 
    };


    template <int A>
    struct Attrib
    {
        enum { vectoronly = ((A & ~AllAliasStatus) <= 0x7) };
        enum { conj = !(!(A & Conj)) };
        enum { fort = !(!(A & FortranStyle)) };
        enum { unit = !(!(A & Unit)) };
        enum { colmajor = !(!(A & ColMajor)) };
        enum { rowmajor = !(!(A & RowMajor)) };
        enum { diagmajor = !(!(A & DiagMajor)) };
        enum { nonmajor = !(colmajor || rowmajor || diagmajor) };
        enum { nonunitdiag = !(!(A & NonUnitDiag)) };
        enum { unitdiag = !(!(A & UnitDiag)) };
        enum { zerodiag = !(!(A & ZeroDiag)) };
        enum { unknowndiag = !(unitdiag || nonunitdiag || zerodiag) };
        enum { packed = !(!(A & Packed)) };
        enum { nodivider = !(!(A & NoDivider)) };
        enum { withdivider = !(!(A & WithDivider)) };
        enum { noalias = !(!(A & NoAlias)) };
        enum { checkalias = !(!(A & CheckAlias)) };
        enum { lower = !(!(A & Lower)) };
        enum { upper = !(!(A & Upper)) };

        enum { stor = (
                rowmajor ? RowMajor : colmajor ? ColMajor :
                diagmajor ? DiagMajor : NonMajor ) };

        static TMV_INLINE void checkzero() 
        { NonZeroMessage<A==0>::call(); }
        static std::string text()
        {
            return 
                std::string() +
                ( (A & ColMajor) ? "ColMajor" : 
                  (A & RowMajor) ? "RowMajor" :
                  (A & DiagMajor) ? "DiagMajor" : "NonMajor") +
                ((A & Conj) ? "|Conj" : "") +
                ((A & FortranStyle) ? "|FortranStyle " : "") +
                ((A & NonUnitDiag) ? "|NonUnitDiag " : "") +
                ((A & UnitDiag) ? "|UnitDiag" : "") +
                ((A & ZeroDiag) ? "|ZeroDiag" : "") +
                ((A & Packed) ? "|Packed" : "") +
                ((A & Lower) ? "|Lower" : "") +
                ((A & Upper) ? "|upper" : "") +
                ((A & NoDivider) ? "|NoDivider" : "") +
                ((A & WithDivider) ? "|WithDivider" : "") +
                ((A & NoAlias) ? "|NoAlias" : "") +
                ((A & CheckAlias) ? "|CheckAlias" : "");
        }
        static std::string vtext()
        {
            return 
                std::string() +
                ((A & Unit) ? "Unit" : "NonUnit") +
                ((A & Conj) ? "|Conj" : "") +
                ((A & FortranStyle) ? "|FortranStyle" : "") +
                ((A & NoAlias) ? "|NoAlias" : "") +
                ((A & CheckAlias) ? "|CheckAlias" : "");
        }
    };

    // TMV_UNKNOWN is the value of _size, _step, etc. whenever it
    // is not known at compile time.
    // We use for this value the maximally negative int.
    // In binary, this is a 1 followed by all zeros.
    const int TMV_UNKNOWN = (1<<(sizeof(int)*8-1));

    enum DivType {
        XX=0, LU=1, CH=2, QR=4, QRP=8, SV=16,
        // We store the divtype in a binary field integer.
        // In addition to the above, we also use the same object to 
        // store the following other flags related to division.
        // So these values must not clash with the above DivType values.
        // These aren't technically DivType's but since they are 
        // stored together, I think this adds to type-safety.
        DivInPlaceFlag = 32,
        SaveDivFlag = 64,
        // And finally shorthand for "one of the real DivType's":
        DivTypeFlags = 31
    };
    TMV_INLINE DivType operator|(DivType a, DivType b) 
    { 
        return static_cast<DivType>(
            static_cast<int>(a) | static_cast<int>(b)); 
    }
    TMV_INLINE DivType operator&(DivType a, DivType b) 
    { 
        return static_cast<DivType>(
            static_cast<int>(a) & static_cast<int>(b)); 
    }
    TMV_INLINE DivType& operator|=(DivType& a, DivType b) 
    { a = (a|b); return a; }
    TMV_INLINE DivType& operator&=(DivType& a, DivType b) 
    { a = (a&b); return a; }
    TMV_INLINE DivType operator~(DivType a) 
    { return static_cast<DivType>(~static_cast<int>(a)); }

    // This is a type to use when there is no valid return type for a 
    // particular function. 
    // It will give a compiler error if it is ever used.
    class InvalidType { private: InvalidType(); };

    template <class T> class ConjRef;
    template <class T, bool C> class TriRef;

    // Declare all the Base types, since sometimes I reference them
    // before they are defined.  e.g. in the *_Funcs.h files.
    template <class V> class BaseVector;
    template <class V> class BaseVector_Calc;
    template <class V> class BaseVector_Mutable;
    template <class M> class BaseMatrix;
    template <class M> class BaseMatrix_Calc;
    template <class M> class BaseMatrix_Mutable;
    template <class M> class BaseMatrix_Rec;
    template <class M> class BaseMatrix_Rec_Mutable;
    template <class M> class BaseMatrix_Diag;
    template <class M> class BaseMatrix_Diag_Mutable;
    template <class M> class BaseMatrix_Tri;
    template <class M> class BaseMatrix_Tri_Mutable;

    // Define which classes get instantiated.
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

        typedef T type;
        typedef T real_type;
        typedef std::complex<T> complex_type;
        typedef typename TypeSelect<isinteger,double,T>::type float_type;

        typedef T& conj_reference;

        static TMV_INLINE const T& convert(const T& x) { return x; }
        template <class T2>
        static TMV_INLINE T convert(const T2& x) { return T(x); }

        static TMV_INLINE T constr_value() { return T(888); }
        static TMV_INLINE T destr_value() { return T(999); }
    };

    template <class T>
    struct Traits<std::complex<T> >
    {
        enum { isreal = false };
        enum { iscomplex = true };
        enum { isinst = InstType<T>::inst };
        enum { isinteger = Traits<T>::isinteger };
        enum { twoifcomplex = 2 };

        typedef std::complex<T> type;
        typedef T real_type;
        typedef std::complex<T> complex_type;
        typedef std::complex<typename Traits<T>::float_type> float_type;

        typedef ConjRef<type>& conj_reference;

        static TMV_INLINE const type& convert(const type& x) { return x; }
        template <class T2>
        static TMV_INLINE type convert(const T2& x) { return type(x); }
        // This last one is the real reason to have a convert function.
        // Because std::complex doesn't have a constructor from a 
        // std::complex of a different base type.
        // So this makes it easier to do the conversion.
        template <class T2>
        static TMV_INLINE std::complex<T> convert(const std::complex<T2>& x) 
        { return std::complex<T>(x.real(),x.imag()); }

        static TMV_INLINE type constr_value() { return type(888,888); }
        static TMV_INLINE type destr_value() { return type(999,999); }
    };

    template <class T> struct Traits<T&> : public Traits<T> {};
    template <class T> struct Traits<const T> : public Traits<T> {};
    template <class T> struct Traits<const T&> : public Traits<T> {};

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
    { static TMV_INLINE T apply(const T& x) { return x; } };
    template <class T>
    struct DoConj_Helper<true,std::complex<T> >
    {
        static TMV_INLINE std::complex<T> apply(std::complex<T> x) 
        { return std::conj(x); } 
    };
    template <bool C, class T>
    TMV_INLINE T DoConj(const T& x) { return DoConj_Helper<C,T>::apply(x); }

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
                TMV_UNKNOWN ) 
        };
        enum { half_roundup = (
                // For very large S, just call it TMV_UNKNOWN to keep from 
                // having big complicated recursive structures.
                S > 128 ? TMV_UNKNOWN :
                S > 16 ? ((((S-1)>>5)+1)<<4) :
                (S>>1)  )
        };
        static inline int text() { return S; }
    };
    template <>
    struct IntTraits<TMV_UNKNOWN>
    {
        enum { negS = TMV_UNKNOWN };
        enum { twoS = TMV_UNKNOWN };
        enum { halfS = TMV_UNKNOWN };
        enum { Sm1 = TMV_UNKNOWN };
        enum { Sp1 = TMV_UNKNOWN };
        enum { Sp2 = TMV_UNKNOWN };
        enum { ispowerof2 = false };
        enum { log = TMV_UNKNOWN };
        enum { half_roundup = TMV_UNKNOWN };
        static inline const char* text() { return "UNKNOWN"; }
    };

    template <int S1, int S2>
    struct IntTraits2
    {
        enum { sum = S1 + S2 };
        enum { diff = S1 - S2 };
        enum { prod = S1 * S2 };
        enum { safeprod = IntTraits2<
            (S1<300 ? S1 : TMV_UNKNOWN), (S2<300 ? S2 : TMV_UNKNOWN)>::prod };
        //enum { quot = S2 == 0 ? 0 : S1 / S2 };
        enum { quot = S1 / S2 };
        enum { min = S1 < S2 ? S1 : S2 };
        enum { max = S1 > S2 ? S1 : S2 };
    };
    template <int S1>
    struct IntTraits2<S1,0>
    { // Specialization just to avoid division by zero warning.
        enum { sum = S1 };
        enum { diff = S1 };
        enum { prod = 0 };
        enum { safeprod = 0 };
        enum { quot = 0 };
        enum { min = S1 < 0 ? S1 : 0 };
        enum { max = S1 > 0 ? S1 : 0 };
    };
    template <int S1>
    struct IntTraits2<S1,TMV_UNKNOWN>
    {
        enum { sum = TMV_UNKNOWN };
        enum { diff = TMV_UNKNOWN };
        enum { prod = TMV_UNKNOWN };
        enum { safeprod = TMV_UNKNOWN };
        enum { quot = TMV_UNKNOWN };
        enum { min = TMV_UNKNOWN };
        enum { max = TMV_UNKNOWN };
    };
    template <int S2>
    struct IntTraits2<TMV_UNKNOWN,S2>
    {
        enum { sum = TMV_UNKNOWN };
        enum { diff = TMV_UNKNOWN };
        enum { prod = TMV_UNKNOWN };
        enum { safeprod = TMV_UNKNOWN };
        enum { quot = TMV_UNKNOWN };
        enum { min = TMV_UNKNOWN };
        enum { max = TMV_UNKNOWN };
    };
    template <>
    struct IntTraits2<TMV_UNKNOWN,TMV_UNKNOWN>
    {
        enum { sum = TMV_UNKNOWN };
        enum { diff = TMV_UNKNOWN };
        enum { prod = TMV_UNKNOWN };
        enum { safeprod = TMV_UNKNOWN };
        enum { quot = TMV_UNKNOWN };
        enum { min = TMV_UNKNOWN };
        enum { max = TMV_UNKNOWN };
    };
    template <>
    struct IntTraits2<TMV_UNKNOWN,0>
    {
        enum { sum = TMV_UNKNOWN };
        enum { diff = TMV_UNKNOWN };
        enum { prod = TMV_UNKNOWN };
        enum { safeprod = TMV_UNKNOWN };
        enum { quot = TMV_UNKNOWN };
        enum { min = TMV_UNKNOWN };
        enum { max = TMV_UNKNOWN };
    };

    template <class T>
    inline T TMV_SQR(const T& x) 
    { return x*x; }

    template <class T>
    inline typename Traits<T>::float_type TMV_SQRT(const T& x) 
    { return std::sqrt(typename Traits<T>::float_type(x)); }

    template <class T>
    inline typename Traits<T>::float_type TMV_EXP(const T& x) 
    { return std::exp(typename Traits<T>::float_type(x)); }

    template <class T>
    inline typename Traits<T>::float_type TMV_LOG(const T& x) 
    { return std::log(typename Traits<T>::float_type(x)); }

    template <class T>
    inline T TMV_NORM(const T& x) 
    { return x*x; }

    template <class T>
    TMV_INLINE T TMV_NORM(const std::complex<T>& x) 
    { return std::norm(x); }

    template <class T>
    TMV_INLINE T TMV_CONJ(const T& x)
    { return x; }

    template <class T>
    TMV_INLINE std::complex<T> TMV_CONJ(const std::complex<T>& x)
    { return std::conj(x); }

    template <class T>
    TMV_INLINE T TMV_REAL(const T& x)
    { return x; }

    template <class T>
    TMV_INLINE T TMV_REAL(const std::complex<T>& x)
    { return x.real(); }

    template <class T>
    TMV_INLINE T TMV_IMAG(const T& )
    { return T(0); }

    template <class T>
    TMV_INLINE T TMV_IMAG(const std::complex<T>& x)
    { return x.imag(); }

    // Many implemenations of complex allow x.real() and x.imag() to return 
    // references to the actual locations of these elements.
    // However, the standard technically says that these return by 
    // value, not reference.  So we need to use a reinterpret_cast to 
    // really make sure we get the right thing.
    template <class T>
    TMV_INLINE T& TMV_REAL_PART(std::complex<T>& x)
    { return reinterpret_cast<T&>(x); }

    template <class T>
    TMV_INLINE T& TMV_REAL_PART(T& x)
    { return x; }

    template <class T>
    TMV_INLINE T& TMV_IMAG_PART(std::complex<T>& x)
    { return *(reinterpret_cast<T*>(&x)+1); }

    template <class T>
    inline T TMV_ARG(const T& x)
    { return x >= T(0) ? T(1) : T(-1); }

    template <class T>
    inline typename Traits<T>::float_type TMV_ARG(const std::complex<T>& x)
    {
        typedef typename Traits<T>::float_type FT;
        return arg(Traits<std::complex<FT> >::convert(x)); 
    }

    template <class T>
    TMV_INLINE T TMV_ABS(const T& x)
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
    TMV_INLINE T TMV_ABS2(const T& x)
    { return std::abs(x); }

    template <class T>
    inline T TMV_ABS2(const std::complex<T>& x)
    { return std::abs(x.real()) + std::abs(x.imag()); }

    template <class T>
    inline T TMV_SIGN(const T& x, const T& )
    { return x > 0 ? T(1) : T(-1); }

    template <class T>
    inline std::complex<T> TMV_SIGN(const std::complex<T>& x, const T& absx)
    { return absx > 0 ? x/absx : std::complex<T>(1); }

    template <class T>
    inline T TMV_MIN(const T& x, const T& y)
    { return x > y ? y : x; }

    template <class T>
    inline T TMV_MAX(const T& x, const T& y)
    { return x > y ? x : y; }

    template <class T>
    inline void TMV_SWAP(T& x, T& y)
    { T z = x; x = y; y = z; }

  // Useful for debugging SSE routines:
#ifdef __SSE__
    inline std::ostream& operator<<(std::ostream& os, const __m128& xf)
    {
        os << &xf<<": ";
        float ff[4];
        _mm_storeu_ps(ff,xf);
        os << ff[0]<<","<<ff[1]<<","<<ff[2]<<","<<ff[3]<<" ";
        return os;
    }
#endif

#ifdef __SSE2__
    inline std::ostream& operator<<(std::ostream& os, const __m128d& xd)
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
        static TMV_INLINE typename Traits2<T1,T2>::type rprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*ry - ix*iy); }
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type iprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*iy + ix*ry); }
    };
    template <>
    struct ZProd2<true,false>
    {
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type rprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*ry + ix*iy); }
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type iprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*iy - ix*ry); }
    };
    template <>
    struct ZProd2<false,true>
    {
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type rprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*ry + ix*iy); }
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type iprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (ry*ix - iy*rx); }
    };
    template <>
    struct ZProd2<true,true>
    {
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type rprod(
            const T1& rx, const T1& ix, const T2& ry, const T2& iy)
        { return (rx*ry - ix*iy); }
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type iprod(
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
        static inline typename Traits2<T1,T2>::type rprod(
            const std::complex<T1>& x, const std::complex<T2>& y)
        {
            return ZProd2<conj1,conj2>::rprod(
                x.real(),x.imag(),y.real(),y.imag()); 
        }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type iprod(
            const std::complex<T1>& x, const std::complex<T2>& y)
        {
            return ZProd2<conj1,conj2>::iprod(
                x.real(),x.imag(),y.real(),y.imag()); 
        }

        // real * complex
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type rprod(
            const T1& x, const std::complex<T2>& y)
        { return x*y.real(); }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type iprod(
            const T1& x, const std::complex<T2>& y)
        { return Maybe<conj2>::neg(x*y.imag()); }

        // complex * real
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type rprod(
            const std::complex<T1>& x, const T2& y)
        { return y*x.real(); }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type iprod(
            const std::complex<T1>& x, const T2& y)
        { return Maybe<conj1>::neg(y*x.imag()); }

        // real * real
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type rprod(
            const T1& x, const T2& y)
        { return x*y; }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type iprod(
            const T1& x, const T2& y)
        { return typename Traits2<T1,T2>::type(0); }

        // scaling * real
        template <class T1, class T2>
        static inline typename Traits2<typename Traits<T1>::real_type,T2>::type rprod(
            const Scaling<0,T1>& x, const T2& y)
        { return rprod(x.x,y); }
        template <class T1, class T2>
        static inline typename Traits2<typename Traits<T1>::real_type,T2>::type iprod(
            const Scaling<0,T1>& x, const T2& y)
        { return iprod(x.x,y); }
        template <class T1, class T2>
        static inline T2 rprod(const Scaling<1,T1>& x, const T2& y)
        { return y; }
        template <class T1, class T2>
        static inline T2 iprod(const Scaling<1,T1>& x, const T2& y)
        { return T2(0); }
        template <class T1, class T2>
        static inline T2 rprod(const Scaling<-1,T1>& x, const T2& y)
        { return -y; }
        template <class T1, class T2>
        static inline T2 iprod(const Scaling<-1,T1>& x, const T2& y)
        { return T2(0); }

        // scaling * complex
        template <class T1, class T2>
        static inline typename Traits2<typename Traits<T1>::real_type,T2>::type rprod(
            const Scaling<0,T1>& x, const std::complex<T2>& y)
        { return rprod(x.x,y); }
        template <class T1, class T2>
        static inline typename Traits2<typename Traits<T1>::real_type,T2>::type iprod(
            const Scaling<0,T1>& x, const std::complex<T2>& y)
        { return iprod(x.x,y); }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type rprod(
            const Scaling<1,T1>& x, const std::complex<T2>& y)
        { return y.real(); }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type iprod(
            const Scaling<1,T1>& x, const std::complex<T2>& y)
        { return Maybe<conj2>::neg(y.imag()); }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type rprod(
            const Scaling<-1,T1>& x, const std::complex<T2>& y)
        { return -y.real(); }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type iprod(
            const Scaling<-1,T1>& x, const std::complex<T2>& y)
        { return Maybe<!conj2>::neg(y.imag()); }

        // full prod
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type prod(
            const T1& x, const T2& y)
        {
            typedef typename Traits2<T1,T2>::type PT;
            const bool iscomplex = Traits<PT>::iscomplex;
            return makeprod<iscomplex,T1,T2>::call(x,y);
        }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type prod(
            const Scaling<1,T1>& x, const T2& y)
        { 
            return tmv::Traits<typename Traits2<T1,T2>::type>::convert(
                Maybe<conj2>::conj(y)); 
        }
        template <int ix, class T1, class T2>
        static inline typename Traits2<T1,T2>::type prod(
            const T1& x, const Scaling<1,T2>& y)
        { 
            return tmv::Traits<typename Traits2<T1,T2>::type>::convert(
                Maybe<conj1>::conj(x)); 
        }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type prod(
            const Scaling<-1,T1>& x, const T2& y)
        { 
            return tmv::Traits<typename Traits2<T1,T2>::type>::convert(
                -Maybe<conj2>::conj(y)); 
        }
        template <int ix, class T1, class T2>
        static inline typename Traits2<T1,T2>::type prod(
            const T1& x, const Scaling<-1,T2>& y)
        { 
            return tmv::Traits<typename Traits2<T1,T2>::type>::convert(
                -Maybe<conj1>::conj(x)); 
        }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type prod(
            const Scaling<0,T1>& x, const T2& y)
        { return prod(x.x , y); }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type prod(
            const T1& x, const Scaling<0,T2>& y)
        { return prod(x , y.x); }

        template <bool iscomplex, class T1, class T2>
        struct makeprod;
        template <class T1, class T2>
        struct makeprod<true,T1,T2>
        {
            typedef typename Traits2<T1,T2>::type PT;
            static inline PT call(const T1& x, const T2& y)
            { return PT(rprod(x,y),iprod(x,y)); }
        };
        template <class T1, class T2>
        struct makeprod<false,T1,T2>
        {
            typedef typename Traits2<T1,T2>::type PT;
            static inline PT call(const T1& x, const T2& y)
            { return x*y; }
        };

        // quot
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type quot(
            const T1& x, const T2& y)
        {
            const bool c2 = Traits<T2>::iscomplex;
            return makequot<c2,T1,T2>::call(x,y);
        }
        template <int ix, class T1, class T2>
        static inline typename Traits2<T1,T2>::type quot(
            const Scaling<ix,T1>& x, const T2& y)
        { return quot(T1(x),y); }
        template <int ix, class T1, class T2>
        static inline typename Traits2<T1,T2>::type quot(
            const T1& x, const Scaling<ix,T2>& y)
        { return prod(x,y); }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type quot(
            const T1& x, const Scaling<0,T2>& y)
        { return quot(x , y.x); }

        template <bool c2, class T1, class T2>
        struct makequot;
        template <class T1, class T2>
        struct makequot<true,T1,T2>
        {
            typedef typename Traits2<T1,T2>::type PT;
            typedef typename Traits<PT>::real_type RT;
            static inline PT call(const T1& x, const T2& y)
            { return ZProd<conj1,!conj2>::prod(x,y) / std::norm(y); }
        };
        template <class T1, class T2>
        struct makequot<false,T1,T2>
        {
            typedef typename Traits2<T1,T2>::type PT;
            static inline PT call(const T1& x, const T2& y)
            { return Maybe<conj1>::conj(x/y); }
        };

    };

    // A simpler structure that is mostly just to handing things like
    // int + complex<double>
    // which don't have a definition normally.
    struct ZSum
    {
        // complex + complex
        template <class T>
        static TMV_INLINE std::complex<T> sum(
            const std::complex<T>& x, const std::complex<T>& y)
        { return x + y; }
        template <class T1, class T2>
        static inline std::complex<typename Traits2<T1,T2>::type> sum(
            const std::complex<T1>& x, const std::complex<T2>& y)
        {
            typedef typename Traits2<T1,T2>::type T12;
            return std::complex<T12>(x.real()+y.real(),x.imag()+y.imag());
        }

        // real + complex
        template <class T>
        static TMV_INLINE std::complex<T> sum(const T& x, const std::complex<T>& y)
        { return x + y; }
        template <class T1, class T2>
        static inline std::complex<typename Traits2<T1,T2>::type> sum(
            const T1& x, const std::complex<T2>& y)
        {
            typedef typename Traits2<T1,T2>::type T12;
            return std::complex<T12>(x+y.real(),y.imag());
        }

        // complex + real
        template <class T>
        static TMV_INLINE std::complex<T> sum(const std::complex<T>& x, const T& y)
        { return x + y; }
        template <class T1, class T2>
        static inline std::complex<typename Traits2<T1,T2>::type> sum(
            const std::complex<T1>& x, const T2& y)
        {
            typedef typename Traits2<T1,T2>::type T12;
            return std::complex<T12>(x.real(),x.imag()+y);
        }

        // real + real
        template <class T>
        static TMV_INLINE T sum(const T& x, const T& y)
        { return x + y; }
        template <class T1, class T2>
        static inline typename Traits2<T1,T2>::type sum(
            const T1& x, const T2& y)
        { return x + y; }
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
        template <class T1, class T2>
        static TMV_INLINE const T1& select(const T1& x, const T2& ) 
        { return x; }
        template <class T1, class T2>
        static TMV_INLINE T1& select_ref(T1& x, const T2& ) { return x; }

        // abs(x) or x
        template <class T>
        static TMV_INLINE typename Traits<T>::real_type abs(const T& x) 
        { return TMV_ABS(x); }

        // real(x) or x
        template <class T>
        static TMV_INLINE typename Traits<T>::real_type real(const T& x) 
        { return TMV_REAL(x); }

        // conj(x) or x
        template <class T>
        static TMV_INLINE T conj(const T& x) { return TMV_CONJ(x); }

        // x = conj(x) or nothing
        template <class T>
        static TMV_INLINE void conjval(T& x) { x = TMV_CONJ(x); }

        // -x or x
        template <class T>
        static TMV_INLINE T neg(const T& x) { return -x; }

        // x<y or x>y
        template <class T>
        static TMV_INLINE bool less(const T& x, const T& y) { return x<y; }

        // 2*x or x
        template <class T>
        static TMV_INLINE T twox(const T& x)
        { return 2*x; }
        static TMV_INLINE int twox(const int& x)
        { return x>>1; }

        // ++x or nothing
        template <class T>
        static TMV_INLINE void increment(T& x) { ++x; }

        // --x or nothing
        template <class T>
        static TMV_INLINE void decrement(T& x) { --x; }

        // 2*x or x  (intended for integer types)
        template <class T>
        static TMV_INLINE T multby2(const T& x) { return x<<1; }
        // x/2 or x  (intended for integer types)
        template <class T>
        static TMV_INLINE T divby2(const T& x) { return x>>1; }

        // x += y or x = y
        template <class T1, class T2>
        static TMV_INLINE void add(T1& x, const T2& y) { x += y; }
        template <class T1, class T2>
        static TMV_INLINE void add(ConjRef<T1> x, const T2& y) { x += y; }
        template <class T1, bool C, class T2>
        static TMV_INLINE void add(TriRef<T1,C> x, const T2& y) { x += y; }

        // x *= y or nothing
        template <class T1, class T2>
        static TMV_INLINE void scale(T1& x, const T2& y) { x *= y; }
        template <class T1, class T2>
        static TMV_INLINE void scale(ConjRef<T1> x, const T2& y) { x *= y; }
        template <class T1, bool C, class T2>
        static TMV_INLINE void scale(TriRef<T1,C> x, const T2& y) { x *= y; }

        // x /= y or nothing
        template <class T1, class T2>
        static TMV_INLINE void invscale(T1& x, const T2& y) { x /= y; }
        template <class T1, class T2>
        static TMV_INLINE void invscale(ConjRef<T1> x, const T2& y) { x /= y; }
        template <class T1, bool C, class T2>
        static TMV_INLINE void invscale(TriRef<T1,C> x, const T2& y) { x /= y; }

        // x + y or y
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type sum(
            const T1& x, const T2& y) 
        { return x + y; }

        template <class T1, int ix2, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type sum(
            const T1& x, const Scaling<ix2,T2>& y) 
        { return x + T2(y); }

        // x * y or y
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type prod(
            const T1& x, const T2& y) 
        { return ZProd<false,false>::prod(x , y); }

        template <class T1, int ix2, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type prod(
            const T1& x, const Scaling<ix2,T2>& y) 
        { return ZProd<false,false>::prod(x , y); }

        // x^-1 * y or y
        template <class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type invprod(
            const T1& x, const T2& y) 
        { return ZProd<false,false>::quot(y , x); }

        template <class T1, int ix2, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type invprod(
            const T1& x, const Scaling<ix2,T2>& y) 
        { return ZProd<false,false>::quot(y , x); }

        // x * y or y
        template <bool c1, bool c2, class T1, class T2>
        static TMV_INLINE typename Traits2<T1,T2>::type zprod(
            const T1& x, const T2& y)
        { return ZProd<c1,c2>::prod(x , y); }

        // x - y or x + y
        template <class T>
        static TMV_INLINE T diff(const T& x, const T& y) { return x - y; }

        // x = y or nothing
        template <class T1, class T2>
        static TMV_INLINE void set(T1& x, const T2& y) { x = y; }
        template <class T1, class T2>
        static TMV_INLINE void set(ConjRef<T1> x, const T2& y) { x = y; }
        template <class T1, bool C, class T2>
        static TMV_INLINE void set(TriRef<T1,C> x, const T2& y) { x = y; }

        // x.assignTo(y) or nothing
        template <class T1, class T2>
        static TMV_INLINE void assignTo(const T1& x, T2& y) { x.assignTo(y); }
        template <class T1, class T2>
        static TMV_INLINE void newAssignTo(const T1& x, T2& y) 
        { x.newAssignTo(y); }

        // m.ref(i,i) = x or nothing
        template <class M, class T>
        static TMV_INLINE void setdiag(M& m, int i, const T& x) 
        { m.ref(i,i) = x; }
        template <class M, class T>
        static TMV_INLINE void setdiag2(M m, int i, const T& x) 
        { m.ref(i,i) = x; }

        // m.diag() = x or nothing
        template <class M, class T>
        static TMV_INLINE void setdiag(M& m, const T& x) 
        { m.diag().setAllTo(x); }
        template <class M, class T>
        static TMV_INLINE void setdiag2(M m, const T& x) 
        { m.diag().setAllTo(x); }

        // m.setZero() or nothing
        template <class M>
        static TMV_INLINE void zero(M& m) { m.setZero(); }
        template <class M>
        static TMV_INLINE void zero2(M m) { m.setZero(); }

        // m.conjugateSelf() or nothing
        template <class M>
        static TMV_INLINE void conjself(M& m) { m.conjugateSelf(); }
        template <class M>
        static TMV_INLINE void conjself2(M m) { m.conjugateSelf(); }

        // m.conjugate() or m
        template <class M>
        static TMV_INLINE typename M::const_conjugate_type conjugate(const M& m) 
        { return m.conjugate(); }
        template <class M>
        static TMV_INLINE typename M::conjugate_type conjugate(M& m) 
        { return m.conjugate(); }

        // m.transpose() or m
        template <class M>
        static TMV_INLINE typename M::const_transpose_type transpose(const M& m) 
        { return m.transpose(); }
        template <class M>
        static TMV_INLINE typename M::transpose_type transpose(M& m) 
        { return m.transpose(); }

        // m.transpose() or m.view()
        template <class M>
        static TMV_INLINE typename M::const_transpose_type transposeview(
            const M& m) 
        { return m.transpose(); }
        template <class M>
        static TMV_INLINE typename M::transpose_type transposeview(M& m) 
        { return m.transpose(); }

        // v.addToAll(x) or v.setAllTo(x)
        template <class V, class T>
        static TMV_INLINE void addtoall(V& v, const T& x) { v.addToAll(x); }
        template <class V, class T>
        static TMV_INLINE void addtoall2(V v, const T& x) { v.addToAll(x); }

        // m.diag().setAllTo(1) or nothing
        template <class M>
        static TMV_INLINE void makeunitdiag(M& m) 
        { m.diag().setAllTo(typename M::value_type(1)); }
        template <class M>
        static TMV_INLINE void makeunitdiag2(M m) 
        { m.diag().setAllTo(typename M::value_type(1)); }

        // m.offDiag().setZero() or nothing
        template <class M>
        static TMV_INLINE void zero_offdiag(M& m) { m.offDiag().setZero(); }
        template <class M>
        static TMV_INLINE void zero_offdiag2(M m) { m.offDiag().setZero(); }

        // m.offDiag() or m.view()
        template <class M>
        static TMV_INLINE typename M::offdiag_type offdiag(M& m) 
        { return m.offDiag(); }
        template <class M>
        static TMV_INLINE typename M::offdiag_type offdiag2(M m) 
        { return m.offDiag(); }

        // m2.viewAsUnitDiag() or m2
        template <class M>
        static TMV_INLINE typename M::const_unitdiag_type unitview(const M& m) 
        { return m.viewAsUnitDiag(); }
        template <class M>
        static TMV_INLINE typename M::unitdiag_type unitview(M& m) 
        { return m.viewAsUnitDiag(); }

        // m2.UnitUpperTri() or m2.UpperTri()
        template <class M>
        static TMV_INLINE typename M::const_unit_uppertri_type unit_uppertri(
            const M& m) 
        { return m.unitUpperTri(); }
        template <class M>
        static TMV_INLINE typename M::unit_uppertri_type unit_uppertri(M& m) 
        { return m.unitUpperTri(); }
        template <class M>
        static TMV_INLINE typename M::unit_uppertri_type unit_uppertri2(M m) 
        { return m.unitUpperTri(); }

        // m2.UnitLowerTri() or m2.LowerTri()
        template <class M>
        static TMV_INLINE typename M::const_unit_lowertri_type unit_lowertri(
            const M& m) 
        { return m.unitLowerTri(); }
        template <class M>
        static TMV_INLINE typename M::unit_lowertri_type unit_lowertri(M& m) 
        { return m.unitLowerTri(); }
        template <class M>
        static TMV_INLINE typename M::unit_lowertri_type unit_lowertri2(M m) 
        { return m.unitLowerTri(); }

        // m2.UpperTri() or m2.LowerTri()
        template <class M>
        static TMV_INLINE typename M::const_uppertri_type uppertri(const M& m) 
        { return m.upperTri(); }
        template <class M>
        static TMV_INLINE typename M::uppertri_type uppertri(M& m) 
        { return m.upperTri(); }
        template <class M>
        static TMV_INLINE typename M::uppertri_type uppertri2(M m) 
        { return m.upperTri(); }
        template <class M>
        static TMV_INLINE typename M::const_unknown_uppertri_type uppertri(
            const M& m, int dt) 
        { return m.upperTri(dt); }
        template <class M>
        static TMV_INLINE typename M::unknown_uppertri_type uppertri(
            M& m, DiagType dt) 
        { return m.upperTri(dt); }
        template <class M>
        static TMV_INLINE typename M::unknown_uppertri_type uppertri2(
            M m, DiagType dt) 
        { return m.upperTri(dt); }

#ifdef __SSE__
        // _mm_load_ps or _mm_set_ps
        static TMV_INLINE void sse_load(
            __m128& m,
            const float* x, const float*, const float*, const float*)
        { m = _mm_load_ps(x); }

        // _mm_loadu_ps or _mm_set_ps
        static TMV_INLINE void sse_loadu(
            __m128& m,
            const float* x, const float*, const float*, const float*)
        { m = _mm_loadu_ps(x); }

        // _mm_store_ps or cast and assign
        static TMV_INLINE void sse_store(
            float* x, float*, float*, float*, const __m128& m)
        { _mm_store_ps(x,m); }

        // _mm_storeu_ps or cast and assign
        static TMV_INLINE void sse_storeu(
            float* x, float*, float*, float*, const __m128& m)
        { _mm_storeu_ps(x,m); }

        // _mm_add_ps or cast and add
        static TMV_INLINE void sse_add(
            float* x, float*, float*, float*, const __m128& m)
        { _mm_store_ps(x,_mm_add_ps(_mm_load_ps(x),m)); }
        static TMV_INLINE void sse_addu(
            float* x, float*, float*, float*, const __m128& m)
        { _mm_storeu_ps(x,_mm_add_ps(_mm_loadu_ps(x),m)); }

        // _mm_load_ps or two _mm_load_pi's
        static TMV_INLINE void sse_load(
            __m128& m,
            const std::complex<float>* x, const std::complex<float>* )
        { m = _mm_load_ps( (const float*) x ); }

        // _mm_loadu_ps or two _mm_load?_pi's
        static TMV_INLINE void sse_loadu(
            __m128& m,
            const std::complex<float>* x, const std::complex<float>* )
        { m = _mm_loadu_ps( (const float*) x ); }

        // _mm_store_ps or two _mm_store?_pi's
        static TMV_INLINE void sse_store(
            std::complex<float>* x, std::complex<float>* , const __m128& m)
        { _mm_store_ps((float*)(x) , m); }

        // _mm_storeu_ps or two _mm_store?_pi's
        static TMV_INLINE void sse_storeu(
            std::complex<float>* x, std::complex<float>* , const __m128& m)
        { _mm_storeu_ps((float*)(x) , m); }

        // _mm_add_ps or cast and add
        static TMV_INLINE void sse_add(
            std::complex<float>* x, std::complex<float>* , const __m128& m)
        { _mm_store_ps((float*)x,_mm_add_ps(_mm_load_ps((float*)x),m)); }
        static TMV_INLINE void sse_addu(
            std::complex<float>* x, std::complex<float>* , const __m128& m)
        { _mm_storeu_ps((float*)x,_mm_add_ps(_mm_loadu_ps((float*)x),m)); }

        // A = _mm_mul_ps(x,A) or nothing
        static TMV_INLINE void sse_mult(const __m128& x, __m128& A)
        { A = _mm_mul_ps(x,A); }
#endif
#ifdef __SSE2__
        static TMV_INLINE void sse_load(
            __m128d& m, const double* x, const double* )
        { m = _mm_load_pd(x); }
        static TMV_INLINE void sse_loadu(
            __m128d& m, const double* x, const double* )
        { m = _mm_loadu_pd(x); }

        static TMV_INLINE void sse_store(double* x, double*, const __m128d& m)
        { _mm_store_pd(x,m); }
        static TMV_INLINE void sse_storeu(double* x, double*, const __m128d& m)
        { _mm_storeu_pd(x,m); }

        static TMV_INLINE void sse_add(double* x, double*, const __m128d& m)
        { _mm_store_pd(x,_mm_add_pd(_mm_load_pd(x),m)); }
        static TMV_INLINE void sse_addu(double* x, double*, const __m128d& m)
        { _mm_storeu_pd(x,_mm_add_pd(_mm_loadu_pd(x),m)); }

        static TMV_INLINE void sse_load(
            __m128d& m, const std::complex<double>* x)
        { m = _mm_load_pd((const double*) x); }
        static TMV_INLINE void sse_loadu(
            __m128d& m, const std::complex<double>* x)
        { m = _mm_loadu_pd((const double*) x); }

        static TMV_INLINE void sse_store(
            std::complex<double>* x, const __m128d& m)
        { _mm_store_pd((double*) x,m); }
        static TMV_INLINE void sse_storeu(
            std::complex<double>* x, const __m128d& m)
        { _mm_storeu_pd((double*) x,m); }

        static TMV_INLINE void sse_add(
            std::complex<double>* x, const __m128d& m)
        { _mm_store_pd((double*)x,_mm_add_pd(_mm_load_pd((double*)x),m)); }
        static TMV_INLINE void sse_addu(
            std::complex<double>* x, const __m128d& m)
        { _mm_storeu_pd((double*)x,_mm_add_pd(_mm_loadu_pd((double*)x),m)); }

        static TMV_INLINE void sse_mult(const __m128d& x, __m128d& A)
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

        template <class T1, class T2>
        static TMV_INLINE const T2& select(const T1& , const T2& y) 
        { return y; }
        template <class T1, class T2>
        static TMV_INLINE const T2& select_ref(T1& , const T2& y) 
        { return y; }

        template <class T>
        static TMV_INLINE T abs(const T& x) { return x; }

        template <class T>
        static TMV_INLINE T real(const T& x) { return x; }

        template <class T>
        static TMV_INLINE T conj(const T& x) { return x; }

        template <class T>
        static TMV_INLINE void conjval(T& ) { }

        template <class T>
        static TMV_INLINE T neg(const T& x) { return x; }

        template <class T>
        static TMV_INLINE bool less(const T& x, const T& y) { return x>y; }

        template <class T>
        static TMV_INLINE T twox(const T& x) { return x; }

        template <class T>
        static TMV_INLINE void increment(const T& ) { }

        template <class T>
        static TMV_INLINE void decrement(const T& ) { }

        template <class T>
        static TMV_INLINE const T& multby2(const T& x) { return x; }
        template <class T>
        static TMV_INLINE const T& divby2(const T& x) { return x; }

        template <class T1, class T2>
        static TMV_INLINE void add(T1& x, const T2& y) { x = y; }
        template <class T1, class T2>
        static TMV_INLINE void add(ConjRef<T1> x, const T2& y) { x = y; }
        template <class T1, bool C, class T2>
        static TMV_INLINE void add(TriRef<T1,C> x, const T2& y) { x = y; }

        template <class T1, class T2>
        static TMV_INLINE void scale(T1& , const T2& ) { }
        template <class T1, class T2>
        static TMV_INLINE void scale(ConjRef<T1> , const T2& ) { }
        template <class T1, bool C, class T2>
        static TMV_INLINE void scale(TriRef<T1,C> , const T2& ) { }

        template <class T1, class T2>
        static TMV_INLINE void invscale(T1& , const T2& ) { }
        template <class T1, class T2>
        static TMV_INLINE void invscale(ConjRef<T1> , const T2& ) { }
        template <class T1, bool C, class T2>
        static TMV_INLINE void invscale(TriRef<T1,C> , const T2& ) { }

        template <class T1, class T2>
        static TMV_INLINE T2 sum(const T1& , const T2& y) { return y; }

        template <class T1, class T2>
        static TMV_INLINE T2 prod(const T1& , const T2& y) { return y; }

        template <class T1, int ix2, class T2>
        static TMV_INLINE Scaling<ix2,T2> prod(
            const T1& , const Scaling<ix2,T2>& y) { return y; }

        template <class T1, class T2>
        static TMV_INLINE T2 invprod(const T1& , const T2& y)  { return y; }

        template <class T1, int ix2, class T2>
        static TMV_INLINE Scaling<ix2,T2> invprod(
            const T1& , const Scaling<ix2,T2>& y) { return y; }

        template <bool c1, bool c2, class T1, class T2>
        static TMV_INLINE T2 zprod(const T1& , const T2& y)
        { return Maybe<c2>::conj(y); }

        template <class T>
        static TMV_INLINE T diff(const T& x, const T& y) { return x + y; }

        template <class T1, class T2>
        static TMV_INLINE void set(T1& , const T2& ) { }
        template <class T1, class T2>
        static TMV_INLINE void set(ConjRef<T1> , const T2& ) { }
        template <class T1, bool C, class T2>
        static TMV_INLINE void set(TriRef<T1,C> , const T2& ) { }

        template <class T1, class T2>
        static TMV_INLINE void assignTo(const T1& , T2& ) { }
        template <class T1, class T2>
        static TMV_INLINE void newAssignTo(const T1& , T2& )  { }

        template <class M, class T>
        static TMV_INLINE void setdiag(M& , int , const T& ) { }
        template <class M, class T>
        static TMV_INLINE void setdiag2(M , int , const T& ) { }

        template <class M, class T>
        static TMV_INLINE void setdiag(M& , const T& ) { }
        template <class M, class T>
        static TMV_INLINE void setdiag2(M , const T& ) { }

        template <class M>
        static TMV_INLINE void zero(M& ) { }
        template <class M>
        static TMV_INLINE void zero2(M ) { }

        template <class M>
        static TMV_INLINE void conjself(M& ) { }
        template <class M>
        static TMV_INLINE void conjself2(M ) { }

        template <class M>
        static TMV_INLINE const M& conjugate(const M& m) { return m; }
        template <class M>
        static TMV_INLINE M& conjugate(M& m) { return m; }

        template <class M>
        static TMV_INLINE const M& transpose(const M& m) { return m; }
        template <class M>
        static TMV_INLINE M& transpose(M& m) { return m; }

        template <class M>
        static TMV_INLINE typename M::const_view_type transposeview(const M& m) 
        { return m.view(); }
        template <class M>
        static TMV_INLINE typename M::view_type transposeview(M& m) 
        { return m.view(); }

        template <class V, class T>
        static TMV_INLINE void addtoall(V& v, const T& x) { v.setAllTo(x); }
        template <class V, class T>
        static TMV_INLINE void addtoall2(V v, const T& x) { v.setAllTo(x); }

        template <class M>
        static TMV_INLINE void makeunitdiag(M& ) { }
        template <class M>
        static TMV_INLINE void makeunitdiag2(M ) { }

        template <class M1, class M2>
        static TMV_INLINE void offdiag_copy(const M1& m1, M2& m2) { m2 = m1; }
        template <class M1, class M2>
        static TMV_INLINE void offdiag_copy2(const M1 m1, M2& m2) { m2 = m1; }

        template <class M>
        static TMV_INLINE void zero_offdiag(M& ) { }
        template <class M>
        static TMV_INLINE void zero_offdiag2(M ) { }

        template <class M>
        static TMV_INLINE typename M::view_type offdiag(M& m) 
        { return m.view(); }
        template <class M>
        static TMV_INLINE typename M::view_type offdiag2(M m) 
        { return m.view(); }

        template <class M>
        static TMV_INLINE const M& unitview(const M& m) 
        { return m; }
        template <class M>
        static TMV_INLINE M& unitview(M& m) 
        { return m; }

        template <class M>
        static TMV_INLINE typename M::const_uppertri_type unit_uppertri(
            const M& m) 
        { return m.upperTri(); }
        template <class M>
        static TMV_INLINE typename M::uppertri_type unit_uppertri(M& m) 
        { return m.upperTri(); }
        template <class M>
        static TMV_INLINE typename M::uppertri_type unit_uppertri2(M m) 
        { return m.upperTri(); }

        template <class M>
        static TMV_INLINE typename M::const_lowertri_type unit_lowertri(
            const M& m) 
        { return m.lowerTri(); }
        template <class M>
        static TMV_INLINE typename M::lowertri_type unit_lowertri(M& m) 
        { return m.lowerTri(); }
        template <class M>
        static TMV_INLINE typename M::lowertri_type unit_lowertri2(M m) 
        { return m.lowerTri(); }

        template <class M>
        static TMV_INLINE typename M::const_lowertri_type uppertri(const M& m) 
        { return m.lowerTri(); }
        template <class M>
        static TMV_INLINE typename M::lowertri_type uppertri(M& m) 
        { return m.lowerTri(); }
        template <class M>
        static TMV_INLINE typename M::lowertri_type uppertri2(M m) 
        { return m.lowerTri(); }
        template <class M>
        static TMV_INLINE typename M::const_unknown_lowertri_type uppertri(
            const M& m, int dt) 
        { return m.lowerTri(dt); }
        template <class M>
        static TMV_INLINE typename M::unknown_lowertri_type uppertri(
            M& m, DiagType dt) 
        { return m.lowerTri(dt); }
        template <class M>
        static TMV_INLINE typename M::unknown_lowertri_type uppertri2(
            M m, DiagType dt) 
        { return m.lowerTri(dt); }

#ifdef __SSE__
        static TMV_INLINE void sse_load(
            __m128& m,
            const float* x1, const float* x2, const float* x3, const float* x4)
        { m = _mm_set_ps(*x4,*x3,*x2,*x1); }
        static TMV_INLINE void sse_loadu(
            __m128& m,
            const float* x1, const float* x2, const float* x3, const float* x4)
        { m = _mm_set_ps(*x4,*x3,*x2,*x1); }

        static TMV_INLINE void sse_store(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 = mf[0]; *x2 = mf[1]; *x3 = mf[2]; *x4 = mf[3];
        }
        static TMV_INLINE void sse_storeu(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 = mf[0]; *x2 = mf[1]; *x3 = mf[2]; *x4 = mf[3];
        }

        static TMV_INLINE void sse_add(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 += mf[0]; *x2 += mf[1]; *x3 += mf[2]; *x4 += mf[3];
        }
        static TMV_INLINE void sse_addu(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 += mf[0]; *x2 += mf[1]; *x3 += mf[2]; *x4 += mf[3];
        }

        static TMV_INLINE void sse_load(
            __m128& m,
            const std::complex<float>* x1, const std::complex<float>* x2)
        {
            // The junk variable is to avoid a warning about m being
            // used uninitialized.
            const __m128 junk = _mm_set1_ps(0.F);
            m = _mm_loadl_pi(junk,(const __m64*)x1);
            m = _mm_loadh_pi(m,(const __m64*)x2);
        }
        static TMV_INLINE void sse_loadu(
            __m128& m,
            const std::complex<float>* x1, const std::complex<float>* x2)
        {
            const __m128 junk = _mm_set1_ps(0.F);
            m = _mm_loadl_pi(junk,(const __m64*)x1);
            m = _mm_loadh_pi(m,(const __m64*)x2);
        }

        static TMV_INLINE void sse_store(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        {
            _mm_storel_pi((__m64*)x1,m);
            _mm_storeh_pi((__m64*)x2,m);
        }
        static TMV_INLINE void sse_storeu(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        {
            _mm_storel_pi((__m64*)x1,m);
            _mm_storeh_pi((__m64*)x2,m);
        }

        static TMV_INLINE void sse_add(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 += std::complex<float>(mf[0],mf[1]);
            *x2 += std::complex<float>(mf[2],mf[3]);
        }
        static TMV_INLINE void sse_addu(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        {
            const float* mf= (const float*)(&m);
            *x1 += std::complex<float>(mf[0],mf[1]);
            *x2 += std::complex<float>(mf[2],mf[3]);
        }

        static TMV_INLINE void sse_mult(const __m128& , __m128& ) {}
#endif
#ifdef __SSE2__
        static TMV_INLINE void sse_load(
            __m128d& m, const double* x1, const double* x2)
        { m = _mm_set_pd(*x2,*x1); }
        static TMV_INLINE void sse_loadu(
            __m128d& m, const double* x1, const double* x2)
        { m = _mm_set_pd(*x2,*x1); }

        static TMV_INLINE void sse_store(
            double* x1, double* x2, const __m128d& m)
        {
            const double* md= (const double*)(&m);
            *x1 = md[0]; *x2 = md[1];
        }
        static TMV_INLINE void sse_storeu(
            double* x1, double* x2, const __m128d& m)
        {
            const double* mf= (const double*)(&m);
            *x1 = mf[0]; *x2 = mf[1];
        }

        static TMV_INLINE void sse_add(
            double* x1, double* x2, const __m128d& m)
        {
            const double* md= (const double*)(&m);
            *x1 += md[0]; *x2 += md[1];
        }
        static TMV_INLINE void sse_addu(
            double* x1, double* x2, const __m128d& m)
        {
            const double* mf= (const double*)(&m);
            *x1 += mf[0]; *x2 += mf[1];
        }

        // These are the same for Maybe<true> or Maybe<false>.
        // I include them to make it easier to write the code that uses
        // these loads and stores, since it makes the syntax more similar
        // between all the different varieties (float/double, real/complex).
        static TMV_INLINE void sse_load(
            __m128d& m, const std::complex<double>* x)
        { m = _mm_load_pd((const double*)x); }
        static TMV_INLINE void sse_loadu(
            __m128d& m, const std::complex<double>* x)
        { m = _mm_loadu_pd((const double*)x); }

        static TMV_INLINE void sse_store(
            std::complex<double>* x, const __m128d& m)
        { _mm_store_pd((double*)x,m); }
        static TMV_INLINE void sse_storeu(
            std::complex<double>* x, const __m128d& m)
        { _mm_storeu_pd((double*)x,m); }

        static TMV_INLINE void sse_add(
            std::complex<double>* x, const __m128d& m)
        { _mm_store_pd((double*)x,_mm_add_pd(_mm_load_pd((double*)x),m)); }
        static TMV_INLINE void sse_addu(
            std::complex<double>* x, const __m128d& m)
        { _mm_storeu_pd((double*)x,_mm_add_pd(_mm_loadu_pd((double*)x),m)); }

        static TMV_INLINE void sse_mult(const __m128d& , __m128d& ) {}
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
        static TMV_INLINE void add(T1& x, const T2& y) 
        { Maybe<yn2>::add(x,y); }
        template <class T1, class T2>
        static TMV_INLINE void add(ConjRef<T1> x, const T2& y) 
        { Maybe<yn2>::add(x,y); }
        template <class T1, bool C, class T2>
        static TMV_INLINE void add(TriRef<T1,C> x, const T2& y) 
        { Maybe<yn2>::add(x,y); }

        // Maybe<unit>::unit_uppertri or Maybe<unit>::unit_lowertri
        template <class M>
        static TMV_INLINE typename TypeSelect<yn2,
                          typename M::const_unit_uppertri_type,
                          typename M::const_uppertri_type>::type uppertri(
                              const M& m) 
        { return Maybe<yn2>::unit_uppertri(m); }
        template <class M>
        static TMV_INLINE typename TypeSelect<yn2,
                          typename M::unit_uppertri_type,
                          typename M::uppertri_type>::type uppertri(M& m) 
        { return Maybe<yn2>::unit_uppertri(m); }

        template <class M>
        static TMV_INLINE typename TypeSelect<yn2,
                          typename M::unit_uppertri_type,
                          typename M::uppertri_type>::type uppertri2(M m) 
        { return Maybe<yn2>::unit_uppertri(m); }

#ifdef __SSE__
        // Maybe<unit>::sse_load or Maybe<unit>::sse_loadu
        static TMV_INLINE void sse_load(
            __m128& m,
            const float* x1, const float* x2, const float* x3, const float* x4)
        { Maybe<yn2>::sse_load(m,x1,x2,x3,x4); }
        static TMV_INLINE void sse_load(
            __m128& m,
            const std::complex<float>* x1, const std::complex<float>* x2)
        { Maybe<yn2>::sse_load(m,x1,x2); }

        // Maybe<unit>::sse_add or Maybe<unit>::sse_store
        static TMV_INLINE void sse_add(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        { Maybe<yn2>::sse_add(x1,x2,x3,x4,m); }
        static TMV_INLINE void sse_add(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        { Maybe<yn2>::sse_add(x1,x2,m); }

        // Maybe<unit>::sse_addu or Maybe<unit>::sse_storeu
        static TMV_INLINE void sse_addu(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        { Maybe<yn2>::sse_addu(x1,x2,x3,x4,m); }
        static TMV_INLINE void sse_addu(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        { Maybe<yn2>::sse_addu(x1,x2,m); }
#endif

#ifdef __SSE2__
        static TMV_INLINE void sse_load(
            __m128d& m, const double* x1, const double* x2)
        { Maybe<yn2>::sse_load(m,x1,x2); }

        static TMV_INLINE void sse_add(
            double* x1, double* x2, const __m128d& m)
        { Maybe<yn2>::sse_add(x1,x2,m); }
        static TMV_INLINE void sse_add(
            std::complex<double>* x, const __m128d& m)
        { Maybe<yn2>::sse_add(x,m); }

        static TMV_INLINE void sse_addu(
            double* x1, double* x2, const __m128d& m)
        { Maybe<yn2>::sse_addu(x1,x2,m); }
        static TMV_INLINE void sse_addu(
            std::complex<double>* x, const __m128d& m)
        { Maybe<yn2>::sse_addu(x,m); }
#endif
    };

    template <bool yn2>
    struct Maybe2<false,yn2>
    {
        template <class T1, class T2>
        static TMV_INLINE void add(T1& , const T2& ) {}
        template <class T1, class T2>
        static TMV_INLINE void add(ConjRef<T1> , const T2& ) {}
        template <class T1, bool C, class T2>
        static TMV_INLINE void add(TriRef<T1,C> , const T2& ) {}

        template <class M>
        static TMV_INLINE 
        typename TypeSelect<yn2,
                          typename M::const_unit_lowertri_type,
                          typename M::const_lowertri_type>::type uppertri(
                              const M& m) 
        { return Maybe<yn2>::unit_lowertri(m); }
        template <class M>
        static TMV_INLINE typename TypeSelect<yn2,
                          typename M::unit_lowertri_type,
                          typename M::lowertri_type>::type uppertri(M& m) 
        { return Maybe<yn2>::unit_lowertri(m); }
        template <class M>
        static TMV_INLINE typename TypeSelect<yn2,
                          typename M::unit_lowertri_type,
                          typename M::lowertri_type>::type uppertri2(M m) 
        { return Maybe<yn2>::unit_lowertri(m); }

#ifdef __SSE2__
        static TMV_INLINE void sse_load(
            __m128& m,
            const float* x1, const float* x2, const float* x3, const float* x4)
        { Maybe<yn2>::sse_loadu(m,x1,x2,x3,x4); }
        static TMV_INLINE void sse_load(
            __m128& m,
            const std::complex<float>* x1, const std::complex<float>* x2)
        { Maybe<yn2>::sse_loadu(m,x1,x2); }

        static TMV_INLINE void sse_add(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        { Maybe<yn2>::sse_store(x1,x2,x3,x4,m); }
        static TMV_INLINE void sse_add(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        { Maybe<yn2>::sse_store(x1,x2,m); }

        static TMV_INLINE void sse_addu(
            float* x1, float* x2, float* x3, float* x4, const __m128& m)
        { Maybe<yn2>::sse_storeu(x1,x2,x3,x4,m); }
        static TMV_INLINE void sse_addu(
            std::complex<float>* x1, std::complex<float>* x2, const __m128& m)
        { Maybe<yn2>::sse_storeu(x1,x2,m); }
#endif

#ifdef __SSE2__
        static TMV_INLINE void sse_load(
            __m128d& m, const double* x1, const double* x2)
        { Maybe<yn2>::sse_loadu(m,x1,x2); }

        static TMV_INLINE void sse_add(
            double* x1, double* x2, const __m128d& m)
        { Maybe<yn2>::sse_store(x1,x2,m); }
        static TMV_INLINE void sse_add(
            std::complex<double>* x, const __m128d& m)
        { Maybe<yn2>::sse_store(x,m); }

        static TMV_INLINE void sse_addu(
            double* x1, double* x2, const __m128d& m)
        { Maybe<yn2>::sse_storeu(x1,x2,m); }
        static TMV_INLINE void sse_addu(
            std::complex<double>* x, const __m128d& m)
        { Maybe<yn2>::sse_storeu(x,m); }
#endif
    };

#define TMV_MAYBE_CREF(T1,T2) \
    typename TypeSelect<Traits2<T1,T2>::sametype,const T1&,T2>::type
#define TMV_MAYBE_REF(T1,T2) \
    typename TypeSelect<Traits2<T1,T2>::sametype,T1&,T2>::type

#ifndef TMV_NO_THROW
    class Error : 
        public std::runtime_error
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

    class FailedAssert : 
        public Error
    {
    public :
        std::string failed_assert;
        unsigned long line;
        std::string file;

        FailedAssert(
            std::string s, unsigned long l, std::string f) throw() :
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
    class ReadError : 
        public Error
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

    class Singular : 
        public Error
    {
    public :
        Singular() throw() :
            Error("Encountered singular matrix.") {}
        Singular(std::string s) throw() :
            Error("Encountered singular ",s) {}
        ~Singular() throw() {}
    };

    class NonPosDef : 
        public Error
    {
    public:
        NonPosDef() throw() :
            Error("Encountered invalid non-positive-definite matrix.") {}
        NonPosDef(std::string s) throw() :
            Error("Encountered invalid non-positive-definite matrix in ",s) {}
        ~NonPosDef() throw() {}
    };
#endif

    inline void ThrowSingular(std::string s)
    {
#ifdef TMV_NO_THROW
        std::cerr<<"Encountered singular "<<s<<std::endl;
        exit(1);
#else
        throw Singular(s);
#endif
    }

    inline void ThrowNonPosDef(std::string s)
    {
#ifdef TMV_NO_THROW
        std::cerr<<"Encountered invalid non-positive-definite matrix in "<<
            s<<std::endl;
        exit(1);
#else
        throw NonPosDef(s);
#endif
    }


#ifdef TMV_WARN
    // defined in TMV_Vector.cpp
    // initialized with &std::cout;
    extern std::ostream* warn_out;
    inline void TMV_Warning(std::string s)
    {
        if (warn_out) {
            *warn_out << "Warning:\n" << s << std::endl;
        }
    }

    inline std::ostream* WriteWarningsTo(std::ostream* newos)
    { std::ostream* temp = warn_out; warn_out = newos; return temp; }

    inline void NoWarnings()
    { warn_out = 0; }
#else
    TMV_INLINE void TMV_Warning(std::string ) {}
    TMV_INLINE std::ostream* WriteWarningsTo(std::ostream* os) { return 0; }
    TMV_INLINE void NoWarnings() {}
#endif

    // A helper structure that acts like an int,
    // but only bothers to make the integer if S == TMV_UNKNOWN.
    // It also checks the constructor if S != TMV_UNKNOWN
    template <int S>
    struct CheckedInt
    {
        TMV_INLINE CheckedInt(int s) { 
#ifdef TMV_DEBUG
            if (s != S) {
                std::cerr<<"Mismatched CheckInt:\n";
                std::cerr<<"template parameter S = "<<S<<std::endl;
                std::cerr<<"argument s = "<<s<<std::endl;
            }
#endif
            TMVAssert(s == S); 
        }
        TMV_INLINE operator int() const { return S; }
    };
    template <>
    struct CheckedInt<TMV_UNKNOWN>
    {
        int step;
        TMV_INLINE CheckedInt(int s) : step(s) {}
        TMV_INLINE operator int() const 
        {
#ifdef TMV_DEBUG
            TMVAssert(step != -987234);
#endif
            return step; 
        }
        TMV_INLINE ~CheckedInt() {
#ifdef TMV_DEBUG
            step = -987234;
#endif
        }
    };

    // A similar helper to account for possibly unknown DiagType
    template <int A>
    struct DiagInt
    {
        TMV_INLINE DiagInt(int d) { TMVAssert(d == A); }
        TMV_INLINE operator int() const { return A; }
    };
    template <>
    struct DiagInt<UnknownDiag>
    {
        int dt;
        TMV_INLINE DiagInt(int d) : dt(d) {}
        TMV_INLINE operator int() const { return dt; }
    };

    template <class T>
    TMV_INLINE typename Traits<T>::real_type TMV_Epsilon() 
    { return std::numeric_limits<typename Traits<T>::real_type>::epsilon(); }

    template <class T>
    inline bool TMV_Underflow(T x)
    {
        typedef typename Traits<T>::real_type RT;
        return TMV_ABS2(x) < 
            std::numeric_limits<RT>::min() * RT(tmv::Traits<T>::twoifcomplex); 
    }

    TMV_INLINE bool TMV_Underflow(int )
    { return false; }

#ifdef TMV_DEBUG
#define TMV_TEXT
#endif

#ifdef TMV_TEXT
    template <class T>
    static inline std::string TMV_Text(const T&)
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
    static inline std::string TMV_Text(std::complex<T>)
    { return std::string("complex<") + TMV_Text(T()) + ">"; }

    static inline std::string TMV_Text(DivType d)
    {
        return 
            d==LU ? "LU" :
            d==CH ? "CH" :
            d==QR ? "QR" :
            d==QRP ? "QRP" :
            d==SV ? "SV" : "XX";
    }

    static inline std::string TMV_Text(ConjType c)
    { return c==Conj ? "Conj" : "NonConj"; }

    static inline std::string TMV_Text(IndexStyle i)
    { return i==CStyle ? "CStyle" : "FortranStyle"; }

    static inline std::string TMV_Text(StepType s)
    { return s==Unit ? "Unit" : "NonUnit"; }

    static inline std::string TMV_Text(StorageType s)
    { 
        return 
            s==ColMajor ? "ColMajor" : 
            s==RowMajor ? "RowMajor" :
            s==DiagMajor ? "DiagMajor" :
            "NonMajor";
    }

    static inline std::string TMV_Text(DiagType d)
    { 
        return 
            d==NonUnitDiag ? "NonUnitDiag" :
            d==UnitDiag ? "UnitDiag" :
            d==ZeroDiag ? "ZeroDiag" :
            "UnknownDiag";
    }

    static inline std::string TMV_Text(PackType p)
    { return p==Packed ? "Packed" : "NonPacked"; }

    static inline std::string TMV_Text(DivStatus d)
    { return d==WithDivider ? "WithDivider" : "NoDivider"; }

    static inline std::string TMV_Text(AliasStatus a)
    { return a==CheckAlias ? "CheckAlias" : "NoAlias"; }

    static inline std::string TMV_Text(UpLoType u)
    { return u==Upper ? "Upper" : "Lower"; }
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

    // This bit is to workaround a bug in pgCC that was fixed in version 7.
    // I don't know if versions earlier than 6.1 had the bug, but 
    // I apply the workaround to all version before 7.
    template <class T>
    TMV_INLINE const T& Value(const T& x) { return x; }
#ifdef PLATFORM_COMPILER_PGI
#if PLATFORM_COMPILER_VERSION < 0x070000
    TMV_INLINE double Value(long double x) { return double(x); }
    TMV_INLINE std::complex<double> Value(std::complex<long double> x)
    { return tmv::Traits<std::complex<double> >::convert(x); }
#endif
#endif


} // namespace tmv

#endif
