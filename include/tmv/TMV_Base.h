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
//  v   SmallVector
//  m   SmallMatrix
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
// TODO  BlockDiagMatrix and Sym varieties
// TODO  SparseMatrix and Sym, Block varieties


#ifndef TMV_Base_H
#define TMV_Base_H

#if ((!defined(NDEBUG) && !defined(TMV_NDEBUG)) || defined(TMV_EXTRA_DEBUG))
#define TMV_DEBUG
#endif

#ifndef TMV_NO_WARN
#define TMV_WARN
#endif

#include <iosfwd>
#include <limits>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <cstddef>
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

#ifdef TMV_MEM_DEBUG
#include "tmv/mmgr.h"
#endif

#ifndef NO_INST_DOUBLE
#define INST_DOUBLE
#endif

#ifndef NO_INST_FLOAT
#define INST_FLOAT
#endif

#ifndef NO_INST_COMPLEX
#define INST_COMPLEX
#endif

namespace tmv {

    inline std::string TMV_Version() { return "0.74"; }
#define TMV_MAJOR_VERSION 0
#define TMV_MINOR_VERSION 74
#define TMV_VERSION_AT_LEAST(major,minor) \
    ( (TMV_MAJOR_VERSION > major) || \
      (TMV_MAJOR_VERSION == major && TMV_MINOR_VERSION >= minor) )

    // Attributes are stored as a binary field, so they can be |'ed
    // together to calculate the full set of attributes.
    //
    // For each object we have a default set of attributes that is
    // implied by A=0.  This way we can write, for example:
    // UpperTriMatrix<T,UnitDiag> U;
    // and this will imply (UnitDiag | ColMajor | CStyle).
    //
    // Attributes that are always the default (whenever they are allowed)
    // can be defined as 0, since the absence of the converse is sufficient.
    // e.g. CStyle, ColMajor, Lower

    enum IndexStyle { CStyle = 0, FortranStyle = 0x1 };
    enum StorageType {
        ColMajor = 0, RowMajor = 0x2, DiagMajor = 0x4,
        AllStorageType = 0x6 };
    enum DiagType { NonUnitDiag = 0, UnitDiag = 0x8 };
    enum UpLoType { Lower = 0, Upper = 0x10 };

    template <int A>
    struct Attrib
    {
        enum { onestoragetype = !( !!(A & RowMajor) && !!(A & DiagMajor) ) };
        enum { viewok = A <= 0x1 };
        enum { vectorok = A <= 0x1 };
        enum { matrixok = A <= 0x3 };
        enum { diagmatrixok = A <= 0x1 };
        enum { trimatrixok = ( A <= 0xb && !(A & DiagMajor) ) };
        enum { bandmatrixok = A <= 0x5 };
        enum { symmatrixok = (
                A <= 0x13 && !(A & DiagMajor) && !(A & UnitDiag) ) };
        enum { symbandmatrixok = (
                A <= 0x15 && onestoragetype && !(A & UnitDiag) ) };

        enum { fort = !!(A & FortranStyle) };
        enum { rowmajor = !!(A & RowMajor) };
        enum { diagmajor = !!(A & DiagMajor) };
        enum { colmajor = !(rowmajor || diagmajor) };
        enum { unitdiag = !!(A & UnitDiag) };
        enum { nonunitdiag = !unitdiag };
        enum { upper = !!(A & Upper) };
        enum { lower = !upper };

        enum { stor = rowmajor ? RowMajor : diagmajor ? DiagMajor : ColMajor };

        static std::string text()
        {
            std::string ret =
                ( colmajor ? "ColMajor" :
                  rowmajor ? "RowMajor" : "DiagMajor" );
            ret += fort ? "|FortranStyle" : "";
            ret += unitdiag ? "|UnitDiag" : "";
            ret += upper ? "|Upper" : "";
            return ret;
        }
        static std::string vtext()
        { return fort ? "FortranStyle" : "CStyle"; }
    };

    enum ADType { Ascend, Descend };
    enum CompType { RealComp, AbsComp, ImagComp, ArgComp };
    enum StepType { Unit, Step };
    enum ConjType { NonConj, Conj };
    enum SymType { Sym, Herm };

#define TMV_UTransOf(U) (U==int(Upper) ? Lower : Upper)

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
    inline DivType operator|(DivType a, DivType b)
    {
        return static_cast<DivType>(
            static_cast<int>(a) | static_cast<int>(b));
    }
    inline DivType operator&(DivType a, DivType b)
    {
        return static_cast<DivType>(
            static_cast<int>(a) & static_cast<int>(b));
    }
    inline DivType& operator|=(DivType& a, DivType b)
    { a = (a|b); return a; }
    inline DivType& operator&=(DivType& a, DivType b)
    { a = (a&b); return a; }
    inline DivType operator~(DivType a)
    { return static_cast<DivType>(~static_cast<int>(a)); }

#ifdef NOALIASCHECK
    template <class M1, class M2>
    inline bool SameStorage(const M1& , const M2& )
    { return false; }
#else
    template <class M1, class M2>
    inline bool SameStorage(const M1& m1, const M2& m2)
    {
        return
            static_cast<const void*>(m1.realPart().cptr()) ==
            static_cast<const void*>(m2.realPart().cptr());
    }
#endif

    template <typename T>
    inline T TMV_SQR(T x)
    { return x*x; }

    template <typename T>
    inline T TMV_SQRT(T x)
    { return T(std::sqrt(x)); }

    template <typename T>
    inline T TMV_EXP(T x)
    { return T(std::exp(x)); }

    template <typename T>
    inline T TMV_LOG(T x)
    { return T(std::log(x)); }

    template <typename T>
    inline T TMV_NORM(T x)
    { return x*x; }

    template <typename T>
    inline T TMV_NORM(std::complex<T> x)
    { return std::norm(x); }

    template <typename T>
    inline T TMV_CONJ(T x)
    { return x; }

    template <typename T>
    inline std::complex<T> TMV_CONJ(std::complex<T> x)
    { return std::conj(x); }

    template <typename T>
    inline T TMV_REAL(T x)
    { return x; }

    template <typename T>
    inline T TMV_REAL(std::complex<T> x)
    { return std::real(x); }

    template <typename T>
    inline T TMV_IMAG(T )
    { return T(0); }

    template <typename T>
    inline T TMV_IMAG(std::complex<T> x)
    { return std::imag(x); }

    template <typename T>
    inline T TMV_ARG(T x)
    { return x >= T(0) ? T(1) : T(-1); }

    template <typename T>
    inline T TMV_ARG(std::complex<T> x)
    { return arg(x); }

    template <typename T>
    inline T TMV_ABS(T x)
    { return std::abs(x); }

    template <typename T>
    inline T TMV_ABS(std::complex<T> x)
    {
        // This is the same as the usual std::abs algorithm.
        // However, I have come across implementations that don't
        // protext against overflow, so I just duplicate it here.

        T xr = x.real();
        T xi = x.imag();
        const T s = std::max(std::abs(xr),std::abs(xi));
        if (s == T(0)) return s; // Check for s == 0
        xr /= s;
        xi /= s;
        return s * std::sqrt(xr*xr + xi*xi);
    }

    template <typename T>
    inline T TMV_ABS2(T x)
    { return std::abs(x); }

    template <typename T>
    inline T TMV_ABS2(std::complex<T> x)
    { return std::abs(std::real(x)) + std::abs(std::imag(x)); }

    template <typename T>
    inline T TMV_MIN(T x, T y)
    { return x > y ? y : x; }

    template <typename T>
    inline T TMV_MAX(T x, T y)
    { return x > y ? x : y; }

    template <typename T>
    inline void TMV_SWAP(T& x, T& y)
    { T z = x; x = y; y = z; }

    template <typename T>
    struct Traits
    {
        enum { isreal = true };
        enum { iscomplex = false };
        enum { isinteger = std::numeric_limits<T>::is_integer };
        enum { twoifcomplex = 1 };

        typedef T real_type;
        typedef std::complex<T> complex_type;
    };

    template <typename T>
    struct Traits<std::complex<T> >
    {
        enum { isreal = false };
        enum { iscomplex = true };
        enum { isinteger = Traits<T>::isinteger };
        enum { twoifcomplex = 2 };

        typedef T real_type;
        typedef std::complex<T> complex_type;
    };

    template <typename T>
    struct Traits<T&> : public Traits<T> {};
    template <typename T>
    struct Traits<const T> : public Traits<T> {};
    template <typename T>
    struct Traits<const T&> : public Traits<T> {};

    template <typename T1, typename T2>
    struct Traits2
    {
        enum { sametype = false };
        enum { samebase = false };
        typedef T1 type;
    };
    template <typename T1, typename T2>
    struct Traits2<T1,std::complex<T2> >
    {
        enum { sametype = false };
        enum { samebase = false };
        typedef std::complex<typename Traits2<T1,T2>::type> type;
    };
    template <typename T1, typename T2>
    struct Traits2<std::complex<T1>,T2>
    {
        enum { sametype = false };
        enum { samebase = false };
        typedef std::complex<typename Traits2<T1,T2>::type> type;
    };
    template <typename T1, typename T2>
    struct Traits2<std::complex<T1>,std::complex<T2> >
    {
        enum { sametype = false };
        enum { samebase = false };
        typedef std::complex<typename Traits2<T1,T2>::type> type;
    };
    template <typename T>
    struct Traits2<T,T>
    {
        enum { sametype = true };
        enum { samebase = true };
        typedef T type;
    };
    template <typename T>
    struct Traits2<std::complex<T>,T>
    {
        enum { sametype = false };
        enum { samebase = true };
        typedef std::complex<T> type;
    };
    template <typename T>
    struct Traits2<T,std::complex<T> >
    {
        enum { sametype = false };
        enum { samebase = true };
        typedef std::complex<T> type;
    };
    template <typename T>
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


#define TMV_RealType(T) typename tmv::Traits<T>::real_type
#define TMV_ComplexType(T) typename tmv::Traits<T>::complex_type

    template <typename T>
    inline bool isReal(T)
    { return true; }
    template <typename T>
    inline bool isReal(std::complex<T>)
    { return false; }
    template <typename T>
    inline bool isComplex(T x)
    { return !isReal(x); }
    template <typename T1, typename T2>
    inline bool isSameType(T1,T2)
    { return false; }
    template <typename T>
    inline bool isSameType(T,T)
    { return true; }

#ifndef NOTHROW
    class Error : public std::runtime_error
    {
    public :
        inline Error(std::string _s1) throw() :
            std::runtime_error("TMV Error: "+_s1) {}
        virtual inline ~Error() throw() {}
        // Typically, this is overridden by subclasses with more detail.
        virtual inline void write(std::ostream& os) const throw()
        { os << this->what() << std::endl; }
    };

    inline std::ostream& operator<<(std::ostream& os, const Error& e) throw()
    { e.write(os); return os; }

    class FailedAssert :
        public Error
    {
    public :
        std::string failed_assert;
        unsigned long line;
        std::string file;

        inline FailedAssert(std::string s, unsigned long l,
                            std::string f) throw() :
            Error("Failed Assert statement "+s),
            failed_assert(s), line(l), file(f) {}
        virtual inline ~FailedAssert() throw() {}
        virtual inline void write(std::ostream& os) const throw()
        {
            os<<"TMV Failed Assert: "<<failed_assert<<std::endl<<
                "on line "<<line<<" in file "<<file<<std::endl;
        }
    };
#endif

#ifndef TMV_DEBUG
#define TMVAssert(x)
#else
#ifdef NOTHROW
#define TMVAssert(x) do { if(!(x)) { \
    std::cerr<<"Failed Assert: "<<#x<<" on line "<<__LINE__<<" in file "<<__FILE__<<std::endl; exit(1); } } while (false)
#else
#define TMVAssert(x) do { if(!(x)) { \
    throw tmv::FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif
#endif
#ifdef NOTHROW
#define TMVAssert2(x) do { if(!(x)) { \
    std::cerr<<"Failed Assert: "<<#x<<" on line "<<__LINE__<<" in file "<<__FILE__<<std::endl; exit(1); } } while (false)
#else
#define TMVAssert2(x) do { if(!(x)) { \
    throw tmv::FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif

#ifdef __GNUC__
#define TMV_DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define TMV_DEPRECATED(func) __declspec(deprecated) func
#elif defined(__PGI)
#define TMV_DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define TMV_DEPRECATED(func) func
#endif

    // Use DEBUGPARAM(x) for parameters that are only used in TMVAssert
    // statements.  So then they don't give warnings when compiled with
    // -DNDEBUG
#ifdef TMV_DEBUG
#define TMV_DEBUG_PARAM(x) x
#else
#define TMV_DEBUG_PARAM(x)
#endif

#ifndef NOTHROW
    class ReadError :
        public Error
    {
    public :
        inline ReadError() throw() :
            Error("Invalid istream input encountered.") {}
        inline ReadError(std::string s) throw() :
            Error("Invalid istream input encountered while reading "+s) {}
        virtual inline ~ReadError() throw() {}
    };

    // Each Vector, Matrix type has it's own derived ReadError class
    // which outputs whatever has been read so far, and possibly
    // what was expected instead of what was actually read in.

    class Singular :
        public Error
    {
    public :
        inline Singular() throw() :
            Error("Encountered singular matrix.") {}
        inline Singular(std::string s) throw() :
            Error("Encountered singular "+s) {}
        virtual inline ~Singular() throw() {}
    };

    class NonPosDef :
        public Error
    {
    public:
        inline NonPosDef() throw() :
            Error("Invalid non-positive-definite matrix found.") {}
        inline NonPosDef(std::string s) throw() :
            Error("Non-positive-definite matrix found in "+s) {}
        virtual inline ~NonPosDef() throw() {}
    };
#endif

#ifdef TMV_WARN
    class TMV_WarnSingleton
    {
        // Note: This is not thread safe.
        // If multiple threads write to the warning stream at the
        // same time, they can clobber each other.
    public:
        static inline std::ostream*& inst() {
            static std::ostream* warn = 0;
            return warn;
        }
    private:
        TMV_WarnSingleton();
    };

    inline void TMV_Warning(std::string s)
    {
        if (TMV_WarnSingleton::inst()) {
            *TMV_WarnSingleton::inst() << "Warning:\n" << s << std::endl;
        }
    }

    inline std::ostream* WriteWarningsTo(std::ostream* newos)
    {
        std::ostream* temp = TMV_WarnSingleton::inst();
        TMV_WarnSingleton::inst() = newos;
        return temp;
    }

    inline void NoWarnings()
    { TMV_WarnSingleton::inst() = 0; }
#else
    inline void TMV_Warning(std::string ) {}
    inline std::ostream* WriteWarningsTo(std::ostream* os) { return 0; }
    inline void NoWarnings() {}
#endif

    // Unknown is the value of _size, _step, etc. whenever it
    // is not known at compile time.
    // We use for this value the maximally negative int.
    // In binary, this is a 1 followed by all zeros.
    // Note: can't use numeric_limits<int>::min, since it is a function,
    // not a compile-time constant.
    // And we need Unknown as a compile-time constant.
    const int Unknown = 1<<(sizeof(int)*8-1);

    // A helper structure that acts like an int,
    // but only bothers to make the integer if S == Unknown.
    // It also checks the constructor if S != Unknown
    template <int S>
    struct CheckedInt
    {
        CheckedInt(ptrdiff_t TMV_DEBUG_PARAM(s)) {
#ifdef TMV_DEBUG
            if (s != S) {
                std::cerr<<"Mismatched CheckInt:\n";
                std::cerr<<"template parameter S = "<<S<<std::endl;
                std::cerr<<"argument s = "<<s<<std::endl;
            }
#endif
            TMVAssert(s == S);
        }
        operator ptrdiff_t() const { return S; }
    };
    template <>
    struct CheckedInt<Unknown>
    {
        ptrdiff_t step;
        CheckedInt(ptrdiff_t s) : step(s) {}
        operator ptrdiff_t() const { return step; }
        ~CheckedInt() {}
    };

    template <typename T>
    inline TMV_RealType(T) TMV_Epsilon()
    { return std::numeric_limits<typename Traits<T>::real_type>::epsilon(); }

    template <typename T>
    inline bool TMV_Underflow(T x)
    {
        return tmv::Traits<T>::isinteger ? false :
            TMV_ABS2(x) < std::numeric_limits<T>::min();
    }

    template <typename T>
    inline bool TMV_Underflow(std::complex<T> x)
    {
        return tmv::Traits<T>::isinteger ? false :
            TMV_ABS2(x) < T(2) * std::numeric_limits<T>::min();
    }

    template <typename T>
    inline T TMV_Divide(T x, T y)
    { return x / y; }

    template <typename T>
    inline std::complex<T> TMV_Divide(std::complex<T> x, T y)
    {
        // return x / y;
        // Amazingly, some implementations get even this one wrong!
        // So I need to do each component explicitly.
        return std::complex<T>(x.real()/y,x.imag()/y);
    }

    template <typename T>
    inline T TMV_InverseOf(T x)
    { return T(1) / x; }

    template <typename T>
    inline std::complex<T> TMV_InverseOf(std::complex<T> x)
    {
        // The standard library's implemenation of complex division
        // does not guard against overflow/underflow.  So we have
        // our own version here.

        T xr = x.real();
        T xi = x.imag();
        if (std::abs(xr) > std::abs(xi)) {
            xi /= xr;
            T denom = xr*(T(1)+xi*xi);
            return std::complex<T>(T(1)/denom,-xi/denom);
        } else if (std::abs(xi) > T(0)) {
            xr /= xi;
            T denom = xi*(T(1)+xr*xr);
            return std::complex<T>(xr/denom,T(-1)/denom);
        } else {
            // (xr,xi) = (0,0)
            // division by zero.  Just go ahead and do it.
            return T(1) / xi;
        }
    }


    template <typename T>
    inline std::complex<T> TMV_Divide(T x, std::complex<T> y)
    { return x * TMV_InverseOf(y); }

    template <typename T>
    inline std::complex<T> TMV_Divide(std::complex<T> x, std::complex<T> y)
    { return x * TMV_InverseOf(y); }

    template <typename T>
    inline T TMV_SIGN(T x, T )
    { return x > 0 ? T(1) : T(-1); }

    template <typename T>
    inline std::complex<T> TMV_SIGN(std::complex<T> x, T absx)
    { return absx > 0 ? TMV_Divide(x,absx) : std::complex<T>(1); }


    extern bool TMV_FALSE;
    // = false (in TMV_Vector.cpp), but without the unreachable returns

    inline int TMV_ABS(const std::complex<int>& x)
    { return int(TMV_ABS(std::complex<double>(std::real(x),std::imag(x)))); }

    inline int TMV_ARG(const std::complex<int>& x)
    { return int(TMV_ARG(std::complex<double>(std::real(x),std::imag(x)))); }

    inline int TMV_SQRT(int x)
    { return int(TMV_SQRT(double(x))); }

    inline std::complex<int> TMV_SQRT(const std::complex<int>& x)
    {
        std::complex<double> temp =
            TMV_SQRT(std::complex<double>(std::real(x),std::imag(x)));
        return std::complex<int>(int(real(temp)),int(imag(temp)));
    }

    inline int TMV_EXP(int x)
    { return int(TMV_EXP(double(x))); }

    inline std::complex<int> TMV_EXP(const std::complex<int>& x)
    {
        std::complex<double> temp =
            TMV_EXP(std::complex<double>(std::real(x),std::imag(x)));
        return std::complex<int>(int(real(temp)),int(imag(temp)));
    }

    inline int TMV_LOG(int x)
    { return int(TMV_LOG(double(x))); }

    inline std::complex<int> TMV_LOG(const std::complex<int>& x)
    {
        std::complex<double> temp =
            TMV_LOG(std::complex<double>(std::real(x),std::imag(x)));
        return std::complex<int>(int(real(temp)),int(imag(temp)));
    }

    inline bool TMV_Underflow(int )
    { return false; }

    template <typename T>
    inline std::string TMV_Text(const T&)
    { return std::string("Unknown (") + typeid(T).name() + ")"; }

    inline std::string TMV_Text(const double&)
    { return "double"; }

    inline std::string TMV_Text(const float&)
    { return "float"; }

    inline std::string TMV_Text(const int&)
    { return "int"; }

    inline std::string TMV_Text(const long double&)
    { return "long double"; }

    template <typename T>
    inline std::string TMV_Text(std::complex<T>)
    { return std::string("complex<") + TMV_Text(T()) + ">"; }

    inline std::string TMV_Text(ConjType c)
    { return c == Conj ? "Conj" : c == NonConj ? "NonConj" : "VarConj"; }

    inline std::string TMV_Text(IndexStyle i)
    { return i == CStyle ? "CStyle" : "FortranStyle"; }

    inline std::string TMV_Text(StepType s)
    { return s == Unit ? "Unit" : "Step"; }

    inline std::string TMV_Text(StorageType s)
    {
        return
            s==ColMajor ? "ColMajor" :
            s==RowMajor ? "RowMajor" : "DiagMajor";
    }

    inline std::string TMV_Text(DiagType u)
    { return u == UnitDiag ? "UnitDiag" : "NonUnitDiag"; }

    inline std::string TMV_Text(UpLoType u)
    { return u == Upper ? "Upper" : "Lower"; }

    inline std::string TMV_Text(SymType s)
    { return s == Sym ? "Sym" : "Herm"; }

    inline std::string TMV_Text(DivType d)
    {
        return
            d==XX ? "XX" :
            d==LU ? "LU" :
            d==CH ? "CH" :
            d==QR ? "QR" :
            d==QRP ? "QRP" :
            d==SV ? "SV" :
            "unkown DivType";
    }

//#define TMVFLDEBUG

#ifdef TMVFLDEBUG
#define TMV_FIRSTLAST ,_first,_last
#define TMV_PARAMFIRSTLAST(T) , const T* _first, const T* _last
#define TMV_DEFFIRSTLAST(f,l) ,_first(f),_last(l)
#define TMV_DEFFIRSTLAST2(f,l) :_first(f),_last(l)
#define TMV_FIRSTLAST1(f,l) ,f,l
#define TMV_SETFIRSTLAST(f,l) _first=(f); _last=(l);
#else
#define TMV_FIRSTLAST
#define TMV_PARAMFIRSTLAST(T)
#define TMV_DEFFIRSTLAST(f,l)
#define TMV_DEFFIRSTLAST2(f,l)
#define TMV_FIRSTLAST1(f,l)
#define TMV_SETFIRSTLAST(f,l)
#endif

#define TMV_ConjOf(T,C) ((isReal(T()) || C==Conj) ? NonConj : Conj)

    // Use DEBUGPARAM(x) for parameters that are only used in TMVAssert
    // statements.  So then they don't give warnings when compiled with
    // -DNDEBUG
#ifdef TMV_DEBUG
#define TMV_DEBUGPARAM(x) x
#else
#define TMV_DEBUGPARAM(x)
#endif

// Workaround for migration from auto_ptr -> unique_ptr.  Use auto_ptr name, but all the
// functionality should work with either class.
#if __cplusplus >= 201103L
    template <typename T>
    using auto_ptr = std::unique_ptr<T>;
#else
    using std::auto_ptr;
#endif

} // namespace tmv

#endif
