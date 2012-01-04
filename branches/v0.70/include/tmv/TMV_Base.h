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
// MJ  BlockDiagMatrix and Sym varieties
// MJ  SparseMatrix and Sym, Block varieties


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

#ifdef MEMTEST
#include "mmgr.h"
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

    inline std::string TMV_Version() { return "0.70"; }
#define TMV_MAJOR_VERSION 0
#define TMV_MINOR_VERSION 70
#define TMV_VERSION_AT_LEAST(major,minor) \
    ( (major > TMV_MAJOR_VERSION) || \
      (major == TMV_MAJOR_VERSION && minor >= TMV_MINOR_VERSION) )

#if 0
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

        static std::string text()
        {
            std::string ret;
            ret += ( (A & ColMajor) ? "ColMajor|" :
                     (A & RowMajor) ? "RowMajor|" :
                     (A & DiagMajor) ? "DiagMajor|" : "");
            ret += (A & Conj) ? "Conj|" : "";
            ret += (A & FortranStyle) ? "FortranStyle|" : "";
            ret += (A & NonUnitDiag) ? "NonUnitDiag|" : "";
            ret += (A & UnitDiag) ? "UnitDiag|" : "";
            ret += (A & ZeroDiag) ? "ZeroDiag|" : "";
            ret += (A & Packed) ? "Packed|" : "";
            ret += (A & Lower) ? "Lower|" : "";
            ret += (A & Upper) ? "Upper|" : "";
            ret += (A & NoDivider) ? "NoDivider|" : "";
            ret += (A & WithDivider) ? "WithDivider|" : "";
            ret += (A & NoAlias) ? "NoAlias|" : "";
            ret += (A & CheckAlias) ? "CheckAlias|" : "";
            if (ret.size() == 0) return "Default";
            else return ret.substr(0,ret.size()-1);
        }
        static std::string vtext()
        {
            std::string ret;
            ret += (A & Unit) ? "Unit|" : "";
            ret += (A & Conj) ? "Conj|" : "";
            ret += (A & FortranStyle) ? "FortranStyle|" : "";
            ret += (A & NoAlias) ? "NoAlias|" : "";
            ret += (A & CheckAlias) ? "CheckAlias|" : "";
            if (ret.size() == 0) return "Default";
            else return ret.substr(0,ret.size()-1);
        }
    };
#endif

    // MJ: Add the packed storage varieties: RowPack, ColPack, DiagPack
    enum StorageType { RowMajor, ColMajor, DiagMajor, NoMajor };
    enum UpLoType { Upper, Lower };
    enum DiagType { UnitDiag, NonUnitDiag };
    // MJ: Any reason to add AntiSym, AntiHerm? Are they useful?
    enum SymType { Sym, Herm /*, AntiSym, AntiHerm*/ };
    enum IndexStyle { CStyle, FortranStyle };
    enum ADType { Ascend, Descend };
    enum CompType { RealComp, AbsComp, ImagComp, ArgComp };
    // Not sure why I did these in a different format, but let's normalize them.
    enum OldADType { ASCEND, DESCEND };
    enum OldCompType { REAL_COMP, ABS_COMP, IMAG_COMP, ARG_COMP };
    enum StepType { Unit, Step };
    enum ConjType { NonConj, Conj };

#define TMV_TransOf(s) \
    (s==RowMajor ? ColMajor : s==ColMajor ? RowMajor : s)
#define TMV_UTransOf(U) (U==Upper ? Lower : Upper)

    template <class M> 
    inline StorageType BaseStorOf(const M& m)
    { return m.stor()==RowMajor ? RowMajor : ColMajor; }

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

    template <class T> 
    inline T TMV_SQR(T x) 
    { return x*x; }

    template <class T> 
    inline T TMV_SQRT(T x) 
    { return T(std::sqrt(x)); }

    template <class T> 
    inline T TMV_EXP(T x) 
    { return T(std::exp(x)); }

    template <class T> 
    inline T TMV_LOG(T x) 
    { return T(std::log(x)); }

    template <class T> 
    inline T TMV_NORM(T x) 
    { return x*x; }

    template <class T> 
    inline T TMV_NORM(std::complex<T> x) 
    { return std::norm(x); }

    template <class T> 
    inline T TMV_CONJ(T x)
    { return x; }

    template <class T> 
    inline std::complex<T> TMV_CONJ(std::complex<T> x)
    { return std::conj(x); }

    template <class T> 
    inline T TMV_REAL(T x)
    { return x; }

    template <class T> 
    inline T TMV_REAL(std::complex<T> x)
    { return std::real(x); }

    template <class T> 
    inline T TMV_IMAG(T )
    { return T(0); }

    template <class T> 
    inline T TMV_IMAG(std::complex<T> x)
    { return std::imag(x); }

    template <class T> 
    inline T TMV_ARG(T x)
    { return x >= T(0) ? T(1) : T(-1); }

    template <class T> 
    inline T TMV_ARG(std::complex<T> x)
    { return arg(x); }

    template <class T> 
    inline T TMV_ABS(T x)
    { return std::abs(x); }

    template <class T> 
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

    template <class T> 
    inline T TMV_ABS2(T x)
    { return std::abs(x); }

    template <class T> 
    inline T TMV_ABS2(std::complex<T> x)
    { return std::abs(std::real(x)) + std::abs(std::imag(x)); }

    template <class T> 
    inline T TMV_SIGN(T x, T )
    { return x > 0 ? T(1) : T(-1); }

    template <class T> 
    inline std::complex<T> TMV_SIGN(std::complex<T> x, T absx)
    { return absx > 0 ? x/absx : std::complex<T>(1); }

    template <class T> 
    inline T TMV_MIN(T x, T y)
    { return x > y ? y : x; }

    template <class T> 
    inline T TMV_MAX(T x, T y)
    { return x > y ? x : y; }

    template <class T> 
    inline void TMV_SWAP(T& x, T& y)
    { T z = x; x = y; y = z; }

    template <class T>
    struct Traits
    {
        enum { isreal = true };
        enum { iscomplex = false };
        enum { isinteger = std::numeric_limits<T>::is_integer };
        enum { twoifcomplex = 1 };

        typedef T real_type;
        typedef std::complex<T> complex_type;
    };

    template <class T>
    struct Traits<std::complex<T> >
    {
        enum { isreal = false };
        enum { iscomplex = true };
        enum { isinteger = Traits<T>::isinteger };
        enum { twoifcomplex = 2 };

        typedef T real_type;
        typedef std::complex<T> complex_type;
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
        enum { sametype = false };
        enum { samebase = false };
        typedef T1 type; 
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


#define TMV_RealType(T) typename tmv::Traits<T>::real_type
#define TMV_ComplexType(T) typename tmv::Traits<T>::complex_type

    template <class T> 
    inline bool isReal(T) 
    { return true; }
    template <class T> 
    inline bool isReal(std::complex<T>) 
    { return false; }
    template <class T> 
    inline bool isComplex(T x) 
    { return !isReal(x); }
    template <class T1, class T2> 
    inline bool isSameType(T1,T2) 
    { return false; }
    template <class T> 
    inline bool isSameType(T,T) 
    { return true; }

#ifndef NOTHROW
    class Error : public std::runtime_error
    {
    public :
        std::string s1;
        std::string s2;
        inline Error(std::string _s1, std::string _s2="") throw() :
            std::runtime_error("TMV Error"), s1(_s1), s2(_s2) {}
        virtual inline ~Error() throw() {}
        virtual inline void write(std::ostream& os) const throw()
        { os << "TMV Error: " << s1 << s2 << std::endl; }
        virtual inline const char* what() const throw()
        { return s1.c_str(); }
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
            Error("Failed Assert statement ",s),
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

#ifndef NOTHROW
    class ReadError :
        public Error
    {
    public :
        inline ReadError() throw() :
            Error("Invalid istream input encountered.") {}
        inline ReadError(std::string s) throw() :
            Error("Invalid istream input encountered while reading ",s) {}
        virtual inline ~ReadError() throw() {}
        virtual inline void write(std::ostream& os) const throw()
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
        inline Singular() throw() :
            Error("Encountered singular matrix.") {}
        inline Singular(std::string s) throw() :
            Error("Encountered singular ",s) {}
        virtual inline ~Singular() throw() {}
        virtual inline void write(std::ostream& os) const throw()
        {
            os << "TMV Singular: " << Error::s1 << ' ' << Error::s2 <<
                std::endl; 
        }
    };

    class NonPosDef :
        public Error
    {
    public:
        inline NonPosDef() throw() :
            Error("Invalid non-positive-definite matrix found.") {}
        inline NonPosDef(std::string s) throw() :
            Error("Non-positive-definite matrix found in ",s) {}
        virtual inline ~NonPosDef() throw() {}
        virtual inline void write(std::ostream& os) const throw()
        {
            os << "TMV NonPosDef: " << Error::s1 << ' ' << Error::s2 << 
                std::endl; 
        }
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

    template <class T> 
    inline TMV_RealType(T) TMV_Epsilon()
    { return std::numeric_limits<typename Traits<T>::real_type>::epsilon(); }

    template <class T>
    inline bool TMV_Underflow(T x)
    {
        typedef typename Traits<T>::real_type RT;
        return TMV_ABS2(x) < 
            std::numeric_limits<RT>::min() * RT(tmv::Traits<T>::twoifcomplex); 
    }

    extern bool TMV_FALSE; 
    // = false (in TMV_Vector.cpp), but without the unreachable returns

    inline double TMV_ABS(const std::complex<int>& x)
    { return TMV_ABS(std::complex<double>(real(x),imag(x))); }

    inline double TMV_ARG(const std::complex<int>& x)
    { return TMV_ARG(std::complex<double>(real(x),imag(x))); }

    inline double TMV_SQRT(int x) 
    { return TMV_SQRT(double(x)); }

    inline std::complex<double> TMV_SQRT(const std::complex<int>& x)
    { return TMV_SQRT(std::complex<double>(real(x),imag(x))); }

    inline double TMV_EXP(int x) 
    { return TMV_EXP(double(x)); }

    inline std::complex<double> TMV_EXP(const std::complex<int>& x)
    { return TMV_EXP(std::complex<double>(real(x),imag(x))); }

    inline double TMV_LOG(int x) 
    { return TMV_LOG(double(x)); }

    inline std::complex<double> TMV_LOG(const std::complex<int>& x)
    { return TMV_LOG(std::complex<double>(real(x),imag(x))); }

    inline bool TMV_Underflow(int )
    { return false; }

    template <class T> 
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

    template <class T> 
    inline std::string TMV_Text(std::complex<T>)
    { return std::string("complex<") + TMV_Text(T()) + ">"; }

    inline std::string TMV_Text(StepType s)
    { return s == Unit ? "Unit" : "Step"; }

    inline std::string TMV_Text(ConjType c)
    { return c == Conj ? "Conj" : c == NonConj ? "NonConj" : "VarConj"; }

    inline std::string TMV_Text(UpLoType u)
    { return u == Upper ? "Upper" : "Lower"; }

    inline std::string TMV_Text(DiagType u)
    { return u == UnitDiag ? "UnitDiag" : "NonUnitDiag"; }

    inline std::string TMV_Text(SymType s)
    { 
        return s == Sym ? "Sym" : "Herm";
        //return s == Sym ? "Sym" : s == Herm ? "Herm" :
        //  s == AntiSym ? "AntiSym" : "AntiHerm";
    }

    inline std::string TMV_Text(DivType d)
    { 
        return 
            d==XX ? "XX" :
            d==LU ? "LU" :
            d==CH ? "CH" :
            d==QR ? "QR" :
            d==QRP ? "QRP" :
            d==SV ? "SV" :
            "unkown dt";
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

    template <class T> 
    inline ConjType conjugateOf(ConjType C)
    {
        return (isReal(T()) || C==Conj) ? NonConj : Conj;
    }
#define TMV_ConjOf(T,C) conjugateOf<T>(C)

    template <class T, StepType S, ConjType C> 
    class VIt;
    template <class T, StepType S, ConjType C> 
    class CVIt;
    template <class T> 
    class VIter;
    template <class T> 
    class CVIter;
    template <class T> 
    class ConjRef;
    template <class T> 
    class VarConjRef;

    template <class T> 
    struct AuxRef
    {
        typedef T& reference;
    };
    template <class T> 
    struct AuxRef<std::complex<T> >
    {
        typedef VarConjRef<std::complex<T> > reference;
    };
#define TMV_RefType(T) typename AuxRef<T>::reference

    // Use DEBUGPARAM(x) for parameters that are only used in TMVAssert
    // statements.  So then they don't give warnings when compiled with 
    // -DNDEBUG
#ifdef TMV_DEBUG
#define TMV_DEBUGPARAM(x) x
#else
#define TMV_DEBUGPARAM(x)
#endif

} // namespace tmv

#endif
