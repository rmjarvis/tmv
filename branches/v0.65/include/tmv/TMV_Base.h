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
// TMVDEBUG defined.  (defined below in this file)
//
// This will turn on some basic debugging in the code and will help 
// you catch programming mistakes.
// For example, if you have two vectors: 
// v1 of length 5 and v2  of lengths 10,
// then v1 = v2 will only produce a runtime error if TMVDEBUG is 
// turned on.  Otherwise some random data will be overwritten and 
// who knows what will happen then.
//
// Once you know your code is working properly, you can recompile without
// the TMVDEBUG flag to speed up the code.
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

#if (!defined(NDEBUG) && !defined(TMVNDEBUG))
#define TMVDEBUG
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

#ifdef TMVDEBUG
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

    inline std::string TMV_Version() { return "0.65"; }
#define TMV_MAJOR_VERSION 0
#define TMV_MINOR_VERSION 65
#define TMV_VERSION_AT_LEAST(major,minor) \
    ( (major > TMV_MAJOR_VERSION) || \
      (major == TMV_MAJOR_VERSION && minor >= TMV_MINOR_VERSION) )

    enum DivType { XXX, LU, CH, QR, QRP, SV };
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

#ifndef TMVDEBUG
#define TMVAssert(x)
#else
#ifdef NOTHROW
#define TMVAssert(x) do { if(!(x)) { \
    std::cerr<<"Failed Assert: "<<#x<<" on line "<<__LINE__<<" in file "<<__FILE__<<std::endl; exit(1); } } while (false)
#else
#define TMVAssert(x) do { if(!(x)) { \
    throw FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif
#endif
#ifdef NOTHROW
#define TMVAssert2(x) do { if(!(x)) { \
    std::cerr<<"Failed Assert: "<<#x<<" on line "<<__LINE__<<" in file "<<__FILE__<<std::endl; exit(1); } } while (false)
#else
#define TMVAssert2(x) do { if(!(x)) { \
    throw FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif

#ifdef __GNUC__
#define TMV_DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define TMV_DEPRECATED(func) __declspec(deprecated) func
#elif defined(__PGI)
#define TMV_DEPRECATED(func) func __attribute__ ((deprecated))
#else
//#warning "WARNING: You should implement TMV_DEPRECATED for this compiler in the file TMV_Base.h"
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

    // defined in TMV_Vector.cpp
    // initialized with &std::cout;
    extern std::ostream* warn_out;
    inline void TMV_Warning(std::string s)
    {
        if (tmv::warn_out) { 
            *tmv::warn_out << "Warning:\n" << s << std::endl;
        }
    }

    inline std::ostream* WriteWarningsTo(std::ostream* os)
    { std::ostream* temp = tmv::warn_out; tmv::warn_out = os; return temp; }

    inline void NoWarnings()
    { tmv::warn_out = 0; }


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

    inline int TMV_ABS(const std::complex<int>& )
    { TMVAssert(TMV_FALSE); return 0; }

    inline int TMV_ARG(const std::complex<int>& )
    { TMVAssert(TMV_FALSE); return 0; }

    inline int TMV_SQRT(int ) 
    { TMVAssert(TMV_FALSE); return 0; }

    inline int TMV_EXP(int ) 
    { TMVAssert(TMV_FALSE); return 0; }

    inline int TMV_LOG(int ) 
    { TMVAssert(TMV_FALSE); return 0; }

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
            d==XXX ? "XXX" :
            d==LU ? "LU" :
            d==CH ? "CH" :
            d==QR ? "QR" :
            d==QRP ? "QRP" :
            d==SV ? "SV" :
            "unkown dt";
    }

#ifdef NOSTLAUTO_PTR
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

    //#define TMVFLDEBUG

#ifdef TMVFLDEBUG
#define TMV_FIRSTLAST ,_first,_last
#define TMV_PARAMFIRSTLAST(T) , const T* _first, const T* _last
#define TMV_DEFFIRSTLAST(f,l) ,_first(f),_last(l)
#define TMV_FIRSTLAST1(f,l) ,f,l
#define TMV_SETFIRSTLAST(f,l) _first=(f); _last=(l);
#else
#define TMV_FIRSTLAST
#define TMV_PARAMFIRSTLAST(T)
#define TMV_DEFFIRSTLAST(f,l)
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
#ifdef TMVDEBUG
#define TMV_DEBUGPARAM(x) x
#else
#define TMV_DEBUGPARAM(x)
#endif

} // namespace tmv

#endif
