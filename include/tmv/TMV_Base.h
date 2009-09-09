///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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
// The WriteCompact routines for sparse Matrix types have an identifying
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
// TODO:  BlockDiagMatrix and Sym varieties
// TODO:  SparseMatrix and Sym, Block varieties


#ifndef TMV_Base_H
#define TMV_Base_H

#ifndef NDEBUG
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

#ifndef TMV_NO_INST_MIX
#define TMV_INST_MIX
#endif

//#define TMV_INST_LONGDOUBLE

//#define TMV_INST_INT


namespace tmv {

  // StorageType defines the order to store the elements of a matrix
  enum StorageType { RowMajor, ColMajor, DiagMajor, NoMajor, 
    RowPacked, ColPacked };

  // IndexStyle defines which kind of indexing to use for a vector
  enum IndexStyle { CStyle, FortranStyle=789234 };

  // UNKNOWN is the value of vsize or vstep whenever it
  // is not known at compile time.
  // We use for this value the maximally negative int.
  // In binary, this is a 1 followed by all zeros.
  const int UNKNOWN = (1<<(sizeof(int)*8-1));

  enum DivType { XXX, LU, CH, QR, QRP, SV };
  //enum UpLoType { Upper, Lower };
  enum DiagType { UnitDiag, NonUnitDiag };
  // TODO:: Any reason to add AntiSym, AntiHerm? Are they useful?
  //enum SymType { Sym, Herm /*, AntiSym, AntiHerm*/ };

//#define TransOf(s) \
  //(s==RowMajor ? ColMajor : s==ColMajor ? RowMajor : s)
//#define UTransOf(U) (U==Upper ? Lower : Upper)

  template <class T> 
  inline T TMV_SQR(const T& x) 
  { return x*x; }

  template <class T> 
  inline T TMV_SQRT(const T& x) 
  { return T(std::sqrt(x)); }

  template <class T> 
  inline T TMV_EXP(const T& x) 
  { return T(std::exp(x)); }

  template <class T> 
  inline T TMV_LOG(const T& x) 
  { return T(std::log(x)); }

  template <class T> 
  inline T TMV_NORM(const T& x) 
  { return x*x; }

  template <class T> 
  inline T TMV_NORM(const std::complex<T>& x) 
  { return std::norm(x); }

  template <class T> 
  inline T TMV_CONJ(const T& x)
  { return x; }

  template <class T> 
  inline std::complex<T> TMV_CONJ(const std::complex<T>& x)
  { return std::conj(x); }

  template <class T> 
  inline T TMV_REAL(const T& x)
  { return x; }

  template <class T> 
  inline T TMV_REAL(const std::complex<T>& x)
  { return std::real(x); }

  template <class T> 
  inline T TMV_IMAG(const T& )
  { return T(0); }

  template <class T> 
  inline T TMV_IMAG(const std::complex<T>& x)
  { return std::imag(x); }

  template <class T> 
  inline T TMV_ARG(const T& x)
  { return x >= T(0) ? T(1) : T(-1); }

  template <class T> 
  inline T TMV_ARG(const std::complex<T>& x)
  { return arg(x); }

  template <class T> 
  inline T TMV_ABS(const T& x)
  { return std::abs(x); }

  template <class T> 
  inline T TMV_ABS(const std::complex<T>& x)
  { return std::abs(x); }

  template <class T> 
  inline T TMV_ABS2(const T& x)
  { return std::abs(x); }

  template <class T> 
  inline T TMV_ABS2(const std::complex<T>& x)
  { return std::abs(std::real(x)) + std::abs(std::imag(x)); }

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

  template <class T>
  struct Traits
  {
    enum { isreal = true };
    enum { iscomplex = false };
    enum { isinst = InstType<T>::inst };

    typedef T real_type;
    typedef std::complex<T> complex_type;

    typedef T& conj_reference;
  };

  template <class T>
  struct Traits<std::complex<T> >
  {
    enum { isreal = false };
    enum { iscomplex = true };
    enum { isinst = InstType<T>::inst };

    typedef T real_type;
    typedef std::complex<T> complex_type;

    typedef ConjRef<std::complex<T> >& conj_reference;
  };

  template <class T>
  struct Traits<T&> : public Traits<T> {};
  template <class T>
  struct Traits<const T> : public Traits<T> {};
  template <class T>
  struct Traits<const T&> : public Traits<T> {};

#define RealType(T) typename tmv::Traits<T>::real_type
#define ComplexType(T) typename tmv::Traits<T>::complex_type

//#define IsReal(T) Traits<T>::isreal
//#define IsComplex(T) Traits<T>::iscomplex

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

  // A helper structure that conjugates a value if necessary:
  template <bool C, class T>
  struct DoConj_Helper // C = false
  { static inline T apply(const T& x) { return x; } };
  template <class T>
  struct DoConj_Helper<true,std::complex<T> >
  { static inline std::complex<T> apply(std::complex<T> x) { return std::conj(x); } };
  template <bool C, class T>
  inline T DoConj(const T& x) { return DoConj_Helper<C,T>::apply(x); }

  template <int S>
  struct IntTraits
  {
    enum { negS = -S };
    enum { twoS = 2*S };
    enum { ispowerof2 = S > 0 && !(S & S-1) };
    static inline int text() { return S; }
  };
  template <>
  struct IntTraits<UNKNOWN>
  {
    enum { negS = UNKNOWN };
    enum { twoS = UNKNOWN };
    enum { ispowerof2 = false };
    static inline const char* text() { return "UNKNOWN"; }
  };

  template <int S> struct IntLog;
  template <> struct IntLog<2> { enum { log = 1 }; };
  template <> struct IntLog<4> { enum { log = 2 }; };
  template <> struct IntLog<8> { enum { log = 3 }; };
  template <> struct IntLog<16> { enum { log = 4 }; };
  template <> struct IntLog<32> { enum { log = 5 }; };
  template <> struct IntLog<64> { enum { log = 6 }; };
  template <> struct IntLog<128> { enum { log = 7 }; };

  template <int Si, int Sj>
  struct IntTraits2
  {
    enum { sum = Si + Sj };
    enum { prod = Si * Sj };
  };
  template <int Si>
  struct IntTraits2<Si,UNKNOWN>
  {
    enum { sum = UNKNOWN };
    enum { prod = UNKNOWN };
  };
  template <int Sj>
  struct IntTraits2<UNKNOWN,Sj>
  {
    enum { sum = UNKNOWN };
    enum { prod = UNKNOWN };
  };
  template <>
  struct IntTraits2<UNKNOWN,UNKNOWN>
  {
    enum { sum = UNKNOWN };
    enum { prod = UNKNOWN };
  };

  template <int ix, class T> struct Scaling;
  template <int algo, int size, class V1, class V2> struct CopyV_Helper;
  template <int algo, int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper;

  template <bool conj1, bool conj2> struct ZProd;
  // A helper struct to pick one of two possibile behaviors
  // according to a template bool argument.  
  template <bool yn> struct Maybe;
  template <> struct Maybe<true>
  {
    // x or y
    template <class T>
    static inline T select(const T& x, const T& /*y*/) { return x; }

    // RealType(T) or T
    template <class T>
    struct RealType { typedef RealType(T) type; };

    // abs(x) or x
    template <class T>
    static inline RealType(T) abs(const T& x) { return TMV_ABS(x); }

    // conj(x) or x
    template <class T>
    static inline T conj(const T& x) { return TMV_CONJ(x); }

    // -x or x
    template <class T>
    static inline T neg(const T& x) { return -x; }

    // ++x or nothing
    template <class T>
    static inline void increment(T& x) { ++x; }

    // --x or nothing
    template <class T>
    static inline void decrement(T& x) { --x; }

    // x += y or x = y
    template <class T1, class T2>
    static inline void add(T1& x, const T2& y) { x += y; }
    template <class T1, class T2>
    static inline void add(ConjRef<T1> x, const T2& y) { x += y; }

    // x -= y or x = -y
    template <class T1, class T2>
    static inline void subtract(T1& x, const T2& y) { x -= y; }

    // x + y or y
    template <class T1, class T2>
    static inline T2 sum(const T1& x, const T2& y) { return x + y; }

    // x * y or y
    template <class T1, class T2>
    static inline T2 prod(const T1& x, const T2& y) { return x * y; }

    // x * y or y
    template <bool c1, bool c2, class T>
    static inline T zprod(const T& x, const T& y) 
    { return ZProd<c1,c2>::prod(x , y); }

    // x - y or x + y
    template <class T>
    static inline T diff(const T& x, const T& y) { return x - y; }

    // x = y or nothing
    template <class T1, class T2>
    static inline void set(T1& x, const T2& y) { x = y; }

    // v.Zero() or nothing
    template <class V>
    static inline void zero(V& v) { v.Zero(); }

    // v.AddToAll(x) or v.SetAllTo(x)
    template <class V, class T>
    static inline void addtoall(V& v, const T& x) 
    { v.AddToAll(x); }

    // m.diag().SetAllTo(1) or nothing
    template <class M>
    static inline void unitdiag(M& m) 
    { m.diag().SetAllTo(typename M::value_type(1)); }

    // x<y or x>y
    template <class T>
    static inline bool less(const T& x, const T& y) { return x<y; }

    // v1 += v2 or v1 = v2
    template <int size, class V1, class V2>
    static inline void addvv(V1& v1, const V2& v2) 
    {
      typedef typename V1::real_type RT;
      typedef typename V1::value_type T1;
      AddVV_Helper<-1,size,true,1,RT,V2,V1>::call(Scaling<1,RT>(1),v2,v1); 
    }

    // 2*x or x
    template <class T>
    static inline T twox(const T& x)
    { return 2*x; }
    static inline int twox(const int& x)
    { return x>>1; }

    // m2.OffDiag() = m1.OffDiag() or m2 = m1
    template <class M1, class M2>
    static inline void offdiag_copy(const M1& m1, M2& m2) 
    { if (m2.size() > 0) m2.OffDiag() = m1.OffDiag(); }

    // m2.UnitUpperTri() or m2.UpperTri()
    template <class M>
    static inline typename M::const_unit_uppertri_type unit_uppertri(
        const M& m) 
    { return m.UnitUpperTri(); }
    template <class M>
    static inline typename M::unit_uppertri_type unit_uppertri( M& m) 
    { return m.UnitUpperTri(); }

    // m2.UnitLowerTri() or m2.LowerTri()
    template <class M>
    static inline typename M::const_unit_lowertri_type unit_lowertri(
        const M& m) 
    { return m.UnitLowerTri(); }
    template <class M>
    static inline typename M::unit_lowertri_type unit_lowertri(M& m) 
    { return m.UnitLowerTri(); }

  };
  template <> struct Maybe<false>
  {
    template <class T>
    static inline T select(const T& , const T& y) { return y; }

    template <class T>
    struct RealType { typedef T type; };

    template <class T>
    static inline T abs(const T& x) { return x; }

    template <class T>
    static inline T conj(const T& x) { return x; }

    template <class T>
    static inline T neg(const T& x) { return x; }

    template <class T>
    static inline void increment(const T& ) { }

    template <class T>
    static inline void decrement(const T& ) { }

    template <class T1, class T2>
    static inline void add(T1& x, const T2& y) { x = y; }
    template <class T1, class T2>
    static inline void add(ConjRef<T1> x, const T2& y) { x = y; }

    template <class T1, class T2>
    static inline void subtract(T1& x, const T2& y) { x = -y; }

    template <class T1, class T2>
    static inline T2 sum(const T1& , const T2& y) { return y; }

    template <class T1, class T2>
    static inline T2 prod(const T1& , const T2& y) { return y; }

    template <bool c1, bool c2, class T>
    static inline T zprod(const T& x, const T& y) 
    { return Maybe<c2>::conj(y); }

    template <class T>
    static inline T diff(const T& x, const T& y) { return x + y; }

    template <class T1, class T2>
    static inline void set(T1& , const T2& ) { }

    template <class V>
    static inline void zero(V& ) { }

    template <class V, class T>
    static inline void addtoall(V& v, const T& x) 
    { v.SetAllTo(x); }

    template <class M>
    static inline void unitdiag(M& ) { }

    template <class T>
    static inline bool less(const T& x, const T& y) { return x>y; }

    template <int size, class V1, class V2>
    static inline void addvv(V1& v1, const V2& v2) 
    {
      typedef typename V1::value_type T1;
      typedef typename V2::value_type T2;
      CopyV_Helper<-1,size,V2,V1>::call(v2,v1);
    }

    template <class T>
    static inline T twox(const T& x)
    { return x; }

    template <class M1, class M2>
    static inline void offdiag_copy(const M1& m1, M2& m2) 
    { m2 = m1; }

    template <class M>
    static inline typename M::const_uppertri_type unit_uppertri(const M& m) 
    { return m.UpperTri(); }
    template <class M>
    static inline typename M::uppertri_type unit_uppertri(M& m) 
    { return m.UpperTri(); }

    template <class M>
    static inline typename M::const_lowertri_type unit_lowertri(const M& m) 
    { return m.LowerTri(); }
    template <class M>
    static inline typename M::lowertri_type unit_lowertri(M& m) 
    { return m.LowerTri(); }

  };

  // The way the STL (typically) implements std::complex means
  // that the exporession z = x * y when x,y,z are complex
  // is not very amenable to optimization by most compilers.
  // It turns out that doing the real and the imaginary parts 
  // separately usually presents the calculation in a form that 
  // the compiler can optimize, so we do this using a helper
  // structure ZProd.

  // First another helper class to do the correct complex multiplication
  // depending on the conj status of x and y:
  template <bool conj1, bool conj2> struct ZProd2;

  template <> struct ZProd2<false,false>
  {
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type rprod(
        const T1& rx, const T1& ix, const T2& ry, const T2& iy)
    { return (rx*ry - ix*iy); }
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type iprod(
        const T1& rx, const T1& ix, const T2& ry, const T2& iy)
    { return (rx*iy + ix*ry); }
  };
  template <> struct ZProd2<true,false>
  {
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type rprod(
        const T1& rx, const T1& ix, const T2& ry, const T2& iy)
    { return (rx*ry + ix*iy); }
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type iprod(
        const T1& rx, const T1& ix, const T2& ry, const T2& iy)
    { return (rx*iy - ix*ry); }
  };
  template <> struct ZProd2<false,true>
  {
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type rprod(
        const T1& rx, const T1& ix, const T2& ry, const T2& iy)
    { return (rx*ry + ix*iy); }
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type iprod(
        const T1& rx, const T1& ix, const T2& ry, const T2& iy)
    { return (ry*ix - iy*rx); }
  };
  template <> struct ZProd2<true,true>
  {
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type rprod(
        const T1& rx, const T1& ix, const T2& ry, const T2& iy)
    { return (rx*ry - ix*iy); }
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type iprod(
        const T1& rx, const T1& ix, const T2& ry, const T2& iy)
    { return -(rx*iy + ix*ry); }
  };

  // Now the real ZProd struct that does the right thing for each real
  // or complex possibility, including the option of x being a Scaling.
  template <bool conj1, bool conj2> struct ZProd
  {
    // complex * complex
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type rprod(
        const std::complex<T1>& x, const std::complex<T2>& y)
    { return ZProd2<conj1,conj2>::rprod(x.real(),x.imag(),y.real(),y.imag()); }
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type iprod(
        const std::complex<T1>& x, const std::complex<T2>& y)
    { return ZProd2<conj1,conj2>::iprod(x.real(),x.imag(),y.real(),y.imag()); }

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
    static inline typename Traits2<RealType(T1),T2>::type rprod(
        const Scaling<0,T1>& x, const T2& y)
    { return rprod(x.x,y); }
    template <class T1, class T2> 
    static inline typename Traits2<RealType(T1),T2>::type iprod(
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
    static inline typename Traits2<RealType(T1),T2>::type rprod(
        const Scaling<0,T1>& x, const std::complex<T2>& y)
    { return rprod(x.x,y); }
    template <class T1, class T2> 
    static inline typename Traits2<RealType(T1),T2>::type iprod(
        const Scaling<0,T1>& x, const std::complex<T2>& y)
    { return iprod(x.x,y); }
    template <class T1, class T2> 
    static inline T2 rprod(const Scaling<1,T1>& x, const std::complex<T2>& y)
    { return std::real(y); }
    template <class T1, class T2> 
    static inline T2 iprod(const Scaling<1,T1>& x, const std::complex<T2>& y)
    { return Maybe<conj2>::neg(y.imag()); }
    template <class T1, class T2> 
    static inline T2 rprod(const Scaling<-1,T1>& x, const std::complex<T2>& y)
    { return -std::real(y); }
    template <class T1, class T2> 
    static inline T2 iprod(const Scaling<-1,T1>& x, const std::complex<T2>& y)
    { return Maybe<!conj2>::neg(y.imag()); }

    // full prod
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type prod(
        const T1& x, const T2& y)
    {
      typedef typename Traits2<T1,T2>::type PT;
      enum { iscomplex = Traits<PT>::iscomplex };
      return makeprod<iscomplex,T1,T2>::call(x,y);
    }
    template <int ix, class T1, class T2> 
    static inline typename Traits2<T1,T2>::type prod(
        const Scaling<ix,T1>& x, const T2& y)
    { return x * Maybe<conj2>::conj(y); }
    template <int ix, class T1, class T2> 
    static inline typename Traits2<T1,T2>::type prod(
        const T1& x, const Scaling<ix,T2>& y)
    { return y * Maybe<conj1>::conj(x); }
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type prod(
        const Scaling<0,T1>& x, const T2& y)
    { return prod(x.x , y); }
    template <class T1, class T2> 
    static inline typename Traits2<T1,T2>::type prod(
        const T1& x, const Scaling<0,T2>& y)
    { return prod(x , y.x); }

    template <bool iscomplex, class T1, class T2> struct makeprod;
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
          
  };


#ifndef TMV_NO_THROW
  class Error : public std::runtime_error
  {
  public :
    std::string s1;
    std::string s2;
    inline Error(std::string _s1, std::string _s2="") throw() :
      std::runtime_error("TMV Error"), s1(_s1), s2(_s2) {}
    virtual inline ~Error() throw() {}
    virtual inline void Write(std::ostream& os) const throw()
    { os << "TMV Error: " << s1 << s2 << std::endl; }
    virtual inline const char* what() const throw()
    { return s1.c_str(); }
  };

  inline std::ostream& operator<<(std::ostream& os, const Error& e) throw()
  { e.Write(os); return os; }

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
    virtual inline void Write(std::ostream& os) const throw()
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
  throw FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif
#endif
#ifdef TMV_NO_THROW
#define TMVAssert2(x) do { if(!(x)) { \
  std::cerr<<"Failed Assert: "<<#x<<" on line " \
  <<__LINE__<<" in file "<<__FILE__<<std::endl; exit(1); } } while (false)
#else
#define TMVAssert2(x) do { if(!(x)) { \
  throw FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif

#if 1
  // This is based on the BOOST implementation BOOST_STATIC_ASSERT:
  template<bool> struct static_assert;
  template<> struct static_assert<true>
  { static inline void call() {} };
#define TMVStaticAssert(e) do { static_assert<(e) != 0>::call(); } while (false)
#else
  // This is copied from the CERT secure coding practices version:
#define TMV_JOIN(x, y) JOIN_AGAIN(x, y)
#define TMV_JOIN_AGAIN(x, y) x ## y

#define TMVStaticAssert(e) \
  typedef char TMV_JOIN(assertion_failed_at_line_, __LINE__) [(e) ? 1 : -1]
#endif

#ifndef TMV_NO_THROW
  class ReadError :
    public Error
  {
  public :
    inline ReadError() throw() :
      Error("Invalid istream input encountered.") {}
    inline ReadError(std::string s) throw() :
      Error("Invalid istream input encountered while reading ",s) {}
    virtual inline ~ReadError() throw() {}
    virtual inline void Write(std::ostream& os) const throw()
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
    virtual inline void Write(std::ostream& os) const throw()
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
    virtual inline void Write(std::ostream& os) const throw()
    {
      os << "TMV NonPosDef: " << Error::s1 << ' ' << Error::s2 << 
      std::endl; 
    }
  };
#endif

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
  inline void TMV_Warning(std::string ) {}
  inline std::ostream* WriteWarningsTo(std::ostream* ) {}
  inline void NoWarnings() {}
#endif

  // A helper structure that acts like 
  // int step;
  // but only bothers to make the integer if S == UNKNOWN.
  // It also checks the constructor if S != UNKNOWN
  template <int S>
  struct StepInt
  {
    StepInt(int s) 
    { TMVAssert(s == S); }
    operator int() const { return S; }
  };
  template <>
  struct StepInt<FortranStyle> 
  // This helps catch programming errors from the change in the 
  // template parameters of VectorView, MatrixView, etc. from
  // version 0.62 to 0.70.
  // It used to be that to make a fortran-style view, you would write:
  // VectorView<T,FortranStyle> v = m.row(i);
  // This is now wrong, since the first template parameter is now the step
  // size, so this would treat FortranStyle as the step size of the view.
  // One thing I did to help catch this was to make FortranStyle a 
  // crazy number that will probably never be used as an actual step
  // (789234).  This specialization of StepInt catches when FortranStyle
  // is used as a step and gives an error.
  {
    StepInt(int s) 
    { 
#if 0 
      TMVStaticAssert(false);
#else
      std::string str = 
      "FortranStyle used as a step size is probably an error...\n";
      str += "Now you should use F at the end of the name (e.g. VectorViewF)\n";
      str += "instead of using I=FotranStyle as a template parameter.\n";
      TMV_Warning(str);
      TMVAssert(s == FortranStyle); 
#endif
    }
    operator int() const { return FortranStyle; }
  };
  template <>
  struct StepInt<UNKNOWN>
  {
    int step;
    StepInt(int s) : step(s) {}
    operator int() const 
    { 
#ifdef TMV_DEBUG
      TMVAssert(step !=  -987234);
#endif
      return step; 
    }
    ~StepInt() {
#ifdef TMV_DEBUG
      step = -987234;
#endif
    }
  };


  template <class T> 
  inline RealType(T) Epsilon() 
  { return std::numeric_limits<RealType(T)>::epsilon(); }

  inline int TMV_ABS(const std::complex<int>& z)
  { 
    return int(floor(std::abs(
            std::complex<double>(std::real(z),std::imag(z))))); 
  }

  template <> 
  inline int TMV_SQRT(const int& x) 
  { return int(floor(std::sqrt(double(x)))); }

  template <> 
  inline int TMV_EXP(const int& x) 
  { return int(floor(std::exp(double(x)))); }

  template <> 
  inline int TMV_LOG(const int& x) 
  { return int(floor(std::log(double(x)))); }

  template <class T> 
  inline std::string TypeText(const T&)
  { return std::string("Unknown (") + typeid(T).name() + ")"; }

  template <> 
  inline std::string TypeText(const double&)
  { return "double"; }

  template <> 
  inline std::string TypeText(const float&)
  { return "float"; }

  template <> 
  inline std::string TypeText(const int&)
  { return "int"; }

  template <> 
  inline std::string TypeText(const long double&)
  { return "long double"; }

  template <class T> 
  inline std::string TypeText(std::complex<T>)
  { return std::string("complex<") + TypeText(T()) + ">"; }

  inline std::string Text(StorageType s)
  {
    return (
        s == RowMajor ? "RowMajor" :
        s == ColMajor ? "ColMajor" :
        s == DiagMajor ? "DiagMajor" :
        s == NoMajor ? "NoMajor" :
        s == RowPacked ? "RowPacked" :
        s == ColPacked ? "ColPacked" :
        "unkown StorageType");
  }

  inline std::string Text(IndexStyle i)
  { return i == CStyle ? "CStyle" : "FortranStyle"; }

  inline std::string Text(DivType d)
  { 
    return (
        d==XXX ? "XXX" :
        d==LU ? "LU" :
        d==CH ? "CH" :
        d==QR ? "QR" :
        d==QRP ? "QRP" :
        d==SV ? "SV" :
        "unkown DivType");
  }

  inline std::string Text(DiagType d)
  { return d == UnitDiag ? "UnitDiag" : "NonUnitDiag"; }

#if 0
  inline std::string Text(UpLoType u)
  { return u == Upper ? "Upper" : "Lower"; }

  inline std::string Text(SymType s)
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
  // Identical, except delete[]
  template <class X> 
  class auto_array {

  private:
    X* ptr;
    template <class Y> 
    struct auto_array_ref {
      Y* ptr;
      explicit auto_array_ref(Y* p) : ptr(p) {}
    };

  public:
    typedef X element_type;
    explicit auto_array(X* p = 0) throw() : ptr(p) {}
    auto_array(auto_array& a) throw() : ptr(a.release()) {}
    auto_array(auto_array_ref<X> ref) throw() : ptr(ref.ptr) {}
    template <class Y> 
    auto_array(auto_array<Y>& a) throw() :
      ptr(a.release()) {}

    auto_array& operator=(auto_array& a) throw() 
    {
      reset(a.release());
      return *this;
    }
    template <class Y> 
    auto_array& operator=(auto_array<Y>& a) throw() 
    {
      reset(a.release());
      return *this;
    }
    auto_array& operator=(auto_array_ref<X> ref) throw() 
    {
      if (ref.ptr != this->get())
      { 
        delete[] ptr;
        ptr = ref.ptr;
      }
      return *this;
    }

    ~auto_array() throw() { delete[] ptr; }

    template <class Y> 
    operator auto_array_ref<Y>() throw()
    { return auto_array_ref<Y>(this->release()); }
    template <class Y> 
    operator auto_array<Y>() throw()
    { return auto_array<Y>(this->release()); }

    X& operator[](const int i) const throw() { return *(ptr+i); }
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
        delete[] ptr;
        ptr = p;
      }
    }
  };

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
