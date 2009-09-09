///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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

#ifndef NDEBUG
#define TMVDEBUG
#endif

#include <iosfwd>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>

#ifdef TMVDEBUG
#include <typeinfo>
#include <iostream>
#endif

#ifdef XDEBUG
#include <string>
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

//#define INST_LONGDOUBLE

//#define INST_INT

namespace tmv {

  enum DivType { XXX, LU, CH, QR, QRP, SV };
  // MJ: Add the packed storage varieties: RowPack, ColPack, DiagPack
  enum StorageType { RowMajor, ColMajor, DiagMajor, NoMajor };
  enum UpLoType { Upper, Lower };
  enum DiagType { UnitDiag, NonUnitDiag };
  // MJ: Any reason to add AntiSym, AntiHerm? Are they useful?
  enum SymType { Sym, Herm /*, AntiSym, AntiHerm*/ };
  enum IndexStyle { CStyle, FortranStyle };
  enum ADType { ASCEND, DESCEND };
  enum COMPType { REAL_COMP, ABS_COMP, IMAG_COMP, ARG_COMP };
  enum StepType { Unit, Step };
  enum ConjType { NonConj, Conj };

#define TransOf(s) \
  (s==RowMajor ? ColMajor : s==ColMajor ? RowMajor : s)
#define UTransOf(U) (U==Upper ? Lower : Upper)

  template <class M> inline StorageType BaseStorOf(const M& m)
  { return m.stor()==RowMajor ? RowMajor : ColMajor; }

#ifdef NOALIASCHECK
  template <class M1, class M2> inline bool SameStorage(const M1& , const M2& )
  { return false; }
#else
  template <class M1, class M2> inline bool SameStorage(
      const M1& m1, const M2& m2)
  { return m1.Real().cptr() == m2.Real().cptr(); }
#endif

  template <class T> inline T SQR(T x) 
  { return x*x; }

  template <class T> inline T SQRT(T x) 
  { return sqrt(x); }

  template <class T> inline T NORM(T x) 
  { return x*x; }

  template <class T> inline T NORM(std::complex<T> x) 
  { return std::norm(x); }

  template <class T> inline T CONJ(T x)
  { return x; }

  template <class T> inline std::complex<T> CONJ(std::complex<T> x)
  { return std::conj(x); }

  template <class T> inline T REAL(T x)
  { return x; }

  template <class T> inline T REAL(std::complex<T> x)
  { return std::real(x); }

  template <class T> inline T IMAG(T )
  { return T(0); }

  template <class T> inline T IMAG(std::complex<T> x)
  { return std::imag(x); }

  template <class T> inline T ARG(T x)
  { return x >= T(0) ? T(1) : T(-1); }

  template <class T> inline T ARG(std::complex<T> x)
  { return arg(x); }

  template <class T> inline T ABS(T x)
  { return std::abs(x); }

  template <class T> inline T ABS(std::complex<T> x)
  { return std::abs(x); }

  template <class T> inline T SIGN(T x, T )
  { return x > 0 ? T(1) : T(-1); }

  template <class T> inline std::complex<T> SIGN(std::complex<T> x, T absx)
  { return absx > 0 ? x/absx : std::complex<T>(1); }

  template <class T> inline T MIN(T x, T y)
  { return x > y ? y : x; }

  template <class T> inline T MAX(T x, T y)
  { return x > y ? x : y; }

  template <class T> inline void SWAP(T& x, T& y)
  { T z = x; x = y; y = z; }

  template <class T> class RCTypeClass 
  {
    public:
      typedef T RealType;
      typedef std::complex<T> ComplexType;
  };

  template <class T> class RCTypeClass<std::complex<T> >
  {
    public:
      typedef T RealType;
      typedef std::complex<T> ComplexType;
  };

#define RealType(T) typename tmv::RCTypeClass<T>::RealType
#define ComplexType(T) typename tmv::RCTypeClass<T>::ComplexType

  template <class T> inline bool IsReal(T) 
  { return true; }
  template <class T> inline bool IsReal(std::complex<T>) 
  { return false; }
  template <class T> inline bool IsComplex(T x) 
  { return !IsReal(x); }
  template <class T1, class T2> inline bool SameType(T1,T2) 
  { return false; }
  template <class T> inline bool SameType(T,T) 
  { return true; }

  class Error : public std::runtime_error
  {
    public :
      const char* s1;
      const char* s2;
      inline Error(const char* _s1, const char* _s2="") throw() :
	std::runtime_error("TMV Error"), s1(_s1), s2(_s2) {}
      virtual inline ~Error() throw() {}
      virtual inline void Write(std::ostream& os) const throw()
      { os << "TMV Error: " << s1 << s2 << std::endl; }
      virtual inline const char* what() const throw()
      { return s1; }
  };

  inline std::ostream& operator<<(std::ostream& os, const Error& e) throw()
  { e.Write(os); return os; }

  class FailedAssert :
    public Error
  {
    public :
      const char* failed_assert;
      unsigned long line;
      const char* file;

      inline FailedAssert(const char* s, unsigned long l, 
	  const char* f) throw() :
	Error("Failed Assert statement ",s),
	failed_assert(s), line(l), file(f) {}
      virtual inline ~FailedAssert() throw() {}
      virtual inline void Write(std::ostream& os) const throw()
      {
	os<<"TMV Failed Assert: "<<failed_assert<<std::endl<< 
	  "on line "<<line<<" in file "<<file<<std::endl;
      }
  };

#ifndef TMVDEBUG
#define TMVAssert(x)
#else
#define TMVAssert(x) do { if(!(x)) { \
  throw FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif
#define TMVAssert2(x) do { if(!(x)) { \
  throw FailedAssert(#x,__LINE__,__FILE__); } } while(false)

  class ReadError :
    public Error
  {
    public :
      inline ReadError() throw() :
	Error("Invalid istream input encountered.") {}
      inline ReadError(const char* s) throw() :
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
      inline Singular(const char* s) throw() :
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
      inline NonPosDef(const char* s) throw() :
	Error("Non-positive-definite matrix found in ",s) {}
      virtual inline ~NonPosDef() throw() {}
      virtual inline void Write(std::ostream& os) const throw()
      {
	os << "TMV NonPosDef: " << Error::s1 << ' ' << Error::s2 << 
	  std::endl; 
      }
  };

  // Put these in Vector.cpp to avoid having to include <limits>.
  // (The fewer includes, the quicker the compile, and it 
  // takes long enough as it is.)
  template <class T> RealType(T) Epsilon();
  template <class T> RealType(T) SqrtEpsilon();

  extern bool FALSE; 
  // = false (in TMV_Vector.cpp), but without the unreachable returns

#ifdef INST_INT
  inline int ABS(const std::complex<int>& z)
  { 
    return int(floor(std::abs(
	    std::complex<double>(std::real(z),std::imag(z))))); 
  }

  template <> inline int SQRT(int x) 
  { return int(floor(std::sqrt(double(x)))); }
#endif

#ifdef INST_LONGDOUBLE
  template <> inline long double SQRT(long double x)
  { return sqrtl(x); }
#endif

#ifdef XDEBUG
  template <class T> inline std::string Type(const T&)
  { return std::string("Unknown (") + typeid(T).name() + ")"; }

  template <> inline std::string Type(const double&)
  { return "double"; }

  template <> inline std::string Type(const float&)
  { return "float"; }

  template <> inline std::string Type(const int&)
  { return "int"; }

  template <> inline std::string Type(const long double&)
  { return "long double"; }

  template <class T> inline std::string Type(std::complex<T>)
  { return std::string("complex<") + Type(T()) + ">"; }

  inline std::string Text(StepType s)
  { return s == Unit ? "Unit" : "Step"; }

  inline std::string Text(ConjType c)
  { return c == Conj ? "Conj" : c == NonConj ? "NonConj" : "VarConj"; }

  inline std::string Text(UpLoType u)
  { return u == Upper ? "Upper" : "Lower"; }

  inline std::string Text(DiagType u)
  { return u == UnitDiag ? "UnitDiag" : "NonUnitDiag"; }

  inline std::string Text(SymType s)
  { 
    return s == Sym ? "Sym" : "Herm";
    //return s == Sym ? "Sym" : s == Herm ? "Herm" :
    //  s == AntiSym ? "AntiSym" : "AntiHerm";
  }

  inline std::string Text(DivType d)
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
#endif

#ifdef NOSTLAUTO_PTR
  // This is copied more or less verbatim from gcc's auto_ptr implementation.
  template <class X> class auto_ptr {

    private:
      X* ptr;
      template <class Y> struct auto_ptr_ref {
	Y* ptr;
	explicit auto_ptr_ref(Y* p) : ptr(p) {}
      };

    public:
      typedef X element_type;
      explicit auto_ptr(X* p = 0) throw() : ptr(p) {}
      auto_ptr(auto_ptr& a) throw() : ptr(a.release()) {}
      auto_ptr(auto_ptr_ref<X> ref) throw() : ptr(ref.ptr) {}
      template <class Y> auto_ptr(auto_ptr<Y>& a) throw() :
	ptr(a.release()) {}

      auto_ptr& operator=(auto_ptr& a) throw() 
      {
	reset(a.release());
	return *this;
      }
      template <class Y> auto_ptr& operator=(auto_ptr<Y>& a) throw() 
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

      template <class Y> operator auto_ptr_ref<Y>() throw()
      { return auto_ptr_ref<Y>(this->release()); }
      template <class Y> operator auto_ptr<Y>() throw()
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
  template <class X> class auto_array {

    private:
      X* ptr;
      template <class Y> struct auto_array_ref {
	Y* ptr;
	explicit auto_array_ref(Y* p) : ptr(p) {}
      };

    public:
      typedef X element_type;
      explicit auto_array(X* p = 0) throw() : ptr(p) {}
      auto_array(auto_array& a) throw() : ptr(a.release()) {}
      auto_array(auto_array_ref<X> ref) throw() : ptr(ref.ptr) {}
      template <class Y> auto_array(auto_array<Y>& a) throw() :
	ptr(a.release()) {}

      auto_array& operator=(auto_array& a) throw() 
      {
	reset(a.release());
	return *this;
      }
      template <class Y> auto_array& operator=(auto_array<Y>& a) throw() 
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

      template <class Y> operator auto_array_ref<Y>() throw()
      { return auto_array_ref<Y>(this->release()); }
      template <class Y> operator auto_array<Y>() throw()
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

  //#define TMVFLDEBUG

#ifdef TMVFLDEBUG
#define FIRSTLAST ,first,last
#define PARAMFIRSTLAST(T) , const T* _first, const T* _last
#define DEFFIRSTLAST(f,l) ,first(f),last(l)
#define FIRSTLAST1(f,l) ,f,l
#define SETFIRSTLAST(f,l) first=(f); last=(l);
#else
#define FIRSTLAST
#define PARAMFIRSTLAST(T)
#define DEFFIRSTLAST(f,l)
#define FIRSTLAST1(f,l)
#define SETFIRSTLAST(f,l)
#endif

  template <class T> inline ConjType ConjugateOf(ConjType C)
  {
    return (IsReal(T()) || C==Conj) ? NonConj : Conj;
  }
#define ConjOf(T,C) ConjugateOf<T>(C)

  template <class T, StepType S, ConjType C> class VIt;
  template <class T, StepType S, ConjType C> class CVIt;
  template <class T> class VIter;
  template <class T> class CVIter;
  template <class T> class ConjRef;
  template <class T> class VarConjRef;

  template <class T> struct AuxRef
  {
    typedef T& reference;
  };
  template <class T> struct AuxRef<std::complex<T> >
  {
    typedef VarConjRef<std::complex<T> > reference;
  };
#define RefType(T) typename AuxRef<T>::reference

  // Use DEBUGPARAM(x) for parameters that are only used in TMVAssert
  // statements.  So then they don't give warnings when compiled with 
  // -DNDEBUG
#ifdef TMVDEBUG
#define DEBUGPARAM(x) x
#else
#define DEBUGPARAM(x)
#endif

} // namespace tmv

#endif
