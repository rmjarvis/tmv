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
// *v   SmallVector
// *m   SmallMatrix
// *sB  SymBandMatrix 
// *hB  HermBandMatrix  
// *kD  BlockDiagMatrix (B was taken, so K = Block) 
// *skD SymBlockDiagMatrix 
// *hkD HermBlockDiagMatrix 
// *Q   SparseMatrix (S was taken, so Q = Sparse) 
// *sQ  SymSparseMatrix
// *hQ  HermSparseMatrix
// *kQ  BlockSparseMatrix 
// *skQ SymBlockSparseMatrix 
// *hkQ HermBlockSparseMatrix 
//
// MJ  SymBandMatrix
// MJ  BlockDiagMatrix and Sym varieties
// MJ  SparseMatrix and Sym, Block varieties
// MJ  SmallMatrix (nrows,ncols are templated)


#ifndef TMV_Base_H
#define TMV_Base_H

#ifndef NDEBUG
#define TMVDEBUG
#endif

#include <iostream>
#include <cmath>
#include <complex>
#include <limits>
#include <valarray>
#include <vector>
#include <algorithm>
#include <string>
#include <typeinfo>
#include <stdexcept>

#ifndef NOFLOAT
#define INST_FLOAT
#endif
#ifdef LONGDOUBLE
#define INST_LONGDOUBLE
#endif
//#define INST_INT

namespace tmv {

  using std::abs;
  using std::real;
  using std::imag;
  using std::complex;
  using std::vector;
  using std::valarray;
  using std::istream;
  using std::ostream;
  using std::endl;
  using std::cerr;
  using std::cout;
  using std::swap;
  using std::min;
  using std::max;
  using std::string;

  enum DivType { XXX, LU, CH, QR, QRP, SV, SVS, SVU, SVV };
  inline bool SV_Type(const DivType dt)
  { return dt >= SV && dt <= SVV; }

  // MJ: Add the packed storage varieties: RowPack, ColPack, DiagPack
  enum StorageType { RowMajor, ColMajor, DiagMajor, NoMajor };

  enum UpLoType { Upper, Lower };

  enum DiagType { UnitDiag, NonUnitDiag };

  // MJ: Any reason to add AntiSym, AntiHerm? Are they useful?
  enum SymType { Sym, Herm /*, AntiSym, AntiHerm*/ };

  template <class T> inline T SQR(T x) 
  { return x*x; }

  template <class T> inline T SQRT(T x) 
  { return sqrt(x); }

  template <class T> inline T NORM(T x) 
  { return x*x; }

  template <class T> inline T NORM(complex<T> x) 
  { return std::norm(x); }

  template <class T> inline T CONJ(T x)
  { return x; }

  template <class T> inline complex<T> CONJ(complex<T> x)
  { return std::conj(x); }

  template <class T> inline T PRINTCONJ(T x)
  { return x; }

  template <class T> inline complex<T> PRINTCONJ(complex<T> x)
  { return imag(x) == T(0) ? complex<T>(real(x),0) : std::conj(x); }

  template <class T> inline T REAL(T x)
  { return x; }

  template <class T> inline T REAL(complex<T> x)
  { return real(x); }

  template <class T> inline T IMAG(T )
  { return T(0); }

  template <class T> inline T IMAG(complex<T> x)
  { return imag(x); }

  template <class T> inline T ARG(T x)
  { return x >= T(0) ? T(1) : T(-1); }

  template <class T> inline T ARG(complex<T> x)
  { return arg(x); }

  template <class T> inline T SIGN(T x, T )
  { return x > 0 ? 1 : -1; }

  template <class T> inline complex<T> SIGN(complex<T> x, T absx)
  { return absx > 0 ? x/absx : complex<T>(1); }

  template <class T> inline T MAXABS(complex<T> x)
  { return max(abs(real(x)),abs(imag(x))); }

  template <class T> inline T MAXABS(T x)
  { return abs(x); }

  template <class T> class RCTypeClass 
  {
    public:
      typedef T RealType;
      typedef complex<T> ComplexType;
  };

  template <class T> class RCTypeClass<complex<T> >
  {
    public:
      typedef T RealType;
      typedef complex<T> ComplexType;
  };

#define RealType(T) typename tmv::RCTypeClass<T>::RealType
#define ComplexType(T) typename tmv::RCTypeClass<T>::ComplexType

  template <class T> inline bool IsReal(const T&) 
  { return true; }
  template <class T> inline bool IsReal(const complex<T>&) 
  { return false; }
  template <class T> inline bool IsComplex(const T&) 
  { return !IsReal(T()); }
  template <class T1, class T2> inline bool SameType(const T1&,const T2&) 
  { return false; }
  template <class T> inline bool SameType(const T&,const T&) 
  { return true; }

  class Error :
    public std::runtime_error
  {
    public :
      string s1, s2;
      Error(string _s1, string _s2="") throw() :
	std::runtime_error("TMV Error"), s1(_s1), s2(_s2) {}
      virtual ~Error() throw() {}
      virtual void Write(ostream& os) const throw()
      { os << "TMV Error: " << s1 << s2 << endl; }
      virtual const char* what() const throw()
      { return (string("TMV Error: ")+s1+s2).c_str(); }
  };

  inline ostream& operator<<(ostream& os, const Error& e) throw()
  { e.Write(os); return os; }

  class FailedAssert :
    public Error
  {
    public :
      string failed_assert;
      unsigned long line;
      string file;

      FailedAssert(string s, unsigned long l, string f) throw() :
	Error("Failed Assert statement ",s.c_str()),
	failed_assert(s), line(l), file(f) {}
      virtual ~FailedAssert() throw() {}
      virtual void Write(ostream& os) const throw()
      {
	os<<"TMV Failed Assert: "<<failed_assert<<endl;
	os<<"on line "<<line<<" in file "<<file<<endl;
      }
  };

#ifndef TMVDEBUG
#define TMVAssert(x)
#else
#define TMVAssert(x) do { if(!(x)) { \
  throw FailedAssert(#x,__LINE__,__FILE__); } } while(false)
#endif

  class ReadError :
    public Error
  {
    public :
      ReadError() throw() :
	Error("Reading istream input") {}
      ReadError(const char* s) throw() :
	Error("Reading istream input for ",s) {}
      virtual ~ReadError() throw() {}
      virtual void Write(ostream& os) const throw()
      { os << "TMV Read Error: " << Error::s1 << ' ' << Error::s2 << endl; }
  };

  // Each Vector, Matrix type has it's own derived ReadError class
  // which outputs whatever has been read so far, and possibly
  // what was expected instead of what was actually read in.

  class Singular :
    public Error
  {
    public :
      Singular() throw() :
	Error("Division by a singular matrix") {}
      Singular(const char* s) throw() :
	Error("Division by a singular ",s) {}
      virtual ~Singular() throw() {}
      virtual void Write(ostream& os) const throw()
      { os << "TMV Singular: " << Error::s1 << ' ' << Error::s2 << endl; }
  };

  class NonPosDef :
    public Error
  {
    public:
      NonPosDef() throw() :
	Error("Cholesky decomposition of non-positive-definite matrix") {}
      NonPosDef(const char* s) throw() :
	Error("Cholesky decomposition of non-positive-definite ",s) {}
      virtual ~NonPosDef() throw() {}
      virtual void Write(ostream& os) const throw()
      { os << "TMV NonPosDef: " << Error::s1 << ' ' << Error::s2 << endl; }
  };

  template <class T> inline RealType(T) Epsilon()
  { return std::numeric_limits<RealType(T)>::epsilon(); }

  template <class T> inline RealType(T) SqrtEpsilon()
  { 
    static const RealType(T) save = SQRT(Epsilon<T>());
    return save;
  }

  template <class T> inline string Type(const T&)
  { return string("Unknown (") + typeid(T()).name() + ")"; }

  template <> inline string Type(const double&)
  { return "double"; }

#ifdef INST_FLOAT
  template <> inline string Type(const float&)
  { return "float"; }
#endif

  extern bool FALSE; 
  // = false (in TMV_Vector.cpp), but without the unreachable returns

#ifdef INST_INT
  inline int Epsilon<int>()
  { TMVAssert(FALSE); return 1; }

  inline int ABS(const complex<int>& z)
  { TMVAssert(FALSE); return int(floor(abs(complex<double>(z)))); }

  template <> inline int SQRT(int x) 
  { TMVAssert(FALSE); return int(floor(sqrt(x))); }

  template <> inline string Type(const int&)
  { return "int"; }
#endif

#ifdef INST_LONGDOUBLE
  template <> inline long double SQRT(long double x)
  { return sqrtl(x); }

  template <> inline string Type(const long double&)
  { return "long double"; }
#endif

  template <class T> inline string Type(complex<T>)
  { return string("complex<") + Type(T()) + ">"; }

  inline string Text(UpLoType u)
  { return u == Upper ? "Upper" : "Lower"; }

  inline string Text(DiagType u)
  { return u == UnitDiag ? "UnitDiag" : "NonUnitDiag"; }

  inline string Text(SymType s)
  { 
    return s == Sym ? "Sym" : "Herm";
    //return s == Sym ? "Sym" : s == Herm ? "Herm" :
    //  s == AntiSym ? "AntiSym" : "AntiHerm";
  }

  inline string Text(DivType d)
  { 
    return 
      d==XXX ? "XXX" :
      d==LU ? "LU" :
      d==CH ? "CH" :
      d==QR ? "QR" :
      d==QRP ? "QRP" :
      d==SV ? "SV" :
      d==SVS ? "SVS" :
      d==SVU ? "SVU" :
      d==SVV ? "SVV" :
      "unkown dt";
  }

#ifdef NOSTLAUTO_PTR
  // This is copied more or less verbatim from gcc's auto_ptr implementation.
  // Although I ignore the conversions, since I don't use them.
  template <class _Tp> class auto_ptr {
    private:
      _Tp* _M_ptr;

    public:
      typedef _Tp element_type;
      explicit auto_ptr(_Tp* __p = 0) throw() : _M_ptr(__p) {}
      auto_ptr(auto_ptr& __a) throw() : _M_ptr(__a.release()) {}
      template <class _Tp1> auto_ptr(auto_ptr<_Tp1>& __a) throw()
	: _M_ptr(__a.release()) {}
      auto_ptr& operator=(auto_ptr& __a) throw() {
	if (&__a != this) {
	  delete _M_ptr;
	  _M_ptr = __a.release();
	}
	return *this;
      }
      template <class _Tp1>
	auto_ptr& operator=(auto_ptr<_Tp1>& __a) throw() {
	  if (__a.get() != this->get()) {
	    delete _M_ptr;
	    _M_ptr = __a.release();
	  }
	  return *this;
	}
      ~auto_ptr() throw() { delete _M_ptr; }

      _Tp& operator*() const throw() {
	return *_M_ptr;
      }
      _Tp* operator->() const throw() {
	return _M_ptr;
      }
      _Tp* get() const throw() {
	return _M_ptr;
      }
      _Tp* release() throw() {
	_Tp* __tmp = _M_ptr;
	_M_ptr = 0;
	return __tmp;
      }
      void reset(_Tp* __p = 0) throw() {
	delete _M_ptr;
	_M_ptr = __p;
      }
  };
#else
  using std::auto_ptr;
#endif
  // Identical, except delete[]
  template <class _Tp> class auto_array {
    private:
      _Tp* _M_array;

    public:
      typedef _Tp element_type;
      explicit auto_array(_Tp* __p = 0) throw() : _M_array(__p) {}
      auto_array(auto_array& __a) throw() : _M_array(__a.release()) {}
      template <class _Tp1> auto_array(auto_array<_Tp1>& __a) throw()
	: _M_array(__a.release()) {}
      auto_array& operator=(auto_array& __a) throw() {
	if (&__a != this) {
	  delete[] _M_array;
	  _M_array = __a.release();
	}
	return *this;
      }
      template <class _Tp1>
	auto_array& operator=(auto_array<_Tp1>& __a) throw() {
	  if (__a.get() != this->get()) {
	    delete[] _M_array;
	    _M_array = __a.release();
	  }
	  return *this;
	}
      ~auto_array() throw() { delete[] _M_array; }

      _Tp& operator*() const throw() {
	return *_M_array;
      }
      _Tp* operator->() const throw() {
	return _M_array;
      }
      _Tp* get() const throw() {
	return _M_array;
      }
      _Tp* release() throw() {
	_Tp* __tmp = _M_array;
	_M_array = 0;
	return __tmp;
      }
      void reset(_Tp* __p = 0) throw() {
	delete[] _M_array;
	_M_array = __p;
      }
  };

} // namespace tmv

#endif
