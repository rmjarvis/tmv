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
// This file defines the BaseVector, BaseVector_Calc, and BaseVector_Mutable
// classes.
//
// BaseVector is the base class for all of the various types of vectors.
// It defines all of the functions that you can use on any vector.
//
// There is a single template argument, which is the type of the most
// derived class for the object.
// So, for example, Vector<double> inherits from BaseVector<Vector<double> >.
//
// This is a pretty clever meta-programming trick that allows us to 
// get most of the functionality of virtual methods, but without the 
// expensive vtable lookups.  
// (I didn't invent it, so I'm allowed to say it is clever.  :)  )
//
// Basically the way it works is that BaseVector has a method called
// vec() that returns a down-cast to V:
//
// const V& vec() const { return *static_cast<const V*>(this); }
//
// Then, for any method foo(), where you would normally use virtual 
// inheritance, you instead call vec().foo(). 
//
// The methods that are defined for BaseVector are described in 
// TMV_Vector.h and include all of the constant access methods,
// the non-modifying functions, write, and arithmetic operators.
//
//
// BaseVector_Calc is the base class for all vectors that have their
// data calculated in memory, possibly with a non-unit step size, and
// possibly being the conjugate of the underlying data.
//
// The methods that are defined for BaseVector_Calc are described in
// TMV_Vector.h and include all of the vector views and iterators.
//
//
// BaseVector_Mutable is the base class for all vectors with calculated
// data (like BaseVector_Calc) whose values are allowed to be modified.
// 
// The methods that are defined for BaseVector_Calc are described in
// TMV_Vector.h and include all of the modifying functions, and the
// non-const access, read, vector views, and iterators.
//

#ifndef TMV_BaseVector_H
#define TMV_BaseVector_H

#include <sstream>
#include "TMV_Base.h"
#include "TMV_ListInit.h"

namespace tmv {

  // BaseVector is the base class for all vector classes.
  // All non-modifying functions are defined for BaseVector, such as
  // Norm, SumElements, MinElement, etc.
  template <class V>
  class BaseVector;

  // BaseVector gets a lot of the relevant information about V 
  // from the Traits class: Traits<V>
  // The types are defined as typedef statements.
  // The integer (or boolean) values are defined as enum statements.
  // See TMV_Vector.h or TMV_SmallVector.h for some concrete examples.
  //
  //
  // value_type = The type of the individual elements
  //
  // type = shorthand for the derived type 
  //
  // calc_type = The type of the calculated version of the vector.
  // i.e. where all of the values are stored in memory somewhere.
  // This is the return type of calc()
  // Use this type when you will access each element
  // of a composite vector multiple times.
  // Also use this when you will use cptr().
  //
  // eval_type = The type of the evaluated version of the vector
  // This is the return type of eval()
  // Use this type when you will access each element
  // of a composite vector only once.
  //
  // copy_type = The type of a new copy of the vector
  // This is the return type of copy()
  // Use this type when you need elements in the original vector
  // after overwriting those elements.
  //
  // vsize = size of vector (Use UNKNOWN if unknown at compile time)
  //
  // vfort = does the indexing use fortran style?
  //
  // vcalc = are the element values already calculated in memory?


  // BaseVector_Calc is derived from BaseVector, and is used
  // for vectors that have their values already calculated somewhere.
  // So composite classes inherit directly from BaseVector, rather
  // than BaseVector_Calc.
  template <class V>
  class BaseVector_Calc;

  // BaseVector_Calc adds some more requirements to the Traits<V> class:
  //
  //  vstep = the step size if known (else UNKNOWN)
  //  vconj = is the vector the conjugate of the underlying data?
  //
  //  const_subvector_type = return type from SubVector(i1,i2) const
  //  const_subvector_step_type = return type from SubVector(i1,i2,istep) const
  //  const_view_type = return type from View() const
  //  const_cview_type = return type from CView() const
  //  const_fview_type = return type from FView() const
  //  const_xview_type = return type from XView() const
  //  const_unitview_type = return type from UnitView() const
  //  const_conjugate_type = return type from Conjugate() const
  //  const_reverse_type = return type from Reverse() const
  //  const_realview_type = return type from Real() const
  //  const_imagview_type = return type from Imag() const
  //  const_flatten_type = return type from Flatten() const
  //  const_nonconj_type = return type from NonConj() const
  //  
  //  const_iterator = return type from begin(), end() const
  //  const_reverse_iterator = return type from rbegin(), rend() const


  // BaseVector_Mutable is derived from BaseVector_Calc, and is used
  // for vectors that are allowed to have their data modified.
  template <class V>
  class BaseVector_Mutable;

  // BaseVector_Mutable adds some more requirements to the Traits<V> class:
  //
  //  reference = return type of v(i)
  //
  //  subvector_type = return type from SubVector(i1,i2) 
  //  subvector_step_type = return type from SubVector(i1,i2,istep) 
  //  view_type = return type from View() 
  //  cview_type = return type from CView()
  //  fview_type = return type from FView()
  //  xview_type = return type from XView()
  //  unitview_type = return type from UnitView()
  //  conjugate_type = return type from Conjugate() 
  //  reverse_type = return type from Reverse() 
  //  realview_type = return type from Real()
  //  imagview_type = return type from Imag() 
  //  flatten_type = return type from Flatten() 
  //  nonconj_type = return type from NonConj()
  //  
  //  iterator = return type from begin(), end() 
  //  reverse_iterator = return type from rbegin(), rend() 

  // The following all derive from BaseVector_Mutable
  // See TMV_Vector.h and TMV_SmallVector.h for their definitions:
  template <class T, IndexStyle I=CStyle>
  class Vector;
  template <class T, int S=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class ConstVectorView;
  template <class T, int S=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class VectorView;
  template <class T, int N, IndexStyle I=CStyle>
  class SmallVector;
  template <class T, int N, int S=1, bool C=false, IndexStyle I=CStyle>
  class ConstSmallVectorView;
  template <class T, int N, int S=1, bool C=false, IndexStyle I=CStyle>
  class SmallVectorView;

  // These are effectively aliases for I = FortranStyle
  // so you don't have to write all the S,C values if you are using
  // the defaults when you want to specify FortranStyle.
  template <class T>
  class VectorF;
  template <class T, int S=UNKNOWN, bool C=false>
  class ConstVectorViewF;
  template <class T, int S=UNKNOWN, bool C=false>
  class VectorViewF;
  template <class T, int N>
  class SmallVectorF;
  template <class T, int N, int S=1, bool C=false>
  class ConstSmallVectorViewF;
  template <class T, int N, int S=1, bool C=false>
  class SmallVectorViewF;

  //
  // Helper functions and values:
  //

  // Values used to define how to sort a vector.
  enum ADType { ASCEND, DESCEND };
  enum COMPType { REAL_COMP, ABS_COMP, ABS2_COMP, IMAG_COMP, ARG_COMP, NORM_COMP, VALUE };

  // A helper class to pick out a component of a (usually) complex value
  template <COMPType comp, class T> struct Component;
  // First the real version:
  template <class T> struct Component<REAL_COMP,T>
  { 
    // return f(x);
    static inline T f(const T& x) { return x; } 
    // apply x = f(x)
    // (usually assigning the result to the real part of x when x is complex
    static inline void applyf(T& x) { }
    // get the applied value of f(x) after doing applyf
    // usually this returns real(x), but sometimes the whole x
    static inline T get(const T& x) { return x; }
  };
  template <class T> struct Component<IMAG_COMP,T>
  { 
    static inline T f(const T& x) { return T(0); } 
    static inline void applyf(T& x) { x = T(0); }
    static inline T get(const T& x) { return x; }
  };
  template <class T> struct Component<ABS_COMP,T>
  { 
    static inline T f(const T& x) { return TMV_ABS(x); } 
    static inline void applyf(T& x) { x = TMV_ABS(x); }
    static inline T get(const T& x) { return x; }
  };
  template <class T> struct Component<ABS2_COMP,T>
  { 
    static inline T f(const T& x) { return TMV_ABS2(x); } 
    static inline void applyf(T& x) { x = TMV_ABS2(x); }
    static inline T get(const T& x) { return x; }
  };
  template <class T> struct Component<ARG_COMP,T>
  { 
    static inline T f(const T& x) { return TMV_ARG(x); } 
    static inline void applyf(T& x) { x = TMV_ARG(x); }
    static inline T get(const T& x) { return x; }
  };
  template <class T> struct Component<NORM_COMP,T>
  { 
    static inline T f(const T& x) { return TMV_NORM(x); } 
    static inline void applyf(T& x) { x = TMV_NORM(x); }
    static inline T get(const T& x) { return x; }
  };
  template <class T> struct Component<VALUE,T>
  { 
    static inline T f(const T& x) { return x; } 
    static inline void applyf(T& x) { }
    static inline T get(const T& x) { return x; }
  };

  // Now the complex version:
#define CT std::complex<T>
  template <class T> struct Component<REAL_COMP,CT>
  { 
    static inline T f(const CT& x) { return real(x); } 
    static inline void applyf(CT& x) { }
    static inline T get(const T& x) { return real(x); }
  };
  template <class T> struct Component<IMAG_COMP,CT>
  { 
    static inline T f(const CT& x) { return imag(x); } 
    static inline void applyf(CT& x) { real(x) = imag(x); }
    static inline T get(const CT& x) { return real(x); }
  };
  template <class T> struct Component<ABS_COMP,CT>
  { 
    static inline T f(const CT& x) { return TMV_ABS(x); } 
    static inline void applyf(CT& x) { real(x) = TMV_ABS(x); }
    static inline T get(const CT& x) { return real(x); }
  };
  template <class T> struct Component<ABS2_COMP,CT>
  { 
    static inline T f(const CT& x) { return TMV_ABS2(x); } 
    static inline void applyf(CT& x) { real(x) = TMV_ABS2(x); }
    static inline T get(const CT& x) { return real(x); }
  };
  template <class T> struct Component<ARG_COMP,CT>
  { 
    static inline T f(const CT& x) { return TMV_ARG(x); } 
    static inline void applyf(CT& x) { real(x) = TMV_ARG(x); }
    static inline T get(const CT& x) { return real(x); }
  };
  template <class T> struct Component<NORM_COMP,CT>
  { 
    static inline T f(const CT& x) { return TMV_NORM(x); } 
    static inline void applyf(CT& x) { real(x) = TMV_NORM(x); }
    static inline T get(const CT& x) { return real(x); }
  };
  template <class T> struct Component<VALUE,CT>
  { 
    static inline CT f(const CT& x) { return x; } 
    static inline void applyf(CT& x) { }
    static inline CT get(const CT& x) { return x; }
  };
#undef CT

  // These helper functions check the validity of indices according
  // to whether the vector uses CStyle or FortranStyle indexing.
  // They also update the indices to be consistent with CStyle.
  template <bool vfort>
  inline void CheckIndex(int& i, int n) 
  { TMVAssert(i>=0 && i<n && "index is not valid"); } // CStyle
  template <>
  inline void CheckIndex<true>(int& i, int n) 
  { TMVAssert(i>=1 && i<=n && "index is not valid"); --i; } // FortranStyle
  template <bool vfort>
  inline void CheckRange(int& i1, int i2, int n)
  { // CStyle
    TMVAssert(i1 >= 0 && "first element must be in range");
    TMVAssert(i2 <= n && "last element must be in range");
    TMVAssert(i2 >= i1 && "range must have a non-negative number of elements");
  }
  template <>
  inline void CheckRange<true>(int& i1, int i2, int n)
  { // FortranStyle
    TMVAssert(i1 >= 1 && "first element must be in range");
    TMVAssert(i2 <= n && "last element must be in range");
    TMVAssert(i2 >= i1 && "range must have a positive number of elements");
    --i1;
  }
  template <bool vfort>
  inline void CheckRange(int& i1, int& i2, int istep, int n)
  { // CStyle
    TMVAssert(istep != 0 && "istep cannot be 0");
    TMVAssert(((i1 >= 0 && i1 < n) || i1==i2) && 
        "first element must be in range");
    TMVAssert(((i2-istep >= 0 && i2-istep < n) || i1==i2) &&
        "last element must be in range");
    TMVAssert((i2-i1) % istep == 0 &&
        "range must be an integral multiple of istep");
    TMVAssert((i2-i1) / istep >= 0 &&
        "must have a non-negative number of elements");
  }
  template <>
  inline void CheckRange<true>(int& i1, int& i2, int istep, int n)
  { // FortranStyle
    TMVAssert(istep != 0 && "istep cannot be 0");
    TMVAssert(i1 >= 1 && i1 <= n && "first element must be in range");
    TMVAssert(i2 >= 1 && i2 <= n && "last element must be in range");
    TMVAssert((i2-i1) % istep == 0 &&
        "subvector range must be an integral multiple of istep");
    TMVAssert((i2-i1) / istep >= 0 &&
        "subvector range must have a positive number of elements");
    --i1; i2 += istep-1;
  }

  // This helper class checks two compile time sizes.
  // same = Are they possibly the same size?
  // equal = Are they definitely the same size?
  // size = If knowable, what is that size?
  template <int S1, int S2>
  struct Sizes
  {
    enum { same = (S1==UNKNOWN || S2==UNKNOWN || S1==S2) };
    enum { equal = (S1!=UNKNOWN && S2!=UNKNOWN && S1==S2) };
    enum { size = (S1==UNKNOWN ? S2 : S1) };
  };

  // This helper class helps decide calc_type for composite classes:
  template <class T, int vsize, bool vfort>
  struct VCopyHelper
  { typedef SmallVector<T,vsize,vfort?FortranStyle:CStyle> type; };
  template <class T, bool vfort>
  struct VCopyHelper<T,UNKNOWN,vfort>
  { typedef Vector<T,vfort?FortranStyle:CStyle> type; };

  // This helper class acts as a ? : operator for a typedef
  template <bool first, class T1, class T2>
  struct TypeSelect // first = true
  { typedef T1 type; };
  template <class T1, class T2>
  struct TypeSelect<false,T1,T2> // first = false
  { typedef T2 type; };

  // This is a type to use when there is no valid return type for a 
  // particular function. 
  // It will give a compiler error if it is ever used.
  class InvalidType { private: InvalidType(); };

  // This helper function checks for aliasis in the memory storage of
  // two objects.  We overload it for specific objects that can be 
  // checked to have the same Real().cptr() values.
  template <class V1, class V2>
  inline bool SameStorage(
      const BaseVector<V1>& v1, const BaseVector<V2>& v2)
  { return false; }
#ifndef TMV_NO_ALIAS_CHECK
  template <class V1, class V2>
  inline bool SameStorage(
      const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2)
  { return v1.Real().cptr() == v2.Real().cptr(); }
#endif

  // This next one is sometime used after SameStorage is found to be true.
  // Here we check whether the elements are exactly the same 
  // by whether the steps are the same.
  // We do not check the vconj values, so that is usually the next 
  // step depending on why we are checking this.
  template <class V1, class V2>
  inline bool ExactSameStorage(
      const BaseVector<V1>& v1, const BaseVector<V2>& v2)
  { return false; }
  template <class V>
  inline bool ExactSameStorage(
      const BaseVector_Calc<V>& v1, const BaseVector_Calc<V>& v2)
  { return v1.step() == v2.step(); }

  // Defined in TMV_CopyV.h
  template <class V1, class V2>
  inline void Copy(const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);

  // Defined in TMV_SwapV.h
  template <class V1, class V2>
  inline void Swap(BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2);
  template <class V>
  inline void ReverseSelf(BaseVector_Mutable<V>& v);
  template <class V>
  inline void ConjugateSelf(BaseVector_Mutable<V>& v);

  // Defined in TMV_NormV.h
  template <class V>
  inline typename V::real_type NormSq(const BaseVector_Calc<V>& v);
  template <class V>
  inline typename V::real_type NormSq(const BaseVector_Calc<V>& v, 
      const typename V::real_type scale);
  template <class V>
  inline typename V::real_type Norm2(const BaseVector_Calc<V>& v);
  template <class V>
  inline typename V::value_type SumElements(const BaseVector_Calc<V>& v);
  template <class V>
  inline typename V::real_type SumAbsElements(const BaseVector_Calc<V>& v);
  template <class V>
  inline typename V::real_type SumAbs2Elements(const BaseVector_Calc<V>& v);

  // Defined in TMV_MinMax.h
  template <class V>
  inline typename V::value_type MaxElement(const BaseVector_Calc<V>& v,
      int*const imax=0);
  template <class V>
  inline typename V::real_type MaxAbsElement(const BaseVector_Calc<V>& v,
      int*const imax=0);
  template <class V>
  inline typename V::real_type MaxAbs2Element(const BaseVector_Calc<V>& v,
      int*const imax=0);
  template <class V>
  inline typename V::value_type MinElement(const BaseVector_Calc<V>& v,
      int*const imin=0);
  template <class V>
  inline typename V::real_type MinAbsElement(const BaseVector_Calc<V>& v,
      int*const imin=0);
  template <class V>
  inline typename V::real_type MinAbs2Element(const BaseVector_Calc<V>& v,
      int*const imin=0);

  // Defined in TMV_VectorIO.h
  template <class V>
  inline void Write(std::ostream& os, const BaseVector_Calc<V>& v);
  template <class V>
  inline void Write(std::ostream& os, const BaseVector_Calc<V>& v,
      typename V::real_type thresh) ;
  template <class V>
  inline void Read(std::istream& is, BaseVector_Mutable<V>& v);

  // Defined in TMV_SortV.h
  template <class V>
  inline void Sort(BaseVector_Mutable<V>& v, ADType ad, COMPType comp);
  template <class V>
  inline void Sort(
      BaseVector_Mutable<V>& v, int*const P, ADType ad, COMPType comp);

  //
  // BaseVector
  //

  template <class V> 
  class BaseVector
  {
  public:

    typedef V type;
    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;

    enum { vsize = Traits<type>::vsize }; 
    enum { vfort = Traits<type>::vfort };
    enum { vcalc = Traits<type>::vcalc };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { visreal = Traits<value_type>::isreal };
    enum { viscomplex = Traits<value_type>::iscomplex };


    //
    // Constructor
    //

    inline BaseVector() {}
    inline BaseVector(const BaseVector<type>&) {}
    inline ~BaseVector() {}

  private :
    inline void operator=(const BaseVector<type>& v2);
  public :

    //
    // Access
    // 

    inline value_type operator[](int i) const 
    { return operator()(i); }
    inline value_type operator()(int i) const 
    {
      CheckIndex<vfort>(i,size());
      return cref(i);
    }

    //
    // Functions
    //

    inline value_type SumElements() const
    { return tmv::SumElements(calc().CView()); }

    inline real_type SumAbsElements() const 
    { return tmv::SumAbsElements(calc().CView()); }

    inline real_type SumAbs2Elements() const // BLAS asum
    { return tmv::SumAbs2Elements(calc().CView()); }

    inline value_type MaxElement(int*const imax=0) const
    { return tmv::MaxElement(calc().CView(),imax); }

    inline real_type MaxAbsElement(int*const imax=0) const 
    { return tmv::MaxAbsElement(calc().CView(),imax); }

    inline real_type MaxAbs2Element(int*const imax=0) const // BLAS amax
    { return tmv::MaxAbs2Element(calc().CView(),imax); }

    inline value_type MinElement(int*const imin=0) const
    { return tmv::MinElement(calc().CView(),imin); }

    inline real_type MinAbsElement(int*const imin=0) const 
    { return tmv::MinAbsElement(calc().CView(),imin); }

    inline real_type MinAbs2Element(int*const imin=0) const // BLAS amin
    { return tmv::MinAbs2Element(calc().CView(),imin); }

    inline real_type Norm1() const
    { return SumAbsElements(); }

    inline real_type NormSq() const
    { return tmv::NormSq(calc().CView()); }

    inline real_type NormSq(const real_type scale) const
    { return tmv::NormSq(calc().CView(),scale); }

    inline real_type Norm2() const // BLAS nrm2
    { return tmv::Norm2(calc().CView()); }

    inline real_type Norm() const
    { return Norm2(); }

    inline real_type NormInf() const
    { return size() > 0 ? MaxAbsElement() : real_type(0); }



    // 
    // I/O
    //

    inline void Write(std::ostream& os) const
    { tmv::Write(os,calc().CView()); }
    inline void Write(std::ostream& os, real_type thresh) const
    { tmv::Write(os,calc().CView(),thresh); }



    //
    // Auxilliary routines
    //

    inline const type& vec() const 
    { return *static_cast<const type*>(this); }

    inline calc_type calc() const 
    { return static_cast<calc_type>(vec()); }

    inline eval_type eval() const 
    { return static_cast<eval_type>(vec()); }

    inline copy_type copy() const 
    { return static_cast<copy_type>(vec()); }

    // Note that these last functions need to be defined in a more derived
    // class than this, or an infinite loop will result when compiling.

    inline size_t size() const { return vec().size(); }

    inline value_type cref(int i) const  { return vec().cref(i); }

    template <class V2>
    inline void AssignTo(BaseVector_Mutable<V2>& v2) const
    { vec().AssignTo(v2); }

  }; // BaseVector

  //
  // BaseVector_Calc
  //

  template <class V> 
  class BaseVector_Calc : 
    public BaseVector<V>
  {
  public:

    typedef V type;

    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { visreal = Traits<value_type>::isreal };
    enum { viscomplex = Traits<value_type>::iscomplex };

    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;

    enum { vsize = Traits<type>::vsize };
    enum { vfort = Traits<type>::vfort };
    enum { vcalc = Traits<type>::vcalc };
    enum { vstep = Traits<type>::vstep };
    enum { vconj = Traits<type>::vconj };

    typedef typename Traits<type>::const_subvector_type const_subvector_type;
    typedef typename Traits<type>::const_subvector_step_type const_subvector_step_type;
    typedef typename Traits<type>::const_view_type const_view_type;
    typedef typename Traits<type>::const_cview_type const_cview_type;
    typedef typename Traits<type>::const_fview_type const_fview_type;
    typedef typename Traits<type>::const_xview_type const_xview_type;
    typedef typename Traits<type>::const_unitview_type const_unitview_type;
    typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
    typedef typename Traits<type>::const_reverse_type const_reverse_type;
    typedef typename Traits<type>::const_realview_type const_realview_type;
    typedef typename Traits<type>::const_imagview_type const_imagview_type;
    typedef typename Traits<type>::const_flatten_type const_flatten_type;
    typedef typename Traits<type>::const_nonconj_type const_nonconj_type;

    typedef typename Traits<type>::const_iterator const_iterator;
    typedef typename Traits<type>::const_reverse_iterator const_reverse_iterator;


    //
    // Constructor
    //

    inline BaseVector_Calc() {}
    inline BaseVector_Calc(const BaseVector_Calc<V>&) {}
    inline ~BaseVector_Calc() {}

  private:
    void operator=(const BaseVector_Calc<V>&);
  public:


    //
    // Access 
    //

    inline const_iterator begin() const
    { return const_iterator(cptr(),step()); }
    inline const_iterator end() const
    { return begin() + size(); }
    inline const_reverse_iterator rbegin() const
    { return const_reverse_iterator(cptr(),step()); }
    inline const_reverse_iterator rend() const
    { return rbegin() + size(); }


    //
    // SubVector
    //

    // CSubVector always uses CStyle
    inline const_subvector_type CSubVector(int i1, int i2) const
    { return const_subvector_type(cptr()+i1*step(),i2-i1,step()); }

    inline const_subvector_step_type CSubVector(
        int i1, int i2, int istep) const
    {
      return const_subvector_step_type(
          cptr()+i1*step(), (i2-i1)/istep, istep*step());
    }

    inline const_subvector_type SubVector(int i1, int i2) const
    {
      CheckRange<vfort>(i1,i2,size());
      return CSubVector(i1,i2);
    }

    inline const_subvector_step_type SubVector(int i1, int i2, int istep) const
    {
      CheckRange<vfort>(i1,i2,istep,size());
      return CSubVector(i1,i2,istep);
    }


    //
    // Views
    //

    inline const_view_type View() const
    { return const_view_type(cptr(),size(),step()); }

    inline const_cview_type CView() const
    { return const_cview_type(cptr(),size(),step()); }

    inline const_fview_type FView() const
    { return const_fview_type(cptr(),size(),step()); }

    inline const_xview_type XView() const
    { return const_xview_type(cptr(),size(),step()); }

    inline const_unitview_type UnitView() const
    {
      TMVAssert(step() == 1);
      return const_unitview_type(cptr(),size(),1); 
    }

    inline const_conjugate_type Conjugate() const
    { return const_conjugate_type(cptr(),size(),step()); }

    inline const_reverse_type Reverse() const
    {
      return const_reverse_type(
          cptr()+(int(size())-1)*step(), size(), -step());
    }

    inline const_realview_type Real() const
    {
      return const_realview_type(
          reinterpret_cast<const real_type*>(cptr()),
          size(), visreal ? step() : 2*step());
    }

    inline const_imagview_type Imag() const
    {
      TMVStaticAssert(viscomplex);
      return const_imagview_type(
          reinterpret_cast<const real_type*>(cptr())+1, size(), 2*step());
    }

    inline const_flatten_type Flatten() const
    {
      TMVStaticAssert(viscomplex);
      TMVStaticAssert(vstep == UNKNOWN || vstep == 1);
      TMVAssert(step() == 1);
      return const_flatten_type(
          reinterpret_cast<const real_type*>(cptr()),
          visreal ? size() : 2*size(), 1);
    }

    inline const_nonconj_type NonConj() const
    { return const_nonconj_type(cptr(),size(),step()); }


    //
    // Auxilliary routines
    //

    template <class V2>
    inline void AssignTo(BaseVector_Mutable<V2>& v2) const
    { tmv::Copy(*this,v2); }

    inline const type& vec() const
    { return *static_cast<const type*>(this); }

    inline bool isconj() const { return vconj; }

    // Note that these last functions need to be defined in a more derived
    // class than this, or an infinite loop will result when compiling.

    inline size_t size() const { return vec().size(); }
    inline int step() const { return vec().step(); }
    inline const value_type* cptr() const { return vec().cptr(); }
    inline value_type cref(int i) const  { return vec().cref(i); }

  }; // BaseVector_Calc

  //
  // BaseVector_Mutable
  //

  template <class V> 
  class BaseVector_Mutable : 
    public BaseVector_Calc<V>
  {
  public:

    typedef V type;
    typedef BaseVector_Calc<V> base_inst;
    typedef BaseVector_Mutable<V> base_mut;

    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { visreal = Traits<value_type>::isreal };
    enum { viscomplex = Traits<value_type>::iscomplex };

    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;

    enum { vsize = Traits<type>::vsize };
    enum { vfort = Traits<type>::vfort };
    enum { vcalc = Traits<type>::vcalc };
    enum { vstep = Traits<type>::vstep };
    enum { vconj = Traits<type>::vconj };

    typedef typename Traits<type>::const_subvector_type const_subvector_type;
    typedef typename Traits<type>::const_subvector_step_type const_subvector_step_type;
    typedef typename Traits<type>::const_view_type const_view_type;
    typedef typename Traits<type>::const_cview_type const_cview_type;
    typedef typename Traits<type>::const_fview_type const_fview_type;
    typedef typename Traits<type>::const_xview_type const_xview_type;
    typedef typename Traits<type>::const_unitview_type const_unitview_type;
    typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
    typedef typename Traits<type>::const_reverse_type const_reverse_type;
    typedef typename Traits<type>::const_realview_type const_realview_type;
    typedef typename Traits<type>::const_imagview_type const_imagview_type;
    typedef typename Traits<type>::const_flatten_type const_flatten_type;
    typedef typename Traits<type>::const_nonconj_type const_nonconj_type;

    typedef typename Traits<type>::const_iterator const_iterator;
    typedef typename Traits<type>::const_reverse_iterator const_reverse_iterator;

    typedef typename Traits<type>::subvector_type subvector_type;
    typedef typename Traits<type>::subvector_step_type subvector_step_type;
    typedef typename Traits<type>::view_type view_type;
    typedef typename Traits<type>::cview_type cview_type;
    typedef typename Traits<type>::fview_type fview_type;
    typedef typename Traits<type>::xview_type xview_type;
    typedef typename Traits<type>::unitview_type unitview_type;
    typedef typename Traits<type>::conjugate_type conjugate_type;
    typedef typename Traits<type>::reverse_type reverse_type;
    typedef typename Traits<type>::realview_type realview_type;
    typedef typename Traits<type>::imagview_type imagview_type;
    typedef typename Traits<type>::flatten_type flatten_type;
    typedef typename Traits<type>::nonconj_type nonconj_type;

    typedef typename Traits<type>::iterator iterator;
    typedef typename Traits<type>::reverse_iterator reverse_iterator;
    typedef typename Traits<type>::reference reference;

    //
    // Constructor
    //

    inline BaseVector_Mutable() {}
    inline BaseVector_Mutable(const base_mut&) {}
    inline ~BaseVector_Mutable() {}


    //
    // Access 
    //

    inline reference operator[](int i)
    { return operator()(i); }
    inline reference operator()(int i)
    {
      CheckIndex<vfort>(i,size());
      return ref(i);
    }

    inline iterator begin()
    { return iterator(ptr(),step()); }
    inline iterator end()
    { return begin() + size(); }
    inline reverse_iterator rbegin()
    { return reverse_iterator(ptr(),step()); }
    inline reverse_iterator rend()
    { return rbegin() + size(); }

    // We need to repeat the const versions so the non-const ones
    // don't clobber them.
    inline value_type operator[](int i) const
    { return base_inst::operator[](i); }
    inline value_type operator()(int i) const
    { return base_inst::operator()(i); }

    inline const_iterator begin() const
    { return base_inst::begin(); }
    inline const_iterator end() const
    { return base_inst::end(); }
    inline const_reverse_iterator rbegin() const
    { return base_inst::rbegin(); }
    inline const_reverse_iterator rend() const
    { return base_inst::rend(); }


    //
    // Op =
    //

    inline base_mut& operator=(base_mut& v2) 
    {
      TMVAssert(size() == v2.size());
      v2.AssignTo(*this);
      return *this; 
    }

    template <class V2>
    inline base_mut& operator=(const BaseVector<V2>& v2) 
    {
      TMVStaticAssert((Sizes<vsize,V2::vsize>::same));
      TMVAssert(size() == v2.size());
      v2.AssignTo(*this);
      return *this; 
    }

    inline ListAssigner<value_type,iterator> operator<<(value_type x)
    { return (ListAssigner<value_type,iterator>(begin(),size()),x); }


    //
    // Modifying Functions
    //

    inline type& Zero() 
    {
      const int n=size();
      for(int i=0;i<n;++i) ref(i) = value_type(0);
      return vec();
    }

    inline type& Clip(real_type thresh) 
    {
      const int n=size();
      for(int i=0;i<n;++i) {
        const real_type temp = TMV_ABS(cref(i));
        if (temp < thresh) ref(i) = value_type(0);
      }
      return vec();
    }

    inline type& SetAllTo(value_type x) 
    {
      const int n=size();
      for(int i=0;i<n;++i) ref(i) = x;
      return vec();
    }

    inline type& AddToAll(value_type x) 
    {
      const int n=size();
      for(int i=0;i<n;++i) ref(i) += x;
      return vec();
    }

    inline type& ConjugateSelf() // LAP lacgv
    { tmv::ConjugateSelf(vec()); return vec(); }

    template <class F>
    inline type& ApplyToAll(const F& f)
    {
      const int n=size();
      for(int i=0;i<n;++i) ref(i) = f(cref(i));
      return vec();
    }

    inline base_mut& CMakeBasis(int i, value_type x=value_type(1)) 
    {
      Zero(); ref(i) = x;
      return *this;
    }
    inline base_mut& MakeBasis(int i, value_type x=value_type(1)) 
    {
      CheckIndex<vfort>(i,size());
      return CMakeBasis(i,x);
    }

    inline base_mut& CSwap(int i1, int i2) 
    {
      TMV_SWAP(ref(i1),ref(i2)); 
      return *this;
    }
    inline base_mut& Swap(int i1, int i2) 
    {
      CheckIndex<vfort>(i1,size());
      CheckIndex<vfort>(i2,size());
      return CSwap(i1,i2);
    }

    inline base_mut& CPermute(const int*const p, int i1, int i2) 
    {
      for(int i=i1;i<i2;++i) CSwap(i,p[i]); 
      return *this;
    }
    inline base_mut& Permute(const int*const p, int i1, int i2) 
    {
      CheckRange<vfort>(i1,i2,size());
      return CPermute(p,i1,i2);
    }
    inline base_mut& Permute(const int*const p) 
    { return CPermute(p,0,size()); }

    inline base_mut& CReversePermute(const int*const p, int i1, int i2) 
    {
      for(int i=i2;i>i1;) { --i; CSwap(i,p[i]); }
      return *this;
    }
    inline base_mut& ReversePermute(const int*const p, int i1, int i2) 
    {
      CheckRange<vfort>(i1,i2,size());
      return CReversePermute(p,i1,i2);
    }
    inline base_mut& ReversePermute(const int*const p) 
    { return CReversePermute(p,0,size()); }

    inline base_mut& ReverseSelf() 
    { tmv::ReverseSelf(*this); return *this; }

    inline base_mut& Sort(ADType ad=ASCEND, COMPType comp=REAL_COMP) 
    { tmv::Sort(*this,ad,comp); return *this; }

    inline base_mut& Sort(
        int*const P, ADType ad=ASCEND, COMPType comp=REAL_COMP)
    { tmv::Sort(*this,P,ad,comp); return *this; }

    //
    // SubVector
    //

    inline subvector_type CSubVector(int i1, int i2) 
    { return subvector_type(ptr()+i1*step(),i2-i1,step()); }

    inline subvector_step_type CSubVector(int i1, int i2, int istep) 
    {
      return subvector_step_type(
          ptr()+i1*step(), (i2-i1)/istep, istep*step());
    }

    inline subvector_type SubVector(int i1, int i2) 
    {
      CheckRange<vfort>(i1,i2,size());
      return CSubVector(i1,i2);
    }

    inline subvector_step_type SubVector(int i1, int i2, int istep) 
    {
      CheckRange<vfort>(i1,i2,istep,size());
      return CSubVector(i1,i2,istep);
    }

    inline reverse_type Reverse() 
    { return reverse_type(ptr()+(int(size())-1)*step(), size(),-step()); }

    inline view_type View() 
    { return view_type(ptr(),size(),step()); }

    inline cview_type CView() 
    { return cview_type(ptr(),size(),step()); }

    inline fview_type FView() 
    { return fview_type(ptr(),size(),step()); }

    inline xview_type XView() 
    { return xview_type(ptr(),size(),step()); }

    inline unitview_type UnitView() 
    {
      TMVAssert(step() == 1);
      return unitview_type(ptr(),size(),1); 
    }

    inline conjugate_type Conjugate() 
    { return conjugate_type(ptr(),size(),step()); }


    inline realview_type Real() 
    {
      return realview_type(reinterpret_cast<real_type*>(ptr()),
          size(), visreal ? step() : 2*step());
    }

    inline imagview_type Imag() 
    {
      TMVStaticAssert(viscomplex);
      return imagview_type(
          reinterpret_cast<real_type*>(ptr())+1, size(), 2*step());
    }

    inline flatten_type Flatten() 
    {
      TMVStaticAssert(viscomplex);
      TMVStaticAssert(vstep == UNKNOWN || vstep == 1);
      TMVAssert(step() == 1);
      return flatten_type(reinterpret_cast<real_type*>(ptr()),
          visreal ? size() : 2*size(), 1);
    }

    inline nonconj_type NonConj()
    { return nonconj_type(ptr(),size(),step()); }



    // Repeat the const versions:
    inline const_subvector_type CSubVector(int i1, int i2) const
    { return base_inst::CSubVector(i1,i2); }
    inline const_subvector_step_type CSubVector(int i1, int i2, int istep) const
    { return base_inst::CSubVector(i1,i2,istep); }
    inline const_subvector_type SubVector(int i1, int i2) const
    { return base_inst::SubVector(i1,i2); }
    inline const_subvector_step_type SubVector(int i1, int i2, int istep) const
    { return base_inst::SubVector(i1,i2,istep); }
    inline const_reverse_type Reverse() const
    { return base_inst::Reverse(); }
    inline const_view_type View() const
    { return base_inst::View(); }
    inline const_cview_type CView() const
    { return base_inst::CView(); }
    inline const_fview_type FView() const
    { return base_inst::FView(); }
    inline const_xview_type XView() const
    { return base_inst::XView(); }
    inline const_unitview_type UnitView() const
    { return base_inst::UnitView(); }
    inline const_conjugate_type Conjugate() const
    { return base_inst::Conjugate(); }
    inline const_realview_type Real() const
    { return base_inst::Real(); }
    inline const_imagview_type Imag() const
    { return base_inst::Imag(); }
    inline const_flatten_type Flatten() const
    { return base_inst::Flatten(); }
    inline const_nonconj_type NonConj() const
    { return base_inst::NonConj(); }



    //
    // I/O
    //

    // Defined in TMV_VectorIO.h
    inline void Read(std::istream& is)
    {
      cview_type vcv = CView();
      tmv::Read(is,vcv); 
    }


    //
    // Arithmetic
    //

    // These operators need to be here, rather than just defining
    // non-member operator*=, etc., since the argument to a non-member
    // function would have to be BaseVector_Mutable& (ie. a reference).
    // But then you couldn't write something like:
    // m.row(i) *= 3.;
    // since the m.row(i) function returns a view by value, which is not
    // castable to the non-const reference argument of operator*=.
    // So we define all these here with unspecified right hand sides.
    // They just have to be valid arguments to a MultEq, AddEq, etc.
    // function defined as non-member functions for various objects.

    template <class X2>
    inline base_mut& operator+=(const X2& x2)
    { AddEq(vec(),x2); return *this; }

    template <class X2>
    inline base_mut& operator-=(const X2& x2)
    { SubtractEq(vec(),x2); return *this; }

    template <class X2>
    inline base_mut& operator*=(const X2& x2)
    { MultEq(vec(),x2); return *this; }

    template <class X2>
    inline base_mut& operator/=(const X2& x2)
    { DivEq(vec(),x2); return *this; }

    template <class X2>
    inline base_mut& operator%=(const X2& x2)
    { RDivEq(vec(),x2); return *this; }



    //
    // Auxilliary routines
    //

    inline const type& vec() const
    { return *static_cast<const type*>(this); }
    inline type& vec()
    { return *static_cast<type*>(this); }

    // Note that these last functionsneed to be defined in a more derived
    // class than this, or an infinite loop will result when compiling.

    inline size_t size() const { return vec().size(); }
    inline int step() const { return vec().step(); }
    inline value_type* ptr() { return vec().ptr(); }
    inline const value_type* cptr() { return vec().cptr(); }
    inline reference ref(int i) { return vec().ref(i); }
    inline value_type cref(int i) { return vec().cref(i); }

  }; // BaseVector_Mutable

  //
  // Vector ==, != Vector
  //

  template <class V1, class V2>
  static bool operator==(const BaseVector<V1>& v1, const BaseVector<V2>& v2)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    enum { size = Sizes<V1::vsize,V2::vsize>::size };
    const int n = (size == UNKNOWN ? int(v1.size()) : size);
    for(int i=0;i<n;++i) if (v1.cref(i) != v2.cref(i)) return false;
    return true;
  }

  template <class V1, class V2>
  inline bool operator!=(
      const BaseVector<V1>& v1, const BaseVector<V2>& v2)
  { return !(v1 == v2); }



  //
  // Other Functions of Vectors
  // (These just call the method version.)
  //

  template <class V>
  inline typename V::real_type Norm(const BaseVector<V>& v)
  { return v.Norm(); }

  template <class V>
  inline typename V::real_type Norm1(const BaseVector<V>& v)
  { return v.Norm1(); }

  template <class V>
  inline typename V::real_type Norm2(const BaseVector<V>& v)
  { return v.Norm2(); }

  template <class V>
  inline typename V::real_type NormSq(const BaseVector<V>& v)
  { return v.NormSq(); }

  template <class V>
  inline typename V::real_type NormInf(const BaseVector<V>& v)
  { return v.NormInf(); }

  template <class V>
  inline typename V::const_conjugate_type Conjugate(const BaseVector_Calc<V>& v)
  { return v.Conjugate(); }

  template <class V>
  inline typename V::conjugate_type Conjugate(BaseVector_Mutable<V>& v)
  { return v.Conjugate(); }

  //
  // TypeText 
  //

  inline std::string Text(ADType ad)
  { return ad == ASCEND ? "Ascend" : "Descend"; }

  inline std::string Text(COMPType comp)
  {
    return comp == REAL_COMP ? "Real" : comp == ABS_COMP ? "Abs" :
    comp == ABS2_COMP ? "Abs2" : comp == IMAG_COMP ? "Imag" : "Arg";
  }

  template <class V>
  static std::string TypeText(const BaseVector<V>& v)
  {
    std::ostringstream s;
    s << "BaseVector< "<<TypeText(v.vec())<<" >";
    return s.str();
  }

  template <class V>
  static std::string TypeText(const BaseVector_Calc<V>& v)
  {
    std::ostringstream s;
    s << "BaseVector_Calc< "<<TypeText(v.vec())<<" >";
    return s.str();
  }

  template <class V>
  static std::string TypeText(const BaseVector_Mutable<V>& v)
  {
    std::ostringstream s;
    s << "BaseVector_Mutable< "<<TypeText(v.vec())<<" >";
    return s.str();
  }

} // namespace tmv

#endif
