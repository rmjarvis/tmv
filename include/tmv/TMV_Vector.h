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
// This file defines the TMV Vector class.
//
// The Vector class and all associated functions are contained
// in the namespace tmv.  Alse, the Vector class is a template, 
// so for a Vector of doubles, one would write 
// tmv::Vector<double>.  
//
// A second optional template argument can tell the Vector to
// use Fortran-style indexing instead of C-style.
// ie. the first element is v(1), and the last is v(n).
// So you would write Vector<T,FortranStyle> v(n), for example.
// If you want to use the Fortran-style indexing, you should
// be aware that it only exists for a regular Vector and the two
// view classes, VectorView and ConstVectorView.  
// GenVector and VectorComposite always use C-style indexing,
// so (v1+v2)(1) would return the second element in the sum, even
// if v1 and v2 both use Fortran-style indexing.
// Also, permutation arrays always use the C-style indexing.
//
// Constructors:
//
//    explicit Vector<T>(size_t n)
//        Makes a Vector of size n with _uninitialized_ values
//
//    Vector<T>(size_t n, T x)
//        Makes a Vector of size n with all values = x
//
//    Vector<T>(size_t n, const T* vv)
//    Vector<T>(const vector<T>& vv)
//        Makes a vector which copies the elements of vv
//        For the second one, n specifies the length of the vector
//
// Special Constructors
//
//    BasisVector(size_t n, int i)
//        Makes a Vector whose elements are all 0, except v(i) = 1
//
//    VectorViewOf(T* vv, size_t n)
//    VectorViewOf(const T* vv, size_t n)
//        Makes a Vector View (see below) which refers to the exact
//        elements of vv, not copying them to new storage.
//        The first one returns a VectorView, the second a ConstVectorView.
//
// Access Functions
//
//    size_t size() const
//        Returns the size of the Vector
//
//    T& operator[](int i)
//    T& operator()(int i)
//    T operator[](int i) const
//    T operator()(int i) const
//        Return the ith element of the Vector
//
//    Vector<T>::iterator begin()
//    Vector<T>::iterator end()
//    Vector<T>::reverse_iterator rbegin()
//    Vector<T>::reverse_iterator rend()
//    Vector<T>::const_iterator begin() const
//    Vector<T>::const_iterator end() const
//    Vector<T>::const_reverse_iterator rbegin() const
//    Vector<T>::const_reverse_iterator rend() const
//        Return iterators that work in the usual way to traverse the Vector
//
// Modifying Functions
//
//    Vector& Zero()
//        Sets all elements to 0
//
//    Vector& Clip(RT thresh)
//        Set to 0 all elements whose absolute values is < thresh
//
//    Vector& SetAllTo(T x)
//        Sets all elements to x
//        We don't call this v = x, since it doesn't really make sense to
//        think of v as being 'equal' to a scalar.
//        Hence the slightly verbose function name SetAllTo.
//
//    Vector& AddToAll(T x)
//        Add x to each element of v
//
//    Vector& ConjugateSelf()
//        Sets all elements to its conjugate
//
//    Vector& MakeBasis(int i)
//        Set all elements to 0, except v(i) = 1
//
//    Vector& Swap(int i1, int i2)
//        Swap elements v(i1) and v(i2)
//
//    Vector& Permute(const int* p)
//    Vector& Permute(const int* p, int i1, int i2)
//        Perform a series of swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (i1,p[i1])...(i2-1,p[i2-1])
//    Vector& ReversePermute(const int* p)
//    Vector& ReversePermute(const int* p, int i1, int i2)
//        The same but perform the swaps in reverse order.
//
//    Swap(v1,v2)
//        Swap vectors v1 and v2
//
//    Vector& ReverseSelf()
//        Reverse the order of the elements of v
//
//    Sort(AD, COMP)
//    Sort(int* P, AD, COMP)
//        Sorts the vector, returning the swaps required in P.
//        If you do not care about P, you may omit the P parameter.
//        AD = ASCEND or DESCEND (ASCEND=default)
//        COMP = REAL_COMP, ABS_COMP, IMAG_COMP, or ARG_COMP (REAL_COMP=default)
//
// VectorViews:
//
//    A VectorView<T> object refers to the elements of a regular Vector
//    or Matrix, so that altering the elements in the view alters the
//    corresponding elements in the original object.  A View can have non-unit
//    steps between elements and can also be a conjugation of the original
//    elements, so that v[3] = z would actually set the original element
//    to z*. 
//
//    Also, since we want to be able to pass these Views around, we have to 
//    keep track of whether we are allowed to alter the original values or
//    just look at them.  Thus, there are two VectorView objects: 
//    ConstVectorView and VectorView.  The first is only allowed to view,
//    not modify, the original elements.  The second is allowed to modify them.
//    This is akin to the const_iterator and iterator types for the
//    standard template library.
//
//    Below, VectorView is written for both of these, and I don't bother to
//    indicate the template modifiers.  In general, a ConstVectorView is 
//    returned from either a ConstVectorView object, or a const Vector object.
//    A (mutable) VectorView is returned from either a VectorView object or
//    a (non-const) Vector object.
//
//    VectorView SubVector(int i1, int i2, int istep=1)
//        This member function will return a subvector view which refers
//        to the same physical elements as the original.
//        i1 is the first element in the subvector.
//        i2 is the usual "one past the end" of the subvector.
//        istep is an optional step size.
//        Thus, if you have a vector v of length 10, and you want to
//        set all the even elements to 0, you could write:
//        v.SubVector(0,10,2).Zero();
//        And then to output the first 4 elements of v, you could write:
//        std::cout << v.SubVector(0,4);
//
//    VectorView Reverse()
//        Returns a subvector whose elements are the same as v but 
//        in the reverse order
//
//    VectorView Conjugate()
//        Returns the conjugate of a Vector (as a view, so it still points
//        to the same physical elements, but modifying this will set the 
//        actual elements to the conjugate of what you set.
//
//    VectorView View()
//        Returns a view of a Vector. 
//
//    VectorView Real()
//    VectorView Imag()
//        Returns a view of the real/imag elements of a complex Vector.
//    VectorView Flatten()
//        For complex vectors with unit step, returns a real vector view to 
//        both the real and imag elements in a single view.
//
// Functions of Vectors:
//        (These are all both member functions and functions of a Vector,
//         so Norm(v) and v.Norm() for example are equivalent)
//
//    Norm(v) or Norm2(v) 
//        Returns the 2-norm of a Vector
//        = sqrt( sum_i |v_i|^2 )
//
//    NormSq(v)
//    v.NormSq(RT scale=1.)
//        Returns the square of Norm()
//        In the method version, you can provide an optional scale, in 
//        which case the output is equal to NormSq(scale*v).
//
//    Norm1(v) 
//        Returns the 1-norm of a Vector
//        = sum_i |v_i|
//
//    NormInf(v) 
//        Returns the infinity-norm of a Vector
//        = max_i |v_i|
//
//    SumElements(v)
//        Returns the sum of all elements in the Vector
//
//    SumAbsElements(v)
//        Returns the sum of absolute values of  elements in the Vector
//
//    MaxElement(v)
//    v.MaxElement(int* imax = 0)
//        Returns the maximum value of any element in the Vector
//        In the method version, you can provide an optional argument imax.
//        On return, *imax holds the index of element with the maximum value.
//        The parameter imax can be omitted if it is not desired.
//        As "max" doesn't make sense for complex values, for these
//        we use just the real components.
//
//    MinElement(v)
//    v.MinElement(int* imin = 0)
//        Returns the minimum value of any element in the Vector
//        On return, (in the method version) *imin holds the index of this 
//        element.  Again, the parameter imin can be omitted.
//
//    MaxAbsElement(v)
//    MaxAbsElement(int* imax=0)
//        The same as MaxElement, except absolute values are used
//
//    MinAbsElement(v)
//    v.MinAbsElement(int* imin=0)
//        The same as MinElement, except absolute values are used
//
//
// Operators:
//        Here we use v for a Vector and x for a Scalar.
//
//        You can also mix real and complex vectors of the same
//        underlying type.  eg. Vector<double> and Vector<complex<double> >
//
//    -v
//
//    v = v
//
//    v += v
//    v + v
//
//    v -= v
//    v - v
//
//    v *= x
//    x * v
//    v * x
//    v * v   (inner product)
//
//    v /= x
//    v / x
//
//    v == v
//    v != v
//
//       These all behave in the logical way for dealing with vectors.
//
//
// I/O: 
//
//    os << v 
//        Writes v to ostream os in the following format:
//          n ( v[0] v[1] v[2] ... v[n] )
//
//    v.Write(ostream& os, Real(T) minnonzero)
//        Write v to os as above, but if |v[i]| < minnonzero, write 0 instead
//
//    is >> v
//        Reads v from istream is in the same format
//        Note: v must already be the correct size
//
//    is >> vptr
//        If you do not know the size of the Vector to be read in, you can
//        use this form where vptr is an auto_ptr to an undefined Vector.
//
//


#ifndef TMV_Vector_H
#define TMV_Vector_H

#include "tmv/TMV_BaseVector.h"
#include <vector>
#include "tmv/TMV_ListInit.h"

#ifdef TMVFLDEBUG
#include "tmv/TMV_VIt.h"
#endif

namespace tmv {

  template <class T, class T1> 
  inline void Copy(const GenVector<T1>& v1, const VectorView<T>& v2);

  template <class T> 
  class GenVector : 
    public AssignableToVector<T>
  {
  public:

    //
    // Constructor
    //

    inline GenVector() {}
    inline GenVector(const GenVector<T>&) {}
    virtual inline ~GenVector() {}


    //
    // Access Functions
    //

    using AssignableToVector<T>::size;
    void AssignToV(const VectorView<RealType(T)>& rhs) const
    {
      TMVAssert(rhs.size() == size());
      TMVAssert(IsReal(T()));
      if (!SameAs(rhs)) Copy(*this,rhs); 
    }
    void AssignToV(const VectorView<ComplexType(T)>& rhs) const
    { 
      TMVAssert(rhs.size() == size());
      if (!SameAs(rhs)) Copy(*this,rhs); 
    }

    inline T operator[](int i) const 
    { 
      TMVAssert(i>=0 && i<int(size()));
      return cref(i); 
    }
    inline T operator()(int i) const 
    { 
      TMVAssert(i>=0 && i<int(size()));
      return cref(i); 
    }

    typedef T value_type;
    typedef CVIter<T> const_iterator;
    typedef CVIter<T> const_reverse_iterator;

    inline const_iterator begin() const
    { return const_iterator(cptr(),step(),ct()); }
    inline const_iterator end() const
    { return begin() + size(); }
    inline const_reverse_iterator rbegin() const
    { return const_reverse_iterator(cptr(),step(),ct()); }
    inline const_reverse_iterator rend() const
    { return rbegin() + size(); }

    template <class T2> 
    inline bool SameAs(const GenVector<T2>&) const
    { return false; }

    inline bool SameAs(const GenVector<T>& v2) const
    {
      return (this == &v2 || (cptr()==v2.cptr() && 
            size()==v2.size() && step()==v2.step() && ct()==v2.ct()));
    }

    //
    // SubVector
    //

    bool OKSubVector(int i1, int i2, int istep) const;

    inline ConstVectorView<T> SubVector(int i1, int i2) const
    {
      TMVAssert(OKSubVector(i1,i2,1));
      return ConstVectorView<T>(cptr()+i1*step(),i2-i1,step(),ct());
    }

    inline ConstVectorView<T> SubVector(int i1, int i2, int istep) const
    {
      TMVAssert(OKSubVector(i1,i2,istep));
      return ConstVectorView<T>(cptr()+i1*step(),
          (i2-i1)/istep,istep*step(),ct());
    }

    inline ConstVectorView<T> Reverse() const
    { 
      return ConstVectorView<T>(cptr()+(int(size())-1)*step(),
          size(),-step(),ct()); 
    }

    inline ConstVectorView<T> View() const
    { return ConstVectorView<T>(cptr(),size(),step(),ct()); }

    inline ConstVectorView<T> Conjugate() const
    { return ConstVectorView<T>(cptr(),size(),step(),ConjOf(T,ct())); }

    inline ConstVectorView<RealType(T)> Real() const
    { 
      return ConstVectorView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(cptr()),
          size(), IsReal(T()) ? step() : 2*step(), NonConj);
    }

    inline ConstVectorView<RealType(T)> Imag() const
    { 
      TMVAssert(IsComplex(T()));
      return ConstVectorView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(cptr())+1,
          size(), 2*step(), NonConj);
    }

    inline ConstVectorView<RealType(T)> Flatten() const
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(step() == 1);
      return ConstVectorView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(cptr()),
          IsReal(T()) ? size() : 2*size(), 1, NonConj);
    }

    //
    // Functions of Vector
    //

    inline RealType(T) Norm() const 
    { return Norm2(); }

    RealType(T) NormSq(const RealType(T) scale = RealType(T)(1)) const; 

    inline RealType(T) Norm1() const 
    { return SumAbsElements(); }

    RealType(T) Norm2() const; 

    inline RealType(T) NormInf() const 
    { return size() > 0 ? MaxAbsElement() : RealType(T)(0); }

    T SumElements() const;

    RealType(T) SumAbsElements() const;

    T MinElement(int* iminout=0) const;

    T MaxElement(int* imaxout=0) const;

    RealType(T) MinAbsElement(int* iminout=0) const;

    RealType(T) MaxAbsElement(int* imaxout=0) const;

    //
    // I/O
    //

    void Write(std::ostream& os) const;
    void Write(std::ostream& os, RealType(T) minnonzero) const;

    virtual const T* cptr() const =0;
    virtual int step() const =0;
    virtual ConjType ct() const =0;
    inline bool isconj() const 
    { return (IsComplex(T()) && ct()==Conj); }

    virtual T cref(int i) const;

  private:

    inline GenVector<T>& operator=(const GenVector<T>&) 
    { TMVAssert(FALSE); return *this; }

  }; // GenVector

  template <class T, IndexStyle I> 
  class ConstVectorView : 
    public GenVector<T>
  {
  public:

    inline ConstVectorView(const ConstVectorView<T,I>& rhs) : 
      itsv(rhs.itsv), itssize(rhs.itssize), itsstep(rhs.itsstep),
      itsct(rhs.itsct) {}
    inline ConstVectorView(const GenVector<T>& rhs) : 
      itsv(rhs.cptr()), itssize(rhs.size()), itsstep(rhs.step()),
      itsct(rhs.ct()) {}
    inline ConstVectorView(const T* inv, size_t insize, int instep, 
        ConjType inct) : 
      itsv(inv), itssize(insize), itsstep(instep), itsct(inct) {}
    virtual inline ~ConstVectorView() 
    {
#ifdef TMVDEBUG
      const_cast<const T*&>(itsv) = 0;
#endif
    }

    virtual inline size_t size() const { return itssize; }
    virtual inline const T* cptr() const { return itsv; }
    virtual inline int step() const { return itsstep; }
    virtual inline ConjType ct() const { return itsct; }

  private:

    const T*const itsv;
    const size_t itssize;
    const int itsstep;
    const ConjType itsct;

    inline ConstVectorView<T,I>& operator=(const ConstVectorView<T,I>&) 
    { TMVAssert(FALSE); return *this; }

  }; // ConstVectorView

  template <class T> 
  class ConstVectorView<T,FortranStyle> : 
  public ConstVectorView<T,CStyle>
  {
  public:

    inline ConstVectorView(const ConstVectorView<T,FortranStyle>& rhs) : 
      ConstVectorView<T,CStyle>(rhs) {}
    inline ConstVectorView(const ConstVectorView<T,CStyle>& rhs) : 
      ConstVectorView<T,CStyle>(rhs) {}
    inline ConstVectorView(const GenVector<T>& rhs) : 
      ConstVectorView<T,CStyle>(rhs) {}
    inline ConstVectorView(const T* inv, size_t insize, int instep, 
        ConjType inct) : 
      ConstVectorView<T,CStyle>(inv,insize,instep,inct) {}
    virtual inline ~ConstVectorView() {}

    //
    // Access Functions
    //

    inline T operator[](int i) const 
    { 
      TMVAssert(i>0 && i<=int(this->size()));
      return GenVector<T>::cref(i-1); 
    }
    inline T operator()(int i) const 
    { 
      TMVAssert(i>0 && i<=int(this->size()));
      return GenVector<T>::cref(i-1); 
    }

    // 
    // SubVector
    //

    bool OKSubVector(int i1, int i2, int istep) const;

    inline ConstVectorView<T,FortranStyle> SubVector(int i1, int i2) const
    {
      TMVAssert(OKSubVector(i1,i2,1));
      return GenVector<T>::SubVector(i1-1,i2);
    }

    inline ConstVectorView<T,FortranStyle> SubVector(
        int i1, int i2, int istep) const
    {
      TMVAssert(OKSubVector(i1,i2,istep));
      return GenVector<T>::SubVector(i1-1,i2-1+istep,istep);
    }

    inline ConstVectorView<T,FortranStyle> Reverse() const
    { return GenVector<T>::Reverse(); }

    inline ConstVectorView<T,FortranStyle> View() const
    { return GenVector<T>::View(); }

    inline ConstVectorView<T,FortranStyle> Conjugate() const
    { return GenVector<T>::Conjugate(); }

    inline ConstVectorView<RealType(T),FortranStyle> Real() const
    { return GenVector<T>::Real(); }

    inline ConstVectorView<RealType(T),FortranStyle> Imag() const
    { return GenVector<T>::Imag(); }

    inline ConstVectorView<RealType(T),FortranStyle> Flatten() const
    { return GenVector<T>::Flatten(); }

    inline T MinElement(int* iminout=0) const
    { 
      T temp = GenVector<T>::MinElement(iminout);
      if (iminout) ++(*iminout);
      return temp;
    }

    inline T MaxElement(int* imaxout=0) const
    { 
      T temp = GenVector<T>::MaxElement(imaxout);
      if (imaxout) ++(*imaxout);
      return temp;
    }

    inline RealType(T) MinAbsElement(int* iminout=0) const
    { 
      RealType(T) temp = GenVector<T>::MinAbsElement(iminout);
      if (iminout) ++(*iminout);
      return temp;
    }

    inline RealType(T) MaxAbsElement(int* imaxout=0) const
    { 
      RealType(T) temp = GenVector<T>::MaxAbsElement(imaxout);
      if (imaxout) ++(*imaxout);
      return temp;
    }

  private:

    inline ConstVectorView<T,FortranStyle>& operator=(
        const ConstVectorView<T,FortranStyle>&) 
    { TMVAssert(FALSE); return *this; }

  }; // FortranStyle ConstVectorView

  template <class T, IndexStyle I> 
  class VectorView : 
    public GenVector<T>
  {
  public:

    //
    // Constructors 
    //

    inline VectorView(const VectorView<T,I>& rhs) : 
      itsv(rhs.itsv), itssize(rhs.itssize), itsstep(rhs.itsstep),
      itsct(rhs.itsct) DEFFIRSTLAST(rhs.first,rhs.last) {}

    inline VectorView(T* inv, size_t insize, int instep, ConjType inct 
        PARAMFIRSTLAST(T) ) :
      itsv(inv), itssize(insize), itsstep(instep),
      itsct(inct) DEFFIRSTLAST(_first,_last) {}

    virtual inline ~VectorView() 
    {
#ifdef TMVDEBUG
      const_cast<T*&>(itsv) = 0;
#endif
    }


    //
    // Op =
    //

    inline const VectorView<T,I>& operator=(const VectorView<T,I>& v2) const
    { 
      TMVAssert(size() == v2.size());
      v2.AssignToV(*this);
      return *this; 
    }

    inline const VectorView<T,I>& operator=(const VectorView<T,I>& v2) 
    { 
      TMVAssert(size() == v2.size());
      v2.AssignToV(*this);
      return *this; 
    }

    inline const VectorView<T,I>& operator=(
        const GenVector<RealType(T)>& v2) const
    { 
      TMVAssert(size() == v2.size());
      v2.AssignToV(*this);
      return *this; 
    }

    inline const VectorView<T,I>& operator=(
        const GenVector<ComplexType(T)>& v2) const
    { 
      TMVAssert(size() == v2.size());
      TMVAssert(IsComplex(T()));
      v2.AssignToV(*this);
      return *this; 
    }

    template <class T2> 
    inline const VectorView<T,I>& operator=(
        const GenVector<T2>& v2) const
    { 
      TMVAssert(IsReal(T2()) || IsComplex(T()));
      TMVAssert(size() == v2.size());
      Copy(v2,*this);
      return *this; 
    }

    inline const VectorView<T,I>& operator=(
        const AssignableToVector<RealType(T)>& v2) const
    { 
      TMVAssert(size() == v2.size());
      v2.AssignToV(*this);
      return *this; 
    }

    inline const VectorView<T,I>& operator=(
        const AssignableToVector<ComplexType(T)>& v2) const
    { 
      TMVAssert(size() == v2.size());
      TMVAssert(IsComplex(T()));
      v2.AssignToV(*this);
      return *this; 
    }

    template <class T2, int N, IndexStyle I2> 
    inline const VectorView<T,I>& operator=(
        const SmallVector<T2,N,I2>& v2) const
    { 
      TMVAssert(size() == v2.size());
      Copy(v2.View(),*this);
      return *this; 
    }

    template <int N> 
    inline const VectorView<T,I>& operator=(
        const SmallVectorComposite<RealType(T),N>& v2) const
    {
      TMVAssert(size() == N);
      v2.AssignToV(*this);
      return *this; 
    }

    template <int N> 
    inline const VectorView<T,I>& operator=(
        const SmallVectorComposite<ComplexType(T),N>& v2) const
    {
      TMVAssert(size() == N);
      TMVAssert(IsComplex(T()));
      v2.AssignToV(*this);
      return *this; 
    }

    //
    // Access Functions
    //

    typedef T value_type;
    typedef VIter<T> iterator;
    typedef CVIter<T> const_iterator;
    typedef VIter<T> reverse_iterator;
    typedef CVIter<T> const_reverse_iterator;
    typedef RefType(T) reference;

    inline reference operator[](int i) const 
    { 
      TMVAssert(i>=0 && i<int(size()));
      return ref(i); 
    }
    inline reference operator()(int i) const 
    { 
      TMVAssert(i>=0 && i<int(size()));
      return ref(i); 
    }

    inline iterator begin() const
    { return iterator(ptr(),step(),ct() FIRSTLAST ); }
    inline iterator end() const
    { return begin() + size(); }
    inline reverse_iterator rbegin() const
    {
      return reverse_iterator(ptr()+step()*(int(size())-1),-step(),ct() 
          FIRSTLAST ); 
    }
    inline reverse_iterator rend() const
    { return rbegin() + size(); }

    ListAssigner<T,iterator> inline operator=(ListInitClass)
    { return ListAssigner<T,iterator>(begin(),size()); }


    //
    // Modifying Functions
    //

    const VectorView<T,I>& Zero() const;

    const VectorView<T,I>& Clip(RealType(T) thresh) const;

    const VectorView<T,I>& SetAllTo(T x) const;

    const VectorView<T,I>& AddToAll(T x) const;

    const VectorView<T,I>& ConjugateSelf() const;

    const VectorView<T,I>& DoMakeBasis(int i) const;
    inline const VectorView<T,I>& MakeBasis(int i) const
    { TMVAssert(i>=0 && i<int(size())); return DoMakeBasis(i); }

    const VectorView<T,I>& DoSwap(int i1, int i2) const;
    inline const VectorView<T,I>& Swap(int i1, int i2) const
    { 
      TMVAssert(i1>=0 && i1<int(size()));
      TMVAssert(i2>=0 && i2<int(size()));
      return DoSwap(i1,i2);
    }

    const VectorView<T,I>& DoPermute(const int* p, 
        int i1, int i2) const;
    inline const VectorView<T,I>& Permute(const int* p, 
        int i1, int i2) const
    {
      TMVAssert(i1>=0 && i1 <= i2 && i2 <= int(size()));
      return DoPermute(p,i1,i2);
    }
    inline const VectorView<T,I>& Permute(const int* p) const
    { return DoPermute(p,0,size()); }

    const VectorView<T,I>& DoReversePermute(const int* p, 
        int i1, int i2) const;
    inline const VectorView<T,I>& ReversePermute(const int* p, 
        int i1, int i2) const
    {
      TMVAssert(i1>=0 && i1 <= i2 && i2 <= int(size()));
      return DoReversePermute(p,i1,i2);
    }
    inline const VectorView<T,I>& ReversePermute(const int* p) const
    { return ReversePermute(p,0,size()); }

    const VectorView<T,I>& ReverseSelf() const;

    const VectorView<T,I>& Sort(int* P, ADType ad=ASCEND, 
        COMPType comp=REAL_COMP) const;

    inline const VectorView<T,I>& Sort(ADType ad=ASCEND, 
        COMPType comp=REAL_COMP) const
    { Sort(0,ad,comp); return *this; }

    //
    // SubVector
    //

    inline VectorView<T,I> SubVector(int i1, int i2) const
    {
      TMVAssert(GenVector<T>::OKSubVector(i1,i2,1));
      return VectorView<T,I>(ptr()+i1*step(),(i2-i1),step(),
          ct() FIRSTLAST );
    }

    inline VectorView<T,I> SubVector(int i1, int i2, int istep) const
    {
      TMVAssert(GenVector<T>::OKSubVector(i1,i2,istep));
      return VectorView<T,I>(ptr()+i1*step(),(i2-i1)/istep,istep*step(),
          ct() FIRSTLAST );
    }

    inline VectorView<T,I> Reverse() const
    { 
      return VectorView<T,I>(ptr()+(int(size())-1)*step(),size(),-step(),
          ct() FIRSTLAST ); 
    }

    inline VectorView<T,I> View() const
    { return *this; }

    inline VectorView<T,I> Conjugate() const
    {
      return VectorView<T,I>(ptr(),size(),step(),
          ConjOf(T,ct()) FIRSTLAST ); 
    }

    inline VectorView<RealType(T),I> Real() const
    { 
      return VectorView<RealType(T),I>(
          reinterpret_cast<RealType(T)*>(ptr()),
          size(), IsReal(T()) ? step() : 2*step(), NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)
          ,reinterpret_cast<const RealType(T)*>(last)
#endif
          );
    }

    inline VectorView<RealType(T),I> Imag() const
    { 
      TMVAssert(IsComplex(T()));
      return VectorView<RealType(T),I>(
          reinterpret_cast<RealType(T)*>(ptr())+1,size(),2*step(), NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)+1
          ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
          );
    }

    inline VectorView<RealType(T),I> Flatten() const
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(step() == 1);
      return VectorView<RealType(T),I>(
          reinterpret_cast<RealType(T)*>(ptr()),
          IsReal(T()) ? size() : 2*size(), 1, NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)
          ,reinterpret_cast<const RealType(T)*>(last)+ (IsReal(T())?0:1)
#endif
          );
    }

    //
    // I/O
    //

    void Read(std::istream& is) const;

    // 
    // Iterator Typedefs
    //

    virtual inline size_t size() const { return itssize; }
    virtual inline const T* cptr() const { return itsv; }
    virtual inline T* ptr() const { return itsv; }
    virtual inline int step() const { return itsstep; }
    virtual inline ConjType ct() const { return itsct; }

    reference ref(int i) const;

  private:

    T*const itsv;
    const size_t itssize;
    const int itsstep;
    const ConjType itsct;

#ifdef TMVFLDEBUG
  public:
    const T*const first;
    const T*const last;
#endif

  }; // VectorView

  template <class T> 
  class VectorView<T,FortranStyle> : 
  public VectorView<T,CStyle>
  {
  public:

    //
    // Constructors 
    //

    inline VectorView(const VectorView<T,CStyle>& rhs) : 
      VectorView<T,CStyle>(rhs) {}

    inline VectorView(const VectorView<T,FortranStyle>& rhs) : 
      VectorView<T,CStyle>(rhs) {}

    inline VectorView(T* inv, size_t insize, int instep, ConjType inct
        PARAMFIRSTLAST(T) ) :
      VectorView<T,CStyle>(inv,insize,instep,inct
          FIRSTLAST1(_first,_last) ) {}

    virtual inline ~VectorView() {}


    //
    // Op =
    //

    inline const VectorView<T,FortranStyle>& operator=(
        const VectorView<T,FortranStyle>& v2) const
    { VectorView<T,CStyle>::operator=(v2); return *this; }

    inline const VectorView<T,FortranStyle>& operator=(
        const GenVector<RealType(T)>& v2) 
    { VectorView<T,CStyle>::operator=(v2); return *this; }

    inline const VectorView<T,FortranStyle>& operator=(
        const GenVector<ComplexType(T)>& v2) const
    { VectorView<T,CStyle>::operator=(v2); return *this; }

    template <class T2> 
    inline const VectorView<T,FortranStyle>& operator=(
        const GenVector<T2>& v2) const
    { VectorView<T,CStyle>::operator=(v2); return *this; }

    inline const VectorView<T,FortranStyle>& operator=(
        const AssignableToVector<RealType(T)>& v2) const
    { VectorView<T,CStyle>::operator=(v2); return *this; }

    inline const VectorView<T,FortranStyle>& operator=(
        const AssignableToVector<ComplexType(T)>& v2) const
    { VectorView<T,CStyle>::operator=(v2); return *this; }

    template <class T2, int N, IndexStyle I2> 
    inline const VectorView<T,FortranStyle>& operator=(
        const SmallVector<T2,N,I2>& v2) const
    { VectorView<T,CStyle>::operator=(v2); return *this; }

    template <int N> 
    inline const VectorView<T,FortranStyle>& operator=(
        const SmallVectorComposite<RealType(T),N>& v2) const
    { VectorView<T,CStyle>::operator=(v2); return *this; }

    template <int N> 
    inline const VectorView<T,FortranStyle>& operator=(
        const SmallVectorComposite<ComplexType(T),N>& v2) const
    { VectorView<T,CStyle>::operator=(v2); return *this; }

    ListAssigner<T,typename VectorView<T,CStyle>::iterator> 
    inline operator=(ListInitClass li)
    { return VectorView<T,CStyle>::operator=(li); }

    //
    // Access Functions
    //

    inline RefType(T) operator[](int i) const 
    { 
      TMVAssert(i>0 && i<=int(this->size()));
      return VectorView<T,CStyle>::ref(i-1); 
    }
    inline RefType(T) operator()(int i) const 
    { 
      TMVAssert(i>0 && i<=int(this->size()));
      return VectorView<T,CStyle>::ref(i-1); 
    }

    //
    // Modifying Functions
    //

    inline const VectorView<T,FortranStyle>& Zero() const 
    { VectorView<T,CStyle>::Zero(); return *this; }

    inline const VectorView<T,FortranStyle>& Clip(RealType(T) thresh) const
    { VectorView<T,CStyle>::Clip(thresh); return *this; }

    inline const VectorView<T,FortranStyle>& SetAllTo(T x) const
    { VectorView<T,CStyle>::SetAllTo(x); return *this; }

    inline const VectorView<T,FortranStyle>& AddToAll(T x) const
    { VectorView<T,CStyle>::AddToAll(x); return *this; }

    inline const VectorView<T,FortranStyle>& ConjugateSelf() const
    { VectorView<T,CStyle>::ConjugateSelf(); return *this; }

    inline const VectorView<T,FortranStyle>& MakeBasis(int i) const
    {
      TMVAssert(i>0 && i<=int(this->size()));
      VectorView<T,CStyle>::MakeBasis(i-1); 
      return *this; 
    }

    inline const VectorView<T,FortranStyle>& Swap(int i1, int i2) const
    {
      TMVAssert(i1>0 && i1<=int(this->size()));
      TMVAssert(i2>0 && i2<=int(this->size()));
      if (i1 != i2) VectorView<T,CStyle>::Swap(i1-1,i2-1);
      return *this;
    }

    inline const VectorView<T,FortranStyle>& Permute(const int* p, 
        int i1, int i2) const
    { 
      TMVAssert(i1>0 && i1 <= i2 && i2 <= int(this->size()));
      VectorView<T,CStyle>::Permute(p,i1-1,i2); 
      return *this; 
    }

    inline const VectorView<T,FortranStyle>& Permute(const int* p) const
    { VectorView<T,CStyle>::Permute(p); return *this; }

    inline const VectorView<T,FortranStyle>& ReversePermute(
        const int* p, int i1, int i2) const
    { 
      TMVAssert(i1>0 && i1 <= i2 && i2 <= int(this->size()));
      VectorView<T,CStyle>::ReversePermute(p,i1-1,i2); 
      return *this; 
    }

    inline const VectorView<T,FortranStyle>& ReversePermute(
        const int* p) const
    { VectorView<T,CStyle>::ReversePermute(p); return *this; }

    inline const VectorView<T,FortranStyle>& ReverseSelf() const
    { VectorView<T,CStyle>::ReverseSelf(); return *this; }

    inline const VectorView<T,FortranStyle>& Sort(
        int* P, ADType ad=ASCEND, COMPType comp=REAL_COMP) const
    { VectorView<T,CStyle>::Sort(P,ad,comp); return *this; }

    inline const VectorView<T,FortranStyle>& Sort(
        ADType ad=ASCEND, COMPType comp=REAL_COMP) const
    { VectorView<T,CStyle>::Sort(0,ad,comp); return *this; }

    //
    // SubVector
    //

    inline bool OKSubVector(int i1, int i2, int istep) const
    { 
      return ConstVectorView<T,FortranStyle>(*this).OKSubVector(
          i1,i2,istep); 
    }

    inline VectorView<T,FortranStyle> SubVector(int i1, int i2) const
    {
      TMVAssert(OKSubVector(i1,i2,1));
      return VectorView<T,CStyle>::SubVector(i1-1,i2);
    }

    inline VectorView<T,FortranStyle> SubVector(
        int i1, int i2, int istep) const
    {
      TMVAssert(OKSubVector(i1,i2,istep));
      return VectorView<T,CStyle>::SubVector(i1-1,i2-1+istep,istep);
    }

    inline VectorView<T,FortranStyle> Reverse() const
    { return VectorView<T,CStyle>::Reverse(); }

    inline VectorView<T,FortranStyle> View() const
    { return VectorView<T,CStyle>::View(); }

    inline VectorView<T,FortranStyle> Conjugate() const
    { return VectorView<T,CStyle>::Conjugate(); }

    inline VectorView<RealType(T),FortranStyle> Real() const
    { return VectorView<T,CStyle>::Real(); }

    inline VectorView<RealType(T),FortranStyle> Imag() const
    { return VectorView<T,CStyle>::Imag(); }

    inline VectorView<RealType(T),FortranStyle> Flatten() const
    { return VectorView<T,CStyle>::Flatten(); }

    inline T MinElement(int* iminout=0) const
    { return ConstVectorView<T,FortranStyle>(*this).MinElement(iminout); }

    inline T MaxElement(int* imaxout=0) const
    { return ConstVectorView<T,FortranStyle>(*this).MaxElement(imaxout); }

    inline RealType(T) MinAbsElement(int* iminout=0) const
    { return ConstVectorView<T,FortranStyle>(*this).MinAbsElement(iminout); }

    inline RealType(T) MaxAbsElement(int* imaxout=0) const
    { return ConstVectorView<T,FortranStyle>(*this).MaxAbsElement(imaxout); }

  }; // FortranStyle VectorView

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, IndexStyle I> 
  class Vector : 
    public GenVector<T>
  {
  public:

    //
    // Constructors
    //

#define NEW_SIZE(n) \
    itsv(new T[(n)]), itssize(n) DEFFIRSTLAST(itsv.get(),itsv.get()+n)

    explicit inline Vector(size_t n) : 
      NEW_SIZE(n)
    {
#ifdef TMVDEBUG
      SetAllTo(T(888));
#endif
    }

    inline Vector(size_t n, T val) : NEW_SIZE(n)
    {
      SetAllTo(val); 
    }

    inline Vector(size_t n, const T* vv) : NEW_SIZE(n)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      std::copy(vv,vv+n,itsv.get());
    }

    inline explicit Vector(const std::vector<T>& vv) : 
      NEW_SIZE(vv.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      std::copy(vv.begin(),vv.end(),itsv.get());
    }

    inline Vector(const Vector<T,I>& rhs) : 
      NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      std::copy(rhs.cptr(),rhs.cptr()+itssize,itsv.get());
    }

    template <IndexStyle I2> 
    inline Vector(const Vector<T,I2>& rhs) :
      NEW_SIZE(rhs.size())
    {
#ifdef TMVDEBUG
      SetAllTo(T(888));
#endif
      std::copy(rhs.cptr(),rhs.cptr()+itssize,itsv.get());
    }

    inline Vector(const GenVector<RealType(T)>& rhs) :
      NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      if (IsReal(T()) && rhs.step() == 1 && !rhs.isconj())
        std::copy(rhs.cptr(),rhs.cptr()+itssize,itsv.get());
      else rhs.AssignToV(View());
    }

    inline Vector(const GenVector<ComplexType(T)>& rhs) :
      NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      if (rhs.step() == 1 && !rhs.isconj())
        std::copy(rhs.cptr(),rhs.cptr()+itssize,itsv.get());
      else rhs.AssignToV(View());
    }

    template <class T2> 
    inline Vector(const GenVector<T2>& rhs) : 
      NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsReal(T2()) || IsComplex(T()));
      Copy(rhs,View()); 
    }

    inline Vector(const AssignableToVector<RealType(T)>& v2) :
      NEW_SIZE(v2.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      v2.AssignToV(View()); 
    }

    inline Vector(const AssignableToVector<ComplexType(T)>& v2) :
      NEW_SIZE(v2.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      v2.AssignToV(View()); 
    }

    template <class T2, int N, IndexStyle I2> 
    inline Vector(const SmallVector<T2,N,I2>& rhs) : 
      NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsReal(T2()) || IsComplex(T()));
      Copy(rhs.View(),View()); 
    }

    template <int N> 
    inline Vector(const SmallVectorComposite<RealType(T),N>& v2) : 
      NEW_SIZE(v2.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      v2.AssignToV(View()); 
    }

    template <int N> 
    inline Vector(const SmallVectorComposite<ComplexType(T),N>& v2) :
      NEW_SIZE(v2.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      v2.AssignToV(View()); 
    }

#undef NEW_SIZE

    virtual inline ~Vector() 
    {
#ifdef TMVDEBUG
      SetAllTo(T(999));
#endif
    }


    //
    // Op =
    //

    inline Vector<T,I>& operator=(Vector<T,I>& v2)
    { 
      TMVAssert(v2.size() == size());
      if (&v2 != this) 
        std::copy(v2.cptr(),v2.cptr()+itssize,itsv.get());
      return *this; 
    }

    template <IndexStyle I2> 
    inline Vector<T,I>& operator=(Vector<T,I2>& v2)
    { 
      TMVAssert(v2.size() == size());
      std::copy(v2.cptr(),v2.cptr()+itssize,itsv.get());
      return *this; 
    }

    inline Vector<T,I>& operator=(const GenVector<RealType(T)>& v2) 
    { 
      TMVAssert(v2.size() == size());
      v2.AssignToV(View());
      return *this; 
    }

    inline Vector<T,I>& operator=(const GenVector<ComplexType(T)>& v2) 
    { 
      TMVAssert(v2.size() == size());
      TMVAssert(IsComplex(T()));
      v2.AssignToV(View());
      return *this; 
    }

    template <class T2> 
    inline Vector<T,I>& operator=(const GenVector<T2>& v2) 
    { 
      TMVAssert(v2.size() == size());
      View() = v2; 
      return *this; 
    }

    inline Vector<T,I>& operator=(
        const AssignableToVector<RealType(T)>& v2)
    { 
      TMVAssert(v2.size() == size());
      v2.AssignToV(View()); 
      return *this; 
    }


    inline Vector<T,I>& operator=(
        const AssignableToVector<ComplexType(T)>& v2)
    { 
      TMVAssert(v2.size() == size());
      TMVAssert(IsComplex(T()));
      v2.AssignToV(View()); 
      return *this; 
    }

    template <class T2, int N, IndexStyle I2> 
    inline Vector<T,I>& operator=(const SmallVector<T2,N,I2>& v2) 
    { 
      TMVAssert(v2.size() == size());
      View() = v2.View(); 
      return *this; 
    }

    template <int N> 
    inline Vector<T,I>& operator=(
        const SmallVectorComposite<RealType(T),N>& v2)
    { 
      TMVAssert(N == size());
      v2.AssignToV(View()); 
      return *this; 
    }


    template <int N> 
    inline Vector<T,I>& operator=(
        const SmallVectorComposite<ComplexType(T),N>& v2)
    { 
      TMVAssert(N == size());
      TMVAssert(IsComplex(T()));
      v2.AssignToV(View()); 
      return *this; 
    }


    //
    // Access Functions
    //

    typedef T value_type;
    typedef VIt<T,Unit,NonConj> iterator;
    typedef CVIt<T,Unit,NonConj> const_iterator;
    typedef VIt<T,Step,NonConj> reverse_iterator;
    typedef CVIt<T,Step,NonConj> const_reverse_iterator;
    typedef T& reference;

    inline const_iterator begin() const
    { return const_iterator(cptr(),1); }
    inline const_iterator end() const
    { return begin() + size(); }
    inline const_reverse_iterator rbegin() const
    { return const_reverse_iterator(cptr()+(int(size())-1),-1); }
    inline const_reverse_iterator rend() const
    { return rbegin() + size(); }
    inline iterator begin()
    { return iterator(ptr(),1); }
    inline iterator end()
    { return begin() + size(); }
    inline reverse_iterator rbegin()
    { return reverse_iterator(ptr()+(int(size())-1),-1); }
    inline reverse_iterator rend()
    { return rbegin() + size(); }

    ListAssigner<T,iterator> inline operator=(ListInitClass)
    { return ListAssigner<T,iterator>(begin(),size()); }

    inline T operator[](int i) const 
    { 
      if (I == CStyle) { TMVAssert(i>=0 && i<int(size())); return cref(i); } 
      else { TMVAssert(i>0 && i<=int(size())); return cref(i-1); }
    }
    inline T operator()(int i) const 
    { 
      if (I == CStyle) { TMVAssert(i>=0 && i<int(size())); return cref(i); } 
      else { TMVAssert(i>0 && i<=int(size())); return cref(i-1); }
    }

    inline T& operator[](int i) 
    { 
      if (I == CStyle) { TMVAssert(i>=0 && i<int(size())); return ref(i); } 
      else { TMVAssert(i>0 && i<=int(size())); return ref(i-1); }
    }
    inline T& operator()(int i) 
    { 
      if (I == CStyle) { TMVAssert(i>=0 && i<int(size())); return ref(i); } 
      else { TMVAssert(i>0 && i<=int(size())); return ref(i-1); }
    }

    //
    // Modifying Functions
    //

    Vector<T,I>& Zero();

    Vector<T,I>& Clip(RealType(T) thresh);

    Vector<T,I>& SetAllTo(T x);

    Vector<T,I>& AddToAll(T x);

    Vector<T,I>& ConjugateSelf();

    Vector<T,I>& DoMakeBasis(int i);
    inline Vector<T,I>& MakeBasis(int i)
    { 
      if (I == CStyle) TMVAssert(i>=0 && i<int(size()));
      else TMVAssert(i>0 && i<=int(size()));
      return DoMakeBasis(i);
    }

    Vector<T,I>& DoSwap(int i1, int i2);
    inline Vector<T,I>& Swap(int i1, int i2)
    {
      if (I == CStyle) 
        TMVAssert(i1>=0 && i1<int(size()) && i2>=0 && i2<int(size()));
      else TMVAssert(i1>0 && i1<=int(size()) && i2>0 && i2<=int(size()));
      return DoSwap(i1,i2);
    }

    inline Vector<T,I>& Permute(const int* p, int i1, int i2)
    {
      if (I==FortranStyle) TMVAssert(i1>0);
      TMVAssert(i1>=0 && i1 <= i2 && i2 <= int(size()));
      View().Permute(p,i1,i2); return *this; 
    }
    inline Vector<T,I>& Permute(const int* p) 
    { View().Permute(p); return *this; }
    inline Vector<T,I>& ReversePermute(const int* p, int i1, int i2)
    {
      if (I==FortranStyle) TMVAssert(i1>0);
      TMVAssert(i1>=0 && i1 <= i2 && i2 <= int(size()));
      View().ReversePermute(p,i1,i2); return *this; 
    }
    inline Vector<T,I>& ReversePermute(const int* p) 
    { View().ReversePermute(p,0,size()); return *this; }

    inline Vector<T,I>& ReverseSelf()
    { View().ReverseSelf(); return *this; }

    inline Vector<T,I>& Sort(int* P, ADType ad=ASCEND, 
        COMPType comp=REAL_COMP) 
    { View().Sort(P,ad,comp); return *this; }

    inline Vector<T,I>& Sort(ADType ad=ASCEND, 
        COMPType comp=REAL_COMP) 
    { View().Sort(0,ad,comp); return *this; }

    //
    // SubVector
    //

    inline ConstVectorView<T,I> SubVector(int i1, int i2) const
    {
      TMVAssert(View().OKSubVector(i1,i2,1));
      if (I==FortranStyle) --i1;
      return ConstVectorView<T,I>(cptr()+i1,i2-i1,1,NonConj);
    }

    inline VectorView<T,I> SubVector(int i1, int i2)
    {
      TMVAssert(View().OKSubVector(i1,i2,1));
      if (I==FortranStyle) --i1;
      return VectorView<T,I>(ptr()+i1,i2-i1,1,NonConj FIRSTLAST );
    }

    inline ConstVectorView<T,I> SubVector(int i1, int i2, int istep) const
    {
      TMVAssert(View().OKSubVector(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2 += istep-1; }
      return ConstVectorView<T,I>(cptr()+i1,(i2-i1)/istep,istep,NonConj);
    }

    inline VectorView<T,I> SubVector(int i1, int i2, int istep)
    {
      TMVAssert(View().OKSubVector(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2 += istep-1; }
      return VectorView<T,I>(ptr()+i1,(i2-i1)/istep,istep,NonConj 
          FIRSTLAST );
    }

    inline ConstVectorView<T,I> Reverse() const
    { return ConstVectorView<T,I>(cptr()+size()-1,size(),-1,NonConj); }

    inline VectorView<T,I> Reverse()
    { return VectorView<T,I>(ptr()+size()-1,size(),-1,NonConj FIRSTLAST ); }

    inline ConstVectorView<T,I> View() const
    { return ConstVectorView<T,I>(cptr(),size(),1,NonConj); }

    inline VectorView<T,I> View()
    { return VectorView<T,I>(ptr(),size(),1,NonConj FIRSTLAST ); }

    inline ConstVectorView<T,I> Conjugate() const
    { return ConstVectorView<T,I>(cptr(),size(),1,ConjOf(T,NonConj)); }

    inline VectorView<T,I> Conjugate()
    { return VectorView<T,I>(ptr(),size(),1,ConjOf(T,NonConj) FIRSTLAST ); }

    inline ConstVectorView<RealType(T),I> Real() const
    { return View().Real(); }

    inline ConstVectorView<RealType(T),I> Imag() const
    { return View().Imag(); }

    inline ConstVectorView<RealType(T),I> Flatten() const
    { return View().Flatten(); }

    inline VectorView<RealType(T),I> Real()
    { return View().Real(); }

    inline VectorView<RealType(T),I> Imag()
    { return View().Imag(); }

    inline VectorView<RealType(T),I> Flatten()
    { return View().Flatten(); }

    inline T MinElement(int* iminout=0) const
    { return View().MinElement(iminout); }

    inline T MaxElement(int* imaxout=0) const
    { return View().MaxElement(imaxout); }

    inline RealType(T) MinAbsElement(int* iminout=0) const
    { return View().MinAbsElement(iminout); }

    inline RealType(T) MaxAbsElement(int* imaxout=0) const
    { return View().MaxAbsElement(imaxout); }

    // 
    // I/O
    //

    inline void Read(std::istream& is)
    { View().Read(is); }

    virtual inline size_t size() const { return itssize; }
    virtual inline const T* cptr() const { return itsv.get(); }
    virtual inline T* ptr() { return itsv.get(); }
    virtual inline int step() const { return 1; }
    virtual inline ConjType ct() const { return NonConj; }
    inline bool isconj() const { return false; }

    inline T cref(int i) const
    { return itsv.get()[i]; }

    inline T& ref(int i)
    { return itsv.get()[i]; }

  private:

    auto_array<T> itsv;
    const size_t itssize;

#ifdef TMVFLDEBUG
  public:
    const T*const first;
    const T*const last;
#endif

  }; // Vector

  //
  // Special Constructors
  //

  template <class T, IndexStyle I> 
  Vector<T,I> DoBasisVector(size_t n, int i);
  template <class T, IndexStyle I> 
  inline Vector<T,I> BasisVector(size_t n, int i)
  { 
    if (I == CStyle) { TMVAssert(i>=0 && i<n); } 
    else { TMVAssert(i>0 && i<=n); }
    return DoBasisVector<T,I>(n,i); 
  }
  template <class T> 
  inline Vector<T,CStyle> BasisVector(size_t n, int i)
  { return DoBasisVector<T,CStyle>(n,i); }

  template <IndexStyle I, class T> 
  inline VectorView<T,I> VectorViewOf(T* v, size_t size)
  { return VectorView<T,I>(v,size,1,NonConj FIRSTLAST1(v,v+size)); }

  template <IndexStyle I, class T> 
  inline ConstVectorView<T,I> VectorViewOf(const T* v, size_t size)
  { return ConstVectorView<T,I>(v,size,1,NonConj); }

  template <class T> 
  inline VectorView<T,CStyle> VectorViewOf(T* v, size_t size)
  {
    return VectorView<T,CStyle>(v,size,1,NonConj 
        FIRSTLAST1(v,v+size)); 
  }

  template <class T> 
  inline ConstVectorView<T,CStyle> VectorViewOf(const T* v, size_t size)
  { return ConstVectorView<T,CStyle>(v,size,1,NonConj); }

  //
  // Copy Vectors
  //

  inline bool ShouldReverse(const int step1, const int step2)
  {
    return ( (step2 < 0 && (step1 != 1 || step2 == -1)) ||
        (step1 == -1 && step2 != 1) );
  }

  template <class T> 
  void DoCopySameType(const GenVector<T>& v1, const VectorView<T>& v2);

  template <class T> 
  inline void DoCopy(const GenVector<T>& v1, const VectorView<T>& v2)
  { if (!v1.SameAs(v2)) DoCopySameType(v1,v2); }

  template <class T, class T1> 
  inline void DoCopyDiffType(const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    TMVAssert(IsReal(T1()) || IsComplex(T()));
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(!v2.SameAs(v1));

    const T1* v1ptr = v1.cptr();
    T* v2ptr = v2.ptr();
    const int step1 = v1.step();
    const int step2 = v2.step();

    if (step1 == 1 && step2 == 1)
      for(int i=v2.size();i>0;--i,++v1ptr,++v2ptr) {
#ifdef TMVFLDEBUG
        TMVAssert(v2ptr >= v2.first);
        TMVAssert(v2ptr < v2.last);
#endif
        *v2ptr = *v1ptr;
      }
    else
      for(int i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) {
#ifdef TMVFLDEBUG
        TMVAssert(v2ptr >= v2.first);
        TMVAssert(v2ptr < v2.last);
#endif
        *v2ptr = *v1ptr;
      }
  }

  template <class T, class T1> 
  inline void DoCopy(const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    TMVAssert(IsReal(T1()) || IsComplex(T()));
    DoCopyDiffType(v1,v2);
  }

  template <class T> 
  inline void DoCopy(const GenVector<std::complex<T> >&, const VectorView<T>&)
  { TMVAssert(FALSE); }

  template <class T, class T1> 
  inline void Copy(const GenVector<T1>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(IsReal(T1()) || IsComplex(T()));
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.step()!=0 || v1.step() == 0);

    if (v1.size() > 0) {
      if (ShouldReverse(v1.step(),v2.step()))
        Copy(v1.Reverse(),v2.Reverse());
      else {
        if (v1.isconj()) {
          if (v2.isconj()) {
            DoCopy(v1.Conjugate(),v2.Conjugate());
          } else {
            DoCopy(v1.Conjugate(),v2);
            v2.ConjugateSelf();
          }
        } else {
          if (v2.isconj()) {
            DoCopy(v1,v2.Conjugate());
            v2.ConjugateSelf();
          }
          else DoCopy(v1,v2);
        }
      }
    }
  }


  //
  // Swap Vectors
  //

  template <class T> 
  void Swap(const VectorView<T>& v1, const VectorView<T>& v2);
  template <class T, IndexStyle I> 
  inline void Swap(const VectorView<T>& v1, Vector<T,I>& v2)
  { Swap(v1.View(),v2.View()); }
  template <class T, IndexStyle I> 
  inline void Swap(Vector<T,I>& v1, const VectorView<T>& v2) 
  { Swap(v1.View(),v2.View()); }
  template <class T, IndexStyle I1, IndexStyle I2> 
  inline void Swap(Vector<T,I1>& v1, Vector<T,I2>& v2)
  { Swap(v1.View(),v2.View()); }

  //
  // Functions of Vectors
  //

  template <class T> 
  inline RealType(T) Norm(const GenVector<T>& v)
  { return v.Norm(); }

  template <class T> 
  inline RealType(T) Norm1(const GenVector<T>& v)
  { return v.Norm1(); }

  template <class T> 
  inline RealType(T) NormSq(const GenVector<T>& v)
  { return v.NormSq(); }

  template <class T> 
  inline RealType(T) Norm2(const GenVector<T>& v)
  { return v.Norm2(); }

  template <class T> 
  inline RealType(T) NormInf(const GenVector<T>& v)
  { return v.NormInf(); }

  template <class T> 
  inline T SumElements(const GenVector<T>& v)
  { return v.SumElements(); }

  template <class T> 
  inline RealType(T) SumAbsElements(const GenVector<T>& v)
  { return v.SumAbsElements(); }

  template <class T> 
  inline T MinElement(const GenVector<T>& v)
  { return v.MinElement(); }

  template <class T, IndexStyle I> 
  inline T MinElement(const ConstVectorView<T,I>& v)
  { return v.MinElement(); }

  template <class T, IndexStyle I> 
  inline T MinElement(const VectorView<T,I>& v)
  { return v.MinElement(); }

  template <class T, IndexStyle I> 
  inline T MinElement(const Vector<T,I>& v)
  { return v.MinElement(); }

  template <class T> 
  inline T MaxElement(const GenVector<T>& v)
  { return v.MaxElement(); }

  template <class T, IndexStyle I> 
  inline T MaxElement(const ConstVectorView<T,I>& v)
  { return v.MaxElement(); }

  template <class T, IndexStyle I> 
  inline T MaxElement(const VectorView<T,I>& v)
  { return v.MaxElement(); }

  template <class T, IndexStyle I> 
  inline T MaxElement(const Vector<T,I>& v)
  { return v.MaxElement(); }

  template <class T> 
  inline RealType(T) MinAbsElement(const GenVector<T>& v)
  { return v.MinAbsElement(); }

  template <class T, IndexStyle I> 
  inline RealType(T) MinAbsElement(const ConstVectorView<T,I>& v)
  { return v.MinAbsElement(); }

  template <class T, IndexStyle I> 
  inline RealType(T) MinAbsElement(const VectorView<T,I>& v)
  { return v.MinAbsElement(); }

  template <class T, IndexStyle I> 
  inline RealType(T) MinAbsElement(const Vector<T,I>& v)
  { return v.MinAbsElement(); }

  template <class T> 
  inline RealType(T) MaxAbsElement(const GenVector<T>& v)
  { return v.MaxAbsElement(); }

  template <class T, IndexStyle I> 
  inline RealType(T) MaxAbsElement(const ConstVectorView<T,I>& v)
  { return v.MaxAbsElement(); }

  template <class T, IndexStyle I> 
  inline RealType(T) MaxAbsElement(const VectorView<T,I>& v)
  { return v.MaxAbsElement(); }

  template <class T, IndexStyle I> 
  inline RealType(T) MaxAbsElement(const Vector<T,I>& v)
  { return v.MaxAbsElement(); }

  template <class T> 
  inline ConstVectorView<T> Conjugate(const GenVector<T>& v)
  { return v.Conjugate(); }

  template <class T, IndexStyle I> 
  inline ConstVectorView<T,I> Conjugate(const ConstVectorView<T,I>& v)
  { return v.Conjugate(); }

  template <class T, IndexStyle I> 
  inline ConstVectorView<T,I> Conjugate(const Vector<T,I>& v)
  { return v.Conjugate(); }

  template <class T, IndexStyle I> 
  inline VectorView<T,I> Conjugate(const VectorView<T,I>& v)
  { return v.Conjugate(); }

  template <class T, IndexStyle I> 
  inline VectorView<T,I> Conjugate(Vector<T,I>& v)
  { return v.Conjugate(); }


  //
  // Vector ==, != Vector
  //

  template <class T1, class T2> 
  bool operator==(const GenVector<T1>& v1, const GenVector<T2>& v2);

  template <class T1, class T2> 
  inline bool operator!=(const GenVector<T1>& v1, const GenVector<T2>& v2)
  { return !(v1 == v2); }

  //
  // I/O
  //

  template <class T> 
  inline std::ostream& operator<<(std::ostream& os, const GenVector<T>& v)
  { v.Write(os); return os;}

  template <class T> 
  std::istream& operator>>(std::istream& is, const VectorView<T>& v);

  template <class T, IndexStyle I> 
  inline std::istream& operator>>(std::istream& is, Vector<T,I>& v)
  { return is >> v.View(); }

  template <class T, IndexStyle I> 
  std::istream& operator>>(std::istream& is, auto_ptr<Vector<T,I> >& v);

} // namespace tmv

#endif
