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
//    Vector<T>(const valarray<T>& vv)
//    Vector<T>(size_t n, const T* vv)
//    Vector<T>(const vector<T>& vv)
//        Makes a vector which copies the elements of vv
//        For the second one, n specifies the length of the vector
//
// Special Constructors
//
//    BasisVector(size_t n, size_t i)
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
//    T& operator[](size_t i)
//    T& operator()(size_t i)
//    T operator[](size_t i) const
//    T operator()(size_t i) const
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
//    Vector& MakeBasis(size_t i)
//        Set all elements to 0, except v(i) = 1
//
//    Vector& Swap(size_t i1, size_t i2)
//        Swap elements v(i1) and v(i2)
//
//    Vector& Permute(const size_t* p)
//    Vector& Permute(const size_t* p, size_t i1, size_t i2)
//        Perform a series of swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (i1,p[i1])...(i2-1,p[i2-1])
//    Vector& ReversePermute(const size_t* p)
//    Vector& ReversePermute(const size_t* p, size_t i1, size_t i2)
//        The same but perform the swaps in reverse order.
//
//    Swap(v1,v2)
//        Swap vectors v1 and v2
//
//    Vector& ReverseSelf()
//        Reverse the order of the elements of v
//
//    Sort(size_t* P, AD, COMP)
//        Sorts the vector, returning the swaps required in P.
//        If you do not care about P, you may set P = 0 on entry.
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
//        Returns the square of Norm()
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
//    MaxElement(v,size_t* imax)
//        Returns the maximum value of any element in the Vector
//        On return, *imax_value holds the index of this element.
//        The parameter max_value can be omitted if it is not desired.
//        As "max" doesn't make sense for complex values, for these
//        we use just the real components.
//
//    MinElement(v,size_t* imin)
//        Returns the minimum value of any element in the Vector
//        On return, *imin_value holds the index of this element.
//        Again, the parameter min_value can be omitted.
//
//    MaxAbsElement(v,size_t* imax)
//        The same as MaxElement, except absolute values are used
//
//    MinAbsElement(v,size_t* imin)
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
//        If you do not know the size of the Vector to be read in,
//        this version will create a new Vector of the correct size.
//        So you should subsequently delete the Vector.
//
//


#ifndef TMV_Vector_H
#define TMV_Vector_H

#include "TMV_VIt.h"

namespace tmv {

  enum IndexStyle { CStyle, FortranStyle };

  template <class T> class GenVector;
  template <class T, IndexStyle I=CStyle> class ConstVectorView;
  template <class T, IndexStyle I=CStyle> class VectorView;
  template <class T, IndexStyle I=CStyle> class Vector;
  template <class T> class VectorComposite;

  template <class T, class T1> inline void Copy(
      const GenVector<T1>& v1, const VectorView<T>& v2);

#define RefType(T) typename VIter<T>::reference
  template <class T> inline T& REF(T* vi, ConjItType ct)
  { return *vi; }
  template <class T> inline VarConjRef<complex<T> > REF(
      complex<T>* vi, ConjItType ct)
  { return VarConjRef<complex<T> >(*vi,ct); }

  enum ADType { ASCEND, DESCEND };
  enum COMPType { REAL_COMP, ABS_COMP, IMAG_COMP, ARG_COMP };

  template <class T> class GenVector
  {
    public:

      //
      // Constructor
      //

      GenVector(ConjItType ct) : itsct(ct) {}
      GenVector(const GenVector<T>& rhs) : itsct(rhs.itsct) {}
      virtual ~GenVector() {}


      //
      // Access Functions
      //

      virtual size_t size() const =0;

      inline T operator[](size_t i) const 
      { 
	TMVAssert(i<size());
	return cref(i); 
      }
      inline T operator()(size_t i) const 
      { 
	TMVAssert(i<size());
	return cref(i); 
      }

      typedef CVIter<T> const_iterator;
      typedef CVIter<T> const_reverse_iterator;

      inline CVIter<T> begin() const 
      { return CVIter<T>(cptr(),step(),itsct); }
      inline CVIter<T> end() const 
      { return begin()+size(); }
      inline CVIter<T> rbegin() const 
      { return CVIter<T>(cptr()+step()*(size()-1),-step(),itsct); }
      inline CVIter<T> rend() const 
      { return rbegin()+size(); }

      template <class T2> inline bool SameStorageAs(
	  const GenVector<T2>& v2) const
      { return false; }

      inline bool SameStorageAs(const GenVector<T>& v2) const
      { return (cptr()==v2.cptr()); }

      template <class T2> inline bool SameAs(const GenVector<T2>& v2) const
      { return false; }

      inline bool SameAs(const GenVector<T>& v2) const
      {
	return (this == &v2 || (cptr()==v2.cptr() && 
	    size()==v2.size() && step()==v2.step() && itsct==v2.itsct));
      }

      //
      // SubVector
      //

      bool OKSubVector(int i1, int i2, int istep) const;

      inline ConstVectorView<T> SubVector(int i1, int i2) const
      {
	TMVAssert(OKSubVector(i1,i2,1));
	return ConstVectorView<T>(cptr()+i1*step(),i2-i1,step(),itsct);
      }

      inline ConstVectorView<T> SubVector(int i1, int i2, int istep) const
      {
	TMVAssert(OKSubVector(i1,i2,istep));
	return ConstVectorView<T>(cptr()+i1*step(),
	    (i2-i1)/istep,istep*step(),itsct);
      }

      inline ConstVectorView<T> Reverse() const
      { 
	return ConstVectorView<T>(cptr()+(size()-1)*step(),
	    size(),-step(),itsct); 
      }

      inline ConstVectorView<T> View() const
      { return ConstVectorView<T>(cptr(),size(),step(),itsct); }

      inline ConstVectorView<T> Conjugate() const
      { return ConstVectorView<T>(cptr(),size(),step(),ConjOf(T,itsct)); }

      inline ConstVectorView<RealType(T)> Real() const
      { 
	return ConstVectorView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),
	    size(), IsReal(T()) ? step() : 2*step(), NonConj);
      }

      inline ConstVectorView<RealType(T)> Imag() const
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(itsct == NonConj);
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

      inline RealType(T) Norm() const // = Norm2
      { return Norm2(); }

      RealType(T) NormSq() const; // Norm()^2

      inline RealType(T) Norm1() const // sum_i |v_i|
      { return SumAbsElements(); }

      RealType(T) Norm2() const; // sqrt( sum_i |v_i|^2 )

      inline RealType(T) NormInf() const // max_i |v_i|
      { return size() > 0 ? MaxAbsElement() : RealType(T)(0); }

      T SumElements() const;

      RealType(T) SumAbsElements() const;

      T MinElement(size_t* iminout=0) const;

      T MaxElement(size_t* imaxout=0) const;

      RealType(T) MinAbsElement(size_t* iminout=0) const;

      RealType(T) MaxAbsElement(size_t* imaxout=0) const;

      //
      // I/O
      //

      void Write(ostream& fout) const;
      void Write(ostream& fout, RealType(T) minnonzero) const;

      virtual const T* cptr() const =0;
      virtual int step() const =0;
      virtual inline ConjItType ct() const { return itsct; }
      inline bool isconj() const 
      { 
	TMVAssert(IsComplex(T()) || itsct==NonConj);
	return (IsComplex(T()) && itsct==Conj);
      }

    protected:

      const ConjItType itsct;

      virtual T cref(size_t i) const;

    private:

      inline void operator=(const GenVector<T>&) { TMVAssert(false); }

  }; // GenVector

  template <class T, IndexStyle I> class ConstVectorView : 
    public GenVector<T>
  {
    public:

      ConstVectorView(const ConstVectorView<T,I>& rhs) : 
	GenVector<T>(rhs),
	itsv(rhs.itsv), itssize(rhs.itssize), itsstep(rhs.itsstep) {}
      ConstVectorView(const GenVector<T>& rhs) : 
	GenVector<T>(rhs),
	itsv(rhs.cptr()), itssize(rhs.size()), itsstep(rhs.step()) {}
      ConstVectorView(const T* inv, size_t insize, int instep, 
	  ConjItType inct) : 
	GenVector<T>(inct), itsv(inv), itssize(insize), itsstep(instep) {}
      ~ConstVectorView() {}

      inline size_t size() const { return itssize; }
      inline const T* cptr() const { return itsv; }
      inline int step() const { return itsstep; }

    private:

      const T*const itsv;
      const size_t itssize;
      const size_t itsstep;

      inline void operator=(const ConstVectorView<T,I>&) { TMVAssert(false); }

  }; // ConstVectorView

  template <class T> class ConstVectorView<T,FortranStyle> : 
    public ConstVectorView<T,CStyle>
  {
    public:

      ConstVectorView(const ConstVectorView<T,FortranStyle>& rhs) : 
	ConstVectorView<T,CStyle>(rhs) {}
      ConstVectorView(const ConstVectorView<T,CStyle>& rhs) : 
	ConstVectorView<T,CStyle>(rhs) {}
      ConstVectorView(const GenVector<T>& rhs) : 
	ConstVectorView<T,CStyle>(rhs) {}
      ConstVectorView(const T* inv, size_t insize, int instep, 
	  ConjItType inct) : 
	ConstVectorView<T,CStyle>(inv,insize,instep,inct) {}
      ~ConstVectorView() {}

      inline T operator[](size_t i) const 
      { 
	TMVAssert(i>0 && i<=size());
	return cref(i-1); 
      }
      inline T operator()(size_t i) const 
      { 
	TMVAssert(i>0 && i<=size());
	return cref(i-1); 
      }

      bool OKSubVector(int i1, int i2, int istep) const;

      inline ConstVectorView<T,FortranStyle> SubVector(int i1, int i2) const
      {
	TMVAssert(OKSubVector(i1,i2,1));
	return GenVector<T>::SubVector(i1-1,i2);
      }

      inline ConstVectorView<T,FortranStyle> SubVector(int i1, int i2, int istep) const
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

      inline T MinElement(size_t* iminout=0) const
      { 
	T temp = GenVector<T>::MinElement(iminout);
	if (iminout) ++(*iminout);
	return temp;
      }

      inline T MaxElement(size_t* imaxout=0) const
      { 
	T temp = GenVector<T>::MaxElement(imaxout);
	if (imaxout) ++(*imaxout);
	return temp;
      }

      inline RealType(T) MinAbsElement(size_t* iminout=0) const
      { 
	RealType(T) temp = GenVector<T>::MinAbsElement(iminout);
	if (iminout) ++(*iminout);
	return temp;
      }

      inline RealType(T) MaxAbsElement(size_t* imaxout=0) const
      { 
	RealType(T) temp = GenVector<T>::MaxAbsElement(imaxout);
	if (imaxout) ++(*imaxout);
	return temp;
      }

      //using ConstVectorView<T,CStyle>::size;
      using GenVector<T>::size;

    protected:

      using GenVector<T>::cref;

    private:

      inline void operator=(const ConstVectorView<T,FortranStyle>&) 
      { TMVAssert(false); }

  }; // FortranStyle ConstVectorView

  template <class T, IndexStyle I> class VectorView : 
    public GenVector<T>
  {
    public:

      //
      // Constructors 
      //

      VectorView(const VectorView<T,I>& rhs) : 
	GenVector<T>(rhs),
	itsv(rhs.itsv), itssize(rhs.itssize), itsstep(rhs.itsstep)
	DEFFIRSTLAST(rhs.first,rhs.last) {}

      VectorView(T* inv, size_t insize, int instep, ConjItType inct 
	  PARAMFIRSTLAST(T) ) :
	GenVector<T>(inct), itsv(inv), itssize(insize), itsstep(instep)
	  DEFFIRSTLAST(_first,_last) {}

      ~VectorView() {}


      //
      // Op =
      //

      inline const VectorView<T,I>& operator=(const VectorView<T,I>& v2) const
      { 
	TMVAssert(size() == v2.size());
	if (!SameAs(v2)) Copy(v2,*this); 
	return *this; 
      }

      inline const VectorView<T,I>& operator=(const GenVector<T>& v2) const
      { 
	TMVAssert(size() == v2.size());
	if (!SameAs(v2)) Copy(v2,*this);
	return *this; 
      }

      template <class T2> inline const VectorView<T,I>& operator=(
	  const GenVector<T2>& v2) const
      { 
	TMVAssert(size() == v2.size());
	Copy(v2,*this);
	return *this; 
      }

      inline const VectorView<T,I>& operator=(
	  const VectorComposite<T>& vcomp) const
      { vcomp.AssignTo(*this); return *this; }


      //
      // Access Functions
      //

      inline RefType(T) operator[](size_t i) const 
      { 
	TMVAssert(i<size());
	return ref(i); 
      }
      inline RefType(T) operator()(size_t i) const 
      { 
	TMVAssert(i<size());
	return ref(i); 
      }

      typedef VIter<T> iterator;
      typedef CVIter<T> const_iterator;
      typedef VIter<T> reverse_iterator;
      typedef CVIter<T> const_reverse_iterator;
      typedef RefType(T) reference;

      inline VIter<T> begin() const 
      { return VIter<T>(ptr(),step(),ct() FIRSTLAST ); }
      inline VIter<T> end() const 
      { return begin() + size(); }
      inline VIter<T> rbegin() const 
      { return VIter<T>(ptr()+step()*(size()-1),-step(),ct() FIRSTLAST ); }
      inline VIter<T> rend() const 
      { return rbegin()+size(); }

      //
      // Modifying Functions
      //

      const VectorView<T,I>& Zero() const 
      {
	if (IsComplex(T()) && step() == 1) Flatten().SetAllTo(0);
	else SetAllTo(0); 
	return *this;
      }
      const VectorView<T,I>& Clip(RealType(T) thresh) const;
      const VectorView<T,I>& SetAllTo(T x) const;
      const VectorView<T,I>& AddToAll(T x) const;
      const VectorView<T,I>& ConjugateSelf() const;
      const VectorView<T,I>& MakeBasis(size_t i) const
      { TMVAssert(i < size()); Zero(); ref(i) = T(1); return *this; }
      const VectorView<T,I>& Swap(size_t i1, size_t i2) const
      { 
	TMVAssert(i1 < size());
	TMVAssert(i2 < size());
	if (i1 != i2) swap(ref(i1),ref(i2));
	return *this;
      }
      const VectorView<T,I>& Permute(const size_t* p, 
	  size_t i1, size_t i2) const;
      inline const VectorView<T,I>& Permute(const size_t* p) const
      { return Permute(p,0,size()); }
      const VectorView<T,I>& ReversePermute(const size_t* p, 
	  size_t i1, size_t i2) const;
      inline const VectorView<T,I>& ReversePermute(const size_t* p) const
      { return ReversePermute(p,0,size()); }
      const VectorView<T,I>& ReverseSelf() const;
      const VectorView<T,I>& Sort(size_t* P, ADType ad=ASCEND, 
	  COMPType comp=REAL_COMP) const;

      //
      // SubVector
      //

      inline VectorView<T,I> SubVector(int i1, int i2) const
      {
	TMVAssert(GenVector<T>::OKSubVector(i1,i2,1));
	return VectorView<T,I>(ptr()+i1*step(),(i2-i1),step(),ct() FIRSTLAST );
      }

      inline VectorView<T,I> SubVector(int i1, int i2, int istep) const
      {
	TMVAssert(GenVector<T>::OKSubVector(i1,i2,istep));
	return VectorView<T,I>(
	    ptr()+i1*step(),(i2-i1)/istep,istep*step(),ct() FIRSTLAST );
      }

      inline VectorView<T,I> Reverse() const
      { 
	return VectorView<T,I>(
	    ptr()+(size()-1)*step(),size(),-step(),ct() FIRSTLAST ); 
      }

      inline VectorView<T,I> View() const
      { return *this; }

      inline VectorView<T,I> Conjugate() const
      {
	return VectorView<T,I>(ptr(),size(),step(),ConjOf(T,ct()) FIRSTLAST ); 
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
	TMVAssert(ct() == NonConj);
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

      void Read(istream& fin) const;

      // 
      // Iterator Typedefs
      //

      inline size_t size() const { return itssize; }
      inline const T* cptr() const { return itsv; }
      inline T* ptr() const { return itsv; }
      inline int step() const { return itsstep; }
      using GenVector<T>::ct;
      using GenVector<T>::isconj;

    protected:

      RefType(T) ref(size_t i) const;

    private:

      T*const itsv;
      const size_t itssize;
      const int itsstep;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
#endif


  }; // VectorView

  template <class T> class VectorView<T,FortranStyle> : 
    public VectorView<T,CStyle>
  {
    public:

      //
      // Constructors 
      //

      VectorView(const VectorView<T,CStyle>& rhs) : 
	VectorView<T,CStyle>(rhs) {}

      VectorView(const VectorView<T,FortranStyle>& rhs) : 
	VectorView<T,CStyle>(rhs) {}

      VectorView(T* inv, size_t insize, int instep, ConjItType inct
	  PARAMFIRSTLAST(T) ) :
	VectorView<T,CStyle>(inv,insize,instep,inct
	    FIRSTLAST1(_first,_last) ) {}

      ~VectorView() {}


      //
      // Op =
      //

      inline const VectorView<T,FortranStyle>& operator=(
	  const VectorView<T,FortranStyle>& v2) const
      { VectorView<T,CStyle>::operator=(v2); return *this; }

      inline const VectorView<T,FortranStyle>& operator=(
	  const GenVector<T>& v2) const
      { VectorView<T,CStyle>::operator=(v2); return *this; }

      template <class T2> inline const VectorView<T,FortranStyle>& operator=(
	  const GenVector<T2>& v2) const
      { VectorView<T,CStyle>::operator=(v2); return *this; }

      inline const VectorView<T,FortranStyle>& operator=(
	  const VectorComposite<T>& vcomp) const
      { VectorView<T,CStyle>::operator=(vcomp); return *this; }


      //
      // Access Functions
      //

      inline RefType(T) operator[](size_t i) const 
      { 
	TMVAssert(i>0 && i<=size());
	return ref(i-1); 
      }
      inline RefType(T) operator()(size_t i) const 
      { 
	TMVAssert(i>0 && i<=size());
	return ref(i-1); 
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
      inline const VectorView<T,FortranStyle>& MakeBasis(size_t i) const
      {
	TMVAssert(i>0 && i<=size());
	VectorView<T,CStyle>::MakeBasis(i-1); return *this; 
      }
      inline const VectorView<T,FortranStyle>& Swap(size_t i1, size_t i2) const
      {
	TMVAssert(i1>0 && i1<=size());
	TMVAssert(i2>0 && i2<=size());
	if (i1 != i2) 
	  VectorView<T,CStyle>::Swap(i1-1,i2-1);
	return *this;
      }
      inline const VectorView<T,FortranStyle>& Permute(const size_t* p, 
	  size_t i1, size_t i2) const
      { 
	TMVAssert(i1 > 0);
	VectorView<T,CStyle>::Permute(p,i1-1,i2); 
	return *this; 
      }
      inline const VectorView<T,FortranStyle>& Permute(const size_t* p) const
      { VectorView<T,CStyle>::Permute(p); return *this; }
      inline const VectorView<T,FortranStyle>& ReversePermute(
	  const size_t* p, size_t i1, size_t i2) const
      { 
	TMVAssert(i1 > 0);
	VectorView<T,CStyle>::ReversePermute(p,i1-1,i2); 
	return *this; 
      }
      inline const VectorView<T,FortranStyle>& ReversePermute(
	  const size_t* p) const
      { VectorView<T,CStyle>::ReversePermute(p); return *this; }
      inline const VectorView<T,FortranStyle>& ReverseSelf() const
      { VectorView<T,CStyle>::ReverseSelf(); return *this; }
      inline const VectorView<T,FortranStyle>& Sort(
	  size_t* P, ADType ad=ASCEND, COMPType comp=REAL_COMP) const
      { VectorView<T,CStyle>::Sort(P,ad,comp); return *this; }

      //
      // SubVector
      //

      inline bool OKSubVector(int i1, int i2, int istep) const
      {
	return ConstVectorView<T,FortranStyle>(*this).OKSubVector(i1,i2,istep);
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

      using VectorView<T,CStyle>::size;

    protected :
   
      using VectorView<T,CStyle>::ref;

  }; // FortranStyle VectorView

  template <class T, IndexStyle I> class Vector : 
    public GenVector<T>
  {
    public:

      //
      // Constructors
      //

#define NEW_SIZE(n) \
      GenVector<T>(NonConj), itsv(new T[(n)]), itssize(n) \
      DEFFIRSTLAST(itsv,itsv+n)

      explicit Vector(size_t n) : NEW_SIZE(n) 
      {
#ifdef TMVDEBUG
	SetAllTo(T(888));
#endif
      }

      Vector(size_t n,T val) : NEW_SIZE(n)
      { SetAllTo(val); }

      explicit Vector(const valarray<T>& vv) : NEW_SIZE(vv.size())
      { 
	T* p = itsv;
	for(size_t i=0;i<size();++i,++p) *p = vv[i]; 
      }

      Vector(size_t n, const T* vv) : NEW_SIZE(n)
      { memmove(itsv,vv,n*sizeof(T)); }

      explicit Vector(const vector<T>& vv) : NEW_SIZE(vv.size())
      {
	T* p = itsv;
	for(size_t i=0;i<size();++i,++p) *p = vv[i]; 
      }

      Vector(const Vector<T,I>& rhs) : NEW_SIZE(rhs.size())
      { 
#ifdef TMVFLDEBUG
	TMVAssert(rhs.cptr() >= rhs.first);
	TMVAssert(rhs.cptr()+itssize <= rhs.last);
#endif
	memmove(itsv,rhs.cptr(),itssize*sizeof(T)); 
      }

      template <IndexStyle I2> Vector(const Vector<T,I2>& rhs) : 
	NEW_SIZE(rhs.size())
      { 
#ifdef TMVFLDEBUG
	TMVAssert(rhs.cptr() >= rhs.first);
	TMVAssert(rhs.cptr()+itssize <= rhs.last);
#endif
	memmove(itsv,rhs.cptr(),itssize*sizeof(T)); 
      }

      Vector(const GenVector<T>& rhs) : NEW_SIZE(rhs.size())
      { 
	if (rhs.step() == 1 && !rhs.isconj()) 
	  memmove(itsv,rhs.cptr(),itssize*sizeof(T)); 
	else Copy(rhs,View());
      }

      template <class T2> Vector(const GenVector<T2>& rhs) : 
	NEW_SIZE(rhs.size())
      { Copy(rhs,View()); }

      Vector(const VectorComposite<T>& vcomp) : NEW_SIZE(vcomp.size())
      { vcomp.AssignTo(View()); }


#undef NEW_SIZE

      ~Vector() { if (itsv) delete[] itsv; }


      //
      // Op =
      //

      inline Vector<T,I>& operator=(Vector<T,I>& v2)
      { 
	TMVAssert(v2.size() == size());
	if (&v2 != this) memmove(itsv,v2.cptr(),itssize*sizeof(T)); 
	return *this; 
      }

      template <IndexStyle I2> inline Vector<T,I>& operator=(Vector<T,I2>& v2)
      { 
	TMVAssert(v2.size() == size());
	memmove(itsv,v2.cptr(),itssize*sizeof(T)); 
	return *this; 
      }

      inline Vector<T,I>& operator=(const GenVector<T>& v2) 
      { View() = v2; return *this; }

      template <class T2> inline Vector<T,I>& operator=(
	  const GenVector<T2>& v2) 
      { View() = v2; return *this; }

      inline Vector<T,I>& operator=(const VectorComposite<T>& vcomp)
      { vcomp.AssignTo(View()); return *this; }


      //
      // Access Functions
      //

      typedef VIt<T,Unit,NonConj> iterator;
      typedef CVIt<T,Unit,NonConj> const_iterator;
      typedef VIt<T,Step,NonConj> reverse_iterator;
      typedef CVIt<T,Step,NonConj> const_reverse_iterator;
      typedef RefType(T) reference;

      inline CVIt<T,Unit,NonConj> begin() const 
      { return CVIt<T,Unit,NonConj>(cptr(),1); }
      inline CVIt<T,Unit,NonConj> end() const 
      { return begin()+size(); }
      
      inline CVIt<T,Step,NonConj> rbegin() const 
      { return const_reverse_iterator(cptr()+(size()-1),-1); }
      inline CVIt<T,Step,NonConj> rend() const 
      { return rbegin()+size(); }

      inline T operator[](size_t i) const 
      { 
	if (I == CStyle) { TMVAssert(i<size()); return cref(i); } 
	else { TMVAssert(i>0 && i<=size()); return cref(i-1); }
      }
      inline T operator()(size_t i) const 
      { 
	if (I == CStyle) { TMVAssert(i<size()); return cref(i); } 
	else { TMVAssert(i>0 && i<=size()); return cref(i-1); }
      }

      inline VIt<T,Unit,NonConj> begin() 
      { return iterator(ptr(),1 FIRSTLAST ); }
      inline VIt<T,Unit,NonConj> end() 
      { return begin() + size(); }

      inline VIt<T,Step,NonConj> rbegin() 
      { return reverse_iterator(ptr()+size()-1,-1 FIRSTLAST ); }
      inline VIt<T,Step,NonConj> rend() 
      { return rbegin()+size(); }

      inline T& operator[](size_t i) 
      { 
	if (I == CStyle) { TMVAssert(i<size()); return ref(i); } 
	else { TMVAssert(i>0 && i<=size()); return ref(i-1); }
      }
      inline T& operator()(size_t i) 
      { 
	if (I == CStyle) { TMVAssert(i<size()); return ref(i); } 
	else { TMVAssert(i>0 && i<=size()); return ref(i-1); }
      }

      //
      // Modifying Functions
      //

      inline Vector<T,I>& Zero() 
      { 
	if (IsComplex(T())) Flatten().SetAllTo(0);
	else SetAllTo(0); 
	return *this;
      }

      inline Vector<T,I>& Clip(RealType(T) thresh)
      { View().Clip(thresh); return *this; }

      inline Vector<T,I>& SetAllTo(T x)
      { 
	T* p = ptr();
	for(size_t i=size();i>0;--i,++p) *p = x;
	return *this;
      }

      inline Vector<T,I>& AddToAll(T x)
      {
	T* p = ptr();
	for(size_t i=size();i>0;--i,++p) *p += x;
	return *this;
      }

      inline Vector<T,I>& ConjugateSelf()
      { View().ConjugateSelf(); return *this; }

      inline Vector<T,I>& MakeBasis(size_t i)
      { 
	if (I == CStyle) { TMVAssert(i<size()); } 
	else { TMVAssert(i>0 && i<=size()); }
	const size_t ix = (I==CStyle ? i : i-1);
	Zero(); 
	ref(ix) = T(1); 
	return *this; 
      }

      inline Vector<T,I>& Swap(size_t i1, size_t i2)
      {
	if (I == CStyle) { TMVAssert(i1<size() && i2<size()); } 
	else { TMVAssert(i1>0 && i1<=size() && i2>0 && i2<=size()); }
	if (i1 != i2) 
	  if (I == CStyle) swap(ref(i1),ref(i2));
	  else swap(ref(i1-1),ref(i2-1));
	return *this;
      }

      inline Vector<T,I>& Permute(const size_t* p, size_t i1, size_t i2)
      {
	if (I==FortranStyle) { TMVAssert(i1>0); }
	TMVAssert(i1 <= i2 && i2 <= size());
	View().Permute(p,i1,i2); return *this; 
      }
      inline Vector<T,I>& Permute(const size_t* p) 
      { View().Permute(p); return *this; }
      inline Vector<T,I>& ReversePermute(const size_t* p, size_t i1, size_t i2)
      {
	if (I==FortranStyle) { TMVAssert(i1>0); }
	TMVAssert(i1 <= i2 && i2 <= size());
	View().ReversePermute(p,i1,i2); return *this; 
      }
      inline Vector<T,I>& ReversePermute(const size_t* p) 
      { View().ReversePermute(p,0,size()); return *this; }

      inline Vector<T,I>& ReverseSelf()
      { View().ReverseSelf(); return *this; }
      inline Vector<T,I>& Sort(size_t* P, ADType ad=ASCEND, 
	  COMPType comp=REAL_COMP) 
      { View().Sort(P,ad,comp); return *this; }

      //
      // SubVector
      //

      inline ConstVectorView<T,I> SubVector(int i1, int i2) const
      {
	TMVAssert(View().OKSubVector(i1,i2,1));
	const int i1x = (I==CStyle ? i1 : i1-1);
	return ConstVectorView<T,I>(cptr()+i1x,i2-i1x,1,NonConj);
      }

      inline VectorView<T,I> SubVector(int i1, int i2)
      {
	TMVAssert(View().OKSubVector(i1,i2,1));
	const int i1x = (I==CStyle ? i1 : i1-1);
	return VectorView<T,I>(ptr()+i1x,i2-i1x,1,NonConj FIRSTLAST );
      }

      inline ConstVectorView<T,I> SubVector(int i1, int i2, int istep) const
      {
	TMVAssert(View().OKSubVector(i1,i2,istep));
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int i2x = (I==CStyle ? i2 : i2-1+istep);
	return ConstVectorView<T,I>(cptr()+i1x,(i2x-i1x)/istep,istep,NonConj);
      }

      inline VectorView<T,I> SubVector(int i1, int i2, int istep)
      {
	TMVAssert(View().OKSubVector(i1,i2,istep));
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int i2x = (I==CStyle ? i2 : i2-1+istep);
	return VectorView<T,I>(ptr()+i1x,(i2x-i1x)/istep,istep,NonConj 
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

      // 
      // I/O
      //
      
      void Read(istream& fin)
      { View().Read(fin); }

      inline size_t size() const { return itssize; }
      inline const T* cptr() const { return itsv; }
      inline T* ptr() { return itsv; }
      inline int step() const { return 1; }
      inline ConjItType ct() const { return NonConj; }
      inline bool isconj() const { return false; }

    private:

      T*const itsv;
      const size_t itssize;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    private:
#endif

      inline T cref(size_t i) const
      { 
	TMVAssert(i<size()); return *(itsv + i);
      }

      inline T& ref(size_t i)
      { 
	TMVAssert(i<size());
	T*const vi = itsv + i;
#ifdef TMVFLDEBUG
	TMVAssert(vi >= first);
	TMVAssert(vi < last);
#endif
	return *vi; 
      }

  }; // Vector

  //
  // Special Constructors
  //

  template <class T, IndexStyle I> inline Vector<T,I> BasisVector(
      size_t n, size_t i)
  { 
    if (I == CStyle) { TMVAssert(i<n); } 
    else { TMVAssert(i>0 && i<=n); }
    Vector<T,I> temp(n,T(0));
    temp(i) = T(1); 
    return temp; 
  }

  template <class T> inline Vector<T,CStyle> BasisVector(
      size_t n, size_t i)
  { 
    TMVAssert(i<n);
    Vector<T,CStyle> temp(n,T(0)); 
    temp(i) = T(1); 
    return temp; 
  }

  template <IndexStyle I, class T> inline VectorView<T,I> VectorViewOf(
      T* v, size_t size, int step=1)
  { return VectorView<T,I>(v,size,step,NonConj FIRSTLAST1(v,v+size*step)); }

  template <IndexStyle I, class T> inline ConstVectorView<T,I> VectorViewOf(
      const T* v, size_t size, int step=1)
  { return ConstVectorView<T,I>(v,size,step,NonConj); }

  template <class T> inline VectorView<T,CStyle> VectorViewOf(
      T* v, size_t size, int step=1)
  {
    return VectorView<T,CStyle>(v,size,step,NonConj 
	FIRSTLAST1(v,v+size*step)); 
  }

  template <class T> inline ConstVectorView<T,CStyle> VectorViewOf(
      const T* v, size_t size, int step=1)
  { return ConstVectorView<T,CStyle>(v,size,step,NonConj); }

  //
  // Copy Vectors
  //

  inline bool ShouldReverse(const int step1, const int step2)
  {
    return (  (step2 < 0 && (step1 != 1 || step2 == -1)) ||
	(step1 == -1 && step2 != 1) );
  }

  // Weirdly, if this and the below DoCopyDiffType have the same name, then
  // DoCopy can't resolve the overload.  So use different names to be sure.
  template <bool c1, class T> void DoCopySameType(
      const GenVector<T>& v1, const VectorView<T>& v2);

  // Also need to repeat this for same-typed vectors
  template <class T> inline void DoCopy(
      const GenVector<T>& v1, const VectorView<T>& v2)
  {
    if (v1.step() == 0) v2.SetAllTo(v1[0]);
    else if (v1.isconj()) DoCopySameType<true>(v1,v2);
    else DoCopySameType<false>(v1,v2);
  }

  template <bool c1, class T, class T1> inline void DoCopyDiffType(
      const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(v1.step()!=0);
    TMVAssert(v2.step()!=0);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(!v2.SameAs(v1));
    TMVAssert(c1 == v1.isconj());

    const T1* v1ptr = v1.cptr();
    T* v2ptr = v2.ptr();
    const int step1 = v1.step();
    const int step2 = v2.step();

    if (step1 == 1 && step2 == 1)
      for(size_t i=v2.size();i>0;--i,++v1ptr,++v2ptr)
	*v2ptr = (c1 ? CONJ(*v1ptr) : (*v1ptr));
    else
      for(size_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2)
	*v2ptr = (c1 ? CONJ(*v1ptr) : (*v1ptr));
  }

  template <class T, class T1> inline void DoCopy(
      const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    if (v1.step() == 0) v2.SetAllTo(v1[0]);
    else if (v1.isconj()) DoCopyDiffType<true>(v1,v2);
    else DoCopyDiffType<false>(v1,v2);
  }

  template <class T> inline void DoCopy(
      const GenVector<complex<T> >& v1, const VectorView<T>& v2)
  { TMVAssert(false); }
  template <class T, class T1> inline void Copy(
      const GenVector<T1>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.step()!=0 || v1.step() == 0);

    if (v1.size() > 0 && !v2.SameAs(v1))
      if (ShouldReverse(v1.step(),v2.step()))
	Copy(v1.Reverse(),v2.Reverse());
      else if (v2.isconj()) DoCopy(v1.Conjugate(),v2.Conjugate());
      else DoCopy(v1,v2);
  }


  //
  // Swap Vectors
  //

  template <class T> void Swap(const VectorView<T>& v1, const VectorView<T>& v2);
  template <class T, IndexStyle I> inline void Swap(const VectorView<T>& v1, 
      Vector<T,I>& v2)
  { Swap(v1.View(),v2.View()); }
  template <class T, IndexStyle I> inline void Swap(Vector<T,I>& v1,
      const VectorView<T>& v2) 
  { Swap(v1.View(),v2.View()); }
  template <class T, IndexStyle I1, IndexStyle I2> void Swap(Vector<T,I1>& v1, 
      Vector<T,I2>& v2)
  { Swap(v1.View(),v2.View()); }

  //
  // Functions of Vectors
  //

  template <class T> inline RealType(T) Norm(const GenVector<T>& v)
  { return v.Norm(); }

  template <class T> inline RealType(T) Norm1(const GenVector<T>& v)
  { return v.Norm1(); }

  template <class T> RealType(T) NormSq(const GenVector<T>& v)
  { return v.NormSq(); }

  template <class T> RealType(T) Norm2(const GenVector<T>& v)
  { return v.Norm2(); }

  template <class T> inline RealType(T) NormInf(const GenVector<T>& v)
  { return v.NormInf(); }

  template <class T> inline T SumElements(const GenVector<T>& v)
  { return v.SumElements(); }

  template <class T> inline RealType(T) SumAbsElements(const GenVector<T>& v)
  { return v.SumAbsElements(); }

  template <class T> inline T MinElement(const GenVector<T>& v, 
      size_t* iminout=0)
  { return v.MinElement(iminout); }

  template <class T> inline T MaxElement(const GenVector<T>& v, 
      size_t* imaxout=0)
  { return v.MaxElement(imaxout); }

  template <class T> inline RealType(T) MinAbsElement(const GenVector<T>& v, 
      size_t* iminout=0)
  { return v.MinAbsElement(iminout); }

  template <class T> inline RealType(T) MaxAbsElement(const GenVector<T>& v, 
      size_t* imaxout=0)
  { return v.MaxAbsElement(imaxout); }

  template <class T> inline ConstVectorView<T> Conjugate(const GenVector<T>& v)
  { return v.Conjugate(); }

  template <class T, IndexStyle I> inline ConstVectorView<T,I> Conjugate(
      const ConstVectorView<T,I>& v)
  { return v.Conjugate(); }

  template <class T, IndexStyle I> inline ConstVectorView<T,I> Conjugate(
      const Vector<T,I>& v)
  { return v.Conjugate(); }

  template <class T, IndexStyle I> inline VectorView<T,I> Conjugate(
      const VectorView<T,I>& v)
  { return v.Conjugate(); }

  template <class T, IndexStyle I> inline VectorView<T,I> Conjugate(Vector<T,I>& v)
  { return v.Conjugate(); }


  //
  // Vector ==, != Vector
  //

  template <class T1, class T2> bool operator==(
      const GenVector<T1>& v1, const GenVector<T2>& v2);

  template <class T1, class T2> inline bool operator!=(
      const GenVector<T1>& v1, const GenVector<T2>& v2)
  { return !(v1 == v2); }

  //
  // I/O
  //

  template <class T> inline ostream& operator<<(
      ostream& fout, const GenVector<T>& v)
  { v.Write(fout); return fout;}

  template <class T> istream& operator>>(
      istream& fin, const VectorView<T>& v);

  template <class T, IndexStyle I> inline istream& operator>>(istream& fin,
      Vector<T,I>& v)
  { return fin >> v.View(); }

  template <class T, IndexStyle I> istream& operator>>(istream& fin, 
      Vector<T,I>*& v);


  inline string Text(ADType ad)
  { return ad == ASCEND ? "Ascend" : "Descend"; }
  inline string Text(COMPType comp)
  { 
    return comp == REAL_COMP ? "Real" : comp == ABS_COMP ? "Abs" :
      comp == IMAG_COMP ? "Imag" : "Arg"; 
  }
  inline string Text(IndexStyle I)
  { return I == CStyle ? "CStyle" : "FortranStyle"; }

  template <class T, IndexStyle I> inline string Type(const Vector<T,I>& v)
  { return string("Vector<")+Type(T())+","+Text(I)+">"; }
  template <class T> inline string Type(const GenVector<T>& v)
  { return string("GenVector<")+Type(T())+","+Text(v.ct())+">"; }
  template <class T, IndexStyle I> inline string Type(const ConstVectorView<T>& v)
  { return string("ConstVectorView<")+Type(T())+","+Text(I)+","+Text(v.ct())+">"; }
  template <class T, IndexStyle I> inline string Type(const VectorView<T>& v)
  { 
    return string("VectorView<")+Type(T())+","+Text(I)+","+Text(v.ct())+">"; 
  }

} // namespace tmv

#endif
