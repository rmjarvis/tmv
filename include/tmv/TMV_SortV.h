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


#ifndef TMV_SortV_H
#define TMV_SortV_H

#include <algorithm>
#include <vector>

namespace tmv {

  //
  // Sort Vector
  //

  template <ADType ad, COMPType comp, class T> struct Compare
  {
    typedef RealType(T) RT;
    inline bool operator()(const T& x, const T& y) const
    {
      return Maybe<ad==ASCEND>::less(
        Component<comp,T>::f(x),Component<comp,T>::f(y)); 
    }
  };

#ifdef TMV_NO_STL
  // We allow for a possibility of the std::sort function not working
  // correctly on someone's system.  This and a couple of other stl
  // issues may be avoided by compiling with -DNOSTL.
  // So, in this case, we use a median-of-three quicksort.
  // This is not quite as good as the STL sort function, since they 
  // use introsort - a strategy that mixes quick sort with other algorithms
  // under particular situations.  So if the std::sort() works on your
  // system, it is recommended.
  template <class IT, class COMP>
  inline void Sort3(IT x1, IT x2, IT x3, const COMP& comp)
  {
    if (comp(*x3,*x1)) {
      if (comp(*x1,*x2)) { // x3 < x1 < x2
        iterator_traits<IT>::value_type temp = *x1;
        *x1 = *x3;
        *x3 = *x2;
        *x2 = temp;
      } else if (comp(*x2,*x3)) { // x2 < x3 < x1
        iterator_traits<IT>::value_type temp = *x1;
        *x1 = *x2;
        *x2 = *x3;
        *x3 = temp;
      } else { // x3 <= x2 <= x1 and x3 < x1
        TMV_SWAP(*x1,*x3);
      }
    } else {
      if (comp(*x2,*x1)) { // x2 < x1 <= x3
        TMV_SWAP(*x1,*x2);
      } else if (comp(*x3,*x2)) { // x1 <= x3 < x2
        TMV_SWAP(*x2,*x3);
      } else { // x1 <= x2 <= x3
        // nothing to do
      }
    }
  }

  template <class IT, class COMP>
  inline void DoSort(IT begin, IT end, const COMP& comp)
  {
    iterator_traits<IT>::difference_type n = end-begin;
    TMVAssert(N >= 0);
    if (N <= 3) {
      if (N == 3) { // 3 elements
        Sort3(begin,begin+1,begin+2,comp);
      } else if (N == 2) { // 2 elements
        if (comp(*(begin+1),*begin)) TMV_SWAP(*begin,*(begin+1));
      } // else 0 or 1 element
      return;
    } else {
      IT mid = begin + N/2;
      Sort3(begin,mid,end-1,comp);
      iterator_traits<IT>::value_type pivot = *mid;
      TMV_SWAP(*mid,*(end-2));
      IT left = begin+1;
      IT right = end-3;
      while (left < right) {
        while (!comp(*right,pivot) && left < right) --right;
        while (!comp(pivot,*left) && left < right) ++left;
        if (left < right) TMV_SWAP(*left,*right);
      }
      TMVAssert(left == right);
      if (comp(*left,pivot)) ++left;
      TMV_SWAP(*left,*(end-2));
      DoSort(begin,left,comp);
      DoSort(left+1,end,comp);
    }
  }

  template <class IT>
  inline void DoSort(const IT& begin, const IT& end)
  { DoSort(begin,end,std::less<iterator_traits<IT>::value_type>()); }
#else
  // Otherwise, just use the standard library sort routine.
  template <class IT, class COMP>
  inline void DoSort(IT begin, IT end, const COMP& comp)
  { std::sort(begin,end,comp); }

  // Overload the case of REAL_COMP with real values, since then
  // we can call the sort without a comp object
  template <class IT, class T>
  inline void DoSort(IT begin, IT end,
      const Compare<ASCEND,REAL_COMP,T>& comp)
  { std::sort(begin,end); }
  template <class IT, class T>
  inline void DoSort(IT begin, IT end,
      const Compare<ASCEND,REAL_COMP,std::complex<T> >& comp)
  { std::sort(begin,end,comp); }

  template <class IT>
  inline void DoSort(const IT& begin, const IT& end)
  { std::sort(begin,end); }
#endif

  //
  // Sort without returning permutation
  //

  // TODO: I should probably change the syntax to have ad and comp be 
  // template parameters.
  template <class V>
  inline void InlineSort(BaseVector_Mutable<V>& v, ADType ad, COMPType comp)
  {
    typedef typename V::value_type T;
    if (ad == ASCEND) {
      if (comp == REAL_COMP) 
        return DoSort(v.ptr(),v.ptr()+v.size(),Compare<ASCEND,REAL_COMP,T>());
      else if (comp == ABS_COMP) 
        return DoSort(v.ptr(),v.ptr()+v.size(),Compare<ASCEND,ABS_COMP,T>());
      else if (comp == IMAG_COMP) 
        return DoSort(v.ptr(),v.ptr()+v.size(),Compare<ASCEND,IMAG_COMP,T>());
      else 
        return DoSort(v.ptr(),v.ptr()+v.size(),Compare<ASCEND,ARG_COMP,T>());
    } else {
      if (comp == REAL_COMP) 
        return DoSort(v.ptr(),v.ptr()+v.size(),Compare<DESCEND,REAL_COMP,T>());
      else if (comp == ABS_COMP) 
        return DoSort(v.ptr(),v.ptr()+v.size(),Compare<DESCEND,ABS_COMP,T>());
      else if (comp == IMAG_COMP) 
        return DoSort(v.ptr(),v.ptr()+v.size(),Compare<DESCEND,IMAG_COMP,T>());
      else 
        return DoSort(v.ptr(),v.ptr()+v.size(),Compare<DESCEND,ARG_COMP,T>());
    }
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  void InstSort(VectorView<T> v, ADType ad, COMPType comp);

  template <bool vconj, bool inst, class V>
  struct CallSort // vconj = true
  {
    static inline void call(V& v, ADType ad, COMPType comp)
    {
      // Swap ad if necessary according to the conj status of the vector:
      if (V::vconj && (comp==IMAG_COMP || comp==ARG_COMP)) {
        if (ad == ASCEND) ad = DESCEND;
        else ad = ASCEND;
      }
      CallSort<false,inst,V>::call(v.Conjugate(),ad,comp); 
    }
  };
  template <class V>
  struct CallSort<false,false,V> // inst = false
  {
    static inline void call(V& v, ADType ad, COMPType comp)
    { InlineSort(v,ad,comp); }
  };
  template <class V>
  struct CallSort<false,true,V> // inst = true
  {
    static inline void call(V& v, ADType ad, COMPType comp)
    { InstSort(v.XView(),ad,comp); }
  };

  template <class V>
  inline void Sort(BaseVector_Mutable<V>& v, ADType ad, COMPType comp)
  {
    typedef typename V::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        V::vsize == UNKNOWN) };
    CallSort<V::vconj,inst,V>::call(v.vec(),ad,comp);
  }



  //
  // Sort with returning permutation
  //

  // A helper class to keep track of indexes in the sort.
  // We sort the "value" according to which component of the actual
  // data we want (REAL, ABS, IMAG, ARG).
  // And if the sorting uses DESCEND, then we store the negative 
  // of the value instead, which makes the comparison just the normal <.
  template <ADType ad, COMPType comp, class T> class VTIndex
  {
  private :

    typedef typename Traits<T>::real_type RT;

    RT itsvalue;
    int itsi;

  public :

    VTIndex() : itsvalue(RT(0)), itsi(0) {}

    VTIndex(T val, int i) : itsvalue(0), itsi(i)
    { itsvalue = Maybe<ad==DESCEND>::neg(Component<comp,T>::f(val)); }

    int GetI() const { return itsi; }
    RT GetVal() const { return itsvalue; }
    bool operator<(const VTIndex& rhs) const
    { return itsvalue < rhs.itsvalue; }
    operator int() const { return itsi; }

  };

  // Input P is in index form, meaning that permutation is
  // v[i] -> v[P[i]]
  // Ouput P is in swap form, meaning that permutation is series of
  // Swap(i,P[i])
  template <class VI> 
  inline void ConvertIndexToPermute(int n, const VI& newindex, int*const P)
  {
    // newindex[i]=j means value at original j location needs to go to i.
    std::vector<int> currindex(n);
    std::vector<int> origindex(n);
    for(int i=0;i<n;++i) {
      currindex[i] = i;
      origindex[i] = i;
    } 
    // currindex[i]=j means value at original i location is currently at j.
    // origindex[j]=i means value at original i location is currently at j.
    for(int i=0;i<n;++i) {
      int ip = currindex[newindex[i]];
      P[i] = ip;
      if (i != ip) { 
        int origi = origindex[i];
        int origip = origindex[ip];
        currindex[origi] = ip;
        currindex[origip] = i;
        origindex[i] = origip;
        origindex[ip] = origi;
      } 
    } 
  }

  template <ADType ad, COMPType comp, class V>
  inline void InlineSort1(
      BaseVector_Mutable<V>& v, int*const P)
  {
    typedef typename V::value_type T;
    const int n = v.size();
    std::vector<VTIndex<ad,comp,T> > newindex(n);
    for(int i=0;i<n;++i) newindex[i] = VTIndex<ad,comp,T>(v.ref(i),i);
    DoSort(newindex.begin(),newindex.end());
    ConvertIndexToPermute(n,newindex,P);
    v.Permute(P);
  }

  template <class V>
  inline void InlineSort(
      BaseVector_Mutable<V>& v, int*const P, ADType ad, COMPType comp)
  {
    if (ad == ASCEND) {
      if (comp == REAL_COMP) InlineSort1<ASCEND,REAL_COMP>(v,P);
      else if (comp == ABS_COMP) InlineSort1<ASCEND,ABS_COMP>(v,P);
      else if (comp == IMAG_COMP) InlineSort1<ASCEND,IMAG_COMP>(v,P);
      else InlineSort1<ASCEND,ARG_COMP>(v,P);
    } else {
      if (comp == REAL_COMP) InlineSort1<DESCEND,REAL_COMP>(v,P);
      else if (comp == ABS_COMP) InlineSort1<DESCEND,ABS_COMP>(v,P);
      else if (comp == IMAG_COMP) InlineSort1<DESCEND,IMAG_COMP>(v,P);
      else InlineSort1<DESCEND,ARG_COMP>(v,P);
    }
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  void InstSort(VectorView<T> v, int*const P, ADType ad, COMPType comp);

  template <bool vconj, bool inst, class V>
  struct CallSortP // vconj = true
  {
    static inline void call(V& v, int*const P, ADType ad, COMPType comp)
    {
      // Swap ad if necessary according to the conj status of the vector:
      if (V::vconj && (comp==IMAG_COMP || comp==ARG_COMP)) {
        if (ad == ASCEND) ad = DESCEND;
        else ad = ASCEND;
      }
      typename V::conjugate_type vc = v.Conjugate();
      CallSortP<false,inst,V>::call(vc,ad,comp); 
    }
  };
  template <class V>
  struct CallSortP<false,false,V> // inst = false
  {
    static inline void call(V& v, int*const P, ADType ad, COMPType comp)
    { InlineSort(v,P,ad,comp); }
  };
  template <class V>
  struct CallSortP<false,true,V> // inst = true
  {
    static inline void call(V& v, int*const P, ADType ad, COMPType comp)
    { InstSort(v.XView(),P,ad,comp); }
  };

  template <class V>
  inline void Sort(BaseVector_Mutable<V>& v,
      int*const P, ADType ad, COMPType comp)
  {
    if (P == 0) Sort(v,ad,comp);
    else {
      typedef typename V::value_type T;
      enum { inst = (
          Traits<T>::isinst &&
          V::vsize == UNKNOWN) };
      CallSortP<V::vconj,inst,V>::call(v.vec(),P,ad,comp);
    }
  }

} // namespace tmv

#endif
