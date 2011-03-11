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


#ifndef TMV_SortV_H
#define TMV_SortV_H

#include <algorithm>
#include <vector>

#include "TMV_BaseVector.h";

namespace tmv {

    //
    // Sort Vector
    //

    template <ADType ad, CompType comp, class T> 
    struct Compare
    {
        bool operator()(const T& x, const T& y) const
        {
            return Maybe<ad==Ascend>::less(
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
    static void Sort3(IT x1, IT x2, IT x3, const COMP& comp)
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
    static void DoSort(IT begin, IT end, const COMP& comp)
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
    static void DoSort(const IT& begin, const IT& end)
    { DoSort(begin,end,std::less<iterator_traits<IT>::value_type>()); }
#else
    // Otherwise, just use the standard library sort routine.
    template <class IT, class COMP>
    static void DoSort(IT begin, IT end, const COMP& comp)
    { std::sort(begin,end,comp); }

    // Overload the case of RealComp with real values, since then
    // we can call the sort without a comp object
    template <class IT, class T>
    static void DoSort(IT begin, IT end,
                       const Compare<Ascend,RealComp,T>& comp)
    { std::sort(begin,end); }
    template <class IT, class T>
    static void DoSort(IT begin, IT end,
                       const Compare<Ascend,RealComp,std::complex<T> >& comp)
    { std::sort(begin,end,comp); }

    template <class IT>
    static void DoSort(const IT& begin, const IT& end)
    { std::sort(begin,end); }
#endif

    // A helper class to keep track of indexes in the sort.
    // We sort the "value" according to which component of the actual
    // data we want (REAL, ABS, IMAG, ARG).
    // And if the sorting uses Descend, then we store the negative 
    // of the value instead, which makes the comparison just the normal <.
    template <ADType ad, CompType comp, class T> class VTIndex
    {
    private :

        typedef typename Traits<T>::real_type RT;
        typedef typename Traits<RT>::float_type FT;

        FT itsvalue;
        int itsi;

    public :

        VTIndex() : itsvalue(FT(0)), itsi(0) {}

        VTIndex(T val, int i) : itsvalue(0), itsi(i)
        { itsvalue = Maybe<ad==Descend>::neg(Component<comp,T>::f(val)); }

        int getI() const { return itsi; }
        FT getVal() const { return itsvalue; }
        bool operator<(const VTIndex& rhs) const
        { return itsvalue < rhs.itsvalue; }
        operator int() const { return itsi; }

    };

    // Input P is in index form, meaning that permutation is
    // v[i] -> v[P[i]]
    // Ouput P is in swap form, meaning that permutation is series of
    // Swap(i,P[i])
    template <class VI> 
    static void ConvertIndexToPermute(int n, const VI& newindex, int* P)
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

    template <ADType ad, CompType comp, class V>
    static void CallSort(V& v, int* P)
    {
        typedef typename V::value_type T;
        const int n = v.size();
        std::vector<VTIndex<ad,comp,T> > newindex(n);
        for(int i=0;i<n;++i) newindex[i] = VTIndex<ad,comp,T>(v.ref(i),i);
        DoSort(newindex.begin(),newindex.end());
        ConvertIndexToPermute(n,newindex,P);
        v.permute(P);
    }

    template <ADType ad, CompType comp, class V>
    static void CallSort(V& v)
    {
        typedef typename V::value_type T;
        DoSort(v.ptr(),v.ptr()+v.size(),Compare<ad,comp,T>()); 
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    void InstSort(VectorView<T> v, ADType ad, CompType comp);
    template <class T>
    void InstSort(VectorView<T> v, int* P, ADType ad, CompType comp);

    template <int algo, class V>
    struct Sort_Helper;

    // algo 1: Normal case.  Call DoSort with the right Compare object.
    template <class V>
    struct Sort_Helper<1,V>
    {
        static void call(V& v, int* P, ADType ad, CompType comp)
        {
            if (ad == Ascend) {
                if (comp == RealComp) CallSort<Ascend,RealComp>(v,P);
                else if (comp == AbsComp) CallSort<Ascend,AbsComp>(v,P);
                else if (comp == ImagComp) CallSort<Ascend,ImagComp>(v,P);
                else CallSort<Ascend,ArgComp>(v,P);
            } else {
                if (comp == RealComp) CallSort<Descend,RealComp>(v,P);
                else if (comp == AbsComp) CallSort<Descend,AbsComp>(v,P);
                else if (comp == ImagComp) CallSort<Descend,ImagComp>(v,P);
                else CallSort<Descend,ArgComp>(v,P);
            }
        }
        static void call(V& v, ADType ad, CompType comp)
        {
            if (ad == Ascend) {
                if (comp == RealComp) CallSort<Ascend,RealComp>(v);
                else if (comp == AbsComp) CallSort<Ascend,AbsComp>(v);
                else if (comp == ImagComp) CallSort<Ascend,ImagComp>(v);
                else CallSort<Ascend,ArgComp>(v);
            } else {
                if (comp == RealComp) CallSort<Descend,RealComp>(v);
                else if (comp == AbsComp) CallSort<Descend,AbsComp>(v);
                else if (comp == ImagComp) CallSort<Descend,ImagComp>(v);
                else CallSort<Descend,ArgComp>(v);
            }
        }
    };

    // algo 97: Conjugate
    template <class V>
    struct Sort_Helper<97,V>
    {
        static void call(V& v, int* P, ADType ad, CompType comp)
        {
            // Swap ad if necessary 
            if (comp==ImagComp || comp==ArgComp) {
                if (ad == Ascend) ad = Descend;
                else ad = Ascend;
            }
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            Sort_Helper<-1,Vc>::call(vc,P,ad,comp); 
        }
        static void call(V& v, ADType ad, CompType comp)
        {
            // Swap ad if necessary 
            if (comp==ImagComp || comp==ArgComp) {
                if (ad == Ascend) ad = Descend;
                else ad = Ascend;
            }
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            Sort_Helper<-1,Vc>::call(vc,ad,comp); 
        }
    };

    // algo 98: Call inst
    template <class V>
    struct Sort_Helper<98,V>
    {
        static void call(V& v, int* P, ADType ad, CompType comp)
        { InstSort(v.xView(),P,ad,comp); }
        static void call(V& v, ADType ad, CompType comp)
        { InstSort(v.xView(),ad,comp); }
    };

    // algo -1: Check for inst
    template <class V>
    struct Sort_Helper<-1,V>
    {
        static void call(V& v, int* P, ADType ad, CompType comp)
        {
            typedef typename V::value_type T;
            const bool inst = 
                (V::_size == UNKNOWN || V::_size > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::_conj ? 97 :
                inst ? 98 :
                1;
            Sort_Helper<algo,V>::call(v,P,ad,comp);
        }
        static void call(V& v, ADType ad, CompType comp)
        {
            typedef typename V::value_type T;
            const bool inst = 
                (V::_size == UNKNOWN || V::_size > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::_conj ? 97 :
                inst ? 98 :
                1;
            Sort_Helper<algo,V>::call(v,ad,comp);
        }
    };

    template <class V>
    static void InlineSort(
        BaseVector_Mutable<V>& v, int* P, ADType ad, CompType comp)
    {
        typedef typename V::cview_type Vv;
        Vv vv = v.cView();
        if (P == 0) Sort_Helper<1,Vv>::call(vv,ad,comp);
        else Sort_Helper<1,Vv>::call(vv,P,ad,comp);
    }

    template <class V>
    static void Sort(
        BaseVector_Mutable<V>& v, int* P, ADType ad, CompType comp)
    {
        typedef typename V::cview_type Vv;
        Vv vv = v.cView();
        if (P == 0) Sort_Helper<-1,Vv>::call(vv,ad,comp);
        else Sort_Helper<-1,Vv>::call(vv,P,ad,comp);
    }

    template <class V>
    static void InlineSort(BaseVector_Mutable<V>& v, ADType ad, CompType comp)
    {
        typedef typename V::cview_type Vv;
        Vv vv = v.cView();
        Sort_Helper<1,Vv>::call(vv,ad,comp);
    }

    template <class V>
    static void Sort(BaseVector_Mutable<V>& v, ADType ad, CompType comp)
    {
        typedef typename V::cview_type Vv;
        Vv vv = v.cView();
        Sort_Helper<-1,Vv>::call(vv,ad,comp);
    }

} // namespace tmv

#endif
