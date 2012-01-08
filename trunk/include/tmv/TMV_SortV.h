

#ifndef TMV_SortV_H
#define TMV_SortV_H

#include "TMV_BaseVector.h"
#include <vector>

namespace tmv {

    //
    // Sort Vector
    //

#if TMV_OPT >= 1
#define COMPILE_TIME_AD_COMP
#endif

    template <class T, ADType ad, CompType comp>
    struct Compare
    {
        bool operator()(const T& x, const T& y) const
        {
            return Maybe<ad==Ascend>::less(
                Component<comp,T>::f(x),Component<comp,T>::f(y)); 
        }
    };
    template <class T>
    struct Compare2 // real T
    {
        ADType ad;
        CompType comp;
        Compare2(ADType _ad, CompType _comp) : ad(_ad), comp(_comp) {}
        bool operator()(const T& x, const T& y) const
        {
            return 
                ad==Ascend ?  (
                    (comp == RealComp) ?  TMV_REAL(x) < TMV_REAL(y) :
                    TMV_ABS(x) < TMV_ABS(y) ) :
                (
                    (comp == RealComp) ?  TMV_REAL(x) > TMV_REAL(y) :
                    TMV_ABS(x) > TMV_ABS(y) );
        }
    };
    template <class T>
    struct Compare2<std::complex<T> >
    {
        typedef std::complex<T> CT;
        ADType ad;
        CompType comp;
        Compare2(ADType _ad, CompType _comp) : ad(_ad), comp(_comp) {}
        bool operator()(const CT& x, const CT& y) const
        {
            return 
                ad==Ascend ?  (
                    (comp == RealComp) ?  TMV_REAL(x) < TMV_REAL(y) :
                    (comp == AbsComp) ? TMV_ABS(x) < TMV_ABS(y) :
                    (comp == ImagComp) ?  TMV_IMAG(x) < TMV_IMAG(y) :
                    TMV_ARG(x) < TMV_ARG(y) ) :
                (
                    (comp == RealComp) ?  TMV_REAL(x) > TMV_REAL(y) :
                    (comp == AbsComp) ? TMV_ABS(x) > TMV_ABS(y) :
                    (comp == ImagComp) ?  TMV_IMAG(x) > TMV_IMAG(y) :
                    TMV_ARG(x) > TMV_ARG(y) );
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
    inline void DoSort(const IT& begin, const IT& end)
    { DoSort(begin,end,std::less<iterator_traits<IT>::value_type>()); }
#else
    // Otherwise, just use the standard library sort routine.
    template <class IT, class COMP>
    inline void DoSort(IT begin, IT end, const COMP& comp)
    { std::sort(begin,end,comp); }

    // Overload the case of RealComp with real values, since then
    // we can call the sort without a comp object
    template <class IT, class T>
    inline void DoSort(
        IT begin, IT end, const Compare<T,Ascend,RealComp>& comp)
    { std::sort(begin,end); }
    template <class IT, class T>
    inline void DoSort(
        IT begin, IT end, const Compare<std::complex<T>,Ascend,RealComp>& comp)
    { std::sort(begin,end,comp); }

    template <class IT>
    inline void DoSort(const IT& begin, const IT& end)
    { std::sort(begin,end); }
#endif

    // A helper class to keep track of indexes in the sort.
    // We sort the "value" according to which component of the actual
    // data we want (REAL, ABS, IMAG, ARG).
    // And if the sorting uses Descend, then we store the negative 
    // of the value instead, which makes the comparison just the normal <.
    template <class T>
    class VTIndex
    {
    private :

        T itsvalue;
        int itsi;

    public :

        void set(T val, int i) { itsvalue = val; itsi = i; }
        int getI() const { return itsi; }
        T getVal() const { return itsvalue; }
        bool operator<(const VTIndex& rhs) const
        { return itsvalue < rhs.itsvalue; }
        operator int() const { return itsi; }

    };

    // Input P is in index form, meaning that permutation is
    // v[i] -> v[P[i]]
    // Ouput P is in swap form, meaning that permutation is series of
    // v.swap(i,P[i])
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

    // Defined in TMV_Vector.cpp
    template <class T>
    void InstSort(VectorView<T> v, ADType ad, CompType comp);
    template <class T>
    void InstSort(VectorView<T> v, int* P, ADType ad, CompType comp);

    template <int algo, class V>
    struct Sort_Helper;

    // algo 11: Normal case.  Call DoSort with the right Compare object.
    template <class V>
    struct Sort_Helper<11,V>
    {
        static void call(V& v, int* P, ADType ad, CompType comp)
        {
            typedef typename V::value_type T;
            typedef typename V::float_type FT;
            const int n = v.size();
            std::vector<VTIndex<FT> > newindex(n);
            if (ad == Ascend) {
                if (comp == RealComp) 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(TMV_REAL(v.cref(i)),i);
                else if (comp == AbsComp) 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(TMV_ABS(v.cref(i)),i);
                else if (comp == ImagComp) 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(TMV_IMAG(v.cref(i)),i);
                else 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(TMV_ARG(v.cref(i)),i);
            } else {
                if (comp == RealComp) 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(-TMV_REAL(v.cref(i)),i);
                else if (comp == AbsComp) 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(-TMV_ABS(v.cref(i)),i);
                else if (comp == ImagComp) 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(-TMV_IMAG(v.cref(i)),i);
                else 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(-TMV_ARG(v.cref(i)),i);
            }
            DoSort(newindex.begin(),newindex.end());
            ConvertIndexToPermute(n,newindex,P);
            v.permute(P);
        }
        static void call(V& v, ADType ad, CompType comp)
        {
            typedef typename V::value_type T;
#ifdef COMPILE_TIME_AD_COMP
            if (ad == Ascend) {
                if (comp == RealComp) DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Ascend,RealComp>());
                else if (comp == AbsComp) DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Ascend,AbsComp>());
                else if (comp == ImagComp) DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Ascend,ImagComp>());
                else DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Ascend,ArgComp>());
            } else {
                if (comp == RealComp) DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Descend,RealComp>());
                else if (comp == AbsComp) DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Descend,AbsComp>());
                else if (comp == ImagComp) DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Descend,ImagComp>());
                else DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Descend,ArgComp>());
            }
#else
            DoSort(v.ptr(),v.ptr()+v.size(),Compare2<T>(ad,comp)); 
#endif
        }
    };

    // algo 12: v is real, so comp can't be Imag or Arg.
    template <class V>
    struct Sort_Helper<12,V>
    {
        static void call(V& v, int* P, ADType ad, CompType comp)
        {
            typedef typename V::value_type T;
            typedef typename V::float_type FT;
            const int n = v.size();
            std::vector<VTIndex<FT> > newindex(n);
            if (ad == Ascend) {
                if (comp == RealComp) 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(TMV_REAL(v.cref(i)),i);
                else 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(TMV_ABS(v.cref(i)),i);
            } else {
                if (comp == RealComp) 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(-TMV_REAL(v.cref(i)),i);
                else 
                    for(int i=0;i<n;++i) 
                        newindex[i].set(-TMV_ABS(v.cref(i)),i);
            }
            DoSort(newindex.begin(),newindex.end());
            ConvertIndexToPermute(n,newindex,P);
            v.permute(P);
        }
        static void call(V& v, ADType ad, CompType comp)
        {
            typedef typename V::value_type T;
#ifdef COMPILE_TIME_AD_COMP
            if (ad == Ascend) {
                if (comp == RealComp) DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Ascend,RealComp>());
                else DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Ascend,AbsComp>());
            } else {
                if (comp == RealComp) DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Descend,RealComp>());
                else DoSort(
                    v.ptr(),v.ptr()+v.size(),Compare<T,Descend,AbsComp>());
            }
#else
            DoSort(v.ptr(),v.ptr()+v.size(),Compare2<T>(ad,comp)); 
#endif
        }
    };

    // algo 90: Call inst
    template <class V>
    struct Sort_Helper<90,V>
    {
        static TMV_INLINE void call(V& v, int* P, ADType ad, CompType comp)
        { InstSort(v.xView(),P,ad,comp); }
        static TMV_INLINE void call(V& v, ADType ad, CompType comp)
        { InstSort(v.xView(),ad,comp); }
    };

    // algo 97: Conjugate
    template <class V>
    struct Sort_Helper<97,V>
    {
        static TMV_INLINE void call(V& v, int* P, ADType ad, CompType comp)
        {
            // Swap ad if necessary 
            if (comp==ImagComp || comp==ArgComp) {
                if (ad == Ascend) ad = Descend;
                else ad = Ascend;
            }
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            Sort_Helper<-2,Vc>::call(vc,P,ad,comp); 
        }
        static TMV_INLINE void call(V& v, ADType ad, CompType comp)
        {
            // Swap ad if necessary 
            if (comp==ImagComp || comp==ArgComp) {
                if (ad == Ascend) ad = Descend;
                else ad = Ascend;
            }
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            Sort_Helper<-2,Vc>::call(vc,ad,comp); 
        }
    };

    // algo -3: Determine which algorithm to use
    // TODO: For complex ABS and ARG, the VTIndex way might be faster,
    // than the direct comparison.  I should do some timings to figure
    // out if that's true and if so where the crossover point is
    // for the vector's size.
    template <class V>
    struct Sort_Helper<-3,V>
    {
        static TMV_INLINE void call(V& v, int* P, ADType ad, CompType comp)
        {
            const int algo = V::isreal ? 12 : 11;
            Sort_Helper<algo,V>::call(v,P,ad,comp);
        }
        static TMV_INLINE void call(V& v, ADType ad, CompType comp)
        {
            const int algo = V::isreal ? 12 : 11;
            Sort_Helper<algo,V>::call(v,ad,comp);
        }
    };

    // algo -2: Check for inst
    template <class V>
    struct Sort_Helper<-2,V>
    {
        static TMV_INLINE void call(V& v, int* P, ADType ad, CompType comp)
        {
            typedef typename V::value_type T;
            const bool inst = 
                (V::_size == Unknown || V::_size > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::_conj ? 97 :
                inst ? 90 :
                -3;
            Sort_Helper<algo,V>::call(v,P,ad,comp);
        }
        static TMV_INLINE void call(V& v, ADType ad, CompType comp)
        {
            typedef typename V::value_type T;
            const bool inst = 
                (V::_size == Unknown || V::_size > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::_conj ? 97 :
                inst ? 90 :
                -3;
            Sort_Helper<algo,V>::call(v,ad,comp);
        }
    };

    template <class V>
    inline void InlineSort(
        BaseVector_Mutable<V>& v, int* P, ADType ad, CompType comp)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        if (P == 0) Sort_Helper<-3,Vv>::call(vv,ad,comp);
        else Sort_Helper<-3,Vv>::call(vv,P,ad,comp);
    }

    template <class V>
    inline void Sort(
        BaseVector_Mutable<V>& v, int* P, ADType ad, CompType comp)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        if (P == 0) Sort_Helper<-2,Vv>::call(vv,ad,comp);
        else Sort_Helper<-2,Vv>::call(vv,P,ad,comp);
    }

    template <class V>
    inline void InlineSort(BaseVector_Mutable<V>& v, ADType ad, CompType comp)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        Sort_Helper<-3,Vv>::call(vv,ad,comp);
    }

    template <class V>
    inline void Sort(BaseVector_Mutable<V>& v, ADType ad, CompType comp)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        Sort_Helper<-2,Vv>::call(vv,ad,comp);
    }

} // namespace tmv

#endif
