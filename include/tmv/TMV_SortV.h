

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

    // Just use the standard library sort routine.
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
        ptrdiff_t itsi;

    public :

        void set(T val, ptrdiff_t i) { itsvalue = val; itsi = i; }
        ptrdiff_t getI() const { return itsi; }
        T getVal() const { return itsvalue; }
        bool operator<(const VTIndex& rhs) const
        { return itsvalue < rhs.itsvalue; }
        operator ptrdiff_t() const { return itsi; }

    };

    // Input P is in index form, meaning that permutation is
    // v[i] -> v[P[i]]
    // Ouput P is in swap form, meaning that permutation is series of
    // v.swap(i,P[i])
    template <class VI>
    static void ConvertIndexToPermute(ptrdiff_t n, const VI& newindex, ptrdiff_t* P)
    {
        // newindex[i]=j means value at original j location needs to go to i.
        std::vector<ptrdiff_t> currindex(n);
        std::vector<ptrdiff_t> origindex(n);
        for(ptrdiff_t i=0;i<n;++i) {
            currindex[i] = i;
            origindex[i] = i;
        } 
        // currindex[i]=j means value at original i location is currently at j.
        // origindex[j]=i means value at original i location is currently at j.
        for(ptrdiff_t i=0;i<n;++i) {
            ptrdiff_t ip = currindex[newindex[i]];
            P[i] = ip;
            if (i != ip) {
                ptrdiff_t origi = origindex[i];
                ptrdiff_t origip = origindex[ip];
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
    void InstSort(VectorView<T> v, ptrdiff_t* P, ADType ad, CompType comp);

    template <int algo, class V>
    struct Sort_Helper;

    // algo 11: Normal case.  Call DoSort with the right Compare object.
    template <class V>
    struct Sort_Helper<11,V>
    {
        static void call(V& v, ptrdiff_t* P, ADType ad, CompType comp)
        {
            typedef typename V::float_type FT;
            const ptrdiff_t n = v.size();
            std::vector<VTIndex<FT> > newindex(n);
            if (ad == Ascend) {
                if (comp == RealComp) 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(TMV_REAL(v.cref(i)),i);
                else if (comp == AbsComp) 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(TMV_ABS(v.cref(i)),i);
                else if (comp == ImagComp) 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(TMV_IMAG(v.cref(i)),i);
                else 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(TMV_ARG(v.cref(i)),i);
            } else {
                if (comp == RealComp) 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(-TMV_REAL(v.cref(i)),i);
                else if (comp == AbsComp) 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(-TMV_ABS(v.cref(i)),i);
                else if (comp == ImagComp) 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(-TMV_IMAG(v.cref(i)),i);
                else 
                    for(ptrdiff_t i=0;i<n;++i) 
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
        static void call(V& v, ptrdiff_t* P, ADType ad, CompType comp)
        {
            typedef typename V::float_type FT;
            const ptrdiff_t n = v.size();
            std::vector<VTIndex<FT> > newindex(n);
            if (ad == Ascend) {
                if (comp == RealComp) 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(TMV_REAL(v.cref(i)),i);
                else 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(TMV_ABS(v.cref(i)),i);
            } else {
                if (comp == RealComp) 
                    for(ptrdiff_t i=0;i<n;++i) 
                        newindex[i].set(-TMV_REAL(v.cref(i)),i);
                else 
                    for(ptrdiff_t i=0;i<n;++i) 
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
        static TMV_INLINE void call(V& v, ptrdiff_t* P, ADType ad, CompType comp)
        { InstSort(v.xView(),P,ad,comp); }
        static TMV_INLINE void call(V& v, ADType ad, CompType comp)
        { InstSort(v.xView(),ad,comp); }
    };

    // algo 97: Conjugate
    template <class V>
    struct Sort_Helper<97,V>
    {
        static TMV_INLINE void call(V& v, ptrdiff_t* P, ADType ad, CompType comp)
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
        static TMV_INLINE void call(V& v, ptrdiff_t* P, ADType ad, CompType comp)
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
        static TMV_INLINE void call(V& v, ptrdiff_t* P, ADType ad, CompType comp)
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
        BaseVector_Mutable<V>& v, ptrdiff_t* P, ADType ad, CompType comp)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        if (P == 0) Sort_Helper<-3,Vv>::call(vv,ad,comp);
        else Sort_Helper<-3,Vv>::call(vv,P,ad,comp);
    }

    template <class V>
    inline void Sort(
        BaseVector_Mutable<V>& v, ptrdiff_t* P, ADType ad, CompType comp)
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
