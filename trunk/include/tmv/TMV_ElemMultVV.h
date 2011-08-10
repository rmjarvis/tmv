

#ifndef TMV_ElemMultVV_H
#define TMV_ElemMultVV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"
#include "TMV_MultXV_Funcs.h"
#include "TMV_MultVV_Funcs.h"

namespace tmv {

    // Defined in TMV_ElemMultVV.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstElemMultVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddElemMultVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasElemMultVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddElemMultVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    //
    // ElementProd functions:
    //

    template <int algo, int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper;

    // algo 0: s == 0, nothing to do
    template <bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<0,0,add,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& , const V1& , const V2& , V3& ) 
        {}
    };

    // algo 11: simple for loop
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<11,s,add,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == TMV_UNKNOWN ? int(v3.size()) : s;
            call2(n,x1,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(
            int n, const Scaling<ix,T>& x1, IT1 it1, IT2 it2, IT3 it3)
        {
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;
            if (n) do {
                Maybe<add>::add(
                    *it3++ , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(*it1++ , *it2++)));
            } while (--n);
        }
    };  

    // algo 12: 2 at a time
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<12,s,add,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == TMV_UNKNOWN ? int(v3.size()) : s;
            call2(n,x1,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(
            const int n, const Scaling<ix,T>& x1, IT1 x, IT2 y, IT3 z)
        {
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;

            if (n_2) do {
                Maybe<add>::add(
                    *z ,
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(*x , *y)) );
                Maybe<add>::add(
                    z[1] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[1] , y[1]) ));
                x += 2; y += 2; z += 2;
            } while (--n_2);
            if (nb) 
                Maybe<add>::add(
                    *z , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(*x , *y)) );
        }
    };

    // algo 13: 4 at a time
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<13,s,add,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == TMV_UNKNOWN ? int(v3.size()) : s;
            call2(n,x1,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(
            const int n, const Scaling<ix,T>& x1, IT1 x, IT2 y, IT3 z)
        {
            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;

            if (n_4) do {
                Maybe<add>::add(
                    *z , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(*x , *y)));
                Maybe<add>::add(
                    z[1] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[1] , y[1])));
                Maybe<add>::add(
                    z[2] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[2] , y[2])));
                Maybe<add>::add(
                    z[3] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[3] , y[3])));
                x += 4; y += 4; z += 4;
            } while (--n_4);
            if (nb) do {
                Maybe<add>::add(
                    *z++ , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(*x++ , *y++)));
            } while (--nb);
        }
    };

    // algo 14: 8 at a time
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<14,s,add,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == TMV_UNKNOWN ? int(v3.size()) : s;
            call2(n,x1,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(
            const int n, const Scaling<ix,T>& x1, IT1 x, IT2 y, IT3 z)
        {
            int n_8 = (n>>3);
            int nb = n-(n_8<<3);
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;

            if (n_8) do {
                Maybe<add>::add(
                    *z , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(*x , *y)));
                Maybe<add>::add(
                    z[1] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[1] , y[1])));
                Maybe<add>::add(
                    z[2] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[2] , y[2])));
                Maybe<add>::add(
                    z[3] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[3] , y[3])));
                Maybe<add>::add(
                    z[4] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[4] , y[4])));
                Maybe<add>::add(
                    z[5] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[5] , y[5])));
                Maybe<add>::add(
                    z[6] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[6] , y[6])));
                Maybe<add>::add(
                    z[7] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[7] , y[7])));
                x += 8; y += 8; z += 8;
            } while (--n_8);
            if (nb) do {
                Maybe<add>::add(
                    *z++ , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(*x++ , *y++)));
            } while (--nb);
        }
    };

    // algo 15: fully unroll
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<15,s,add,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        template <int I, int N>
        struct Unroller
        {
            static inline void unroll(
                const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
            {
                Unroller<I,N/2>::unroll(x1,v1,v2,v3);
                Unroller<I+N/2,N-N/2>::unroll(x1,v1,v2,v3);
            }
            static inline void unroll2(
                const Scaling<ix,T>& x1,
                const IT1& x, const IT2& y, const IT3& z)
            {
                Unroller<I,N/2>::unroll2(x1,x,y,z);
                Unroller<I+N/2,N-N/2>::unroll2(x1,x,y,z);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
            {
                Maybe<add>::add( 
                    v3.ref(I) , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<false,false>(v1.cref(I) , v2.cref(I)) )); 
            }
            static inline void unroll2(
                const Scaling<ix,T>& x1,
                const IT1& x, const IT2& y, const IT3& z)
            {
                const bool c1 = V1::_conj;
                const bool c2 = V2::_conj;
                Maybe<add>::add( 
                    z[I] , 
                    ZProd<false,false>::prod(
                        x1 , ZProd<c1,c2>::prod(x[I],y[I]))); 
            }
        };
        template <int I>
        struct Unroller<I,0>
        {
            static inline void unroll(
                const Scaling<ix,T>& , const V1& , const V2& , V3& ) {}
            static inline void unroll2(
                const Scaling<ix,T>& , const IT1& , const IT2& , const IT3& )
            {}
        };
        static inline void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        { Unroller<0,s>::unroll(x1,v1,v2,v3); }
        static inline void call2(
            const int , const Scaling<ix,T>& x1,
            const IT1& x, const IT2& y, const IT3& z)
        { Unroller<0,s>::unroll2(x1,x,y,z); }
    };

    // algo 90: Call inst
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<90,s,true,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x1));
            InstAddElemMultVV(xx,v1.xView(),v2.xView(),v3.xView()); 
        }
    };
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<90,s,false,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x1));
            InstElemMultVV(xx,v1.xView(),v2.xView(),v3.xView()); 
        }
    };

    // algo 91: Call inst alias
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<91,s,true,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x1));
            InstAliasAddElemMultVV(xx,v1.xView(),v2.xView(),v3.xView()); 
        }
    };
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<91,s,false,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x1));
            InstAliasElemMultVV(xx,v1.xView(),v2.xView(),v3.xView()); 
        }
    };

    // algo 97: Conjugate
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<97,s,add,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            ElemMultVV_Helper<-2,s,add,ix,T,V1c,V2c,V3c>::call(
                TMV_CONJ(x1),v1c,v2c,v3c);
        }
    };

    // algo 197: Conjugate
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<197,s,add,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            ElemMultVV_Helper<99,s,add,ix,T,V1c,V2c,V3c>::call(
                TMV_CONJ(x1),v1c,v2c,v3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<98,s,add,ix,T,V1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            const bool noclobber1 =
                !SameStorage(v1,v3) ||
                ExactSameStorage(v1,v3) || 
                v1.step()*v3.step() < 0 || 
                std::abs(v3.step()) < std::abs(v1.step());
            const bool noclobber2 =
                !SameStorage(v2,v3) ||
                ExactSameStorage(v2,v3) || 
                v2.step()*v3.step() < 0 || 
                std::abs(v3.step()) < std::abs(v2.step());
            if (noclobber1) {
                if (noclobber2) {
                    // No aliasing (or no clobering)
                    ElemMultVV_Helper<-2,s,add,ix,T,V1,V2,V3>::call(
                        x1,v1,v2,v3); 
                } else {
                    // Need a temporary for v2
                    NoAliasElemMultVV<add>(x1,v1,v2.copy(),v3);
                }
            } else {
                if (noclobber2) {
                    // Need a temporary for v1
                    NoAliasElemMultVV<add>(x1,v1.copy(),v2,v3);
                } else {
                    // Need a temporary for v3
                    typename V3::copy_type v3c(v3.size());
                    Scaling<1,typename Traits<T>::real_type> one;
                    NoAliasElemMultVV<false>(one,v1,v2,v3c);
                    NoAliasMultXV<add>(x1,v3c,v3);
                }
            }
        }
    };

    // algo 99: Check for aliases
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<99,s,add,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst =
                (s == TMV_UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                V3::_conj ? 197 : 
                inst ? 91 : 
                98;
            ElemMultVV_Helper<algo,s,add,ix,T,V1,V2,V3>::call(x1,v1,v2,v3); 
        }
    };

    // algo -4: No branches or copies
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<-4,s,add,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        typedef typename V3::real_type RT;
        typedef typename V3::value_type VT;
        enum { allunit =
            V1::_step == 1 && V2::_step == 1 && V3::_step == 1 };
        enum { algo = (
                s == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                // Unrolling doesn't ever seem to be faster for this function.
                //(s != TMV_UNKNOWN && s <= int(128/sizeof(VT))) ? 15 :
                (allunit && sizeof(RT) == 8) ? V3::iscomplex ? 12 : 13 :
                (allunit && sizeof(RT) == 4) ? V3::iscomplex ? 13 : 14 :
                11 ) };
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            ElemMultVV_Helper<algo,s,add,ix,T,V1,V2,V3>::call(x1,v1,v2,v3); 
        }
        static TMV_INLINE void call2(
            const int n,
            const Scaling<ix,T>& x1, const IT1& x, const IT2& y, IT3& z)
        {
            TMVStaticAssert(!V3::_conj);
            ElemMultVV_Helper<algo,s,add,ix,T,V1,V2,V3>::call2(n,x1,x,y,z); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<-3,s,add,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        { ElemMultVV_Helper<-4,s,add,ix,T,V1,V2,V3>::call(x1,v1,v2,v3); }
    };

    // algo -2: Check for inst
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<-2,s,add,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst =
                (s == TMV_UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                V3::_conj ? 97 : 
                inst ? 90 : 
                -4;
            ElemMultVV_Helper<algo,s,add,ix,T,V1,V2,V3>::call(
                x1,v1,v2,v3); 
        }
    };

    // algo -1: Check for aliases?
    template <int s, bool add, int ix, class T, class V1, class V2, class V3>
    struct ElemMultVV_Helper<-1,s,add,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x1, const V1& v1, const V2& v2, V3& v3)
        {
            const bool noclobber =
                VStepHelper<V1,V3>::noclobber &&
                VStepHelper<V2,V3>::noclobber;
            const bool checkalias =
                V3::_checkalias && !noclobber;
            const int algo =
                checkalias ? 99 : 
                -2;
            ElemMultVV_Helper<algo,s,add,ix,T,V1,V2,V3>::call(
                x1,v1,v2,v3); 
        }
    };

    template <bool add, int ix, class T, class V1, class V2, class V3>
    static inline void ElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::_size,V2::_size>::size;
        const int s = Sizes<s12,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        ElemMultVV_Helper<-1,s,add,ix,T,V1v,V2v,V3v>::call(
            x1,v1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class V1, class V2, class V3>
    static inline void NoAliasElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::_size,V2::_size>::size;
        const int s = Sizes<s12,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        ElemMultVV_Helper<-2,s,add,ix,T,V1v,V2v,V3v>::call(
            x1,v1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class V1, class V2, class V3>
    static inline void InlineElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::_size,V2::_size>::size;
        const int s = Sizes<s12,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        ElemMultVV_Helper<-4,s,add,ix,T,V1v,V2v,V3v>::call(
            x1,v1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class V1, class V2, class V3>
    static inline void InlineAliasElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::_size,V2::_size>::size;
        const int s = Sizes<s12,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        ElemMultVV_Helper<98,s,add,ix,T,V1v,V2v,V3v>::call(
            x1,v1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class V1, class V2, class V3>
    static inline void AliasElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::_size,V2::_size>::size;
        const int s = Sizes<s12,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        ElemMultVV_Helper<99,s,add,ix,T,V1v,V2v,V3v>::call(
            x1,v1v,v2v,v3v);
    }

} // namespace tmv

#endif 
