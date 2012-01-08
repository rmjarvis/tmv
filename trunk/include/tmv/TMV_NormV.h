
#ifndef TMV_NormV_H
#define TMV_NormV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"

namespace tmv {

    // TODO: Write SSE algo's.
    // The real version of SumAbsElements is one routine where the 
    // BLAS function (dasum) on my computer is significantly 
    // (~30%) faster than TMV.
    // I haven't figured out any way to make the function faster though.
    // Maybe the SSE algos will be the difference.
    
    // 
    // SumElements
    // SumAbsElements
    // SumAbs2Elements
    // NormSq
    //

    // Defined in TMV_Vector.cpp
    template <class T>
    T InstSumElements(const ConstVectorView<T>& v); 
    template <class T>
    typename ConstVectorView<T>::float_type InstSumAbsElements(
        const ConstVectorView<T>& v);
    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstVectorView<T>& v);
    template <class T>
    typename Traits<T>::real_type InstNormSq(const ConstVectorView<T>& v); 
    template <class T>
    typename ConstVectorView<T>::float_type InstNormSq(
        const ConstVectorView<T>& v, 
        typename ConstVectorView<T>::float_type scale); 

    // UNROLL is the maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_NORMV_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_NORMV_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_NORMV_UNROLL 9
#else
#define TMV_NORMV_UNROLL 0
#endif

    template <int algo, int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper;

    // algo 11: simple for loop
    template <int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper<11,s,comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const V& v, const Scaling<ix,RT>& x)
        {
            typedef typename V::value_type VT;
            typedef typename TypeSelect<
                V::isreal || Traits<ret>::iscomplex , ret ,
                std::complex<ret> >::type BT;
            const int n = s == Unknown ? v.size() : s;
            ret sum(0);
            for(int i=0;i<n;++i) 
                sum += Component<comp,BT>::f(
                    x * Traits<BT>::convert(v.cref(i)));
            return sum;
        }
    };

    // algo 12: 2 at a time
    template <int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper<12,s,comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const V& v, const Scaling<ix,RT>& x)
        {
            typedef typename V::value_type VT;
            typedef typename TypeSelect<
                V::isreal || Traits<ret>::iscomplex , ret ,
                std::complex<ret> >::type BT;
            const int n = s == Unknown ? v.size() : s;
            ret sum0(0), sum1(0);
            BT v0, v1;
            typedef typename V::const_iterator IT;
            IT it = v.begin();
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);

            if (n_2) {
                do {
                    v0 = x * Traits<BT>::convert(it[0]); 
                    v1 = x * Traits<BT>::convert(it[1]); it += 2;
                    Component<comp,BT>::applyf(v0);
                    Component<comp,BT>::applyf(v1);
                    sum0 += Component<comp,BT>::get(v0);
                    sum1 += Component<comp,BT>::get(v1);
                } while (--n_2);
                sum0 += sum1;
            }
            if (nb) {
                v0 = x * Traits<BT>::convert(it[0]); 
                Component<comp,BT>::applyf(v0);
                sum0 += Component<comp,BT>::get(v0);
            }
            return sum0;
        }
    };

    // algo 13: 4 at a time
    template <int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper<13,s,comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const V& v, const Scaling<ix,RT>& x)
        {
            typedef typename V::value_type VT;
            typedef typename TypeSelect<
                V::isreal || Traits<ret>::iscomplex , ret ,
                std::complex<ret> >::type BT;
            const int n = s == Unknown ? v.size() : s;
            ret sum0(0), sum1(0);
            BT v0, v1;
            typedef typename V::const_iterator IT;
            IT it = v.begin();

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);

            if (n_4) {
                do {
                    v0 = x * Traits<BT>::convert(it[0]); 
                    v1 = x * Traits<BT>::convert(it[1]); 
                    Component<comp,BT>::applyf(v0);
                    Component<comp,BT>::applyf(v1);
                    sum0 += Component<comp,BT>::get(v0);
                    sum1 += Component<comp,BT>::get(v1);

                    v0 = x * Traits<BT>::convert(it[2]); 
                    v1 = x * Traits<BT>::convert(it[3]); it += 4;
                    Component<comp,BT>::applyf(v0);
                    Component<comp,BT>::applyf(v1);
                    sum0 += Component<comp,BT>::get(v0);
                    sum1 += Component<comp,BT>::get(v1);
                } while (--n_4);
                sum0 += sum1;
            }
            if (nb) do {
                v0 = x * Traits<BT>::convert(*it++); 
                Component<comp,BT>::applyf(v0);
                sum0 += Component<comp,BT>::get(v0);
            } while (--nb);
            return sum0;
        }
    };

    // algo 14: 8 at a time
    template <int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper<14,s,comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const V& v, const Scaling<ix,RT>& x)
        {
            typedef typename V::value_type VT;
            typedef typename TypeSelect<
                V::isreal || Traits<ret>::iscomplex , ret ,
                std::complex<ret> >::type BT;
            const int n = s == Unknown ? v.size() : s;
            ret sum0(0), sum1(0);
            BT v0, v1;
            typedef typename V::const_iterator IT;
            IT it = v.begin();

            int n_8 = (n>>3);
            int nb = n-(n_8<<3);

            if (n_8) {
                do {
                    v0 = x * Traits<BT>::convert(it[0]); 
                    v1 = x * Traits<BT>::convert(it[1]); 
                    Component<comp,BT>::applyf(v0);
                    Component<comp,BT>::applyf(v1);
                    sum0 += Component<comp,BT>::get(v0);
                    sum1 += Component<comp,BT>::get(v1);

                    v0 = x * Traits<BT>::convert(it[2]); 
                    v1 = x * Traits<BT>::convert(it[3]);
                    Component<comp,BT>::applyf(v0);
                    Component<comp,BT>::applyf(v1);
                    sum0 += Component<comp,BT>::get(v0);
                    sum1 += Component<comp,BT>::get(v1);

                    v0 = x * Traits<BT>::convert(it[4]); 
                    v1 = x * Traits<BT>::convert(it[5]); 
                    Component<comp,BT>::applyf(v0);
                    Component<comp,BT>::applyf(v1);
                    sum0 += Component<comp,BT>::get(v0);
                    sum1 += Component<comp,BT>::get(v1);

                    v0 = x * Traits<BT>::convert(it[6]); 
                    v1 = x * Traits<BT>::convert(it[7]); it += 8;
                    Component<comp,BT>::applyf(v0);
                    Component<comp,BT>::applyf(v1);
                    sum0 += Component<comp,BT>::get(v0);
                    sum1 += Component<comp,BT>::get(v1);
                } while (--n_8);
                sum0 += sum1;
            }
            if (nb) do {
                v0 = x * Traits<BT>::convert(*it++); 
                Component<comp,BT>::applyf(v0);
                sum0 += Component<comp,BT>::get(v0);
            } while (--nb);
            return sum0;
        }
    };

    // algo 15: fully unroll
    template <int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper<15,s,comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        typedef typename V::value_type VT;
        typedef typename TypeSelect<
            V::isreal || Traits<ret>::iscomplex , ret ,
            std::complex<ret> >::type BT;
        template <int I, int N>
        struct Unroller
        {
            static TMV_INLINE ret unroll(const V& v, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,N/2>::unroll(v,x) +
                    Unroller<I+N/2,N-N/2>::unroll(v,x));
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static TMV_INLINE ret unroll(const V& v, const Scaling<ix,RT>& x)
            {
                return Component<comp,BT>::f(
                    x * Traits<BT>::convert(v.cref(I)));
            }
        };
        template <int I>
        struct Unroller<I,0>
        {
            static TMV_INLINE ret unroll(const V& v, const Scaling<ix,RT>& x)
            { return ret(0); }
        };
        static inline ret call(const V& v, const Scaling<ix,RT>& x)
        { return Unroller<0,s>::unroll(v,x); }
    };

    // algo 90: Call inst
    template <int s, int ix, class ret, class V>
    struct SumElementsV_Helper<90,s,ValueComp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static TMV_INLINE ret call(const V& v, const Scaling<ix,RT>& x)
        { return InstSumElements(v.xView()); }
    };
    template <int s, int ix, class ret, class V>
    struct SumElementsV_Helper<90,s,AbsComp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static TMV_INLINE ret call(const V& v, const Scaling<ix,RT>& x)
        { return InstSumAbsElements(v.xView()); }
    };
    template <int s, int ix, class ret, class V>
    struct SumElementsV_Helper<90,s,Abs2Comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static TMV_INLINE ret call(const V& v, const Scaling<ix,RT>& x)
        { return InstSumAbs2Elements(v.xView()); }
    };
    template <int s, class ret, class V>
    struct SumElementsV_Helper<90,s,NormComp,1,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static TMV_INLINE ret call(const V& v, const Scaling<1,RT>& x)
        { return InstNormSq(v.xView()); }
    };
    template <int s, class ret, class V>
    struct SumElementsV_Helper<90,s,NormComp,0,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static TMV_INLINE ret call(const V& v, const Scaling<0,RT>& x)
        { return InstNormSq(v.xView(),x); }
    };

    // algo 97: Conjugate
    template <int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper<97,s,comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static TMV_INLINE ret call(const V& v, const Scaling<ix,RT>& x)
        {
            typedef typename V::const_nonconj_type Vnc;
            Vnc vnc = v.nonConj();
            return SumElementsV_Helper<-2,s,comp,ix,ret,Vnc>::call(vnc,x);
        }
    };
    template <int s, int ix, class ret, class V>
    struct SumElementsV_Helper<97,s,ValueComp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static TMV_INLINE ret call(const V& v, const Scaling<ix,RT>& x)
        {
            typedef typename V::const_conjugate_type Vc;
            Vc vc = v.conjugate();
            return TMV_CONJ(
                SumElementsV_Helper<-2,s,ValueComp,ix,ret,Vc>::call(vc,x));
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper<-3,s,comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        typedef typename V::value_type VT;
        static TMV_INLINE ret call(const V& v, const Scaling<ix,RT>& x)
        {
            typedef typename V::real_type RT;
            const int maxunroll = 80;
            const int algo = 
                TMV_OPT == 0 ? 11 :
                ( s != Unknown && s <= maxunroll ) ? 15 :
                (sizeof(RT) == 8 && V::_step == 1) ? (V::iscomplex ? 12 : 13) :
                (sizeof(RT) == 4 && V::_step == 1) ? (V::iscomplex ? 13 : 14) :
                11;
            return SumElementsV_Helper<algo,s,comp,ix,ret,V>::call(v,x);
        }
    };

    // algo -2: Check for inst
    template <int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper<-2,s,comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static TMV_INLINE ret call(const V& v, const Scaling<ix,RT>& x)
        {
            typedef typename V::value_type VT;
            const bool inst = 
                (s == Unknown || s > 16) &&
                Traits<VT>::isinst;
            const int algo = 
                V::_conj ? 97 :
                inst ? 90 :
                -3;
            return SumElementsV_Helper<algo,s,comp,ix,ret,V>::call(v,x);
        }
    };

    template <int s, CompType comp, int ix, class ret, class V>
    struct SumElementsV_Helper<-1,s,comp,ix,ret,V>
    {
        typedef typename Traits<ret>::real_type RT;
        static TMV_INLINE ret call(const V& v, const Scaling<ix,RT>& x)
        { return SumElementsV_Helper<-2,s,comp,ix,ret,V>::call(v,x); }
    };

    template <class V>
    inline typename V::value_type InlineSumElements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-3,V::_size,ValueComp,1,VT,Vv>::call(
            vv,Scaling<1,RT>());
    }

    template <class V>
    inline typename V::value_type DoSumElements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-2,V::_size,ValueComp,1,VT,Vv>::call(
            vv,Scaling<1,RT>());
    }

    template <class V>
    inline typename V::float_type InlineSumAbsElements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::float_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-3,V::_size,AbsComp,1,RT,Vv>::call(
            vv,Scaling<1,RT>());
    }

    template <class V>
    inline typename V::float_type DoSumAbsElements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::float_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-2,V::_size,AbsComp,1,RT,Vv>::call(
            vv,Scaling<1,RT>());
    }

    template <class V>
    inline typename V::real_type InlineSumAbs2Elements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-3,V::_size,Abs2Comp,1,RT,Vv>::call(
            vv,Scaling<1,RT>());
    }

    template <class V>
    inline typename V::real_type DoSumAbs2Elements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-2,V::_size,Abs2Comp,1,RT,Vv>::call(
            vv,Scaling<1,RT>());
    }

    template <class V>
    inline typename V::real_type InlineNormSq(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-3,V::_size,NormComp,1,RT,Vv>::call(
            vv,Scaling<1,RT>());
    }

    template <class V>
    inline typename V::real_type DoNormSq(const BaseVector_Calc<V>& v)
    {
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-2,V::_size,NormComp,1,RT,Vv>::call(
            vv,Scaling<1,RT>());
    }

    template <class V>
    inline typename V::float_type InlineNormSq(
        const BaseVector_Calc<V>& v, typename V::float_type scale)
    {
        typedef typename V::float_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-3,V::_size,NormComp,0,RT,Vv>::call(
            vv,Scaling<0,RT>(scale));
    }

    template <class V>
    inline typename V::float_type DoNormSq(
        const BaseVector_Calc<V>& v, const typename V::float_type scale)
    {
        typedef typename V::float_type RT;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return SumElementsV_Helper<-2,V::_size,NormComp,0,RT,Vv>::call(
            vv,Scaling<0,RT>(scale));
    }

} // namespace tmv

#endif
