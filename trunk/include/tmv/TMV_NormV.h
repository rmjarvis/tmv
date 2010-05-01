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

#ifndef TMV_NormV_H
#define TMV_NormV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"

namespace tmv {

    // TODO: Convert to new algo numbering scheme
    // TODO: Write SSE algo's.
    
    // 
    // SumElements
    //

    template <int algo, int size, CompType comp, int ix, class V>
    struct SumElementsV_Helper;

    // algo 1: simple for loop
    template <int size, CompType comp, int ix, class V>
    struct SumElementsV_Helper<1,size,comp,ix,V> 
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<VT>::type ret;

        static inline ret call(const V& v, const Scaling<ix,RT>& x)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            ret sum(0);
            for(int i=0;i<n;++i) sum += Component<comp,VT>::f(x * v.cref(i));
            return sum;
        }
    };

    // algo 2: 2 at a time
    template <int size, CompType comp, int ix, class V>
    struct SumElementsV_Helper<2,size,comp,ix,V> 
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<VT>::type ret;

        static inline ret call(const V& v, const Scaling<ix,RT>& x)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            ret sum0(0), sum1(0);
            VT v0, v1;
            typedef typename V::const_iterator IT;
            IT it = v.begin();
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);

            if (n_2) {
                do {
                    v0 = x * *it; v1 = x*it[1]; it += 2;
                    Component<comp,VT>::applyf(v0);
                    Component<comp,VT>::applyf(v1);
                    sum0 += Component<comp,VT>::get(v0);
                    sum1 += Component<comp,VT>::get(v1);
                } while (--n_2);
                sum0 += sum1;
            }
            if (nb) {
                v0 = x * *it;
                Component<comp,VT>::applyf(v0);
                sum0 += Component<comp,VT>::get(v0);
            }
            return sum0;
        }
    };

    // algo 3: 4 at a time
    template <int size, CompType comp, int ix, class V>
    struct SumElementsV_Helper<3,size,comp,ix,V> 
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<VT>::type ret;

        static inline ret call(const V& v, const Scaling<ix,RT>& x)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            ret sum0(0), sum1(0);
            VT v0, v1;
            typedef typename V::const_iterator IT;
            IT it = v.begin();

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);

            if (n_4) {
                do {
                    v0 = x*it[0]; v1 = x*it[1];
                    Component<comp,VT>::applyf(v0);
                    Component<comp,VT>::applyf(v1);
                    sum0 += Component<comp,VT>::get(v0);
                    sum1 += Component<comp,VT>::get(v1);

                    v0 = x*it[2]; v1 = x*it[3]; it += 4;
                    Component<comp,VT>::applyf(v0);
                    Component<comp,VT>::applyf(v1);
                    sum0 += Component<comp,VT>::get(v0);
                    sum1 += Component<comp,VT>::get(v1);
                } while (--n_4);
                sum0 += sum1;
            }
            if (nb) do {
                v0 = x*(*it++);
                Component<comp,VT>::applyf(v0);
                sum0 += Component<comp,VT>::get(v0);
            } while (--nb);
            return sum0;
        }
    };

    // algo 4: 8 at a time
    template <int size, CompType comp, int ix, class V>
    struct SumElementsV_Helper<4,size,comp,ix,V> 
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<VT>::type ret;

        static inline ret call(const V& v, const Scaling<ix,RT>& x)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            ret sum0(0), sum1(0);
            VT v0, v1;
            typedef typename V::const_iterator IT;
            IT it = v.begin();

            int n_8 = (n>>3);
            int nb = n-(n_8<<3);

            if (n_8) {
                do {
                    v0 = x*it[0]; v1 = x*it[1];
                    Component<comp,VT>::applyf(v0);
                    Component<comp,VT>::applyf(v1);
                    sum0 += Component<comp,VT>::get(v0);
                    sum1 += Component<comp,VT>::get(v1);

                    v0 = x*it[2]; v1 = x*it[3];
                    Component<comp,VT>::applyf(v0);
                    Component<comp,VT>::applyf(v1);
                    sum0 += Component<comp,VT>::get(v0);
                    sum1 += Component<comp,VT>::get(v1);

                    v0 = x*it[4]; v1 = x*it[5];
                    Component<comp,VT>::applyf(v0);
                    Component<comp,VT>::applyf(v1);
                    sum0 += Component<comp,VT>::get(v0);
                    sum1 += Component<comp,VT>::get(v1);

                    v0 = x*it[6]; v1 = x*it[7]; it += 8;
                    Component<comp,VT>::applyf(v0);
                    Component<comp,VT>::applyf(v1);
                    sum0 += Component<comp,VT>::get(v0);
                    sum1 += Component<comp,VT>::get(v1);
                } while (--n_8);
                sum0 += sum1;
            }
            if (nb) do {
                v0 = x*(*it++);
                Component<comp,VT>::applyf(v0);
                sum0 += Component<comp,VT>::get(v0);
            } while (--nb);
            return sum0;
        }
    };

    // algo 5: fully unroll
    template <int size, CompType comp, int ix, class V>
    struct SumElementsV_Helper<5,size,comp,ix,V> // known size, unroll
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<VT>::type ret;

        template <int I, int N>
        struct Unroller
        {
            static inline ret unroll(const V& v, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,N/2>::unroll(v,x) +
                    Unroller<I+N/2,N-N/2>::unroll(v,x));
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline ret unroll(const V& v, const Scaling<ix,RT>& x)
            { return Component<comp,VT>::f(x * v.cref(I)); }
        };
        template <int I>
        struct Unroller<I,0>
        { 
            static inline ret unroll(const V& v, const Scaling<ix,RT>& x)
            { return ret(0); }
        };
        static inline ret call(const V& v, const Scaling<ix,RT>& x)
        { return Unroller<0,size>::unroll(v,x); }
    };

    // algo 7: complex with unit step, convert to real
    // (This only makes sense for comp = NORM and ABS2.)
    template <int size, CompType comp, int ix, class V>
    struct SumElementsV_Helper<7,size,comp,ix,V> 
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<VT>::type ret;

        static inline ret call(const V& v, const Scaling<ix,RT>& x)
        {
            typedef typename V::const_flatten_type Vf;
            Vf vf = v.flatten();
            const int size2 = IntTraits<size>::twoS;
            const int algo2 = 
#if TMV_OPT >= 1
                (size2 != UNKNOWN && size2 <= 64) ? 5 :
                sizeof(RT) == 8 ? 3 :
                sizeof(RT) == 4 ? 4 :
#endif
                1;
            return SumElementsV_Helper<algo2,size2,comp,ix,Vf>::call(vf,x);
        }
    };

    // algo -1: Determine which algorithm to use
    template <int size, CompType comp, int ix, class V>
    struct SumElementsV_Helper<-1,size,comp,ix,V> 
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<VT>::type ret;

        static inline ret call(const V& v, const Scaling<ix,RT>& x)
        {
#if TMV_OPT >= 1
            const int maxunroll = 
                comp == ValueComp ? 80 :
                comp == AbsComp ? 64 :
                comp == Abs2Comp ? 64 :
                comp == NormComp ? 64 : 
                /* no other options currently */ 0;
            const int algo = 
                ( ( comp == Abs2Comp || comp == NormComp ) && 
                  ( V::iscomplex && V::_step == 1) ) ? 7 :
                ( V::_size != UNKNOWN && V::_size <= int(maxunroll) ) ? 5 :
                (sizeof(RT) == 8 && V::_step == 1) ? (V::iscomplex ? 2 : 3) :
                (sizeof(RT) == 4 && V::_step == 1) ? (V::iscomplex ? 3 : 4) :
                1;
#else 
            const int algo = 1;
#endif
            return SumElementsV_Helper<algo,size,comp,ix,V>::call(v,x);
        }
    };

    template <class V>
    inline typename V::value_type InlineSumElements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return SumElementsV_Helper<-1,V::_size,ValueComp,1,Vv>::call(
            vv,Scaling<1,RT>());
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    T InstSumElements(const ConstVectorView<T>& v); 

    template <bool conj, bool inst, class V>
    struct CallSumElementsv // inst = false
    {
        static inline typename V::value_type call(const V& v)
        {
            typedef typename V::const_conjugate_type Vc;
            return TMV_CONJ(
                CallSumElementsv<false,inst,Vc>::call(v.conjugate()));
        }
    };
    template <class V>
    struct CallSumElementsv<false,false,V>
    {
        static inline typename V::value_type call(const V& v)
        { return InlineSumElements(v); }
    };
    template <class V>
    struct CallSumElementsv<false,true,V>
    {
        static inline typename V::value_type call(const V& v)
        { return InstSumElements(v.xView()); }
    };

    template <class V>
    inline typename V::value_type SumElements(const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type T;
        const bool inst = 
            V::unknownsizes &&
            Traits<T>::isinst;
        return CallSumElementsv<V::_conj,inst,V>::call(v.vec());
    }


    // 
    // SumAbsElements
    //

    // TODO: This (the real version) is one routine where the 
    // BLAS function (dasum) on my computer is significantly 
    // (~30%) faster than TMV.
    // I haven't figured out any way to make the above functions faster though.
    template <class V>
    inline typename V::real_type InlineSumAbsElements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return SumElementsV_Helper<-1,V::_size,AbsComp,1,Vv>::call(
            vv,Scaling<1,RT>());
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    typename Traits<T>::real_type InstSumAbsElements(
        const ConstVectorView<T>& v);

    template <bool inst, class V>
    struct CallSumAbsElementsv // inst = false
    {
        static inline typename V::real_type call(const V& v)
        { return InlineSumAbsElements(v); }
    };
    template <class V>
    struct CallSumAbsElementsv<true,V>
    {
        static inline typename V::real_type call(const V& v)
        { return InstSumAbsElements(v.xView()); }
    };

    template <class V>
    inline typename V::real_type SumAbsElements(const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type T;
        typedef typename V::const_nonconj_type Vn;
        const bool inst = 
            V::unknownsizes &&
            Traits<T>::isinst;
        return CallSumAbsElementsv<inst,Vn>::call(v.nonConj());
    }


    // 
    // SumAbs2Elements
    //

    template <class V>
    inline typename V::real_type InlineSumAbs2Elements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return SumElementsV_Helper<-1,V::_size,Abs2Comp,1,Vv>::call(
            vv,Scaling<1,RT>());
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstVectorView<T>& v);

    template <bool inst, class V>
    struct CallSumAbs2Elementsv // inst = false
    {
        static inline typename V::real_type call(const V& v)
        { return InlineSumAbs2Elements(v); }
    };
    template <class V>
    struct CallSumAbs2Elementsv<true,V>
    {
        static inline typename V::real_type call(const V& v)
        { return InstSumAbs2Elements(v.xView()); }
    };

    template <class V>
    inline typename V::real_type SumAbs2Elements(const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type T;
        typedef typename V::const_nonconj_type Vn;
        const bool inst = 
            V::unknownsizes &&
            Traits<T>::isinst;
        return CallSumAbs2Elementsv<inst,Vn>::call(v.nonConj());
    }


    // 
    // NormSq
    //

    template <class V>
    inline typename V::real_type InlineNormSq(const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return SumElementsV_Helper<-1,V::_size,NormComp,1,Vv>::call(
            vv,Scaling<1,RT>());
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    typename Traits<T>::real_type InstNormSq(const ConstVectorView<T>& v); 

    template <bool inst, class V>
    struct CallNormSqv // inst = false
    {
        static inline typename V::real_type call(const V& v)
        { return InlineNormSq(v); }
    };
    template <class V>
    struct CallNormSqv<true,V>
    {
        static inline typename V::real_type call(const V& v)
        { return InstNormSq(v.xView()); }
    };

    template <class V>
    inline typename V::real_type NormSq(const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type T;
        typedef typename V::const_nonconj_type Vn;
        const bool inst = 
            V::unknownsizes &&
            Traits<T>::isinst;
        return CallNormSqv<inst,Vn>::call(v.nonConj());
    }


    // 
    // NormSq with scaling
    //

    template <class V>
    inline typename V::real_type InlineNormSq(
        const BaseVector_Calc<V>& v, typename V::real_type scale)
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return SumElementsV_Helper<-1,V::_size,NormComp,0,Vv>::call(
            vv,Scaling<0,RT>(scale));
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    typename Traits<T>::real_type InstNormSq(
        const ConstVectorView<T>& v, typename Traits<T>::real_type scale); 

    template <bool inst, class V>
    struct CallNormSq_scalev // inst = false
    {
        static inline typename V::real_type call(
            const V& v, const typename V::real_type scale)
        { return InlineNormSq(v,scale); }
    };
    template <class V>
    struct CallNormSq_scalev<true,V>
    {
        static inline typename V::real_type call(
            const V& v, const typename V::real_type scale)
        { return InstNormSq(v.xView(),scale); }
    };

    template <class V>
    inline typename V::real_type NormSq(
        const BaseVector_Calc<V>& v, const typename V::real_type scale)
    {
        typedef typename V::value_type T;
        typedef typename V::const_nonconj_type Vn;
        const bool inst = 
            V::unknownsizes &&
            Traits<T>::isinst;
        return CallNormSq_scalev<inst,Vn>::call(v.nonConj(),scale);
    }


    //
    // Norm2
    //

    // This helper struct works for either Vector or Matrix "V"
    template <int algo, class V> 
    struct Norm_Helper;

    // algo 1: simple: sqrt(NormSq(v))
    template <class V> 
    struct Norm_Helper<1,V>
    {
        typedef typename V::real_type RT;
        static inline RT call(const V& v)
        { return TMV_SQRT(InlineNormSq(v)); }
    };

    // algo 2: Robust algorithm with checks for overflow and underflow.
    // This one always calls MaxAbsElement and then NormSq.
    // This is inefficient if there are no problems.
    // Since no problems is the usual case, I switched to the below
    // version (algo 3) that calls NormSq first and then redoes it if there
    // are problems.
    template <class V> 
    struct Norm_Helper<2,V>
    {
        typedef typename V::real_type RT;
        static RT call(const V& v)
        {
            const RT eps = TMV_Epsilon<RT>();

            // Start with the maximum |v(i)|.  It will tell us how (and if)
            // we need to use a scaling for NormSq().
            RT vmax = v.maxAbsElement();

            // If vmax = 0, then norm2 = 0:
            if (vmax == RT(0)) return RT(0);

            // If vmax^2 * eps = 0, but vmax != 0, then a naive NormSq()
            // will produce underflow rounding errors.  Find a better scaling.
            // eps is a pure power of 2, so no rounding errors from
            // rescaling by a power of eps.
            else if (vmax * vmax * eps == RT(0)) {
                const RT inveps = RT(1)/eps;
                RT scale = inveps;
                vmax *= scale;
                RT eps2 = eps*eps;
                while (vmax < eps2) { scale *= inveps; vmax *= inveps; }
                return TMV_SQRT(v.normSq(scale))/scale;
            }

            // If 1/vmax == 0, then vmax is already inf, so no hope of
            // making it more accurate.  (And need to check, since otherwise
            // the next section would cause an infinite loop.)
            else if (RT(1)/vmax == RT(0)) {
                return vmax;
            }

            // If 1/(vmax^2) == 0, then a naive NormSq() will produce 
            // overflow.  Find a better scaling.
            else if (RT(1)/(vmax*vmax) == RT(0)) {
                RT scale = eps;
                vmax *= scale;
                while (vmax > RT(1)) { scale *= eps; vmax *= eps; }
                return TMV_SQRT(v.normSq(scale))/scale;
            }

            // No problems with overflow or underflow.
            else return TMV_SQRT(v.normSq());
        }
    };

    // algo 3: Robust algorithm with checks for overflow and underflow.
    // This version is slower if there is a problem, but since
    // there usually isn't a problem, it is generally faster.
    template <class V> 
    struct Norm_Helper<3,V>
    {
        typedef typename V::real_type RT;
        static RT call(const V& v)
        {
            const RT eps = TMV_Epsilon<RT>();
            const RT vnormsq = v.normSq();

            if (vnormsq * eps == RT(0)) {
                // Possible underflow errors:

                // If vmax = 0, then norm2 = 0:
                RT vmax = v.maxAbsElement();
                if (vmax == RT(0)) return RT(0);

                // If vmax^2 * eps = 0, but vmax != 0, then vnormsq has
                // underflow rounding errors.  Find a better scaling.
                // eps is a pure power of 2, so no rounding errors from
                // rescaling by a power of eps.
                else if (vmax * vmax * eps == RT(0)) {
                    const RT inveps = RT(1)/eps;
                    RT scale = inveps;
                    vmax *= scale;
                    RT eps2 = eps*eps;
                    while (vmax < eps2) { scale *= inveps; vmax *= inveps; }
                    return TMV_SQRT(v.normSq(scale))/scale;
                }

                else return TMV_SQRT(vnormsq);
            }

            else if (RT(1)/vnormsq == RT(0)) {
                // Possible overflow errors:

                // If 1/vmax == 0, then vmax is already inf, so no hope of
                // making it more accurate.  (And need to check, since 
                // otherwise the next section would cause an infinite loop.)
                RT vmax = v.maxAbsElement();
                if (RT(1)/vmax == RT(0)) {
                    return vmax;
                }

                // If 1/(vmax^2) == 0, then vnormsq has overflow errors. 
                // Find a better scaling.
                else if (RT(1)/(vmax*vmax) == RT(0)) {
                    RT scale = eps;
                    vmax *= scale;
                    while (vmax > RT(1)) { scale *= eps; vmax *= eps; }
                    return TMV_SQRT(v.normSq(scale))/scale;
                }

                else return TMV_SQRT(vnormsq);
            }

            // No problems with overflow or underflow.
            else return TMV_SQRT(vnormsq);
        }
    };

    template <class V>
    inline typename V::real_type InlineNorm2(const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
#if TMV_OPT == 0
        return Norm_Helper<1,Vv>::call(vv);
#else
        if (v.size() == 0) return typename V::real_type(0);
        else return Norm_Helper<3,Vv>::call(vv);
#endif
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    typename Traits<T>::real_type InstNorm2(const ConstVectorView<T>& v);

    template <bool inst, class V>
    struct CallNorm2v // inst = false
    {
        static inline typename V::real_type call(const V& v)
        { return InlineNorm2(v); }
    };
    template <class V>
    struct CallNorm2v<true,V>
    {
        static inline typename V::real_type call(const V& v)
        { return InstNorm2(v.xView()); }
    };

    template <class V>
    inline typename V::real_type Norm2(const BaseVector_Calc<V>& v)
    {
        typedef typename V::value_type T;
        typedef typename V::const_nonconj_type Vn;
        const bool inst = 
            V::unknownsizes &&
            Traits<T>::isinst;
        return CallNorm2v<inst,Vn>::call(v.nonConj());
    }

} // namespace tmv

#endif
