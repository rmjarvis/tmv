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

#ifndef TMV_NormU_H
#define TMV_NormU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_NormV.h"
#include "TMV_MinMax.h"

namespace tmv {

    // TODO: Convert to new algo numbering scheme
    
    // 
    // SumElements
    //

    // Q1 is the maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_Q1 200 
#elif TMV_OPT >= 2
#define TMV_Q1 25
#elif TMV_OPT >= 1
#define TMV_Q1 9
#else
#define TMV_Q1 0
#endif

    template <int algo, int s, CompType comp, int ix, class M1>
    struct SumElementsU_Helper;

    // algo 1: m1 is unitdiag
    template <int s, CompType comp, int ix, class M1>
    struct SumElementsU_Helper<1,s,comp,ix,M1> 
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<MT>::type ret;

        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int N = (s == UNKNOWN ? m.size() : s);
            typedef typename M1::const_offdiag_type Mo;
            Mo mo = m.offDiag();
            const int sm1 = IntTraits2<s,-1>::sum;
            return SumElementsU_Helper<-1,sm1,comp,ix,Mo>::call(mo,x) 
                + Component<comp,RT>::f(RT(x)) * RT(N);
        }
    };

    // algo 2: unknown diag, figure out which it is.
    template <int s, CompType comp, int ix, class M1>
    struct SumElementsU_Helper<2,s,comp,ix,M1>
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<MT>::type ret;
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            if (m.isunit()) 
                return SumElementsU_Helper<1,s,comp,ix,M1>::call(m,x);
            else 
                return SumElementsU_Helper<-4,s,comp,ix,M1>::call(m,x);
        }
    };

    // algo 11: loop over columns
    template <int s, CompType comp, int ix, class M1>
    struct SumElementsU_Helper<11,s,comp,ix,M1> 
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<MT>::type ret;

        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int N = (s == UNKNOWN ? m.size() : s);
            typedef typename M1::const_col_sub_type Mc;
            ret sum(0);
            for(int j=0;j<N;++j) {
                Mc mc = m.get_col(j,0,j+1);
                sum += SumElementsV_Helper<-1,UNKNOWN,comp,ix,Mc>::call(mc,x);
            }
            return sum;
        }
    };

    // algo 12: loop over rows
    template <int s, CompType comp, int ix, class M1>
    struct SumElementsU_Helper<12,s,comp,ix,M1> 
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<MT>::type ret;

        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int N = (s == UNKNOWN ? m.size() : s);
            typedef typename M1::const_row_sub_type Mr;
            ret sum(0);
            for(int i=0;i<N;++i) {
                Mr mr = m.get_row(i,i,N);
                sum += SumElementsV_Helper<-1,UNKNOWN,comp,ix,Mr>::call(mr,x);
            }
            return sum;
        }
    };

    // algo 15: Fully unroll by rows
    template <int s, CompType comp, int ix, class M1>
    struct SumElementsU_Helper<15,s,comp,ix,M1>
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<MT>::type ret;

        template <int I, int M, int J, int N>
        struct Unroller
        {
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,M/2,J,N>::unroll(m,x) +
                    Unroller<I+M/2,M-M/2,J,N>::unroll(m,x));
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,1,J,N/2>::unroll(m,x) +
                    Unroller<I,1,J+N/2,N-N/2>::unroll(m,x));
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        { 
            static inline ret unroll(const M1& , const Scaling<ix,RT>& ) 
            { return ret(0); } 
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return Component<comp,MT>::f(x * m.cref(I,J)); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        { 
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return ret(0); }
        };
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        { return Unroller<0,s,0,s>::unroll(m,x); }
    };

    // algo 16: Fully unroll by columns
    template <int s, CompType comp, int ix, class M1>
    struct SumElementsU_Helper<16,s,comp,ix,M1>
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<MT>::type ret;

        template <int I, int M, int J, int N>
        struct Unroller
        {
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,M-(N-N/2),J,N/2>::unroll(m,x) +
                    Unroller<I,M,J+N/2,N-N/2>::unroll(m,x));
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,1>
        {
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,M/2,J,1>::unroll(m,x) +
                    Unroller<I+M/2,M-M/2,J,1>::unroll(m,x));
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,0>
        { 
            static inline ret unroll(const M1& , const Scaling<ix,RT>& ) 
            { return ret(0); } 
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return Component<comp,MT>::f(x * m.cref(I,J)); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        { 
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return ret(0); }
        };
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        { return Unroller<0,s,0,s>::unroll(m,x); }
    };

    // algo -4: No branches or copies
    // And m is NonUnitDiag
    template <int s, CompType comp, int ix, class M1>
    struct SumElementsU_Helper<-4,s,comp,ix,M1>
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<MT>::type ret;
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            TMVStaticAssert(M1::_upper);
            TMVStaticAssert(!M1::_unit);
            TMVAssert(!m.isunit());
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s == UNKNOWN ? false :
                // Norm is faster with the regular algorithm except for 
                // very small matrices.
                (s > 3 && comp == NormComp) ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const int algo = 
                unroll ? ( M1::_rowmajor ? 15 : 16 ) :
                M1::_colmajor ? 11 : 12;
            return SumElementsU_Helper<algo,s,comp,ix,M1>::call(m,x);
        }
    };

    // algo -1: Determine which algorithm to use
    template <int s, CompType comp, int ix, class M1>
    struct SumElementsU_Helper<-1,s,comp,ix,M1>
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        typedef typename Maybe<comp!=ValueComp>::
            template RealType<MT>::type ret;
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            TMVStaticAssert(M1::_upper);
            const int algo = 
                M1::_unit ? 1 :
                M1::_unknowndiag ? 2 :
                -4;
            return SumElementsU_Helper<algo,s,comp,ix,M1>::call(m,x);
        }
    };

    template <class M>
    inline typename M::value_type InlineSumElements(const BaseMatrix_Tri<M>& m)
    {
        TMVStaticAssert(M::_upper);
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return SumElementsU_Helper<-1,M::_size,ValueComp,1,Mv>::call(
            mv,Scaling<1,RT>());
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    T InstSumElements(
        const ConstUpperTriMatrixView<T>& m); 

    template <bool lower, bool conj, bool inst, class M>
    struct CallSumElementsU;

    template <bool conj, bool inst, class M>
    struct CallSumElementsU<true,conj,inst,M> // lower = true
    {
        static inline typename M::value_type call(const M& m)
        {
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return CallSumElementsU<false,conj,inst,Mt>::call(mt);
        }
    };
    template <bool inst, class M>
    struct CallSumElementsU<false,true,inst,M> // conj = true
    {
        static inline typename M::value_type call(const M& m)
        {
            typedef typename M::const_conjugate_type Mc;
            return TMV_CONJ(
                CallSumElementsU<false,false,inst,Mc>::call(m.conjugate()));
        }
    };
    template <class M>
    struct CallSumElementsU<false,false,true,M> // inst = true
    {
        static inline typename M::value_type call(const M& m)
        { return InstSumElements(m.xdView()); }
    };
    template <class M>
    struct CallSumElementsU<false,false,false,M> // inst = false
    {
        static inline typename M::value_type call(const M& m)
        { return InlineSumElements(m); }
    };

    template <class M>
    inline typename M::value_type SumElements(const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::value_type T;
        const bool inst = 
            M::unknownsizes &&
            (M::_rowmajor || M::_colmajor) &&
            Traits<T>::isinst;
        const bool lower = M::_lower;
        return CallSumElementsU<lower,M::_conj,inst,M>::call(m.mat());
    }


    // 
    // SumAbsElements
    //

    template <class M>
    inline typename M::real_type InlineSumAbsElements(
        const BaseMatrix_Tri<M>& m)
    {
        TMVStaticAssert(M::_upper);
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return SumElementsU_Helper<-1,M::_size,AbsComp,1,Mv>::call(
            mv,Scaling<1,RT>());
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename Traits<T>::real_type InstSumAbsElements(
        const ConstUpperTriMatrixView<T>& m); 

    template <bool lower, bool inst, class M>
    struct CallSumAbsElementsU;

    template <bool inst, class M>
    struct CallSumAbsElementsU<true,inst,M> // lower = true
    {
        static inline typename M::real_type call(const M& m)
        {
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return CallSumAbsElementsU<false,inst,Mt>::call(mt);
        }
    };
    template <class M>
    struct CallSumAbsElementsU<false,true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstSumAbsElements(m.xdView()); }
    };
    template <class M>
    struct CallSumAbsElementsU<false,false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineSumAbsElements(m); }
    };

    template <class M>
    inline typename M::real_type SumAbsElements(const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            (M::_rowmajor || M::_colmajor) &&
            Traits<T>::isinst;
        const bool lower = M::_lower;
        return CallSumAbsElementsU<lower,inst,Mn>::call(m.nonConj());
    }


    // 
    // SumAbs2Elements
    //

    template <class M>
    inline typename M::real_type InlineSumAbs2Elements(
        const BaseMatrix_Tri<M>& m)
    {
        TMVStaticAssert(M::_upper);
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return SumElementsU_Helper<-1,M::_size,Abs2Comp,1,Mv>::call(
            mv,Scaling<1,RT>());
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstUpperTriMatrixView<T>& m); 

    template <bool lower, bool inst, class M>
    struct CallSumAbs2ElementsU;

    template <bool inst, class M>
    struct CallSumAbs2ElementsU<true,inst,M> // lower = true
    {
        static inline typename M::real_type call(const M& m)
        {
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return CallSumAbs2ElementsU<false,inst,Mt>::call(mt);
        }
    };
    template <class M>
    struct CallSumAbs2ElementsU<false,true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstSumAbs2Elements(m.xdView()); }
    };
    template <class M>
    struct CallSumAbs2ElementsU<false,false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineSumAbs2Elements(m); }
    };

    template <class M>
    inline typename M::real_type SumAbs2Elements(const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            (M::_rowmajor || M::_colmajor) &&
            Traits<T>::isinst;
        const bool lower = M::_lower;
        return CallSumAbs2ElementsU<lower,inst,Mn>::call(m.nonConj());
    }


    // 
    // NormSq
    //

    template <class M>
    inline typename M::real_type InlineNormSq(const BaseMatrix_Tri<M>& m)
    {
        TMVStaticAssert(M::_upper);
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return SumElementsU_Helper<-1,M::_size,NormComp,1,Mv>::call(
            mv,Scaling<1,RT>());
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNormSq(
        const ConstUpperTriMatrixView<T>& m); 

    template <bool lower, bool inst, class M>
    struct CallNormSqU;

    template <bool inst, class M>
    struct CallNormSqU<true,inst,M> // lower = true
    {
        static inline typename M::real_type call(const M& m)
        {
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return CallNormSqU<false,inst,Mt>::call(mt);
        }
    };
    template <class M>
    struct CallNormSqU<false,true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstNormSq(m.xdView()); }
    };
    template <class M>
    struct CallNormSqU<false,false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineNormSq(m); }
    };

    template <class M>
    inline typename M::real_type NormSq(const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            (M::_rowmajor || M::_colmajor) &&
            Traits<T>::isinst;
        const bool lower = M::_lower;
        return CallNormSqU<lower,inst,Mn>::call(m.nonConj());
    }


    // 
    // NormSq with scaling
    //

    template <class M>
    inline typename M::real_type InlineNormSq(
        const BaseMatrix_Tri<M>& m, typename M::real_type scale)
    {
        TMVStaticAssert(M::_upper);
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return SumElementsU_Helper<-1,M::_size,NormComp,0,Mv>::call(
            mv,Scaling<0,RT>(scale));
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNormSq(
        const ConstUpperTriMatrixView<T>& m, 
        typename Traits<T>::real_type scale); 

    template <bool lower, bool inst, class M>
    struct CallNormSq_scaleU;

    template <bool inst, class M>
    struct CallNormSq_scaleU<true,inst,M> // lower = true
    {
        static inline typename M::real_type call(
            const M& m, const typename M::real_type scale)
        {
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return CallNormSq_scaleU<false,inst,Mt>::call(mt,scale);
        }
    };
    template <class M>
    struct CallNormSq_scaleU<false,true,M> // inst = true
    {
        static inline typename M::real_type call(
            const M& m, const typename M::real_type scale)
        { return InstNormSq(m.xdView(),scale); }
    };
    template <class M>
    struct CallNormSq_scaleU<false,false,M> // inst = false
    {
        static inline typename M::real_type call(
            const M& m, const typename M::real_type scale)
        { return InlineNormSq(m,scale); }
    };

    template <class M>
    inline typename M::real_type NormSq(
        const BaseMatrix_Tri<M>& m, const typename M::real_type scale)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            (M::_rowmajor || M::_colmajor) &&
            Traits<T>::isinst;
        const bool lower = M::_lower;
        return CallNormSq_scaleU<lower,inst,Mn>::call(m.nonConj(),scale);
    }



    //
    // NormF
    //

    // Norm_Helper in TMV_NormV.h works for UpperTriMatrix as well, so no
    // need to repeat that here.
    template <class M>
    inline typename M::real_type InlineNormF(const BaseMatrix_Tri<M>& m)
    {
        TMVStaticAssert(M::_upper);
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
#if TMV_OPT == 0
        return Norm_Helper<1,Mv>::call(mv);
#else
        if (m.size() == 0) return typename M::real_type(0);
        else return Norm_Helper<3,Mv>::call(mv);
#endif
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNormF(
        const ConstUpperTriMatrixView<T>& m); 

    template <bool lower, bool inst, class M>
    struct CallNormFU;

    template <bool inst, class M>
    struct CallNormFU<true,inst,M> // lower = true
    {
        static inline typename M::real_type call(const M& m)
        {
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return CallNormFU<false,inst,Mt>::call(mt);
        }
    };
    template <class M>
    struct CallNormFU<false,true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstNormF(m.xdView()); }
    };
    template <class M>
    struct CallNormFU<false,false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineNormF(m); }
    };

    template <class M>
    inline typename M::real_type NormF(const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            (M::_rowmajor || M::_colmajor) &&
            Traits<T>::isinst;
        const bool lower = M::_lower;
        return CallNormFU<lower,inst,Mn>::call(m.nonConj());
    }


    // 
    // MaxAbsElement
    //

    template <int algo, int s, class M1>
    struct MaxAbsElementU_Helper;

    // algo 1: m1 is unitdiag
    template <int s, class M1>
    struct MaxAbsElementU_Helper<1,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            typedef typename M1::const_offdiag_type Mo;
            Mo mo = m.offDiag();
            const int sm1 = IntTraits2<s,-1>::sum;
            const int algo2 = M1::_rowmajor ? 2 : 3;
            RT temp = MaxAbsElementU_Helper<algo2,sm1,Mo>::call(mo);
            return (temp > RT(1)) ? temp : RT(1);
        }
    };

    // algo 2: loop over rows
    template <int s, class M1>
    struct MaxAbsElementU_Helper<2,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int N = (s == UNKNOWN ? m.size() : s);
            RT max(0);
            for(int i=0;i<N;++i) {
                RT temp = InlineMaxAbsElement(m.get_row(i,i,N),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 3: loop over columns
    template <int s, class M1>
    struct MaxAbsElementU_Helper<3,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int N = (s == UNKNOWN ? m.size() : s);
            RT max(0);
            for(int j=0;j<N;++j) {
                RT temp = InlineMaxAbsElement(m.get_col(j,0,j+1),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 4: loop over rows with temp storage
    template <int s, class M1>
    struct MaxAbsElementU_Helper<4,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int N = (s == UNKNOWN ? m.size() : s);
            tmv::Vector<RT> temp(N);
            for(int i=0;i<N;++i) {
                temp(i) = InlineMaxAbsElement(m.get_row(i,i,N),0);
            }
            return temp.maxAbsElement();
        }
    };

    // algo 5: loop over columns with temp storage
    template <int s, class M1>
    struct MaxAbsElementU_Helper<5,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int N = (s == UNKNOWN ? m.size() : s);
            tmv::Vector<RT> temp(N);
            for(int j=0;j<N;++j) {
                temp(j) = InlineMaxAbsElement(m.get_col(j,0,j+1),0);
            }
            return temp.maxAbsElement();
        }
    };

    // algo 90: unknown diag, figure out which it is.
    template <int s, class M1>
    struct MaxAbsElementU_Helper<90,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int algo2 = M1::_rowmajor ? 4 : 5;
            if (m.isunit()) 
                return MaxAbsElementU_Helper<1,s,M1>::call(m);
            else 
                return MaxAbsElementU_Helper<algo2,s,M1>::call(m);
        }
    };

    template <class M>
    inline typename M::real_type InlineMaxAbsElement(
        const BaseMatrix_Tri<M>& m)
    {
        TMVStaticAssert(M::_upper);
        const int s = M::_size;
        const int algo = 
            M::_unit ? 1 : 
            M::_unknowndiag ? 90 :
            M::_rowmajor ? 4 : 
            5;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return MaxAbsElementU_Helper<algo,s,M>::call(mv);
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename Traits<T>::real_type InstMaxAbsElement(
        const ConstUpperTriMatrixView<T>& m); 

    template <bool lower, bool inst, class M>
    struct CallMaxAbsElementU;

    template <bool inst, class M>
    struct CallMaxAbsElementU<true,inst,M> // lower = true
    {
        static inline typename M::real_type call(const M& m)
        {
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return CallMaxAbsElementU<false,inst,Mt>::call(mt);
        }
    };
    template <class M>
    struct CallMaxAbsElementU<false,true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstMaxAbsElement(m.xdView()); }
    };
    template <class M>
    struct CallMaxAbsElementU<false,false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineMaxAbsElement(m); }
    };

    template <class M>
    inline typename M::real_type MaxAbsElement(const BaseMatrix_Tri<M>& m)
    {
        TMVAssert(m.size() > 0);
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            (M::_rowmajor || M::_colmajor) &&
            Traits<T>::isinst;
        const bool lower = M::_lower;
        return CallMaxAbsElementU<lower,inst,Mn>::call(m.nonConj());
    }


    // 
    // Norm1
    //

    // TODO: Norm1 and NormInf would benefit from an unroller.
    // Unlike with a regular matrix, where small matrices can unroll
    // each column (or row), with a triangle matrix, we lose the knowledge
    // of the length of each column in the for loop, so a full unroller
    // would be able to keep that.
    template <int algo, int s, class M1>
    struct Norm1U_Helper;

    // algo 1: loop over columns
    template <int s, class M1>
    struct Norm1U_Helper<1,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int N = (s == UNKNOWN ? m.size() : s);
            RT max(0);
            for(int j=0;j<N;++j) {
                // If unit,     temp = 1 + SumAbsElements(m.col(j,0,j)
                // If non-unit, temp =     SumAbsElements(m.col(j,0,j+1)
                RT temp = Maybe<M1::_unit>::sum( 
                    RT(1) , InlineSumAbsElements(
                        m.get_col(j,0,Maybe<M1::_unit>::select(j,j+1))));
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 2: loop over rows
    template <int s, class M1>
    struct Norm1U_Helper<2,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            int N = (s == UNKNOWN ? m.size() : s);
            if (N <= 8) return Norm1U_Helper<1,s,M1>::call(m);

            typedef typename M1::value_type MT;
            typedef typename tmv::Vector<RT>::iterator IT1;
            typedef typename M1::const_row_sub_type::const_iterator IT2;

            MT value;

            // If unit,     start with all 1's.
            // If non-unit, start with all 0's.
            tmv::Vector<RT> temp(N,
                                 Maybe<M1::_unit>::select(RT(1),RT(0)) );
            const IT1 begin1 = temp.begin();
            const IT1 end1 = temp.end();

            IT1 it1 = begin1;
            IT2 it2 = m.get_row(0,Maybe<M1::_unit>::select(1,0),N).begin();
            // If unit, then we only need to add the offdiag to temp.
            // This is effected by: --N, ++it1, and the above select for it2.
            Maybe<M1::_unit>::decrement(N);
            Maybe<M1::_unit>::increment(it1);
            int end_step = m.diagstep() - N; // back to the start of next row

            do {
                do {
                    value = *it2++;
                    Component<AbsComp,MT>::applyf(value);
                    *it1++ += TMV_REAL(value);
                } while (it1 != end1);
                it2 += end_step++;
                it1 -= (--N);
            } while (N);
            return InlineMaxElement(temp,0);
        }
    };

    template <class M>
    inline typename M::real_type InlineNorm1(const BaseMatrix_Tri<M>& m);

    // algo 90: unknown diag, figure out which it is.
    template <int s, class M1>
    struct Norm1U_Helper<90,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            if (m.isunit()) 
                return InlineNorm1(m.viewAsUnitDiag());
            else 
                return InlineNorm1(m.viewAsNonUnitDiag());
        }
    };

    template <class M>
    inline typename M::real_type InlineNorm1(const BaseMatrix_Tri<M>& m)
    {
        TMVStaticAssert(M::_upper);
        const int s = M::_size;
        const int algo = 
            M::_unknowndiag ? 90 :
#if TMV_OPT >= 1
            ( M::_rowmajor && (s == UNKNOWN || s > 8) ) ? 2 :
#endif
            1;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return Norm1U_Helper<algo,s,Mv>::call(mv);
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNorm1(
        const ConstUpperTriMatrixView<T>& m); 

    template <bool lower, bool inst, class M>
    struct CallNorm1U;
    template <bool lower, bool inst, class M>
    struct CallNormInfU;

    template <bool inst, class M>
    struct CallNorm1U<true,inst,M> // lower = true
    {
        static inline typename M::real_type call(const M& m)
        {
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return CallNormInfU<false,inst,Mt>::call(mt);
        }
    };
    template <class M>
    struct CallNorm1U<false,true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstNorm1(m.xdView()); }
    };
    template <class M>
    struct CallNorm1U<false,false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineNorm1(m); }
    };

    template <class M>
    inline typename M::real_type Norm1(const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            (M::_rowmajor || M::_colmajor) &&
            Traits<T>::isinst;
        const bool lower = M::_lower;
        return CallNorm1U<lower,inst,Mn>::call(m.nonConj());
    }

    // 
    // NormInf
    //

    template <int algo, int s, class M1>
    struct NormInfU_Helper;

    // algo 1: loop over rows
    template <int s, class M1>
    struct NormInfU_Helper<1,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int N = (s == UNKNOWN ? m.size() : s);
            RT max(0);
            for(int i=0;i<N;++i) {
                // If unit,     temp = 1 + SumAbsElements(m.row(i,0,i)
                // If non-unit, temp =     SumAbsElements(m.row(i,0,i+1)
                RT temp = Maybe<M1::_unit>::sum( 
                    RT(1) , InlineSumAbsElements(
                        m.get_row(i,Maybe<M1::_unit>::select(i+1,i),N)));
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 2: loop over cols
    template <int s, class M1>
    struct NormInfU_Helper<2,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            int N = (s == UNKNOWN ? m.size() : s);
            if (N <= 8) return NormInfU_Helper<1,s,M1>::call(m);

            typedef typename M1::value_type MT;
            typedef typename tmv::Vector<RT>::iterator IT1;
            typedef typename M1::const_col_sub_type::const_iterator IT2;

            MT value;

            // If unit,     start with all 1's.
            // If non-unit, start with all 0's.
            tmv::Vector<RT> temp(N, Maybe<M1::_unit>::select(RT(1),RT(0)) );
            const IT1 begin1 = temp.begin();

            IT1 it1 = begin1;
            IT2 it2 = m.get_col(Maybe<M1::_unit>::select(1,0),0,1).begin();
            // If unit, then we only need to add the offdiag to temp.
            // This is effected by: --N, and the above select for it2.
            Maybe<M1::_unit>::decrement(N);
            int end_step = m.stepj()-1; // back to the start of next column

            int M=1, i;
            do {
                i=M; do {
                    value = *it2++;
                    Component<AbsComp,MT>::applyf(value);
                    *it1++ += TMV_REAL(value);
                } while (--i);
                it1 -= M++;
                it2 += end_step--;
            } while (--N);
            return InlineMaxElement(temp,0);
        }
    };

    template <class M>
    inline typename M::real_type InlineNormInf(const BaseMatrix_Tri<M>& m);

    // algo 90: unknown diag, figure out which it is.
    template <int s, class M1>
    struct NormInfU_Helper<90,s,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            if (m.isunit()) 
                return InlineNormInf(m.viewAsUnitDiag());
            else 
                return InlineNormInf(m.viewAsNonUnitDiag());
        }
    };

    template <class M>
    inline typename M::real_type InlineNormInf(const BaseMatrix_Tri<M>& m)
    {
        TMVStaticAssert(M::_upper);
        const int s = M::_size;
        const int algo = 
            M::_unknowndiag ? 90 :
#if TMV_OPT >= 1
            ( M::_colmajor && (s == UNKNOWN || s > 8) ) ? 2 :
#endif
            1;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return NormInfU_Helper<algo,s,Mv>::call(mv);
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNormInf(
        const ConstUpperTriMatrixView<T>& m); 

    template <bool inst, class M>
    struct CallNormInfU<true,inst,M> // lower = true
    {
        static inline typename M::real_type call(const M& m)
        {
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return CallNorm1U<false,inst,Mt>::call(mt);
        }
    };
    template <class M>
    struct CallNormInfU<false,true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstNormInf(m.xdView()); }
    };
    template <class M>
    struct CallNormInfU<false,false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineNormInf(m); }
    };

    template <class M>
    inline typename M::real_type NormInf(const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst =
            M::unknownsizes &&
            (M::_rowmajor || M::_colmajor) &&
            Traits<T>::isinst;
        const bool lower = M::_lower;
        return CallNormInfU<lower,inst,Mn>::call(m.nonConj());
    }

#undef TMV_Q1

} // namespace tmv

#endif
