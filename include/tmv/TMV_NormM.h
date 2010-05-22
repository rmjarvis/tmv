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

#ifndef TMV_NormM_H
#define TMV_NormM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_NormV.h"
#include "TMV_MinMax.h"

namespace tmv {

    // TODO: Convert to new algo numbering scheme

    // 
    // SumElements
    //

    template <int algo, int cs, int rs, CompType comp, int ix, class M1>
    struct SumElementsM_Helper;

    // algo 1: linearize to vector version
    template <int cs, int rs, CompType comp, int ix, class M1>
    struct SumElementsM_Helper<1,cs,rs,comp,ix,M1> 
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        enum { rt = comp != ValueComp };
        typedef typename Maybe<rt>::template RealType<MT>::type ret;

        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::const_linearview_type Ml;
            Ml ml = m.linearView();
            return SumElementsV_Helper<-1,Ml::_size,comp,ix,Ml>::call(ml,x);
        }
    };

    // algo 2: loop over rows
    template <int cs, int rs, CompType comp, int ix, class M1>
    struct SumElementsM_Helper<2,cs,rs,comp,ix,M1> 
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        enum { rt = comp != ValueComp };
        typedef typename Maybe<rt>::template RealType<MT>::type ret;

        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int M = (cs == UNKNOWN ? int(m.colsize()) : cs);
            typedef typename M1::const_row_type Mr;
            ret sum(0);
            for(int i=0;i<M;++i) {
                Mr mr = m.get_row(i);
                sum += SumElementsV_Helper<-1,rs,comp,ix,Mr>::call(mr,x);
            }
            return sum;
        }
    };

    // algo 3: loop over columns
    template <int cs, int rs, CompType comp, int ix, class M1>
    struct SumElementsM_Helper<3,cs,rs,comp,ix,M1> 
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        enum { rt = comp != ValueComp };
        typedef typename Maybe<rt>::template RealType<MT>::type ret;

        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int N = (rs == UNKNOWN ? int(m.rowsize()) : rs);
            typedef typename M1::const_col_type Mc;
            ret sum(0);
            for(int j=0;j<N;++j) {
                Mc mc = m.get_col(j);
                sum += SumElementsV_Helper<-1,cs,comp,ix,Mc>::call(mc,x);
            }
            return sum;
        }
    };

    // algo 4: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, CompType comp, int ix, class M1>
    struct SumElementsM_Helper<4,cs,rs,comp,ix,M1>
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        enum { rt = comp != ValueComp };
        typedef typename Maybe<rt>::template RealType<MT>::type ret;

        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#if TMV_OPT >= 2
            if (m.canLinearize())
                return SumElementsM_Helper<1,cs,rs,comp,ix,M1>::call(m,x);
            else if ( m.isrm() ||
                      (!m.iscm() && (m.colsize() < m.rowsize()) ) )
                return SumElementsM_Helper<2,cs,rs,comp,ix,M1>::call(m,x);
            else
                return SumElementsM_Helper<3,cs,rs,comp,ix,M1>::call(m,x);
#else
            const int algo2 = 
#if TMV_OPT >= 1
                M1::_rowmajor ? 2 :
                M1::_colmajor ? 3 :
                ( cs == UNKNOWN || rs == UNKNOWN ) ? 3 :
                ( cs < rs ) ? 2 :
#endif
                3;
            return SumElementsM_Helper<algo2,cs,rs,comp,ix,M1>::call(m,x);
#endif
        }
    };

    // algo 5: Fully unroll by rows
    template <int cs, int rs, CompType comp, int ix, class M1>
    struct SumElementsM_Helper<5,cs,rs,comp,ix,M1>
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        enum { rt = comp != ValueComp };
        typedef typename Maybe<rt>::template RealType<MT>::type ret;

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
        { return Unroller<0,cs,0,rs>::unroll(m,x); }
    };

    // algo 6: Fully unroll by columns
    template <int cs, int rs, CompType comp, int ix, class M1>
    struct SumElementsM_Helper<6,cs,rs,comp,ix,M1>
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        enum { rt = comp != ValueComp };
        typedef typename Maybe<rt>::template RealType<MT>::type ret;

        template <int I, int M, int J, int N>
        struct Unroller
        {
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,M,J,N/2>::unroll(m,x) +
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
        { return Unroller<0,cs,0,rs>::unroll(m,x); }
    };

    // algo -1: Determine which algorithm to use
    template <int cs, int rs, CompType comp, int ix, class M1>
    struct SumElementsM_Helper<-1,cs,rs,comp,ix,M1>
    {
        typedef typename M1::value_type MT;
        typedef typename M1::real_type RT;
        enum { rt = comp != ValueComp };
        typedef typename Maybe<rt>::template RealType<MT>::type ret;

        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int algo = 
#if TMV_OPT >= 1
                M1::_canlin ? 1 :
                ( cs != UNKNOWN && rs != UNKNOWN ) ? (
                    ( IntTraits2<cs,rs>::prod <= int(128/sizeof(MT)) ) ? (
                        M1::_rowmajor ? 5 : 6 ) :
                    M1::_rowmajor ? 2 : 
                    M1::_colmajor ? 3 :
                    cs > rs ? 2 : 3 ) :
#endif
                4;
            return SumElementsM_Helper<algo,cs,rs,comp,ix,M1>::call(m,x);
        }
    };

    template <class M>
    inline typename M::value_type InlineSumElements(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        return SumElementsM_Helper<-1,cs,rs,ValueComp,1,Mv>::call(
            mv,Scaling<1,RT>());
    }

    // Defined in TMV_Matrix.cpp
    template <class T>
    T InstSumElements(const ConstMatrixView<T>& m); 

    template <class T>
    inline T InstSumElements(const ConstMatrixView<T,UNKNOWN,1>& m)
    { return InstSumElements(m.transpose()); }

    template <bool conj, bool inst, class M>
    struct CallSumElementsm;

    template <bool inst, class M>
    struct CallSumElementsm<true,inst,M> // conj = true
    {
        static inline typename M::value_type call(const M& m)
        {
            typedef typename M::const_conjugate_type Mc;
            return TMV_CONJ(CallSumElementsm<false,inst,Mc>::call(m.conjugate()));
        }
    };
    template <class M>
    struct CallSumElementsm<false,false,M> // inst = false
    {
        static inline typename M::value_type call(const M& m)
        { return InlineSumElements(m); }
    };
    template <class M>
    struct CallSumElementsm<false,true,M> // inst = true
    {
        static inline typename M::value_type call(const M& m)
        { return InstSumElements(m.xView()); }
    };

    template <class M>
    inline typename M::value_type SumElements(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type T;
        const bool inst = 
            M::unknownsizes &&
            Traits<T>::isinst;
        return CallSumElementsm<M::_conj,inst,M>::call(m.mat());
    }


    // 
    // SumAbsElements
    //

    template <class M>
    inline typename M::real_type InlineSumAbsElements(
        const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return SumElementsM_Helper<
            -1,M::_colsize,M::_rowsize,AbsComp,1,Mv>::call(
                mv,Scaling<1,RT>());
    }

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename Traits<T>::real_type InstSumAbsElements(
        const ConstMatrixView<T>& m); 

    template <class T>
    inline typename Traits<T>::real_type InstSumAbsElements(
        const ConstMatrixView<T,UNKNOWN,1>& m)
    { return InstSumAbsElements(m.transpose()); }

    template <bool inst, class M>
    struct CallSumAbsElementsm;

    template <class M>
    struct CallSumAbsElementsm<true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstSumAbsElements(m.xView()); }
    };
    template <class M>
    struct CallSumAbsElementsm<false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineSumAbsElements(m); }
    };

    template <class M>
    inline typename M::real_type SumAbsElements(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            Traits<T>::isinst;
        return CallSumAbsElementsm<inst,Mn>::call(m.nonConj());
    }


    // 
    // SumAbs2Elements
    //

    template <class M>
    inline typename M::real_type InlineSumAbs2Elements(
        const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return SumElementsM_Helper<
            -1,M::_colsize,M::_rowsize,Abs2Comp,1,Mv>::call(
                mv,Scaling<1,RT>());
    }

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstMatrixView<T>& m); 

    template <class T>
    inline typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstMatrixView<T,UNKNOWN,1>& m)
    { return InstSumAbs2Elements(m.transpose()); }

    template <bool inst, class M>
    struct CallSumAbs2Elementsm;

    template <class M>
    struct CallSumAbs2Elementsm<true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstSumAbs2Elements(m.xView()); }
    };
    template <class M>
    struct CallSumAbs2Elementsm<false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineSumAbs2Elements(m); }
    };

    template <class M>
    inline typename M::real_type SumAbs2Elements(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            Traits<T>::isinst;
        return CallSumAbs2Elementsm<inst,Mn>::call(m.nonConj());
    }


    // 
    // NormSq
    //

    template <class M>
    inline typename M::real_type InlineNormSq(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return SumElementsM_Helper<
            -1,M::_colsize,M::_rowsize,NormComp,1,Mv>::call(
                mv,Scaling<1,RT>());
    }

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNormSq(const ConstMatrixView<T>& m); 

    template <class T>
    inline typename Traits<T>::real_type InstNormSq(
        const ConstMatrixView<T,UNKNOWN,1>& m)
    { return InstNormSq(m.transpose()); }

    template <bool inst, class M>
    struct CallNormSqm;

    template <class M>
    struct CallNormSqm<true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstNormSq(m.xView()); }
    };
    template <class M>
    struct CallNormSqm<false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineNormSq(m); }
    };

    template <class M>
    inline typename M::real_type NormSq(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            Traits<T>::isinst;
        return CallNormSqm<inst,Mn>::call(m.nonConj());
    }


    // 
    // NormSq with scaling
    //

    template <class M>
    inline typename M::real_type InlineNormSq(
        const BaseMatrix_Rec<M>& m, typename M::real_type scale)
    {
        typedef typename M::value_type MT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return SumElementsM_Helper<
            -1,M::_colsize,M::_rowsize,NormComp,0,Mv>::call(
                mv,Scaling<0,RT>(scale));
    }

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNormSq(
        const ConstMatrixView<T>& m, typename Traits<T>::real_type scale); 

    template <class T>
    inline typename Traits<T>::real_type InstNormSq(
        const ConstMatrixView<T,UNKNOWN,1>& m,
        typename Traits<T>::real_type scale)
    { return InstNormSq(m.transpose(),scale); }

    template <bool inst, class M>
    struct CallNormSq_scalem;

    template <class M>
    struct CallNormSq_scalem<true,M> // inst = true
    {
        static inline typename M::real_type call(
            const M& m, const typename M::real_type scale)
        { return InstNormSq(m.xView(),scale); }
    };
    template <class M>
    struct CallNormSq_scalem<false,M> // inst = false
    {
        static inline typename M::real_type call(
            const M& m, const typename M::real_type scale)
        { return InlineNormSq(m,scale); }
    };

    template <class M>
    inline typename M::real_type NormSq(
        const BaseMatrix_Rec<M>& m, const typename M::real_type scale)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            Traits<T>::isinst;
        return CallNormSq_scalem<inst,Mn>::call(m.nonConj(),scale);
    }



    //
    // NormF
    //

    // Norm_Helper in TMV_NormV.h works for Matrix as well, so no
    // need to repeat that here.
    template <class M>
    inline typename M::real_type InlineNormF(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
#if TMV_OPT == 0
        return Norm_Helper<1,Mv>::call(mv);
#else
        typedef typename M::real_type RT;
        const bool isint = std::numeric_limits<RT>::is_integer;
        const int algo = isint ? 1 : 3;
        if (m.colsize() == 0 || m.rowsize() == 0) 
            return typename M::real_type(0);
        else 
            return Norm_Helper<algo,Mv>::call(mv);
#endif
    }

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNormF(const ConstMatrixView<T>& m); 

    template <class T>
    inline typename Traits<T>::real_type InstNormF(
        const ConstMatrixView<T,UNKNOWN,1>& m)
    { return InstNormF(m.transpose()); }

    template <bool inst, class M>
    struct CallNormFm;

    template <class M>
    struct CallNormFm<true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstNormF(m.xView()); }
    };
    template <class M>
    struct CallNormFm<false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineNormF(m); }
    };

    template <class M>
    inline typename M::real_type NormF(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            Traits<T>::isinst;
        return CallNormFm<inst,Mn>::call(m.nonConj());
    }


    // 
    // MaxAbsElement
    //

    template <int algo, int cs, int rs, class M1>
    struct MaxAbsElementM_Helper;

    // algo 1: linearize to vector version
    template <int cs, int rs, class M1>
    struct MaxAbsElementM_Helper<1,cs,rs,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        { return InlineMaxAbsElement(m.linearView(),0); }
    };

    // algo 2: loop over rows
    template <int cs, int rs, class M1>
    struct MaxAbsElementM_Helper<2,cs,rs,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int M = (cs == UNKNOWN ? int(m.colsize()) : cs);
            RT max(0);
            for(int i=0;i<M;++i) {
                RT temp = InlineMaxAbsElement(m.get_row(i),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 3: loop over columns
    template <int cs, int rs, class M1>
    struct MaxAbsElementM_Helper<3,cs,rs,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int N = (rs == UNKNOWN ? int(m.rowsize()) : rs);
            RT max(0);
            for(int j=0;j<N;++j) {
                RT temp = InlineMaxAbsElement(m.get_col(j),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 4: loop over rows with temp storage
    template <int cs, int rs, class M1>
    struct MaxAbsElementM_Helper<4,cs,rs,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int M = (cs == UNKNOWN ? int(m.colsize()) : cs);
            Vector<RT> temp(M);
            for(int i=0;i<M;++i) {
                temp(i) = InlineMaxAbsElement(m.get_row(i),0);
            }
            return temp.maxAbsElement();
        }
    };

    // algo 5: loop over columns
    template <int cs, int rs, class M1>
    struct MaxAbsElementM_Helper<5,cs,rs,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int N = (rs == UNKNOWN ? int(m.rowsize()) : rs);
            Vector<RT> temp(N);
            for(int j=0;j<N;++j) {
                temp(j) = InlineMaxAbsElement(m.get_col(j),0);
            }
            return temp.maxAbsElement();
        }
    };

    // algo 10: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, class M1>
    struct MaxAbsElementM_Helper<10,cs,rs,M1>
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            if (m.canLinearize())
                return MaxAbsElementM_Helper<1,cs,rs,M1>::call(m);
            else if ( m.isrm() ||
                      (!m.iscm() && (m.colsize() < m.rowsize()) ) )
                return MaxAbsElementM_Helper<4,cs,rs,M1>::call(m);
            else
                return MaxAbsElementM_Helper<5,cs,rs,M1>::call(m);
        }
    };

    // TODO: I don't have a full unroller here.
    // Usually, this will be run on a matrix that can be linearized, so
    // it will unroll there.  But for the rare case that we have a small matrix
    // that isn't the full SmallMatrix (and hence can't linearize), then
    // there might be a small speed increase available by unrolling.

    template <class M>
    inline typename M::real_type InlineMaxAbsElement(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
#if TMV_OPT == 0
        const int algo = 
            M::_canlin ? 1 : M::_rowmajor ? 2 : 3;
#else
        const int algo = 
            M::_canlin ? 1 :
#if TMV_OPT >= 2
            ( cs == UNKNOWN || rs == UNKNOWN ) ? 10 :
#endif
            M::_rowmajor ? 4 :
            M::_colmajor ? 5 :
            cs < rs ? 4 : 5;
#endif
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return MaxAbsElementM_Helper<algo,cs,rs,Mv>::call(mv);
    }

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename Traits<T>::real_type InstMaxAbsElement(
        const ConstMatrixView<T>& m);

    template <class T>
    inline typename Traits<T>::real_type InstMaxAbsElement(
        const ConstMatrixView<T,UNKNOWN,1>& m)
    { return InstMaxAbsElement(m.transpose()); }

    template <bool inst, class M>
    struct CallMaxAbsElementm;

    template <class M>
    struct CallMaxAbsElementm<true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstMaxAbsElement(m.xView()); }
    };
    template <class M>
    struct CallMaxAbsElementm<false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineMaxAbsElement(m); }
    };

    template <class M>
    inline typename M::real_type MaxAbsElement(const BaseMatrix_Rec<M>& m)
    {
        TMVAssert(m.colsize() > 0);
        TMVAssert(m.rowsize() > 0);
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            Traits<T>::isinst;
        return CallMaxAbsElementm<inst,Mn>::call(m.nonConj());
    }


    // 
    // Norm1
    //

    template <int algo, int cs, int rs, class M1>
    struct Norm1M_Helper;

    // algo 1: loop over columns
    template <int cs, int rs, class M1>
    struct Norm1M_Helper<1,cs,rs,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            const int N = (rs == UNKNOWN ? int(m.rowsize()) : rs);
            RT max(0);
            for(int j=0;j<N;++j) {
                RT temp = InlineSumAbsElements(m.get_col(j));
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 2: loop over rows
    template <int cs, int rs, class M1>
    struct Norm1M_Helper<2,cs,rs,M1> 
    {
        typedef typename M1::real_type RT;
        static inline RT call(const M1& m)
        {
            int M = (cs == UNKNOWN ? int(m.colsize()) : cs);
            if (M <= 8) return Norm1M_Helper<1,cs,rs,M1>::call(m);

            const int N = (rs == UNKNOWN ? int(m.rowsize()) : rs);
            if (N == 0) return RT(0);

            typedef typename M1::value_type MT;
            typedef typename tmv::Vector<RT>::iterator IT1;
            typedef typename M1::const_row_type::const_iterator IT2;

            MT value;

            tmv::Vector<RT> temp(N,RT(0));
            const IT1 begin1 = temp.begin();
            const IT1 end1 = temp.end();

            IT1 it1 = begin1;
            IT2 it2 = m.get_row(0).begin();
            const int end_step = m.stepi()-N; // back to the start of next row

            do {
                do {
                    value = *it2++;
                    Component<AbsComp,MT>::applyf(value);
                    *it1++ += TMV_REAL(value);
                } while (it1 != end1);
                it1 -= N;
                it2 += end_step;
            } while (--M);
            return InlineMaxElement(temp,0);
        }
    };

    template <class M>
    inline typename M::real_type InlineNorm1(const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        const int algo = 
#if TMV_OPT >= 1
            ( M::_rowmajor && 
              ( cs == UNKNOWN || rs == UNKNOWN || 
                (cs > 16 && rs < 512) || (cs > 8 && rs >= 512) ) ) ? 2 :
#endif
            1;
        typedef typename M::const_cview_type Mv;
        Mv mv = m.cView();
        return Norm1M_Helper<algo,cs,rs,Mv>::call(mv);
    }

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNorm1(const ConstMatrixView<T>& m); 

    template <class T>
    inline typename Traits<T>::real_type InstNorm1(
        const ConstMatrixView<T,UNKNOWN,1>& m)
    { return InstNormInf(m.transpose()); }

    template <bool inst, class M>
    struct CallNorm1M;

    template <class M>
    struct CallNorm1M<true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstNorm1(m.xView()); }
    };
    template <class M>
    struct CallNorm1M<false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineNorm1(m); }
    };

    template <class M>
    inline typename M::real_type Norm1(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            Traits<T>::isinst;
        return CallNorm1M<inst,Mn>::call(m.nonConj());
    }

    // 
    // NormInf
    //

    template <class M>
    inline typename M::real_type InlineNormInf(const BaseMatrix_Rec<M>& m)
    { return InlineNorm1(m.transpose()); }

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename Traits<T>::real_type InstNormInf(const ConstMatrixView<T>& m); 

    template <class T>
    inline typename Traits<T>::real_type InstNormInf(
        const ConstMatrixView<T,UNKNOWN,1>& m)
    { return InstNorm1(m.transpose()); }

    template <bool inst, class M>
    struct CallNormInfM;

    template <class M>
    struct CallNormInfM<true,M> // inst = true
    {
        static inline typename M::real_type call(const M& m)
        { return InstNormInf(m.xView()); }
    };
    template <class M>
    struct CallNormInfM<false,M> // inst = false
    {
        static inline typename M::real_type call(const M& m)
        { return InlineNormInf(m); }
    };

    template <class M>
    inline typename M::real_type NormInf(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type T;
        typedef typename M::const_nonconj_type Mn;
        const bool inst = 
            M::unknownsizes &&
            Traits<T>::isinst;
        return CallNormInfM<inst,Mn>::call(m.nonConj());
    }

} // namespace tmv

#endif
