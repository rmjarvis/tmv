
#ifndef TMV_NormM_H
#define TMV_NormM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_NormV.h"
#include "TMV_MinMax.h"

#ifdef PRINTALGO_NormM
#include <iostream>
#include "TMV_MatrixIO.h"
#endif

namespace tmv {

    // 
    // SumElements
    // SumAbsElements
    // SumAbs2Elements
    // NormSq
    //

    // Defined in TMV_Matrix.cpp
    template <class T>
    T InstSumElements(const ConstMatrixView<T>& m); 
    template <class T>
    typename ConstMatrixView<T>::float_type InstSumAbsElements(
        const ConstMatrixView<T>& m); 
    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstMatrixView<T>& m); 
    template <class T>
    typename Traits<T>::real_type InstNormSq(
        const ConstMatrixView<T>& m); 
    template <class T>
    typename ConstMatrixView<T>::float_type InstNormSq(
        const ConstMatrixView<T>& m, 
        typename ConstMatrixView<T>::float_type scale); 

    // The maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_NORMM_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_NORMM_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_NORMM_UNROLL 9
#else
#define TMV_NORMM_UNROLL 0
#endif

    template <int algo, int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper;

    // algo 1: linearize to vector version
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<1,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"SumElementsM algo 1: "<<TMV_Text(comp)<<std::endl;
#endif
            typedef typename M1::const_linearview_type Ml;
            Ml ml = m.linearView();
            return SumElementsV_Helper<-3,Ml::_size,comp,ix,ret,Ml>::call(ml,x);
        }
    };

    // algo 201: linearize to vector version
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<201,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"SumElementsM algo 201: "<<TMV_Text(comp)<<std::endl;
#endif
            typedef typename M1::const_linearview_type Ml;
            Ml ml = m.linearView();
            return SumElementsV_Helper<-2,Ml::_size,comp,ix,ret,Ml>::call(ml,x);
        }
    };

    // algo 11: loop over columns
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<11,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"SumElementsM algo 11: "<<TMV_Text(comp)<<std::endl;
#endif
            const int N = (rs == TMV_UNKNOWN ? m.rowsize() : rs);
            typedef typename M1::const_col_type Mc;
            ret sum(0);
            for(int j=0;j<N;++j) {
                Mc mc = m.get_col(j);
                sum += SumElementsV_Helper<-3,cs,comp,ix,ret,Mc>::call(mc,x);
            }
            return sum;
        }
    };

    // algo 12: loop over rows
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<12,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"SumElementsM algo 12: "<<TMV_Text(comp)<<std::endl;
#endif
            const int M = (cs == TMV_UNKNOWN ? m.colsize() : cs);
            typedef typename M1::const_row_type Mr;
            ret sum(0);
            for(int i=0;i<M;++i) {
                Mr mr = m.get_row(i);
                sum += SumElementsV_Helper<-3,rs,comp,ix,ret,Mr>::call(mr,x);
            }
            return sum;
        }
    };

    // algo 15: Fully unroll by columns
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<15,cs,rs,comp,ix,ret,M1>
    {
        typedef typename M1::value_type VT;
        typedef typename Traits<ret>::real_type RT;

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
            { return Component<comp,VT>::f(x * m.cref(I,J)); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        {
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return ret(0); }
        };
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"SumElementsM algo 15: "<<TMV_Text(comp)<<std::endl;
#endif
            return Unroller<0,cs,0,rs>::unroll(m,x); 
        }
    };

    // algo 16: Fully unroll by rows
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<16,cs,rs,comp,ix,ret,M1>
    {
        typedef typename M1::value_type VT;
        typedef typename Traits<ret>::real_type RT;
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
            { return Component<comp,VT>::f(x * m.cref(I,J)); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        {
            static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return ret(0); }
        };
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"SumElementsM algo 16: "<<TMV_Text(comp)<<std::endl;
#endif
            return Unroller<0,cs,0,rs>::unroll(m,x); 
        }
    };

    // algo 30: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<30,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"SumElementsM algo 30: "<<TMV_Text(comp)<<std::endl;
#endif
            if (m.canLinearize())
                return SumElementsM_Helper<1,cs,rs,comp,ix,ret,M1>::call(m,x);
            else 
                return SumElementsM_Helper<-4,cs,rs,comp,ix,ret,M1>::call(m,x);
        }
    };
    
    // algo 90: Call inst
    template <int cs, int rs, int ix, class ret, class M1>
    struct SumElementsM_Helper<90,cs,rs,ValueComp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return InstSumElements(m.xView()); }
    };
    template <int cs, int rs, int ix, class ret, class M1>
    struct SumElementsM_Helper<90,cs,rs,AbsComp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return InstSumAbsElements(m.xView()); }
    };
    template <int cs, int rs, int ix, class ret, class M1>
    struct SumElementsM_Helper<90,cs,rs,Abs2Comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return InstSumAbs2Elements(m.xView()); }
    };
    template <int cs, int rs, class ret, class M1>
    struct SumElementsM_Helper<90,cs,rs,NormComp,1,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<1,RT>& x)
        { return InstNormSq(m.xView()); }
    };
    template <int cs, int rs, class ret, class M1>
    struct SumElementsM_Helper<90,cs,rs,NormComp,0,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<0,RT>& x)
        { return InstNormSq(m.xView(),x); }
    };
    
    // algo 97: Conjugate
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<97,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return SumElementsM_Helper<-2,cs,rs,comp,ix,ret,Mnc>::call(mnc,x);
        }
    };
    template <int cs, int rs, int ix, class ret, class M1>
    struct SumElementsM_Helper<97,cs,rs,ValueComp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::const_conjugate_type Mc;
            Mc mc = m.conjugate();
            return TMV_CONJ(
                SumElementsM_Helper<-2,cs,rs,ValueComp,ix,ret,Mc>::call(mc,x));
        }
    };
    
    // algo -4: No branches or copies
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<-4,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int cs2 = cs > 20 ? TMV_UNKNOWN : cs;
            const int rs2 = rs > 20 ? TMV_UNKNOWN : rs;
            // nops = m*n
            const int nops = IntTraits2<cs2,rs2>::safeprod;
            const bool unroll = 
                ( cs > 10 && rs > 10 ) ? false :
                ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? false :
                nops == TMV_UNKNOWN ? false :
                // Norm is faster with the regular algorithm except for 
                // very small matrices.
                (nops > 9 && comp == NormComp) ? false :
                nops <= TMV_NORMM_UNROLL;
            const int algo = 
                M1::_canlin ? 1 : 
                unroll ? ( M1::_colmajor ? 15 : 16 ) :
                M1::_rowmajor ? 12 :
                M1::_colmajor ? 11 :
                rs == TMV_UNKNOWN ? 11 :
                cs == TMV_UNKNOWN ? 12 :
                ( cs < rs ) ? 11 : 12;
#ifdef PRINTALGO_NormM
            std::cout<<"SumElementsM algo -4: "<<TMV_Text(comp)<<std::endl;
            std::cout<<"nops = "<<nops<<", unroll = "<<unroll<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return SumElementsM_Helper<algo,cs,rs,comp,ix,ret,M1>::call(m,x);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<-3,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int algo = 
                M1::_canlin ? 1 : 
                TMV_OPT >= 2 && (cs == TMV_UNKNOWN || rs == TMV_UNKNOWN) ? 30 :
                -4;
#ifdef PRINTALGO_NormM
            std::cout<<"Inline SumElementsM: "<<TMV_Text(comp)<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return SumElementsM_Helper<algo,cs,rs,comp,ix,ret,M1>::call(m,x);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<-2,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::value_type VT;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits<VT>::isinst;
            const int algo =
                M1::_canlin ? 201 : 
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            return SumElementsM_Helper<algo,cs,rs,comp,ix,ret,M1>::call(m,x);
        }
    };
    
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsM_Helper<-1,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return SumElementsM_Helper<-2,cs,rs,comp,ix,ret,M1>::call(m,x); }
    };
    
    template <class M>
    static inline typename M::value_type InlineSumElements(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::value_type VT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-3,cs,rs,ValueComp,1,VT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    static inline typename M::value_type DoSumElements(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::value_type VT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-2,cs,rs,ValueComp,1,VT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    static inline typename M::float_type InlineSumAbsElements(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-3,cs,rs,AbsComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    static inline typename M::float_type DoSumAbsElements(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-2,cs,rs,AbsComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    static inline typename M::real_type InlineSumAbs2Elements(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-3,cs,rs,Abs2Comp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    static inline typename M::real_type DoSumAbs2Elements(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-2,cs,rs,Abs2Comp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    static inline typename M::real_type InlineNormSq(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-3,cs,rs,NormComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    static inline typename M::real_type DoNormSq(const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-2,cs,rs,NormComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    static inline typename M::float_type InlineNormSq(
        const BaseMatrix_Rec<M>& m, typename M::float_type scale)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-3,cs,rs,NormComp,0,RT,Mv>::call(
            mv,Scaling<0,RT>(scale));
    }

    template <class M>
    static inline typename M::float_type DoNormSq(
        const BaseMatrix_Rec<M>& m, typename M::float_type scale)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsM_Helper<-2,cs,rs,NormComp,0,RT,Mv>::call(
            mv,Scaling<0,RT>(scale));
    }


    // 
    // MaxAbsElement
    // MaxAbs2Element
    //

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename ConstMatrixView<T>::float_type InstMaxAbsElement(
        const ConstMatrixView<T>& m); 
    template <class T>
    typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstMatrixView<T>& m); 

    template <int algo, int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper;

    // algo 1: linearize to vector version
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<1,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static inline ret call(const M1& m)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"MaxAbsElementM algo 1: "<<TMV_Text(comp)<<std::endl;
#endif
            typedef typename M1::const_linearview_type Ml;
            Ml ml = m.linearView();
            return MinMaxElement_Helper<-3,comp,true,Ml>::call(ml,0);
        }
    };

    // algo 201: linearize to vector version
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<201,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static inline ret call(const M1& m)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"MaxAbsElementM algo 201: "<<TMV_Text(comp)<<std::endl;
#endif
            typedef typename M1::const_linearview_type Ml;
            Ml ml = m.linearView();
            return MinMaxElement_Helper<-2,comp,true,Ml>::call(ml,0);
        }
    };

    // algo 11: loop over columns
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<11,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"MaxAbsElementM algo 11: "<<TMV_Text(comp)<<std::endl;
#endif
            typedef typename M1::const_col_type Mc;
            const int N = (rs == TMV_UNKNOWN ? m.rowsize() : rs);
            ret max(0);
            for(int j=0;j<N;++j) {
                ret temp = MinMaxElement_Helper<-3,comp,true,Mc>::call(
                    m.get_col(j),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 12: loop over rows
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<12,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"MaxAbsElementM algo 12: "<<TMV_Text(comp)<<std::endl;
#endif
            typedef typename M1::const_row_type Mr;
            const int M = (cs == TMV_UNKNOWN ? m.colsize() : cs);
            ret max(0);
            for(int i=0;i<M;++i) {
                ret temp = MinMaxElement_Helper<-3,comp,true,Mr>::call(
                    m.get_row(i),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 13: loop over columns with temp storage
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<13,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"MaxAbsElementM algo 13: "<<TMV_Text(comp)<<std::endl;
#endif
            typedef typename M1::const_col_type Mc;
            typedef tmv::Vector<ret> V;
            const int N = (rs == TMV_UNKNOWN ? m.rowsize() : rs);
            if (N == 0) return ret(0);
            else {
                V temp(N);
                for(int j=0;j<N;++j) {
                    temp(j) = MinMaxElement_Helper<-3,comp,true,Mc>::call(
                        m.get_col(j),0);
                }
                return MinMaxElement_Helper<-3,comp,true,V>::call(temp,0);
            }
        }
    };

    // algo 14: loop over rows with temp storage
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<14,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"MaxAbsElementM algo 14: "<<TMV_Text(comp)<<std::endl;
#endif
            typedef typename M1::const_row_type Mr;
            typedef tmv::Vector<ret> V;
            const int M = (cs == TMV_UNKNOWN ? m.colsize() : cs);
            if (M == 0) return ret(0);
            else {
                V temp(M);
                for(int i=0;i<M;++i) {
                    temp(i) = MinMaxElement_Helper<-3,comp,true,Mr>::call(
                        m.get_row(i),0);
                }
                return MinMaxElement_Helper<-3,comp,true,V>::call(temp,0);
            }
        }
    };

    // TODO: I don't have a full unroller here.
    // Usually, this will be run on a matrix that can be linearized, so
    // it will unroll there.  But for the rare case that we have a small matrix
    // that isn't the full SmallMatrix (and hence can't linearize), then
    // there might be a small speed increase available by unrolling.

    // algo 30: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<30,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"MaxAbsElementM algo 30: "<<TMV_Text(comp)<<std::endl;
#endif
            if (m.colsize() == 0 || m.rowsize() == 0) return ret(0);
            else if (m.canLinearize())
                return MaxAbsElementM_Helper<1,cs,rs,comp,M1>::call(m);
            else
                return MaxAbsElementM_Helper<-4,cs,rs,comp,M1>::call(m);
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, class M1>
    struct MaxAbsElementM_Helper<90,cs,rs,AbsComp,M1>
    {
        typedef typename M1::float_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return InstMaxAbsElement(m.xView()); }
    };
    template <int cs, int rs, class M1>
    struct MaxAbsElementM_Helper<90,cs,rs,Abs2Comp,M1>
    {
        typedef typename M1::real_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return InstMaxAbs2Element(m.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<97,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            typedef typename M1::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return MaxAbsElementM_Helper<-2,cs,rs,comp,Mnc>::call(mnc);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<-4,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            const int algo =
                M1::_canlin ? 1 :
                TMV_OPT == 0 ? ( M1::_rowmajor ? 12 : 11 ) :
                M1::_rowmajor ? 14 :
                M1::_colmajor ? 13 :
                cs < rs ? 14 : 13;
#ifdef PRINTALGO_NormM
            std::cout<<"MaxAbsElementM algo -4: "<<TMV_Text(comp)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return MaxAbsElementM_Helper<algo,cs,rs,comp,M1>::call(m);
        }
    };

    // algo -3: Determine which algo to use
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<-3,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            const int algo =
                M1::_canlin ? 1 :
                TMV_OPT >= 2 && ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? 30 :
                -4;
#ifdef PRINTALGO_NormM
            std::cout<<"Inline MaxAbsElementM: "<<TMV_Text(comp)<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return MaxAbsElementM_Helper<algo,cs,rs,comp,M1>::call(m);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<-2,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            typedef typename M1::value_type VT;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits<VT>::isinst;
            const int algo =
                M1::_canlin ? 201 : 
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            return MaxAbsElementM_Helper<algo,cs,rs,comp,M1>::call(m);
        }
    };

    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementM_Helper<-1,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return MaxAbsElementM_Helper<-2,cs,rs,comp,M1>::call(m); }
    };

    template <class M>
    static inline typename M::float_type InlineMaxAbsElement(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementM_Helper<-3,cs,rs,AbsComp,Mv>::call(mv);
    }

    template <class M>
    static inline typename M::float_type DoMaxAbsElement(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementM_Helper<-2,cs,rs,AbsComp,Mv>::call(mv);
    }

    template <class M>
    static inline typename M::real_type InlineMaxAbs2Element(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementM_Helper<-3,cs,rs,Abs2Comp,Mv>::call(mv);
    }

    template <class M>
    static inline typename M::real_type DoMaxAbs2Element(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementM_Helper<-2,cs,rs,Abs2Comp,Mv>::call(mv);
    }


    // 
    // Norm1
    // NormInf
    //

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename ConstMatrixView<T>::float_type InstNorm1(
        const ConstMatrixView<T>& m); 

    template <int algo, int cs, int rs, class M1>
    struct Norm1M_Helper;

    // algo 11: loop over columns
    template <int cs, int rs, class M1>
    struct Norm1M_Helper<11,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static RT call(const M1& m)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"Norm1M algo 1: "<<std::endl;
#endif
            const int N = (rs == TMV_UNKNOWN ? m.rowsize() : rs);
            RT max(0);
            for(int j=0;j<N;++j) {
                RT temp = InlineSumAbsElements(m.get_col(j));
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 12: loop over rows
    template <int cs, int rs, class M1>
    struct Norm1M_Helper<12,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static RT call(const M1& m)
        {
#ifdef PRINTALGO_NormM
            std::cout<<"Norm1M algo 12: "<<std::endl;
#endif
            int M = (cs == TMV_UNKNOWN ? m.colsize() : cs);
            int N = (rs == TMV_UNKNOWN ? m.rowsize() : rs);
            if (M == 0 || N == 0) return RT(0);

            if (M <= 8) return Norm1M_Helper<11,cs,rs,M1>::call(m);

            typedef typename M1::value_type VT;
            typedef typename tmv::Vector<RT>::iterator IT1;
            typedef typename M1::const_row_type::const_iterator IT2;

            VT value;

            tmv::Vector<RT> temp(N,RT(0));
            const IT1 begin1 = temp.begin();
            const IT1 end1 = temp.end();

            IT1 it1 = begin1;
            IT2 it2 = m.get_row(0).begin();
            int end_step = m.stepi() - N; // back to the start of next row

            do {
                do {
                    value = *it2++;
                    Component<AbsComp,VT>::applyf(value);
                    *it1++ += TMV_REAL(value);
                } while (it1 != end1);
                it1 -= N;
                it2 += end_step;
            } while (--M);
            return InlineMaxElement(temp,0);
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, class M1>
    struct Norm1M_Helper<90,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        { return InstNorm1(m.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1>
    struct Norm1M_Helper<97,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        {
            typedef typename M1::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return Norm1M_Helper<-2,cs,rs,Mnc>::call(mnc);
        }
    };

    // algo -3: Determine which algo to use
    template <int cs, int rs, class M1>
    struct Norm1M_Helper<-3,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        {
            const int algo =
                TMV_OPT == 0 ? 11 :
                ( M1::_rowmajor &&
                  ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ||
                    (cs > 16 && rs < 512) || (cs > 8 && rs >= 512) ) ) ? 12 :
                11;
#ifdef PRINTALGO_NormM
            std::cout<<"Inline Norm1M: "<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return Norm1M_Helper<algo,cs,rs,M1>::call(m);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1>
    struct Norm1M_Helper<-2,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        {
            typedef typename M1::value_type VT;
            const bool inst =
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits<VT>::isinst;
            const int algo =
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            return Norm1M_Helper<algo,cs,rs,M1>::call(m);
        }
    };

    template <int cs, int rs, class M1>
    struct Norm1M_Helper<-1,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        { return Norm1M_Helper<-2,cs,rs,M1>::call(m); }
    };

    template <class M>
    static inline typename M::float_type InlineNorm1(
        const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm1M_Helper<-3,cs,rs,Mv>::call(mv);
    }

    template <class M>
    static inline typename M::float_type DoNorm1(const BaseMatrix_Rec<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm1M_Helper<-2,cs,rs,Mv>::call(mv);
    }

    template <class M>
    static inline typename M::float_type InlineNormInf(
        const BaseMatrix_Rec<M>& m)
    { return InlineNorm1(m.transpose()); }

    template <class M>
    static inline typename M::float_type DoNormInf(const BaseMatrix_Rec<M>& m)
    { return Norm1(m.transpose()); }


} // namespace tmv

#endif
