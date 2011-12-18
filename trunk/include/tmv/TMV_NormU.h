
#ifndef TMV_NormU_H
#define TMV_NormU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_NormV.h"
#include "TMV_MinMax.h"

namespace tmv {

    // 
    // SumElements
    // SumAbsElements
    // SumAbs2Elements
    // NormSq
    //

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    T InstSumElements(const ConstUpperTriMatrixView<T>& m); 
    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstSumAbsElements(
        const ConstUpperTriMatrixView<T>& m); 
    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstUpperTriMatrixView<T>& m); 
    template <class T>
    typename Traits<T>::real_type InstNormSq(
        const ConstUpperTriMatrixView<T>& m); 
    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstNormSq(
        const ConstUpperTriMatrixView<T>& m, 
        typename ConstUpperTriMatrixView<T>::float_type scale); 

    // The maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_NORMU_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_NORMU_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_NORMU_UNROLL 9
#else
#define TMV_NORMU_UNROLL 0
#endif

    template <int algo, int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper;

    // algo 1: m1 is unitdiag
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<1,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormU
            std::cout<<"SumElementsU algo 1: "<<
                TMV_Text(comp)<<"  "<<RT(x)<<std::endl;
#endif

            const int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            typedef typename M1::const_offdiag_type Mo;
            Mo mo = m.offDiag();
            const int sm1 = IntTraits2<s,-1>::sum;
            return SumElementsU_Helper<-4,sm1,comp,ix,ret,Mo>::call(mo,x) 
                + Component<comp,RT>::f(RT(x)) * RT(N);
        }
    };

    // algo 2: unknown diag, figure out which it is.
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<2,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormU
            std::cout<<"SumElementsU algo 2: "<<
                TMV_Text(comp)<<"  "<<RT(x)<<std::endl;
#endif
            if (m.isunit()) 
                return SumElementsU_Helper<1,s,comp,ix,ret,M1>::call(m,x);
            else 
                return SumElementsU_Helper<-4,s,comp,ix,ret,M1>::call(m,x);
        }
    };

    // algo 11: loop over columns
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<11,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormU
            std::cout<<"SumElementsU algo 11: "<<
                TMV_Text(comp)<<"  "<<RT(x)<<std::endl;
#endif
            const int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            typedef typename M1::const_col_sub_type Mc;
            ret sum(0);
            for(int j=0;j<N;++j) {
                Mc mc = m.get_col(j,0,j+1);
                sum += SumElementsV_Helper<-3,TMV_UNKNOWN,comp,ix,ret,Mc>::call(mc,x);
            }
            return sum;
        }
    };

    // algo 12: loop over rows
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<12,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormU
            std::cout<<"SumElementsU algo 12: "<<
                TMV_Text(comp)<<"  "<<RT(x)<<std::endl;
#endif
            const int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            typedef typename M1::const_row_sub_type Mr;
            ret sum(0);
            for(int i=0;i<N;++i) {
                Mr mr = m.get_row(i,i,N);
                sum += SumElementsV_Helper<-3,TMV_UNKNOWN,comp,ix,ret,Mr>::call(mr,x);
            }
            return sum;
        }
    };

    // algo 15: Fully unroll by columns
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<15,s,comp,ix,ret,M1>
    {
        typedef typename M1::value_type VT;
        typedef typename Traits<ret>::real_type RT;

        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE ret unroll(const M1& m, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,M-(N-N/2),J,N/2>::unroll(m,x) +
                    Unroller<I,M,J+N/2,N-N/2>::unroll(m,x));
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,1>
        {
            static TMV_INLINE ret unroll(const M1& m, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,M/2,J,1>::unroll(m,x) +
                    Unroller<I+M/2,M-M/2,J,1>::unroll(m,x));
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,0>
        {
            static TMV_INLINE ret unroll(const M1& , const Scaling<ix,RT>& ) 
            { return ret(0); } 
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return Component<comp,VT>::f(x * m.cref(I,J)); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        {
            static TMV_INLINE ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return ret(0); }
        };
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormU
            std::cout<<"SumElementsU algo 15: "<<
                TMV_Text(comp)<<"  "<<RT(x)<<std::endl;
#endif
            return Unroller<0,s,0,s>::unroll(m,x); 
        }
    };

    // algo 16: Fully unroll by rows
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<16,s,comp,ix,ret,M1>
    {
        typedef typename M1::value_type VT;
        typedef typename Traits<ret>::real_type RT;
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE ret unroll(const M1& m, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,M/2,J,N-(M-M/2)>::unroll(m,x) +
                    Unroller<I+M/2,M-M/2,J,N>::unroll(m,x));
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static TMV_INLINE ret unroll(const M1& m, const Scaling<ix,RT>& x)
            {
                return (
                    Unroller<I,1,J,N/2>::unroll(m,x) +
                    Unroller<I,1,J+N/2,N-N/2>::unroll(m,x));
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        {
            static TMV_INLINE ret unroll(const M1& , const Scaling<ix,RT>& ) 
            { return ret(0); } 
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return Component<comp,VT>::f(x * m.cref(I,J)); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        {
            static TMV_INLINE ret unroll(const M1& m, const Scaling<ix,RT>& x)
            { return ret(0); }
        };
        static inline ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormU
            std::cout<<"SumElementsU algo 16: "<<
                TMV_Text(comp)<<"  "<<RT(x)<<std::endl;
#endif
            return Unroller<0,s,0,s>::unroll(m,x); 
        }
    };

    // algo 90: Call inst
    template <int s, int ix, class ret, class M1>
    struct SumElementsU_Helper<90,s,ValueComp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return InstSumElements(m.xView()); }
    };
    template <int s, int ix, class ret, class M1>
    struct SumElementsU_Helper<90,s,AbsComp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return InstSumAbsElements(m.xView()); }
    };
    template <int s, int ix, class ret, class M1>
    struct SumElementsU_Helper<90,s,Abs2Comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return InstSumAbs2Elements(m.xView()); }
    };
    template <int s, class ret, class M1>
    struct SumElementsU_Helper<90,s,NormComp,1,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<1,RT>& x)
        { return InstNormSq(m.xView()); }
    };
    template <int s, class ret, class M1>
    struct SumElementsU_Helper<90,s,NormComp,0,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<0,RT>& x)
        { return InstNormSq(m.xView(),x); }
    };
    
    // algo 96: Transpose
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<96,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::const_transpose_type Mt;
            Mt mt = m.transpose();
            return SumElementsU_Helper<-2,s,comp,ix,ret,Mt>::call(mt,x);
        }
    };
    
    // algo 396: Transpose 
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<396,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::const_transpose_type Mt;
            Mt mt = m.transpose();
            return SumElementsU_Helper<-3,s,comp,ix,ret,Mt>::call(mt,x);
        }
    };

    // algo 97: Conjugate
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<97,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return SumElementsU_Helper<-2,s,comp,ix,ret,Mnc>::call(mnc,x);
        }
    };
    template <int s, int ix, class ret, class M1>
    struct SumElementsU_Helper<97,s,ValueComp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::const_conjugate_type Mc;
            Mc mc = m.conjugate();
            return TMV_CONJ(
                SumElementsU_Helper<-2,s,ValueComp,ix,ret,Mc>::call(mc,x));
        }
    };
    
    // algo -4: No branches or copies, and m is Upper, NonUnitDiag
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<-4,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            TMVStaticAssert(M1::_upper);
            TMVStaticAssert(!M1::_unit);
            TMVAssert(!m.isunit());
            const int s2 = s > 20 ? TMV_UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits<IntTraits2<s2,s2p1>::safeprod>::halfS;
            const bool unroll = 
                s > 10 ? false :
                s == TMV_UNKNOWN ? false :
                nops == TMV_UNKNOWN ? false :
                // Norm is faster with the regular algorithm except for 
                // very small matrices.
                (s > 3 && comp == NormComp) ? false :
                nops <= TMV_NORMU_UNROLL;
            const int algo = 
                unroll ? ( M1::_colmajor ? 15 : 16 ) :
                M1::_colmajor ? 11 : 12;
#ifdef PRINTALGO_NormU
            std::cout<<"Inline SumElementsU\n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"comp = "<<TMV_Text(comp)<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<RT(x)<<std::endl;
            std::cout<<"unroll = "<<unroll<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return SumElementsU_Helper<algo,s,comp,ix,ret,M1>::call(m,x);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<-3,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int algo = 
                M1::_lower ? 396 :
                M1::_unit ? 1 :
                M1::_unknowndiag ? 2 :
                -4;
#ifdef PRINTALGO_NormU
            std::cout<<"SumElementsU algo -3: "<<
                TMV_Text(comp)<<"  "<<RT(x)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return SumElementsU_Helper<algo,s,comp,ix,ret,M1>::call(m,x);
        }
    };

    // algo -2: Check for inst
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<-2,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::value_type VT;
            const bool inst = 
                (s == TMV_UNKNOWN || s > 16) &&
                Traits<VT>::isinst;
            const int algo =
                M1::_lower ? 96 :
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            return SumElementsU_Helper<algo,s,comp,ix,ret,M1>::call(m,x);
        }
    };
    
    template <int s, CompType comp, int ix, class ret, class M1>
    struct SumElementsU_Helper<-1,s,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return SumElementsU_Helper<-2,s,comp,ix,ret,M1>::call(m,x); }
    };
    
    template <class M>
    inline typename M::value_type InlineSumElements(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::value_type VT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-3,s,ValueComp,1,VT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::value_type DoSumElements(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::value_type VT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-2,s,ValueComp,1,VT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::float_type InlineSumAbsElements(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-3,s,AbsComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::float_type DoSumAbsElements(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-2,s,AbsComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::real_type InlineSumAbs2Elements(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-3,s,Abs2Comp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::real_type DoSumAbs2Elements(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-2,s,Abs2Comp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::real_type InlineNormSq(const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-3,s,NormComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::real_type DoNormSq(const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-2,s,NormComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::float_type InlineNormSq(
        const BaseMatrix_Tri<M>& m, typename M::float_type scale)
    {
        const int s = M::_size;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-3,s,NormComp,0,RT,Mv>::call(
            mv,Scaling<0,RT>(scale));
    }

    template <class M>
    inline typename M::float_type DoNormSq(
        const BaseMatrix_Tri<M>& m, typename M::float_type scale)
    {
        const int s = M::_size;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsU_Helper<-2,s,NormComp,0,RT,Mv>::call(
            mv,Scaling<0,RT>(scale));
    }


    // 
    // MaxAbsElement
    // MaxAbs2Element
    //

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstMaxAbsElement(
        const ConstUpperTriMatrixView<T>& m); 
    template <class T>
    typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstUpperTriMatrixView<T>& m); 

    template <int algo, int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper;

    // algo 1: m1 is unitdiag
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<1,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
            if (s == 0) return ret(0);
            else if (s == 1) return ret(1);
            else {
                typedef typename M1::const_offdiag_type Mo;
                Mo mo = m.offDiag();
                const int sm1 = IntTraits2<s,-1>::sum;
                ret temp = MaxAbsElementU_Helper<-4,sm1,comp,Mo>::call(mo);
                return (temp > ret(1)) ? temp : ret(1);
            }
        }
    };

    // algo 2: unknown diag, figure out which it is.
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<2,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
            const int algo2 = M1::_rowmajor ? 12 : 11;
            if (m.isunit()) 
                return MaxAbsElementU_Helper<1,s,comp,M1>::call(m);
            else 
                return MaxAbsElementU_Helper<algo2,s,comp,M1>::call(m);
        }
    };

    // algo 11: loop over columns
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<11,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
            typedef typename M1::const_col_sub_type Mc;
            const int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            ret max(0);
            for(int j=0;j<N;++j) {
                ret temp = MinMaxElement_Helper<-3,comp,true,Mc>::call(
                    m.get_col(j,0,j+1),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 12: loop over rows
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<12,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
            typedef typename M1::const_row_sub_type Mr;
            const int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            ret max(0);
            for(int i=0;i<N;++i) {
                ret temp = MinMaxElement_Helper<-3,comp,true,Mr>::call(
                    m.get_row(i,i,N),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 13: loop over columns with temp storage
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<13,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
            typedef typename M1::const_col_sub_type Mc;
            typedef tmv::Vector<ret> V;
            const int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            if (N == 0) return ret(0);
            else {
                V temp(N);
                for(int j=0;j<N;++j) {
                    temp(j) = MinMaxElement_Helper<-3,comp,true,Mc>::call(
                        m.get_col(j,0,j+1),0);
                }
                return MinMaxElement_Helper<-3,comp,true,V>::call(temp,0);
            }
        }
    };

    // algo 14: loop over rows with temp storage
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<14,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
            typedef typename M1::const_row_sub_type Mr;
            typedef tmv::Vector<ret> V;
            const int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            if (N == 0) return ret(0);
            else {
                V temp(N);
                for(int i=0;i<N;++i) {
                    temp(i) = MinMaxElement_Helper<-3,comp,true,Mr>::call(
                        m.get_row(i,i,N),0);
                }
                return MinMaxElement_Helper<-3,comp,true,V>::call(temp,0);
            }
        }
    };

    // algo 90: Call inst
    template <int s, class M1>
    struct MaxAbsElementU_Helper<90,s,AbsComp,M1>
    {
        typedef typename M1::float_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return InstMaxAbsElement(m.xView()); }
    };
    template <int s, class M1>
    struct MaxAbsElementU_Helper<90,s,Abs2Comp,M1>
    {
        typedef typename M1::real_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return InstMaxAbs2Element(m.xView()); }
    };

    // algo 96: Transpose
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<96,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            typedef typename M1::const_transpose_type Mt;
            Mt mt = m.transpose();
            return MaxAbsElementU_Helper<-2,s,comp,Mt>::call(mt);
        }
    };

    // algo 396: Transpose
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<396,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            typedef typename M1::const_transpose_type Mt;
            Mt mt = m.transpose();
            return MaxAbsElementU_Helper<-3,s,comp,Mt>::call(mt);
        }
    };

    // algo 97: Conjugate
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<97,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            typedef typename M1::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return MaxAbsElementU_Helper<-2,s,comp,Mnc>::call(mnc);
        }
    };

    // algo -4: No branches or copies, and m is Upper, NonUnitDiag
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<-4,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            TMVStaticAssert(M1::_upper);
            TMVStaticAssert(!M1::_unit);
            TMVAssert(!m.isunit());
            const int algo = 
                M1::_rowmajor ? 14 : 13;
            return MaxAbsElementU_Helper<algo,s,comp,M1>::call(m);
        }
    };

    // algo -3: Determine which algo to use
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<-3,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            const int algo = 
                M1::_lower ? 396 : 
                M1::_unit ? 1 : 
                M1::_unknowndiag ? 2 :
                -4;
            return MaxAbsElementU_Helper<algo,s,comp,M1>::call(m);
        }
    };

    // algo -2: Check for inst
    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<-2,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            typedef typename M1::value_type VT;
            const bool inst = 
                (s == TMV_UNKNOWN || s > 16) &&
                Traits<VT>::isinst;
            const int algo =
                M1::_lower ? 96 :
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            return MaxAbsElementU_Helper<algo,s,comp,M1>::call(m);
        }
    };

    template <int s, CompType comp, class M1>
    struct MaxAbsElementU_Helper<-1,s,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return MaxAbsElementU_Helper<-2,s,comp,M1>::call(m); }
    };

    template <class M>
    inline typename M::float_type InlineMaxAbsElement(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementU_Helper<-3,s,AbsComp,Mv>::call(mv);
    }

    template <class M>
    inline typename M::float_type DoMaxAbsElement(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementU_Helper<-2,s,AbsComp,Mv>::call(mv);
    }

    template <class M>
    inline typename M::real_type InlineMaxAbs2Element(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementU_Helper<-3,s,Abs2Comp,Mv>::call(mv);
    }

    template <class M>
    inline typename M::real_type DoMaxAbs2Element(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementU_Helper<-2,s,Abs2Comp,Mv>::call(mv);
    }


    // 
    // Norm1
    // NormInf
    //

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstNorm1(
        const ConstUpperTriMatrixView<T>& m); 
    template <class T>
    typename ConstLowerTriMatrixView<T>::float_type InstNorm1(
        const ConstLowerTriMatrixView<T>& m); 

    // TODO: Norm1 would benefit from an unroller.
    // Unlike with a regular matrix, where small matrices can unroll
    // each column (or row), with a triangle matrix, we lose the knowledge
    // of the length of each column in the for loop, so a full unroller
    // would be able to keep that.
    
    template <int algo, int s, class M1>
    struct Norm1U_Helper;

    // algo 1: unknown diag, figure out which it is.
    template <int s, class M1>
    struct Norm1U_Helper<1,s,M1>
    {
        typedef typename M1::float_type RT;
        static RT call(const M1& m)
        {
            if (m.isunit()) {
                typedef typename M1::const_unitdiag_type Mu;
                return Norm1U_Helper<-4,s,Mu>::call(m.viewAsUnitDiag());
            } else {
                typedef typename M1::const_nonunitdiag_type Mnu;
                return Norm1U_Helper<-4,s,Mnu>::call(m.viewAsNonUnitDiag());
            }
        }
    };

    // algo 11: loop over columns, uppertri
    template <int s, class M1>
    struct Norm1U_Helper<11,s,M1>
    {
        typedef typename M1::float_type RT;
        static RT call(const M1& m)
        {
            const int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
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

    // algo 12: loop over rows, uppertri
    template <int s, class M1>
    struct Norm1U_Helper<12,s,M1>
    {
        typedef typename M1::float_type RT;
        static RT call(const M1& m)
        {
            int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            if (N <= 8) return Norm1U_Helper<11,s,M1>::call(m);

            typedef typename M1::value_type VT;
            typedef typename tmv::Vector<RT>::iterator IT1;
            typedef typename M1::const_row_sub_type::const_iterator IT2;

            VT value;

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
                    Component<AbsComp,VT>::applyf(value);
                    *it1++ += TMV_REAL(value);
                } while (it1 != end1);
                it2 += end_step++;
                it1 -= (--N);
            } while (N);
            return InlineMaxElement(temp,0);
        }
    };

    // algo 21: loop over columns, lowertri
    template <int s, class M1>
    struct Norm1U_Helper<21,s,M1>
    {
        typedef typename M1::float_type RT;
        static RT call(const M1& m)
        {
            const int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            RT max(0);
            for(int i=0;i<N;++i) {
                // If unit,     temp = 1 + SumAbsElements(m.col(i,i+1,N)
                // If non-unit, temp =     SumAbsElements(m.col(i,i,N)
                RT temp = Maybe<M1::_unit>::sum( 
                    RT(1) , InlineSumAbsElements(
                        m.get_col(i,Maybe<M1::_unit>::select(i+1,i),N)));
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 22: loop over rows, lowertri
    template <int s, class M1>
    struct Norm1U_Helper<22,s,M1>
    {
        typedef typename M1::float_type RT;
        static RT call(const M1& m)
        {
            int N = (s == TMV_UNKNOWN ? int(m.size()) : s);
            if (N <= 8) return Norm1U_Helper<21,s,M1>::call(m);

            typedef typename M1::value_type VT;
            typedef typename tmv::Vector<RT>::iterator IT1;
            typedef typename M1::const_col_sub_type::const_iterator IT2;

            VT value;

            // If unit,     start with all 1's.
            // If non-unit, start with all 0's.
            tmv::Vector<RT> temp(N, Maybe<M1::_unit>::select(RT(1),RT(0)) );
            const IT1 begin1 = temp.begin();

            IT1 it1 = begin1;
            IT2 it2 = m.get_row(Maybe<M1::_unit>::select(1,0),0,1).begin();
            // If unit, then we only need to add the offdiag to temp.
            // This is effected by: --N, and the above select for it2.
            Maybe<M1::_unit>::decrement(N);
            int end_step = m.stepi()-1; // back to the start of next column

            int M=1, i;
            do {
                i=M; do {
                    value = *it2++;
                    Component<AbsComp,VT>::applyf(value);
                    *it1++ += TMV_REAL(value);
                } while (--i);
                it1 -= M++;
                it2 += end_step--;
            } while (--N);
            return InlineMaxElement(temp,0);
        }
    };

    // algo 90: Call inst
    template <int s, class M1>
    struct Norm1U_Helper<90,s,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        { return InstNorm1(m.xView()); }
    };

    // algo 97: Conjugate
    template <int s, class M1>
    struct Norm1U_Helper<97,s,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        {
            typedef typename M1::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return Norm1U_Helper<-2,s,Mnc>::call(mnc);
        }
    };

    // algo -4: No branches or copies
    template <int s, class M1>
    struct Norm1U_Helper<-4,s,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        {
            TMVStaticAssert(!M1::_unknowndiag);
            const int algo = 
                M1::_upper ? 
                TMV_OPT == 0 ? 11 :
                ( M1::_rowmajor && (s == TMV_UNKNOWN || s > 8) ) ? 12 :
                11 :

                // lower
                TMV_OPT == 0 ? 21 :
                ( M1::_rowmajor && (s == TMV_UNKNOWN || s > 8) ) ? 22 :
                21 ;
            return Norm1U_Helper<algo,s,M1>::call(m);
        }
    };

    // algo -3: Determine which algo to use
    template <int s, class M1>
    struct Norm1U_Helper<-3,s,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        {
            const int algo = M1::_unknowndiag ? 1 : -4;
            return Norm1U_Helper<algo,s,M1>::call(m);
        }
    };

    // algo -2: Check for inst
    template <int s, class M1>
    struct Norm1U_Helper<-2,s,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        {
            typedef typename M1::value_type VT;
            const bool inst =
                (s == TMV_UNKNOWN || s > 16) &&
                Traits<VT>::isinst;
            const int algo =
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            return Norm1U_Helper<algo,s,M1>::call(m);
        }
    };

    template <int s, class M1>
    struct Norm1U_Helper<-1,s,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        { return Norm1U_Helper<-2,s,M1>::call(m); }
    };

    template <class M>
    inline typename M::float_type InlineNorm1(
        const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm1U_Helper<-3,s,Mv>::call(mv);
    }

    template <class M>
    inline typename M::float_type DoNorm1(const BaseMatrix_Tri<M>& m)
    {
        const int s = M::_size;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm1U_Helper<-2,s,Mv>::call(mv);
    }

    template <class M>
    inline typename M::float_type InlineNormInf(
        const BaseMatrix_Tri<M>& m)
    { return InlineNorm1(m.transpose()); }

    template <class M>
    inline typename M::float_type DoNormInf(const BaseMatrix_Tri<M>& m)
    { return Norm1(m.transpose()); }

} // namespace tmv

#endif
