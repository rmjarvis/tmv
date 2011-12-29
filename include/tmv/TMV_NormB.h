
#ifndef TMV_NormB_H
#define TMV_NormB_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_NormV.h"
#include "TMV_MinMax.h"

#ifdef PRINTALGO_NormB
#include <iostream>
#include "TMV_BandMatrixIO.h"
#endif

namespace tmv {

    // 
    // SumElements
    // SumAbsElements
    // SumAbs2Elements
    // NormSq
    //

    // Defined in TMV_BandMatrix.cpp
    template <class T>
    T InstSumElements(const ConstBandMatrixView<T>& m); 
    template <class T>
    typename ConstBandMatrixView<T>::float_type InstSumAbsElements(
        const ConstBandMatrixView<T>& m); 
    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstBandMatrixView<T>& m); 
    template <class T>
    typename Traits<T>::real_type InstNormSq(
        const ConstBandMatrixView<T>& m); 
    template <class T>
    typename ConstBandMatrixView<T>::float_type InstNormSq(
        const ConstBandMatrixView<T>& m, 
        typename ConstBandMatrixView<T>::float_type scale); 

    template <int algo, int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsB_Helper;

    // algo 11: loop over columns
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsB_Helper<11,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormB
            std::cout<<"SumElementsM algo 11: "<<TMV_Text(comp)<<std::endl;
#endif
            const int M = cs == TMV_UNKNOWN ? m.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m.rowsize() : rs;
            const int xx = TMV_UNKNOWN;
            typedef typename M1::const_col_sub_type M1c;

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int j1 = m.nhi();
            const int j2 = TMV_MIN(N,M-m.nlo());
            const int j3 = TMV_MIN(N,M+m.nhi());
            int i1 = 0;
            int i2 = m.nlo()+1;
            ret sum(0);
            int j=0;
            for(;j<j1;++j) {
                sum += SumElementsV_Helper<-3,xx,comp,ix,ret,M1c>::call(
                    m.get_col(j,i1,i2),x);
                if (i2 < M) ++i2;
            }
            for(;j<j2;++j) {
                sum += SumElementsV_Helper<-3,lh,comp,ix,ret,M1c>::call(
                    m.get_col(j,i1++,i2++),x);
            }
            for(;j<j3;++j) {
                sum += SumElementsV_Helper<-3,lh,comp,ix,ret,M1c>::call(
                    m.get_col(j,i1++,M),x);
            }
            return sum;
        }
    };

    // algo 12: loop over rows
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsB_Helper<12,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormB
            std::cout<<"SumElementsM algo 12: "<<TMV_Text(comp)<<std::endl;
#endif
            const int M = cs == TMV_UNKNOWN ? m.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m.rowsize() : rs;
            const int xx = TMV_UNKNOWN;
            typedef typename M1::const_row_sub_type M1r;

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int i1 = m.nlo();
            const int i2 = TMV_MIN(M,N-m.nhi());
            const int i3 = TMV_MIN(M,N+m.nlo());
            int j1 = 0;
            int j2 = m.nhi()+1;
            ret sum(0);
            int i=0;
            for(;i<i1;++i) {
                sum += SumElementsV_Helper<-3,xx,comp,ix,ret,M1r>::call(
                    m.get_row(i,j1,j2),x);
                if (j2 < N) ++j2;
            }
            for(;i<i2;++i) {
                sum += SumElementsV_Helper<-3,lh,comp,ix,ret,M1r>::call(
                    m.get_row(i,j1++,j2++),x);
            }
            for(;i<i3;++i) {
                sum += SumElementsV_Helper<-3,lh,comp,ix,ret,M1r>::call(
                    m.get_row(i,j1++,N),x);
            }
            return sum;
        }
    };

    // algo 13: loop over diagonals
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsB_Helper<13,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;
        static ret call(const M1& m, const Scaling<ix,RT>& x)
        {
#ifdef PRINTALGO_NormB
            std::cout<<"SumElementsM algo 13: "<<TMV_Text(comp)<<std::endl;
#endif
            const int xx = TMV_UNKNOWN;
            typedef typename M1::const_diag_sub_type M1d;
            ret sum(0);
            for(int k=-m.nlo();k<=m.nhi();++k) {
                sum += SumElementsV_Helper<-3,xx,comp,ix,ret,M1d>::call(
                    m.get_diag(k),x);
            }
            return sum;
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, int ix, class ret, class M1>
    struct SumElementsB_Helper<90,cs,rs,ValueComp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return InstSumElements(m.xView()); }
    };
    template <int cs, int rs, int ix, class ret, class M1>
    struct SumElementsB_Helper<90,cs,rs,AbsComp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return InstSumAbsElements(m.xView()); }
    };
    template <int cs, int rs, int ix, class ret, class M1>
    struct SumElementsB_Helper<90,cs,rs,Abs2Comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return InstSumAbs2Elements(m.xView()); }
    };
    template <int cs, int rs, class ret, class M1>
    struct SumElementsB_Helper<90,cs,rs,NormComp,1,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<1,RT>& x)
        { return InstNormSq(m.xView()); }
    };
    template <int cs, int rs, class ret, class M1>
    struct SumElementsB_Helper<90,cs,rs,NormComp,0,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<0,RT>& x)
        { return InstNormSq(m.xView(),x); }
    };
    
    // algo 97: Conjugate
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsB_Helper<97,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return SumElementsB_Helper<-2,cs,rs,comp,ix,ret,Mnc>::call(mnc,x);
        }
    };
    template <int cs, int rs, int ix, class ret, class M1>
    struct SumElementsB_Helper<97,cs,rs,ValueComp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            typedef typename M1::const_conjugate_type Mc;
            Mc mc = m.conjugate();
            return TMV_CONJ(
                SumElementsB_Helper<-2,cs,rs,ValueComp,ix,ret,Mc>::call(mc,x));
        }
    };
    
    // algo -4: No branches or copies
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsB_Helper<-4,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        {
            const int algo = 
                TMV_OPT == 0 ? 13 :
                M1::_colmajor ? 11 :
                M1::_rowmajor ? 12 :
                13;
#ifdef PRINTALGO_NormB
            std::cout<<"SumElementsM algo -4: "<<TMV_Text(comp)<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            std::cout<<"x = "<<ix<<" "<<RT(x)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return SumElementsB_Helper<algo,cs,rs,comp,ix,ret,M1>::call(m,x);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsB_Helper<-3,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return SumElementsB_Helper<-4,cs,rs,comp,ix,ret,M1>::call(m,x); }
    };

    // algo -2: Check for inst
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsB_Helper<-2,cs,rs,comp,ix,ret,M1>
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
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            return SumElementsB_Helper<algo,cs,rs,comp,ix,ret,M1>::call(m,x);
        }
    };
    
    template <int cs, int rs, CompType comp, int ix, class ret, class M1>
    struct SumElementsB_Helper<-1,cs,rs,comp,ix,ret,M1>
    {
        typedef typename Traits<ret>::real_type RT;

        static TMV_INLINE ret call(const M1& m, const Scaling<ix,RT>& x)
        { return SumElementsB_Helper<-2,cs,rs,comp,ix,ret,M1>::call(m,x); }
    };
    
    template <class M>
    inline typename M::value_type InlineSumElements(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::value_type VT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-3,cs,rs,ValueComp,1,VT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::value_type DoSumElements(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::value_type VT;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-2,cs,rs,ValueComp,1,VT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::float_type InlineSumAbsElements(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-3,cs,rs,AbsComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::float_type DoSumAbsElements(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-2,cs,rs,AbsComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::real_type InlineSumAbs2Elements(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-3,cs,rs,Abs2Comp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::real_type DoSumAbs2Elements(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-2,cs,rs,Abs2Comp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::real_type InlineNormSq(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-3,cs,rs,NormComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::real_type DoNormSq(const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::real_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-2,cs,rs,NormComp,1,RT,Mv>::call(
            mv,Scaling<1,RT>());
    }

    template <class M>
    inline typename M::float_type InlineNormSq(
        const BaseMatrix_Band<M>& m, typename M::float_type scale)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-3,cs,rs,NormComp,0,RT,Mv>::call(
            mv,Scaling<0,RT>(scale));
    }

    template <class M>
    inline typename M::float_type DoNormSq(
        const BaseMatrix_Band<M>& m, typename M::float_type scale)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::float_type RT;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return SumElementsB_Helper<-2,cs,rs,NormComp,0,RT,Mv>::call(
            mv,Scaling<0,RT>(scale));
    }


    // 
    // MaxAbsElement
    // MaxAbs2Element
    //

    // Defined in TMV_BandMatrix.cpp
    template <class T>
    typename ConstBandMatrixView<T>::float_type InstMaxAbsElement(
        const ConstBandMatrixView<T>& m); 
    template <class T>
    typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstBandMatrixView<T>& m); 

    template <int algo, int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementB_Helper;

    // algo 11: loop over columns
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementB_Helper<11,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
#ifdef PRINTALGO_NormB
            std::cout<<"MaxAbsElementM algo 11: "<<TMV_Text(comp)<<std::endl;
#endif
            const int M = cs == TMV_UNKNOWN ? m.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m.rowsize() : rs;
            typedef typename M1::const_col_sub_type M1c;

            const int j1 = m.nhi();
            const int j2 = TMV_MIN(N,M-m.nlo());
            const int j3 = TMV_MIN(N,M+m.nhi());
            int i1 = 0;
            int i2 = m.nlo()+1;
            ret max(0);
            int j=0;
            for(;j<j1;++j) {
                ret temp = MinMaxElement_Helper<-3,comp,true,M1c>::call(
                    m.get_col(j,i1,i2),0);
                if (temp > max) max = temp;
                if (i2 < M) ++i2;
            }
            for(;j<j2;++j) {
                ret temp = MinMaxElement_Helper<-3,comp,true,M1c>::call(
                    m.get_col(j,i1++,i2++),0);
                if (temp > max) max = temp;
            }
            for(;j<j3;++j) {
                ret temp = MinMaxElement_Helper<-3,comp,true,M1c>::call(
                    m.get_col(j,i1++,M),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 12: loop over rows
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementB_Helper<12,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
#ifdef PRINTALGO_NormB
            std::cout<<"MaxAbsElementM algo 11: "<<TMV_Text(comp)<<std::endl;
#endif
            const int M = cs == TMV_UNKNOWN ? m.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m.rowsize() : rs;
            typedef typename M1::const_row_sub_type M1r;

            const int i1 = m.nlo();
            const int i2 = TMV_MIN(M,N-m.nhi());
            const int i3 = TMV_MIN(M,N+m.nlo());
            int j1 = 0;
            int j2 = m.nhi()+1;
            ret max(0);
            int i=0;
            for(;i<i1;++i) {
                ret temp = MinMaxElement_Helper<-3,comp,true,M1r>::call(
                    m.get_row(i,j1,j2),0);
                if (temp > max) max = temp;
                if (j2 < N) ++j2;
            }
            for(;i<i2;++i) {
                ret temp = MinMaxElement_Helper<-3,comp,true,M1r>::call(
                    m.get_row(i,j1++,j2++),0);
                if (temp > max) max = temp;
            }
            for(;i<i3;++i) {
                ret temp = MinMaxElement_Helper<-3,comp,true,M1r>::call(
                    m.get_row(i,j1++,N),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 13: loop over diagonals
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementB_Helper<13,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static ret call(const M1& m)
        {
#ifdef PRINTALGO_NormB
            std::cout<<"MaxAbsElementM algo 13: "<<TMV_Text(comp)<<std::endl;
#endif
            typedef typename M1::const_diag_sub_type M1d;
            ret max(0);
            for(int k=-m.nlo();k<=m.nhi();++k) {
                ret temp = MinMaxElement_Helper<-3,comp,true,M1d>::call(
                    m.get_diag(k),0);
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, class M1>
    struct MaxAbsElementB_Helper<90,cs,rs,AbsComp,M1>
    {
        typedef typename M1::float_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return InstMaxAbsElement(m.xView()); }
    };
    template <int cs, int rs, class M1>
    struct MaxAbsElementB_Helper<90,cs,rs,Abs2Comp,M1>
    {
        typedef typename M1::real_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return InstMaxAbs2Element(m.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementB_Helper<97,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            typedef typename M1::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return MaxAbsElementB_Helper<-2,cs,rs,comp,Mnc>::call(mnc);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementB_Helper<-4,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        {
            const int algo =
                TMV_OPT == 0 ? 13 :
                M1::_colmajor ? 11 :
                M1::_rowmajor ? 12 :
                13;
#ifdef PRINTALGO_NormB
            std::cout<<"MaxAbsElementM algo -4: "<<TMV_Text(comp)<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return MaxAbsElementB_Helper<algo,cs,rs,comp,M1>::call(m);
        }
    };

    // algo -3: Determine which algo to use
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementB_Helper<-3,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return MaxAbsElementB_Helper<-4,cs,rs,comp,M1>::call(m); }
    };

    // algo -2: Check for inst
    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementB_Helper<-2,cs,rs,comp,M1>
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
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            return MaxAbsElementB_Helper<algo,cs,rs,comp,M1>::call(m);
        }
    };

    template <int cs, int rs, CompType comp, class M1>
    struct MaxAbsElementB_Helper<-1,cs,rs,comp,M1>
    {
        typedef typename Component<comp,typename M1::value_type>::ret_type ret;
        static TMV_INLINE ret call(const M1& m)
        { return MaxAbsElementB_Helper<-2,cs,rs,comp,M1>::call(m); }
    };

    template <class M>
    inline typename M::float_type InlineMaxAbsElement(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementB_Helper<-3,cs,rs,AbsComp,Mv>::call(mv);
    }

    template <class M>
    inline typename M::float_type DoMaxAbsElement(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementB_Helper<-2,cs,rs,AbsComp,Mv>::call(mv);
    }

    template <class M>
    inline typename M::real_type InlineMaxAbs2Element(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementB_Helper<-3,cs,rs,Abs2Comp,Mv>::call(mv);
    }

    template <class M>
    inline typename M::real_type DoMaxAbs2Element(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return MaxAbsElementB_Helper<-2,cs,rs,Abs2Comp,Mv>::call(mv);
    }


    // 
    // Norm1
    // NormInf
    //

    // Defined in TMV_BandMatrix.cpp
    template <class T>
    typename ConstBandMatrixView<T>::float_type InstNorm1(
        const ConstBandMatrixView<T>& m); 

    template <int algo, int cs, int rs, class M1>
    struct Norm1B_Helper;

    // algo 11: loop over columns
    template <int cs, int rs, class M1>
    struct Norm1B_Helper<11,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static RT call(const M1& m)
        {
#ifdef PRINTALGO_NormB
            std::cout<<"Norm1M algo 11: "<<std::endl;
#endif
            const int M = cs == TMV_UNKNOWN ? m.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m.rowsize() : rs;
            typedef typename M1::const_col_sub_type M1c;

            const int j1 = m.nhi();
            const int j2 = TMV_MIN(N,M-m.nlo());
            const int j3 = TMV_MIN(N,M+m.nhi());
            int i1 = 0;
            int i2 = m.nlo()+1;
            RT max(0);
            int j=0;
            for(;j<j1;++j) {
                RT temp = InlineSumAbsElements(m.get_col(j,i1,i2));
                if (temp > max) max = temp;
                if (i2 < M) ++i2;
            }
            for(;j<j2;++j) {
                RT temp = InlineSumAbsElements(m.get_col(j,i1++,i2++));
                if (temp > max) max = temp;
            }
            for(;j<j3;++j) {
                RT temp = InlineSumAbsElements(m.get_col(j,i1++,M));
                if (temp > max) max = temp;
            }
            return max;
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, class M1>
    struct Norm1B_Helper<90,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        { return InstNorm1(m.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1>
    struct Norm1B_Helper<97,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        {
            typedef typename M1::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return Norm1B_Helper<-2,cs,rs,Mnc>::call(mnc);
        }
    };

    // algo -3: Determine which algo to use
    template <int cs, int rs, class M1>
    struct Norm1B_Helper<-3,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        {
            const int algo = 11;
#ifdef PRINTALGO_NormB
            std::cout<<"Inline Norm1M: "<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return Norm1B_Helper<algo,cs,rs,M1>::call(m);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1>
    struct Norm1B_Helper<-2,cs,rs,M1>
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
            return Norm1B_Helper<algo,cs,rs,M1>::call(m);
        }
    };

    template <int cs, int rs, class M1>
    struct Norm1B_Helper<-1,cs,rs,M1>
    {
        typedef typename M1::float_type RT;
        static TMV_INLINE RT call(const M1& m)
        { return Norm1B_Helper<-2,cs,rs,M1>::call(m); }
    };

    template <class M>
    inline typename M::float_type InlineNorm1(
        const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm1B_Helper<-3,cs,rs,Mv>::call(mv);
    }

    template <class M>
    inline typename M::float_type DoNorm1(const BaseMatrix_Band<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm1B_Helper<-2,cs,rs,Mv>::call(mv);
    }

    template <class M>
    inline typename M::float_type InlineNormInf(
        const BaseMatrix_Band<M>& m)
    { return InlineNorm1(m.transpose()); }

    template <class M>
    inline typename M::float_type DoNormInf(const BaseMatrix_Band<M>& m)
    { return Norm1(m.transpose()); }


} // namespace tmv

#endif
