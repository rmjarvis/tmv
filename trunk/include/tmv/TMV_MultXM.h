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

#ifndef TMV_MultXM_H
#define TMV_MultXM_H

#include "TMV_MultXV.h"
#include "TMV_CopyM.h"
#include "TMV_ScaleM.h"
#include "TMV_Scaling.h"

namespace tmv {

    // Defined in TMV_MultXM.cpp
    template <class T1, bool C1, class T2>
    void InstMultXM(
        const T2 x, const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        MatrixView<T2> m2);
    template <class T1, bool C1, class T2>
    void InstAddMultXM(
        const T2 x, const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        MatrixView<T2> m2);

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2);
    //
    // Matrix += x * Matrix
    //

    template <int algo, int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper;

    // algo 1: trivial: ix == 1, !add, so call Copy (with alias check)
    template <int cs, int rs, class T, class M1, class M2>
    struct MultXM_Helper<1,cs,rs,false,1,T,M1,M2>
    {
        static void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        { CopyM_Helper<-3,cs,rs,M1,M2>::call(m1,m2); }
    };

    // algo 101: Same as algo 1, but use algo -3 for Copy
    template <int cs, int rs, class T, class M1, class M2>
    struct MultXM_Helper<101,cs,rs,false,1,T,M1,M2>
    {
        static void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        { CopyM_Helper<-1,cs,rs,M1,M2>::call(m1,m2); }
    };

    // algo 201: Same as algo 1, but use algo -2 for Copy
    template <int cs, int rs, class T, class M1, class M2>
    struct MultXM_Helper<201,cs,rs,false,1,T,M1,M2>
    {
        static void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        {
            // Need a quick alias check here, since ExactSameStorage
            // is allowed for MultXM, but not Copy
            if (!SameStorage(m1,m2)) {
                CopyM_Helper<-2,cs,rs,M1,M2>::call(m1,m2); 
            } else {
                TMVAssert(ExactSameStorage(m1,m2));
                Maybe<(M2::_conj != int(M1::_conj))>::conjself(m2);
            }
        }
    };

    // algo 2: Linearize to vector version
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<2,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_linearview_type M1l;
            typedef typename M2::linearview_type M2l;
            M1l m1l = m1.linearView();
            M2l m2l = m2.linearView();
            const int prod = IntTraits2<cs,rs>::prod;
            MultXV_Helper<-4,prod,add,ix,T,M1l,M2l>::call(x,m1l,m2l);
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<11,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
            int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
            typedef typename M1::const_col_type M1c;
            typedef typename M2::col_type M2c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M2c::iterator IT2;
            const int step1 = m1.stepj();
            const int step2 = m2.stepj();
            IT1 it1 = m1.get_col(0).begin().nonConj();
            IT2 it2 = m2.get_col(0).begin();
            for(;N;--N) {
                MultXV_Helper<-4,cs,add,ix,T,M1c,M2c>::call2(M,x,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 12: Loop over rows
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<12,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
            typedef typename M1::const_row_type M1r;
            typedef typename M2::row_type M2r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2r::iterator IT2;
            const int step1 = m1.stepi();
            const int step2 = m2.stepi();
            IT1 it1 = m1.get_row(0).begin().nonConj();
            IT2 it2 = m2.get_row(0).begin();
            for(;M;--M) {
                MultXV_Helper<-4,rs,add,ix,T,M1r,M2r>::call2(N,x,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 15: Fully unroll by columns
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<15,cs,rs,add,ix,T,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,M,J,N/2>::unroll(x,m1,m2);
                Unroller<I,M,J+N/2,N-N/2>::unroll(x,m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,1>
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,1>::unroll(x,m1,m2);
                Unroller<I+M/2,M-M/2,J,1>::unroll(x,m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,0>
        {
            static void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Maybe<add>::add( 
                    m2.ref(I,J), ZProd<false,false>::prod(x,m1.cref(I,J)) ); 
            }
        };
        template <int I, int J>
        struct Unroller<I,0,J,1>
        {
            static void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };

        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { Unroller<0,cs,0,rs>::unroll(x,m1,m2); }
    };

    // algo 16: Fully unroll by rows
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<16,cs,rs,add,ix,T,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,N>::unroll(x,m1,m2);
                Unroller<I+M/2,M-M/2,J,N>::unroll(x,m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,1,J,N/2>::unroll(x,m1,m2);
                Unroller<I,1,J+N/2,N-N/2>::unroll(x,m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        {
            static void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Maybe<add>::add( 
                    m2.ref(I,J), ZProd<false,false>::prod(x,m1.cref(I,J)) ); 
            }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        {
            static void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };

        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { Unroller<0,cs,0,rs>::unroll(x,m1,m2); }
    };

    // algo 30: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<30,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            if (m1.canLinearize() && m2.canLinearize() &&
                m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) 
                MultXM_Helper<2,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
            else if ((m1.isrm() && m2.isrm()) || 
                     ( !(m1.iscm() && m2.iscm()) &&
                       (m1.colsize() > m1.rowsize()) ) )
                MultXM_Helper<12,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
            else 
                MultXM_Helper<11,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<90,cs,rs,true,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultXM(xx,m1.xView(),m2.xView()); 
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<90,cs,rs,false,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultXM(xx,m1.xView(),m2.xView()); 
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<97,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            MultXM_Helper<-2,cs,rs,add,ix,T,M1c,M2c>::call(
                TMV_CONJ(x),m1c,m2c);
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<99,cs,rs,true,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            // More than just these are non-clobbering.
            // TODO: Add in the steps that work.
            if ( !SameStorage(m1,m2) ||
                 ExactSameStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                MultXM_Helper<-2,cs,rs,true,ix,T,M1,M2>::call(x,m1,m2);
            } else {
                // Need a temporary
                NoAliasMultXM<true>(x,m1.copy(),m2);
            }
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<99,cs,rs,false,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            if ( !SameStorage(m1,m2) ||
                 ExactSameStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                MultXM_Helper<-2,cs,rs,false,ix,T,M1,M2>::call(x,m1,m2);
            } else {
                // Let Copy handle the aliasing
                AliasCopy(m1,m2);
                Scale(x,m2);
            }
        }
    };

    // algo -4: No branches or copies 
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<-4,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type T2;
            const bool canlin = 
                M1::_canlin && M2::_canlin &&
                ( (M1::_rowmajor && M2::_rowmajor) ||
                  (M1::_colmajor && M2::_colmajor) );
            const int algo = 
                ( ix == 1 && !add ) ? 1 :
                canlin ? 2 :
                TMV_OPT == 0 ? 11 :
                ( cs != UNKNOWN && rs != UNKNOWN ) ? (
                    ( IntTraits2<cs,rs>::prod <= int(128/sizeof(T2)) ) ? (
                        ( M1::_rowmajor && M2::_rowmajor ) ? 16 : 15 ) :
                    ( M1::_rowmajor && M2::_rowmajor ) ? 12 : 
                    ( M1::_colmajor && M2::_colmajor ) ? 11 :
                    ( cs > rs ) ? 12 : 11 ) :
                ( M1::_rowmajor && M2::_rowmajor ) ? 12 : 
                11;
#ifdef PRINTALGO_XM
            std::cout<<"XM: algo = "<<algo<<std::endl;
#endif
            MultXM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<-3,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool canlin = 
                M1::_canlin && M2::_canlin &&
                ( (M1::_rowmajor && M2::_rowmajor) ||
                  (M1::_colmajor && M2::_colmajor) );
            const int algo = 
                ( ix == 1 && !add ) ? 1 :
                canlin ? 2 :
                TMV_OPT >= 2 && ( cs == UNKNOWN || rs == UNKNOWN ) ? 30 :
                -4;
#ifdef PRINTALGO_XM
            std::cout<<"InlineMultXM: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"cs,rs = "<<cs<<','<<rs<<
                " = "<<m2.colsize()<<','<<m2.rowsize()<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
            MultXM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<-2,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst =
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                ( ix == 1 && !add ) ? 201 :
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            MultXM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<-1,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool noclobber = MStepHelper<M1,M2>::same;
            const bool checkalias =
                M1::_colsize == UNKNOWN && M1::_rowsize == UNKNOWN &&
                M2::_colsize == UNKNOWN && M2::_rowsize == UNKNOWN &&
                !noclobber;
            const int algo = 
                ( ix == 1 && !add ) ? 101 :
                checkalias ? 99 : 
                -2;
            MultXM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXM_Helper<-1,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXM_Helper<-2,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void InlineMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXM_Helper<-3,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXM_Helper<99,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

} // namespace tmv

#endif 
