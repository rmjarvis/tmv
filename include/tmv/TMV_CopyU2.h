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

#ifndef TMV_CopyU_H
#define TMV_CopyU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_CopyV.h"

namespace tmv {

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);
    template <class M1, class M2>
    inline void NoAliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);
    template <class M1, class M2>
    inline void InlineCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);
    template <class M1, class M2>
    inline void AliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);

    // Defined in TMV_TriMatrix.cpp
    template <class T1, bool C1, class T2>
    void InstCopy(
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        UpperTriMatrixView<T2,UnknownDiag> m2);

    //
    // Copy Matrices
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

    template <int algo, int s, class M1, class M2>
    struct CopyU_Helper;

    // algo 1: m1 is unitdiag
    template <int s, class M1, class M2>
    struct CopyU_Helper<1,s,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_offdiag_type M1o;
            typedef typename M2::offdiag_type M2o;
            M1o m1o = m1.offDiag();
            M2o m2o = m2.offDiag();
            const int sm1 = IntTraits2<s,-1>::sum;
            CopyU_Helper<-3,sm1,M1o,M2o>::call(m1o,m2o);
            if (!m2.isunit()) m2.diag().setAllTo(1);
        }
    };

    // algo 2: Loop over rows
    template <int s, class M1, class M2>
    struct CopyU_Helper<2,s,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            int N = (s == UNKNOWN ? m2.size() : s);
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::row_sub_type M2r;
            typedef typename M1r::const_iterator IT1;
            typedef typename M2r::iterator IT2;
            const int step1 = m1.diagstep();
            const int step2 = m2.diagstep();
            IT1 it1 = m1.get_row(0,0,N).begin();
            IT2 it2 = m2.get_row(0,0,N).begin();
            std::cout<<"Start algo2: N = "<<N<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m1.row = "<<TMV_Text(m1.get_row(0,0,N))<<std::endl;
            std::cout<<"m2.row = "<<TMV_Text(m2.get_row(0,0,N))<<std::endl;
            std::cout<<"it1.step = "<<it1.step()<<std::endl;
            std::cout<<"it2.step = "<<it2.step()<<std::endl;
            for(;N;--N) {
                std::cout<<"N = "<<N<<std::endl;
                std::cout<<"it1 = "<<(it1.getP()-m1.cptr())<<std::endl;
                std::cout<<"it2 = "<<(it2.getP()-m2.cptr())<<std::endl;
                CopyV_Helper<-3,UNKNOWN,M1r,M2r>::call2(N,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 3: Loop over columns
    template <int s, class M1, class M2>
    struct CopyU_Helper<3,s,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            int N = (s == UNKNOWN ? m2.size() : s);
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::col_sub_type M2c;
            typedef typename M1c::const_iterator IT1;
            typedef typename M2c::iterator IT2;
            const int step1 = m1.stepj();
            const int step2 = m2.stepj();
            IT1 it1 = m1.get_col(0,0,1).begin();
            IT2 it2 = m2.get_col(0,0,1).begin();
            int M=1;
            for(;N;--N) {
                CopyV_Helper<-3,UNKNOWN,M1c,M2c>::call2(M++,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            } 
        }
    };

    // algo 5: Fully unroll by rows
    template <int s, class M1, class M2>
    struct CopyU_Helper<5,s,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static inline void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,N>::unroll(m1,m2);
                Unroller<I+M/2,M-M/2,J+M/2,N-M/2>::unroll(m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static inline void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,1,J,N/2>::unroll(m1,m2);
                Unroller<I,1,J+N/2,N-N/2>::unroll(m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        { static inline void unroll(const M1& , M2& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static inline void unroll(const M1& m1, M2& m2)
            { m2.ref(I,J) = m1.cref(I,J); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        { static inline void unroll(const M1& , M2& ) {} };

        static inline void call(const M1& m1, M2& m2)
        { Unroller<0,s,0,s>::unroll(m1,m2); }
    };

    // algo 6: Fully unroll by columns
    template <int s, class M1, class M2>
    struct CopyU_Helper<6,s,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static inline void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,M-(N-N/2),J,N/2>::unroll(m1,m2);
                Unroller<I,M,J+N/2,N-N/2>::unroll(m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,1>
        {
            static inline void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,1>::unroll(m1,m2);
                Unroller<I+M/2,M-M/2,J,1>::unroll(m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,0>
        { static inline void unroll(const M1& , M2& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static inline void unroll(const M1& m1, M2& m2)
            { m2.ref(I,J) = m1.cref(I,J); }
        };
        template <int I, int J>
        struct Unroller<I,0,J,1>
        { static inline void unroll(const M1& , M2& ) {} };

        static inline void call(const M1& m1, M2& m2)
        { Unroller<0,s,0,s>::unroll(m1,m2); }
    };

    // algo 90: UnknownDiag
    template <int s, class M1, class M2>
    struct CopyU_Helper<90,s,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            const int nops = IntTraits2<s2,s2p1>::prod / 2;
            const bool unroll = 
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const int algo2 = 
                unroll ? ( M2::mrowmajor ? 5 : 6 ) :
                M2::mrowmajor ? 2 : 3;
            std::cout<<"algo 90: m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"algo2 = "<<algo2<<std::endl;
            if (m1.isunit())
                CopyU_Helper<1,s,M1,M2>::call(m1,m2);
            else
                CopyU_Helper<algo2,s,M1,M2>::call(m1,m2);
        }
    };

    // algo 95: Transpose (and go back to -3, rather than -2)
    template <int s, class M1, class M2>
    struct CopyU_Helper<95,s,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            std::cout<<"algo 95: m1 = "<<TMV_Text(m1)<<std::endl;
            CopyU_Helper<-3,s,M1t,M2t>::call(m1t,m2t);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class M1, class M2>
    struct CopyU_Helper<-3,s,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::mupper == int(M2::mupper));
            TMVStaticAssert(M1::munit || M1::munknowndiag || !M2::munit);
            TMVAssert(m1.isunit() || !m2.isunit());
            typedef typename M2::value_type T2;
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::prod / 2;
            const bool unroll = 
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const int algo = 
                M2::mlower ? 95 :
                M1::munknowndiag ? 90 :
                M1::munit ? 1 :
                unroll ? ( M2::mrowmajor ? 5 : 6 ) :
                M2::mrowmajor ? 2 : 3;
            std::cout<<"Start CopyU -3:\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            CopyU_Helper<algo,s,M1,M2>::call(m1,m2);
            std::cout<<"m2 => "<<m2<<std::endl;
        }
    };

    // algo 96: Transpose
    template <int size, class M1, class M2>
    struct CopyU_Helper<96,size,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            std::cout<<"algo 96: m1 = "<<TMV_Text(m1)<<std::endl;
            CopyU_Helper<-2,size,M1t,M2t>::call(m1t,m2t);
        }
    };

    // algo 97: Conjugate
    template <int size, class M1, class M2>
    struct CopyU_Helper<97,size,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            std::cout<<"algo 97: m1 = "<<TMV_Text(m1)<<std::endl;
            CopyU_Helper<-2,size,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 98: Call inst
    template <int size, class M1, class M2>
    struct CopyU_Helper<98,size,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        { InstCopy(m1.xdView(),m2.xdView()); }
    };

    // algo -2: Check for inst
    template <int s, class M1, class M2>
    struct CopyU_Helper<-2,s,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                M1::msize == UNKNOWN && 
                M2::msize == UNKNOWN &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const bool conj = M2::mconj;
            const bool lower = M2::mlower;
            const int algo = 
                lower ? 96 :
                conj ? 97 :
                inst ? 98 :
                -3;
            std::cout<<"algo -2: m1 = "<<TMV_Text(m1)<<std::endl;
            CopyU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo 99: Check for aliases
    template <int size, class M1, class M2>
    struct CopyU_Helper<99,size,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            std::cout<<"algo 99: m1 = "<<TMV_Text(m1)<<std::endl;
            if ( !SameStorage(m1,m2) ||
                 OppositeStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                CopyU_Helper<-2,size,M1,M2>::call(m1,m2);
            } else if (ExactSameStorage(m1,m2)) {
                // They are already equal modulo a conjugation
                Maybe<M1::mconj != int(M2::mconj)>::conjself(m2);
                // Except possibly for the diagonal...
                if (m1.isunit() && !m2.isunit()) m2.diag().setAllTo(1);
            } else {
                // Need a temporary
                NoAliasCopy(m1.copy(),m2);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int s, class M1, class M2>
    struct CopyU_Helper<-1,s,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool noclobber = MStepHelper<M1,M2>::opp;
            const bool checkalias = 
                M1::msize == UNKNOWN && 
                M2::msize == UNKNOWN &&
                !noclobber;
            const int algo = 
                checkalias ? 99 : 
                -2;
            std::cout<<"algo -1: m1 = "<<TMV_Text(m1)<<std::endl;
            CopyU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::mupper == int(M2::mupper));
        TMVStaticAssert(M1::munit || M1::munknowndiag || !M2::munit);
        TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());
        const int s = Sizes<M1::msize,M2::msize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        std::cout<<"Copy: m1 = "<<TMV_Text(m1)<<std::endl;
        CopyU_Helper<-1,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void NoAliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::mupper == int(M2::mupper));
        TMVStaticAssert(M1::munit || M1::munknowndiag || !M2::munit);
        TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());
        const int s = Sizes<M1::msize,M2::msize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        std::cout<<"NoAliasCopy: m1 = "<<TMV_Text(m1)<<std::endl;
        CopyU_Helper<-2,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::mupper == int(M2::mupper));
        TMVStaticAssert(M1::munit || M1::munknowndiag || !M2::munit);
        TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());
        const int s = Sizes<M1::msize,M2::msize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        std::cout<<"InlineCopy: m1 = "<<TMV_Text(m1)<<std::endl;
        CopyU_Helper<-3,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void AliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::mupper == int(M2::mupper));
        TMVStaticAssert(M1::munit || M1::munknowndiag || !M2::munit);
        TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());
        const int s = Sizes<M1::msize,M2::msize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        std::cout<<"AliasCopy: m1 = "<<TMV_Text(m1)<<std::endl;
        CopyU_Helper<99,s,M1v,M2v>::call(m1v,m2v);
    }


    // 
    // M = U
    //

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const bool upper = M1::mupper;
        typedef typename TypeSelect<upper,
                typename M2::uppertri_type,
                typename M2::lowertri_type>::type M2u;
        M2u m2u = Maybe<upper>::uppertri(m2);
        Copy(m1,m2u);
        if (m1.size() > 0) 
            Maybe<!upper>::uppertri(m2).offDiag().setZero();
    }

    template <class M1, class M2>
    inline void NoAliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const bool upper = M1::mupper;
        typedef typename TypeSelect<upper,
                typename M2::uppertri_type,
                typename M2::lowertri_type>::type M2u;
        M2u m2u = Maybe<upper>::uppertri(m2);
        m2.setZero();
        NoAliasCopy(m1,m2u);
    }

    template <class M1, class M2>
    inline void AliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const bool upper = M1::mupper;
        typedef typename TypeSelect<upper,
                typename M2::uppertri_type,
                typename M2::lowertri_type>::type M2u;
        M2u m2u = Maybe<upper>::uppertri(m2);
        AliasCopy(m1,m2u);
        if (m1.size() > 0) 
            Maybe<!upper>::uppertri(m2).offDiag().setZero();
    }


    //
    // U = D
    //

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        Copy(m1d,m2d);
        m2.offDiag().setZero();
    }

    template <class M1, class M2>
    inline void NoAliasCopy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        m2.setZero();
        NoAliasCopy(m1d,m2d);
    }

    template <class M1, class M2>
    inline void AliasCopy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        AliasCopy(m1d,m2d);
        m2.offDiag().setZero();
    }

    template <class T1, bool C1, class T2>
    void InstCopy(
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        LowerTriMatrixView<T2,UnknownDiag> m2)
    { InstCopy(m1.transpose(),m2.transpose()); }


#undef TMV_Q1

} // namespace tmv

#endif
