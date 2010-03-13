//////////////////////////////////////////////////////////////////////////////
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

#ifndef TMV_ScaleU_H
#define TMV_ScaleU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_ScaleV.h"

namespace tmv {

    // Defined below:
    template <int ix, class T, class M>
    inline void Scale(
        const Scaling<ix,T>& x, BaseMatrix_Tri_Mutable<M>& m);
    template <int ix, class T, class M>
    inline void InlineScale(
        const Scaling<ix,T>& x, BaseMatrix_Tri_Mutable<M>& m);

    // Defined in TMV_ScaleU.cpp
    template <class T>
    void InstScale(const T x, UpperTriMatrixView<T,NonUnitDiag> m);

    //
    // Matrix *= Scalar 
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

    template <int algo, int s, int ix, class T, class M1>
    struct ScaleU_Helper;

    // algo 0: trivial: s == 0 or ix == 1, so nothing to do
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<0,s,ix,T,M1>
    { static inline void call(const Scaling<1,T>& , M1& ) {} };

    // algo 1: transpose
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<1,s,ix,T,M1>
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::transpose_type Mt;
            Mt mt = m.transpose();
            ScaleU_Helper<-4,s,ix,T,Mt>::call(x,mt);
        }
    };

    // algo 11: Loop over columns
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<11,s,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            int N = (s == UNKNOWN ? m.size() : s);
            typedef typename M1::col_sub_type Mc;
            typedef typename Mc::iterator IT;
            const int step = m.stepj();
            IT it = m.get_col(0,0,1).begin();
            int M=1;
            for(;N;--N) {
                ScaleV_Helper<-4,UNKNOWN,ix,T,Mc>::call2(M++,x,it);
                it.shiftP(step);
            }
        }
    };

    // algo 15: Fully unroll by columns
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<15,s,ix,T,M1>
    {
        template <int I, int M, int J, int N, bool iscomplex>
        struct Unroller
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,M-(N-N/2),J,N/2,iscomplex>::unroll(x,m);
                Unroller<I,M,J+N/2,N-N/2,iscomplex>::unroll(x,m);
            }
        };
        template <int I, int M, int J, bool iscomplex>
        struct Unroller<I,M,J,1,iscomplex>
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,M/2,J,1,iscomplex>::unroll(x,m);
                Unroller<I+M/2,M-M/2,J,1,iscomplex>::unroll(x,m);
            }
        };
        template <int I, int M, int J, bool iscomplex>
        struct Unroller<I,M,J,0,iscomplex>
        { static inline void unroll(const Scaling<ix,T>& , M1& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1,false>
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            { m.ref(I,J) *= x; }
        };
        template <int I, int J>
        struct Unroller<I,1,J,1,true>
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            {
                typedef typename M1::real_type RT;
                typedef typename M1::value_type VT;
                const RT rm =
                    ZProd<false,M1::mconj>::rprod(x,m.nonConj().cref(I,J));
                const RT im =
                    ZProd<false,M1::mconj>::iprod(x,m.nonConj().cref(I,J));
                m.ref(I,J) = VT(rm,im);
            }
        };
        template <int I, int J, bool iscomplex>
        struct Unroller<I,0,J,1,iscomplex>
        { static inline void unroll(const Scaling<ix,T>& , M1& ) {} };

        static inline void call(const Scaling<ix,T>& x, M1& m)
        { Unroller<0,s,0,s,M1::miscomplex>::unroll(x,m); }
    };

    // algo 21: Loop over rows
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<21,s,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            int N = (s == UNKNOWN ? m.size() : s);
            typedef typename M1::row_sub_type Mr;
            typedef typename Mr::iterator IT;
            const int step = m.diagstep();
            IT it = m.get_row(0,0,N).begin();
            for(;N;--N) {
                ScaleV_Helper<-4,UNKNOWN,ix,T,Mr>::call2(N,x,it);
                it.shiftP(step);
            }
        }
    };

    // algo 25: Fully unroll by rows
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<25,s,ix,T,M1>
    {
        template <int I, int M, int J, int N, bool iscomplex>
        struct Unroller
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,M/2,J,N,iscomplex>::unroll(x,m);
                Unroller<I+M/2,M-M/2,J+M/2,N-M/2,iscomplex>::unroll(x,m);
            }
        };
        template <int I, int J, int N, bool iscomplex>
        struct Unroller<I,1,J,N,iscomplex>
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,1,J,N/2,iscomplex>::unroll(x,m);
                Unroller<I,1,J+N/2,N-N/2,iscomplex>::unroll(x,m);
            }
        };
        template <int I, int J, int N, bool iscomplex>
        struct Unroller<I,0,J,N,iscomplex>
        { static inline void unroll(const Scaling<ix,T>& , M1& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1,false>
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            { m.ref(I,J) *= x; }
        };
        template <int I, int J>
        struct Unroller<I,1,J,1,true>
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            {
                typedef typename M1::real_type RT;
                typedef typename M1::value_type VT;
                const RT rm =
                    ZProd<false,M1::mconj>::rprod(x,m.nonConj().cref(I,J));
                const RT im =
                    ZProd<false,M1::mconj>::iprod(x,m.nonConj().cref(I,J));
                m.ref(I,J) = VT(rm,im);
            }
        };
        template <int I, int J, bool iscomplex>
        struct Unroller<I,1,J,0,iscomplex>
        { static inline void unroll(const Scaling<ix,T>& , M1& ) {} };

        static inline void call(const Scaling<ix,T>& x, M1& m)
        { Unroller<0,s,0,s,M1::miscomplex>::unroll(x,m); }
    };

    // algo -4: No branches or copies
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<-4,s,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            TMVStaticAssert(!M1::munit || ix == 1);
            typedef typename M1::value_type T1;
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const int algo = 
                (s == 0 || ix == 1) ? 0 :
                M1::mlower ? 1 :
                unroll ? ( M1::mrowmajor ? 25 : 15 ) :
                M1::mrowmajor ? 21 : 11;
            ScaleU_Helper<algo,s,ix,T,M1>::call(x,m);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<-3,s,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        { ScaleU_Helper<-4,s,ix,T,M1>::call(x,m); }
    };

    // algo 96: Transpose
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<96,s,ix,T,M1>
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::transpose_type Mt;
            Mt mt = m.transpose();
            ScaleU_Helper<-2,s,ix,T,Mt>::call(x,mt);
        }
    };

    // algo 97: Conjugate
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<97,s,ix,T,M1>
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::conjugate_type Mc;
            Mc mc = m.conjugate();
            ScaleU_Helper<-2,s,ix,T,Mc>::call(TMV_CONJ(x),mc);
        }
    };

    // algo 98: Call inst
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<98,s,ix,T,M1>
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            TMVAssert(!m.isunit());
            typename M1::value_type xx(x);
            InstScale(xx,m.xView().viewAsNonUnitDiag());
        }
    };

    // algo -2: Check for inst
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<-2,s,ix,T,M1>
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::value_type T1;
            const bool inst =
                M1::msize == UNKNOWN &&
                Traits<T1>::isinst;
            const int algo =
                ix == 1 ? 0 :
                M1::mlower ? 96 :
                M1::mconj ? 97 :
                inst ? 98 :
                -4;
            ScaleU_Helper<algo,s,ix,T,M1>::call(x,m);
        }
    };

    // algo -1: Check for alias? No.
    template <int s, int ix, class T, class M1>
    struct ScaleU_Helper<-1,s,ix,T,M1>
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        { ScaleU_Helper<-2,s,ix,T,M1>::call(x,m); }
    };

    template <int ix, class T, class M>
    inline void Scale(
        const Scaling<ix,T>& x, BaseMatrix_Tri_Mutable<M>& m)
    {
        TMVStaticAssert(!M::munit || ix == 1);
        TMVAssert(!m.isunit() || ix == 1);
        typedef typename M::cview_type Mv;
        Mv mv = m.cView();
        ScaleU_Helper<-2,M::msize,ix,T,Mv>::call(x,mv); 
    }

    template <int ix, class T, class M>
    inline void InlineScale(
        const Scaling<ix,T>& x, BaseMatrix_Tri_Mutable<M>& m)
    {
        TMVStaticAssert(!M::munit || ix == 1);
        TMVAssert(!m.isunit() || ix == 1);
        typedef typename M::cview_type Mv;
        Mv mv = m.cView();
        ScaleU_Helper<-4,M::msize,ix,T,Mv>::call(x,mv); 
    }

    template <class T>
    inline void InstScale(const T x, LowerTriMatrixView<T,NonUnitDiag> m)
    { InstScale(x,m.transpose()); }


#undef TMV_Q1

} // namespace tmv

#endif
