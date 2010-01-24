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

//
// This file defines the basic composite type for a product of a 
// matrix and a scalar.  It also implements the calculation for 
// dense rectangular matrices.

#ifndef TMV_ScaleM_H
#define TMV_ScaleM_H

#include "TMV_ScaleV.h"
#include "TMV_MultXV.h"
#include "TMV_BaseMatrix_Rec.h"

namespace tmv {

    // Defined below:
    template <int ix, class T, class M>
    inline void Scale(
        const Scaling<ix,T>& x, BaseMatrix_Rec_Mutable<M>& m);
    template <int ix, class T, class M>
    inline void InlineScale(
        const Scaling<ix,T>& x, BaseMatrix_Rec_Mutable<M>& m);

    // Defined in TMV_ScaleM.cpp
    template <class T>
    void InstScale(const T x, MatrixView<T> m);

    //
    // Matrix *= x
    //

    template <int algo, int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper;

    // algo 0: trivial: cs == 0, rs == 0 or ix == 1, so nothing to do
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<0,cs,rs,ix,T,M1>
    { static inline void call(const Scaling<ix,T>& , M1& ) {} };

    // algo 1: Linearize to vector version
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<1,cs,rs,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::linearview_type Ml;
            Ml ml = m.linearView();
            Scale(x,ml);
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<11,cs,rs,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            const int M = cs == UNKNOWN ? int(m.colsize()) : cs;
            int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
            typedef typename M1::col_type Mc;
            typedef typename Mc::iterator IT;
            const int step = m.stepj();
            IT it = m.get_col(0).begin();
            for(;N;--N) {
                ScaleV_Helper<-4,cs,ix,T,Mc>::call2(M,x,it);
                it.shiftP(step);
            }
        }
    };

    // algo 21: Loop over rows
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<21,cs,rs,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            int M = cs == UNKNOWN ? int(m.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
            typedef typename M1::row_type Mr;
            typedef typename Mr::iterator IT;
            const int step = m.stepi();
            IT it = m.get_row(0).begin();
            for(;M;--M) {
                ScaleV_Helper<-4,rs,ix,T,Mr>::call2(N,x,it);
                it.shiftP(step);
            }
        }
    };

    // algo 31: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<31,cs,rs,ix,T,M1>
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            if (m.canLinearize()) 
                ScaleM_Helper<1,cs,rs,ix,T,M1>::call(x,m);
            else if ( m.isrm() ||
                      (!m.iscm() && (m.colsize() > m.rowsize()) ) )
                ScaleM_Helper<21,cs,rs,ix,T,M1>::call(x,m);
            else
                ScaleM_Helper<11,cs,rs,ix,T,M1>::call(x,m);
        }
    };

    // algo 25: Fully unroll by rows
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<25,cs,rs,ix,T,M1>
    {
        template <int I, int M, int J, int N, bool iscomplex>
        struct Unroller
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,M/2,J,N,iscomplex>::unroll(x,m);
                Unroller<I+M/2,M-M/2,J,N,iscomplex>::unroll(x,m);
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
        { Unroller<0,cs,0,rs,M1::miscomplex>::unroll(x,m); }
    };

    // algo 15: Fully unroll by columns
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<15,cs,rs,ix,T,M1>
    {
        template <int I, int M, int J, int N, bool iscomplex>
        struct Unroller
        {
            static inline void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,M,J,N/2,iscomplex>::unroll(x,m);
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
        { Unroller<0,cs,0,rs,M1::miscomplex>::unroll(x,m); }
    };

    // algo -4: No copies or branches
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<-4,cs,rs,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::value_type T1;
            const int algo = 
                (cs == 0 || rs == 0) ? 0 :
                (ix == 1) ? 1 :
                M1::mcanlin ? 2 :
#if TMV_OPT >= 1
                ( cs != UNKNOWN && rs != UNKNOWN ) ? (
                    ( IntTraits2<cs,rs>::prod <= int(128/sizeof(T1)) ) ? (
                        ( M1::mrowmajor ? 25 : 15 ) ) :
                    M1::mrowmajor ? 21 :
                    M1::mcolmajor ? 11 :
                    ( cs > rs ) ? 21 : 11 ) :
                M1::mrowmajor ? 21 :
#endif
                11;
            ScaleM_Helper<algo,cs,rs,ix,T,M1>::call(x,m);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<-3,cs,rs,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            const int algo = 
                (cs == 0 || rs == 0 || ix == 1) ? 0 :
                M1::mcanlin ? 1 :
#if TMV_OPT >= 2
                cs == UNKNOWN || rs == UNKNOWN ? 31 :
#endif
                -4;
            ScaleM_Helper<algo,cs,rs,ix,T,M1>::call(x,m);
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<97,cs,rs,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::conjugate_type Mc;
            Mc mc = m.conjugate();
            ScaleM_Helper<-2,cs,rs,ix,T,Mc>::call(TMV_CONJ(x),mc);
        }
    };

    // algo 98: Call inst
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<98,cs,rs,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            typename M1::value_type xx(x);
            InstScale(xx,m.xView());
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<-2,cs,rs,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::value_type T1;
            const bool inst =
                M1::mcolsize == UNKNOWN && M1::mrowsize == UNKNOWN &&
                Traits<T1>::isinst;
            const bool conj = M1::mconj;
            const int algo = 
                ix == 1 ? 0 :
                conj ? 97 :
                inst ? 98 : 
                -3;
            ScaleM_Helper<algo,cs,rs,ix,T,M1>::call(x,m);
        }
    };

    // algo -1: Check alias? No.
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<-1,cs,rs,ix,T,M1> 
    {
        static inline void call(const Scaling<ix,T>& x, M1& m)
        { ScaleM_Helper<-2,cs,rs,ix,T,M1>::call(x,m); }
    };

    template <int ix, class T, class M>
    inline void Scale(
        const Scaling<ix,T>& x, BaseMatrix_Rec_Mutable<M>& m)
    {
        const int cs = M::mcolsize;
        const int rs = M::mrowsize;
        typedef typename M::cview_type Mv;
        Mv mv = m.cView();
        ScaleM_Helper<-2,cs,rs,ix,T,Mv>::call(x,mv);
    }

    template <int ix, class T, class M>
    inline void InlineScale(
        const Scaling<ix,T>& x, BaseMatrix_Rec_Mutable<M>& m)
    {
        const int cs = M::mcolsize;
        const int rs = M::mrowsize;
        typedef typename M::cview_type Mv;
        Mv mv = m.cView();
        ScaleM_Helper<-3,cs,rs,ix,T,Mv>::call(x,mv);
    }

} // namespace tmv

#endif
