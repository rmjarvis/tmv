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

#ifndef TMV_TransposeM_H
#define TMV_TransposeM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_SwapV.h"

namespace tmv {

    // Defined below:
    template <class M1>
    inline void TransposeSelf(BaseMatrix_Rec_Mutable<M1>& m);
    template <class M1>
    inline void InlineTransposeSelf(BaseMatrix_Rec_Mutable<M1>& m);

    // Defined in TMV_Matrix.cpp
    template <class T>
    void InstTransposeSelf(MatrixView<T> m);

    //
    // TransposeSelf
    //

    template <int algo, int size, class M1>
    struct TransposeSelf_Helper;

    // algo 1: Simple for loop
    template <int size, class M1>
    struct TransposeSelf_Helper<1,size,M1>
    {
        static inline void call(M1& m)
        {
            const int n = (size == UNKNOWN ? int(m.colsize()) : size);
            for(int i=1;i<n;++i) {
                typename M1::row_sub_type v1 = m.get_row(i,0,i);
                typename M1::col_sub_type v2 = m.get_col(i,0,i);
                NoAliasSwap(v1,v2);
            }
        }
    };

    // algo 2: The same thing, but with iterators.
    template <int size, class M1>
    struct TransposeSelf_Helper<2,size,M1>
    {
        static inline void call(M1& m)
        {
            const int n = size == UNKNOWN ? int(m.colsize()) : size;
            if (n <= 1) return;
            typedef typename M1::row_type Mr;
            typedef typename M1::col_type Mc;
            typedef typename Mr::iterator IT1;
            typedef typename Mc::iterator IT2;
            IT1 it1 = m.row(1).begin();
            IT2 it2 = m.col(1).begin();
            const int step1 = m.stepi();
            const int step2 = m.stepj();
            for(int i=1;i<n;++i) {
                SwapV_Helper<-3,UNKNOWN,Mr,Mc>::call2(i,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 3: The other way to do the loop.  
    // This way seems to be a little bit slower.
    template <int size, class M1>
    struct TransposeSelf_Helper<3,size,M1>
    {
        static inline void call(M1& m)
        {
            int n = size == UNKNOWN ? int(m.colsize()) : size;
            if (n <= 1) return;
            typedef typename M1::row_type Mr;
            typedef typename M1::col_type Mc;
            typedef typename Mr::iterator IT1;
            typedef typename Mc::iterator IT2;
            IT1 it1 = m.row(0).begin(); ++it1;
            IT2 it2 = m.col(0).begin(); ++it2;
            const int step = m.stepi() + m.stepj();
            for(--n;n;--n) {
                SwapV_Helper<-3,UNKNOWN,Mr,Mc>::call2(n,it1,it2);
                it1.shiftP(step);
                it2.shiftP(step);
            }
        }
    };

    // algo 5: Fully unroll
    template <int size, class M1>
    struct TransposeSelf_Helper<5,size,M1>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static inline void unroll(M1& m)
            {
                Unroller<I,M/2,0,I>::unroll(m);
                Unroller<I+M/2,M-M/2,0,I+M/2>::unroll(m);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static inline void unroll(M1& m)
            {
                Unroller<I,1,J,N/2>::unroll(m);
                Unroller<I,1,J+N/2,N-N/2>::unroll(m);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        { static inline void unroll(M1& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static inline void unroll(M1& m)
            { TMV_SWAP(m.ref(I,J) , m.ref(J,I) ); }
        };
        template <int I>
        struct Unroller<I,1,I,1>
        { static inline void unroll(M1& m) {} };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        { static inline void unroll(M1& ) {} };

        static inline void call(M1& m)
        { Unroller<0,size,0,size>::unroll(m); }
    };

    // algo -3: Determine which algorithm to use
    template <int size, class M1>
    struct TransposeSelf_Helper<-3,size,M1>
    {
        static inline void call(M1& m)
        {
            TMVStaticAssert(!M1::_conj);
            const int algo = 
#if TMV_OPT >= 1
                ( size != UNKNOWN && size < 8 ) ? 5 :
#endif
                2;
            TransposeSelf_Helper<algo,size,M1>::call(m);
        }
    };

    // algo 97: Conjugate
    template <int size, class M1>
    struct TransposeSelf_Helper<97,size,M1>
    {
        static inline void call(M1& m)
        {
            typedef typename M1::conjugate_type Mc;
            Mc mc = m.conjugate();
            TransposeSelf_Helper<-1,size,Mc>::call(mc);
        }
    };

    // algo 98: Call inst
    template <int size, class M1>
    struct TransposeSelf_Helper<98,size,M1>
    {
        static inline void call(M1& m)
        { InstTransposeSelf(m.xView()); }
    };

    // algo -1: Check for inst
    template <int size, class M1>
    struct TransposeSelf_Helper<-1,size,M1>
    {
        static inline void call(M1& m)
        {
            typedef typename M1::value_type T;
            const bool inst = 
                M1::unknownsizes &&
                Traits<T>::isinst;
            const bool conj = M1::_conj;
            const int algo = 
                conj ? 97 :
                inst ? 98 :
                -3;
            TransposeSelf_Helper<algo,size,M1>::call(m);
        }
    };

    template <class M1>
    inline void TransposeSelf(BaseMatrix_Rec_Mutable<M1>& m)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same)); 
        TMVAssert(m.colsize() == m.rowsize());
        const int size = Sizes<M1::_colsize,M1::_rowsize>::size;
        typedef typename M1::cview_type M1v;
        M1v m1v = m.cView();
        TransposeSelf_Helper<-1,size,M1v>::call(m1v);
    }

    template <class M1>
    inline void InlineTransposeSelf(BaseMatrix_Rec_Mutable<M1>& m)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same)); 
        TMVAssert(m.colsize() == m.rowsize());
        const int size = Sizes<M1::_colsize,M1::_rowsize>::size;
        typedef typename M1::cview_type M1v;
        M1v m1v = m.cView();
        TransposeSelf_Helper<-3,size,M1v>::call(m1v);
    }

} // namespace tmv

#endif
