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


#ifndef TMV_InvertU_H
#define TMV_InvertU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_MultUV.h"
#include "TMV_MultUM.h"
#include "TMV_InvertD.h"
#include "tmv/TMV_Det.h"
#include "tmv/TMV_TriMatrixIO.h"
#include "tmv/TMV_MatrixIO.h"

#ifdef PRINTALGO_InvU
#include <iostream>
#endif

// Use the specialized 1,2,3,4 sized algorithms for the end of the
// recursive algorithm.
#if TMV_OPT >= 2
#define TMV_OPT_CLEANUP
#endif

// Q1 is the maximum nops to unroll
#if TMV_OPT >= 3
#define TMV_Q1 200 
#elif TMV_OPT >= 2
#define TMV_Q1 25
#elif TMV_OPT >= 1
#define TMV_Q1 9
#else
#define TMV_Q1 0
#endif

// Q2 is the minimum size to keep recursing.
#define TMV_Q2 8

namespace tmv {

    // Defined below:
    template <class M>
    inline void InvertSelf(BaseMatrix_Tri_Mutable<M>& m);
    template <class M>
    inline void InlineInvertSelf(BaseMatrix_Tri_Mutable<M>& m);

    // Defined in TMV_InvertU.cpp
    template <class T>
    void InstInvertSelf(UpperTriMatrixView<T> m);
    template <class T>
    void InstInvertSelf(LowerTriMatrixView<T> m);

    //
    // U = U^-1
    //

    template <int algo, int s, class M>
    struct InvertU_Helper;

    // algo 0: nothing to do (either s==0, or invert diag of unitdiag m.)
    template <int s, class M>
    struct InvertU_Helper<0,s,M>
    { static void call(M& ) {} };

    // algo 1: s == 1, so simplifies to a scalar quotient
    template <class M>
    struct InvertU_Helper<1,1,M>
    {
        // m.diag is already inverted, so nothing to do.
        static void call(M& ) {}
    };

    // algo 2: transpose
    template <int s, class M>
    struct InvertU_Helper<2,s,M>
    {
        static void call(M& m)
        { 
#ifdef PRINTALGO_InvU
            int N = (s == UNKNOWN ? m.size() : s);
            std::cout<<"InvU algo 2: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            typedef typename M::transpose_type Mt;
            Mt mt = m.transpose();
            InvertU_Helper<-3,s,Mt>::call(mt);
        }
    };

    // algo 11: column major loop 
    template <int s, class M>
    struct InvertU_Helper<11,s,M>
    {
        static void call(M& m)
        {
            int N = (s == UNKNOWN ? m.size() : s);
#ifdef PRINTALGO_InvU
            std::cout<<"InvU algo 11: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            const bool u = M::munit;
            typedef typename M::col_sub_type Mc;
            typedef typename M::const_col_sub_type Mcc;
            typedef typename M::const_subtrimatrix_type Mst;
            const int ix = u?-1:0;
            typedef typename M::value_type T;
            typedef typename M::real_type RT;
            typedef typename TypeSelect<u,RT,T>::type XT;
            const int xx = UNKNOWN;

            for(int j=0;j<N;++j) {
                // m.col(j,0,j) = -m.subTraiMatrix(0,j)^-1 *
                //                   m.col(j,0,j) * m(j,j)^-1
                Mc mj = m.get_col(j,0,j);
                Mst mst = m.cSubTriMatrix(0,j); 
                Scaling<ix,XT> invAjj(-m.cref(j,j));
                MultUV_Helper<-4,xx,false,ix,XT,Mst,Mcc,Mc>::call(
                    invAjj,mst,mj,mj);
            }
        }
    };

    // algo 12: row major loop
    template <int s, class M>
    struct InvertU_Helper<12,s,M>
    {
        static void call(M& m)
        {
            int N = (s == UNKNOWN ? m.size() : s);
#ifdef PRINTALGO_InvU
            std::cout<<"InvU algo 11: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            const bool u = M::munit;
            typedef typename M::row_sub_type Mr;
            typedef typename M::const_row_sub_type Mcr;
            typedef typename M::const_subtrimatrix_type Mst;
            typedef typename Mst::const_transpose_type Mstct;
            const int ix = u?-1:0;
            typedef typename M::value_type T;
            typedef typename M::real_type RT;
            typedef typename TypeSelect<u,RT,T>::type XT;
            const int xx = UNKNOWN;

            for(int i=N-1;i>=0;--i) {
                // m.row(i,i+1,N) = -(m(i,i)^-1) * m.row(i,i+1,N) * 
                //                      m.subTraiMatrix(i+1,N)^-1
                Mr mi = m.get_row(i,i+1,N);
                Mstct mstt = m.cSubTriMatrix(i+1,N).transpose();
                Scaling<ix,XT> invAii(-m.cref(i,i));
                MultUV_Helper<-4,xx,false,ix,XT,Mstct,Mcr,Mr>::call(
                    invAii,mstt,mi,mi);
            }
        }
    };

    // algo 16: UpperTri: Unroll small case
    template <int s, class M>
    struct InvertU_Helper<16,s,M>
    {
        template <int I, int N>
        struct Unroller
        {
            static void unroll(M& m)
            {
                // [ B00 B01 ] * [ A00 A01 ] = [ 1 0 ]
                // [  0  B11 ]   [  0  A11 ]   [ 0 1 ]
                // B00 A00 = 1
                // B00 A01 + B01 A11 = 0
                // B11 A11 = 1
                //
                // B11 = A11^-1
                // B00 = A00^-1
                // B01 = -B00 A01 B11

                const int Nx = N/2;
                const int Ny = N-Nx;
                const int I1 = I+Nx;
                const int I2 = I+N;
                typedef typename M::submatrix_type Msm;
                typedef typename M::subtrimatrix_type Mst;
                typedef typename M::const_submatrix_type Msmc;
                typedef typename M::const_subtrimatrix_type Mstc;
                typedef typename Msm::transpose_type Msmt;
                typedef typename Mst::const_transpose_type Mstct;
                typedef typename Msm::const_transpose_type Msmct;
                typedef typename M::real_type RT;

                Mst A00 = m.cSubTriMatrix(I,I1);
                Msm A01 = m.cSubMatrix(I,I1,I1,I2);
                Msmt A01t = A01.transpose();
                Msmct A01ct = A01.transpose();
                Mstct A11t = m.cSubTriMatrix(I1,I2).transpose();

                // B00 = A00^-1
                Unroller<I,Nx>::unroll(m);

                // B11 = A11^-1
                Unroller<I1,Ny>::unroll(m);

                // B01 = -B00 A01
                Scaling<-1,RT> mone;
                MultUM_Helper<-4,Nx,Ny,false,-1,RT,Mstc,Msmc,Msm>::call(
                    mone,A00,A01,A01);

                // B01 = B01 B11
                Scaling<1,RT> one;
                MultUM_Helper<-4,Ny,Nx,false,1,RT,Mstct,Msmct,Msmt>::call(
                    one,A11t,A01ct,A01t);
            }
        };
        template <int I>
        struct Unroller<I,1> // diagonal is already done, so nothing to do.
        { static inline void unroll(M&) {} };
        template <int I>
        struct Unroller<I,0>
        { static inline void unroll(M&) {} };
        static void call(M& m)
        {
#ifdef PRINTALGO_InvU
            const int N = m.size();
            std::cout<<"InvU algo 16: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            Unroller<0,s>::unroll(m); 
        }
    };
    template <class M>
    struct InvertU_Helper<16,UNKNOWN,M>
    {
        static void call(M& m)
        {
            const int N = m.size();
#ifdef PRINTALGO_InvU
            std::cout<<"InvU algo 16: N,s,x = "<<N<<','<<UNKNOWN<<std::endl;
#endif
            const int algo2 = M::mrowmajor ? 12 : 11;
            switch (N) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   InvertU_Helper<1,1,M>::call(m);
                   break;
              case 2 :
                   InvertU_Helper<16,2,M>::call(m);
                   break;
              case 3 :
                   InvertU_Helper<16,3,M>::call(m);
                   break;
              default :
                   InvertU_Helper<algo2,UNKNOWN,M>::call(m);
            }
        }
    };

    // algo 17: Split U into 3 sections and recurse 
    // TODO: Combine these the way I do for MultUU, etc.
    template <int s, class M>
    struct InvertU_Helper<17,s,M>
    {
        template <int which, int dummy>
        struct Helper2;

        template <int dummy>
        struct Helper2<0,dummy> // s == UNKNOWN
        {
            static void call(M& m)
            {
                const int N = m.size();
#ifdef TMV_OPT_CLEANUP
                const int algo2 = 16;
#else
                const int algo2 = M::mrowmajor ? 12 : 11;
#endif
                const int xx = UNKNOWN;
                if (N < TMV_Q2) {
                    InvertU_Helper<algo2,xx,M>::call(m);
                } else {
                    // [ B00 B01 ] * [ A00 A01 ] = [ 1 0 ]
                    // [  0  B11 ]   [  0  A11 ]   [ 0 1 ]
                    // B00 A00 = 1
                    // B00 A01 + B01 A11 = 0
                    // B11 A11 = 1
                    //
                    // B11 = A11^-1
                    // B00 = A00^-1
                    // B01 = -B00 A01 B11

                    const int Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                    // (If N > 16, round N/2 up to a multiple of 16.)
                    
                    typedef typename M::submatrix_type Msm;
                    typedef typename M::subtrimatrix_type Mst;
                    typedef typename Msm::transpose_type Msmt;
                    typedef typename M::const_submatrix_type Msmc;
                    typedef typename M::const_subtrimatrix_type Mstc;
                    typedef typename Mst::const_transpose_type Mstct;
                    typedef typename Msm::const_transpose_type Msmct;
                    typedef typename M::real_type RT;

                    Mst A00 = m.cSubTriMatrix(0,Nx);
                    Msm A01 = m.cSubMatrix(0,Nx,Nx,N);
                    Msmt A01t = A01.transpose();
                    Mst A11 = m.cSubTriMatrix(Nx,N);
                    Mstc A00c = A00;
                    Msmct A01ct = A01.transpose();
                    Mstct A11ct = A11.transpose();

                    // B00 = A00^-1
                    InvertU_Helper<17,xx,Mst>::call(A00);

                    // B11 = A11^-1
                    InvertU_Helper<17,xx,Mst>::call(A11);

                    // B01 = -B00 A01
                    Scaling<-1,RT> mone;
                    MultUM_Helper<-4,xx,xx,false,-1,RT,Mstc,Msmc,Msm>::call(
                        mone,A00c,A01,A01);

                    // B01 = B01 B11
                    Scaling<1,RT> one;
                    MultUM_Helper<-4,xx,xx,false,1,RT,Mstct,Msmct,Msmt>::call(
                        one,A11ct,A01ct,A01t);
                }
            }
        };
        template <int dummy>
        struct Helper2<1,dummy> // s < TMV_Q2
        {
            static void call(M& m)
            {
                const int s2 = s > 20 ? UNKNOWN : s;
                const int s2p1 = IntTraits<s2>::Sp1;
                const int s2p2 = IntTraits<s2p1>::Sp1;
                // nops = 1/6 n(n+1)(n+2)
                const int nops =
                    IntTraits2<IntTraits2<s2,s2p1>::safeprod,s2p2>::safeprod/6;
                const bool unroll =
                    s == UNKNOWN ? false :
                    nops > TMV_Q1 ? false :
                    s <= 10;
                const int algo2 =
                    s == 0 ? 0 :
                    s == 1 ? ( M::munit ? 0 : 1 ) :
                    unroll ? 16 :
                    M::mrowmajor ? 12 : 11;
                InvertU_Helper<algo2,s,M>::call(m);
            }
        };
        template <int dummy>
        struct Helper2<2,dummy> // s >= TMV_Q2
        {
            static void call(M& m)
            {
                const int N = s;
                const int Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round N/2 up to a multiple of 16.)
                const int Ny = N-Nx;

                typedef typename M::submatrix_type Msm;
                typedef typename M::subtrimatrix_type Mst;
                typedef typename Msm::transpose_type Msmt;
                typedef typename M::const_submatrix_type Msmc;
                typedef typename M::const_subtrimatrix_type Mstc;
                typedef typename Mst::const_transpose_type Mstct;
                typedef typename Msm::const_transpose_type Msmct;
                typedef typename M::real_type RT;

                Mst A00 = m.cSubTriMatrix(0,Nx);
                Msm A01 = m.cSubMatrix(0,Nx,Nx,N);
                Msmt A01t = A01.transpose();
                Mst A11 = m.cSubTriMatrix(Nx,N);
                Mstc A00c = A00;
                Msmct A01ct = A01.transpose();
                Mstct A11ct = A11.transpose();

                // B00 = A00^-1
                // F(n)
                InvertU_Helper<17,Nx,Mst>::call(A00);

                // B11 = A11^-1
                // F(n)
                InvertU_Helper<17,Ny,Mst>::call(A11);

                // B01 = -B00 A01
                // n*n*(n+1)/2
                Scaling<-1,RT> mone;
                MultUM_Helper<-4,Nx,Ny,false,-1,RT,Mstc,Msmc,Msm>::call(
                    mone,A00c,A01,A01);

                // B01 = B01 B11
                // n*n*(n+1)/2
                Scaling<1,RT> one;
                MultUM_Helper<-4,Ny,Nx,false,1,RT,Mstct,Msmct,Msmt>::call(
                    one,A11ct,A01ct,A01t);
            }
        };
        static void call(M& m)
        {
#ifdef PRINTALGO_InvU
            const int N = m.size();
            std::cout<<"InvU algo 17: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            const int s2 = (s >= 0 && s <= 64) ? s : UNKNOWN;
            const int which = s2 == UNKNOWN ? 0 : s2 < TMV_Q2 ? 1 : 2;
            Helper2<which,1>::call(m);
        }
    };

    // algo 50: Invert the diagonal.
    // I invert the diagonal first, since it is more efficient to 
    // use the sse commands if possible.
    template <int s, class M>
    struct InvertU_Helper<50,s,M>
    {
        static void call(M& m)
        {
            TMVStaticAssert(M::mupper);
            TMVStaticAssert(!M::munit);
            typedef typename M::diag_type Md;
            Md md = m.diag();
            ElemInvert_Helper<-4,s,Md>::call(md);
        }
    };

    // algo 51: Maybe invert the diagonal.  (m is UnknownDiag)
    template <int s, class M>
    struct InvertU_Helper<51,s,M>
    {
        static void call(M& m)
        {
            if (!m.isunit()) 
                InvertU_Helper<50,s,M>::call(m); 
        }
    };

    // algo -4: No branches or copies
    template <int s, class M>
    struct InvertU_Helper<-4,s,M>
    {
        static void call(M& m)
        {
            TMVStaticAssert(M::mupper);
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            const int s2p2 = IntTraits<s2p1>::Sp1;
            // nops = 1/6 n(n+1)(n+2)
            const int nops =
                IntTraits2<IntTraits2<s2,s2p1>::safeprod,s2p2>::safeprod / 6;
            const bool unroll = 
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const int algo1 = 
                M::munit ? 0 : M::munknowndiag ? 51 : 50;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? ( M::munit ? 0 : 1 ) :
                unroll ? 16 :
                17;
#ifdef PRINTALGO_InvU
            std::cout<<"InlineInvertU (algo -4): \n";
            std::cout<<"m "<<TMV_Text(m)<<std::endl;
            std::cout<<"s,algo = "<<s<<"  "<<algo1<<"  "<<algo2<<std::endl;
#endif
            InvertU_Helper<algo1,s,M>::call(m);
            InvertU_Helper<algo2,s,M>::call(m);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class M>
    struct InvertU_Helper<-3,s,M>
    {
        static void call(M& m)
        {
            // Possible algorithms to choose from:
            //
            // Trivial:
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar quotient
            //  2 = transpose
            //
            // UpperTri:
            // 11 = column major
            // 12 = row major
            // 16 = unroll
            // 17 = recurse

            const int algo = 
                ( s == 0 ) ? 0 :
                M::mlower ? 2 :
                -4;
#ifdef PRINTALGO_InvU
            std::cout<<"InlineInvertU: \n";
            std::cout<<"m "<<TMV_Text(m)<<std::endl;
            std::cout<<"s,algo = "<<s<<"  "<<algo<<std::endl;
#endif
            InvertU_Helper<algo,s,M>::call(m);
        }
    };

    // algo 97: Conjugate
    template <int s, class M>
    struct InvertU_Helper<97,s,M>
    {
        static void call(M& m)
        { 
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            InvertU_Helper<-2,s,Mc>::call(mc);
        }
    };

    // algo 98: call InstInvertSelf
    template <int s, class M>
    struct InvertU_Helper<98,s,M>
    {
        static void call(M& m)
        { InstInvertSelf(m.xdView()); }
    };

    // algo -2: Check for inst
    template <int s, class M>
    struct InvertU_Helper<-2,s,M>
    {
        static void call(M& m)
        {
            typedef typename M::value_type T;
            const bool inst = 
                M::msize == UNKNOWN && 
                Traits<T>::isinst;
            const int algo = 
                ( s == 0 ) ? 0 :
                M::mconj ? 97 :
                inst ? 98 : 
                -3;
            InvertU_Helper<algo,s,M>::call(m);
        }
    };

    // algo -1: Check for aliases? No.
    template <int s, class M>
    struct InvertU_Helper<-1,s,M>
    {
        static void call(M& m)
        { InvertU_Helper<-2,s,M>::call(m); }
    };

    template <int algo, class M>
    inline void DoInvertSelf(BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        Mv mv = m.cView();
        if (m.isSingular()) {
#ifdef TMV_NO_THROW
            std::cerr<<"Singular TriMatrix found in InvertSelf\n";
            exit(1);
#else
            throw SingularMatrix<M>(m.mat());
#endif
        }
        InvertU_Helper<algo,M::msize,Mv>::call(mv);
    }

    template <class M>
    inline void InvertSelf(BaseMatrix_Tri_Mutable<M>& m)
    { DoInvertSelf<-2>(m); }

    template <class M>
    inline void InlineInvertSelf(BaseMatrix_Tri_Mutable<M>& m)
    { DoInvertSelf<-3>(m); }

#undef TMV_Q1
#undef TMV_Q2
#undef TMV_Q3

} // namespace tmv

#endif 
