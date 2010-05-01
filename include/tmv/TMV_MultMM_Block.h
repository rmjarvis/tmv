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

#ifndef TMV_MultMM_BLOCK_H
#define TMV_MultMM_BLOCK_H

#include "TMV_MultMM_Kernel.h"

//#define TEST_POINTERS

namespace tmv {

    // Defined in TMV_MultMM_Block.cpp
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM_Block(
        const T3 x, 
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM_Block(
        const T3 x, 
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM_RecursiveBlock(
        const T3 x, 
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM_RecursiveBlock(
        const T3 x, 
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3);

    // Defined below:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_Block(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM_Block(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_RecursiveBlock(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM_RecursiveBlock(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);


#ifdef TEST_POINTERS
    void* glob_m1p_end;
    void* glob_m2p_end;
    void* glob_m3p_end;
#endif

    // This next helper function determines which K_known function to use 
    // for cleanup of the parts that don't fully fit into a K block.

    // If xs is known, then this is easy, since we know the right 
    // function to use at compile time.
    template <int K2, int KB, class T> 
    struct get_Kcleanup
    {
        typedef void Kcleanup(
            const int M, const int N, const int K,
            const Scaling<1,T>& x, const T* A, const T* B, T* C);
        static inline Kcleanup* call(const int K)
        {
            TMVStaticAssert(K2 != UNKNOWN);
            return &call_multmm_16_16_K_known<K2,1,T>;
        }
    };

#ifndef TMV_OPT_CLEANUP
    template <int KB, class T>
    struct get_Kcleanup<UNKNOWN,KB,T>
    {
        typedef void Kcleanup(
            const int M, const int N, const int K,
            const Scaling<1,T>& x, const T* A, const T* B, T* C);
        static inline Kcleanup* call(const int K)
        { return &call_multmm_16_16_K<1,T>; }
    };
#else
    // For unknown xs, we need to figure out K at runtime.
    // We need a version for each different value of KB, so that 
    // the switch statement has the right number of entries.
    // I probably don't need all of these, but I have versions for each
    // power of 2 from 2 to 64, so it is easy to switch KB values without
    // worrying about whether I have the corresponding get_Kcleanup function.
    template <class T>
    struct get_Kcleanup<UNKNOWN,64,T>
    {
        typedef void Kcleanup(
            const int M, const int N, const int K,
            const Scaling<1,T>& x, const T* A, const T* B, T* C);
        static inline Kcleanup* call(const int K)
        {
            TMVAssert(K > 0);
            TMVAssert(K < 64);
            switch(K) 
            {
              case 1 : return &call_multmm_16_16_K_known<1,1,T>; 
                       break;
              case 2 : return &call_multmm_16_16_K_known<2,1,T>;
                       break;
              case 3 : return &call_multmm_16_16_K_known<3,1,T>;
                       break;
              case 4 : return &call_multmm_16_16_K_known<4,1,T>;
                       break;
              case 5 : return &call_multmm_16_16_K_known<5,1,T>;
                       break;
              case 6 : return &call_multmm_16_16_K_known<6,1,T>;
                       break;
              case 7 : return &call_multmm_16_16_K_known<7,1,T>;
                       break;
              case 8 : return &call_multmm_16_16_K_known<8,1,T>;
                       break;
              case 9 : return &call_multmm_16_16_K_known<9,1,T>;
                       break;
              case 10 : return &call_multmm_16_16_K_known<10,1,T>;
                        break;
              case 11 : return &call_multmm_16_16_K_known<11,1,T>;
                        break;
              case 12 : return &call_multmm_16_16_K_known<12,1,T>;
                        break;
              case 13 : return &call_multmm_16_16_K_known<13,1,T>;
                        break;
              case 14 : return &call_multmm_16_16_K_known<14,1,T>;
                        break;
              case 15 : return &call_multmm_16_16_K_known<15,1,T>;
                        break;
              case 16 : return &call_multmm_16_16_K_known<16,1,T>;
                        break;
              case 17 : return &call_multmm_16_16_K_known<17,1,T>;
                        break;
              case 18 : return &call_multmm_16_16_K_known<18,1,T>;
                        break;
              case 19 : return &call_multmm_16_16_K_known<19,1,T>;
                        break;
              case 20 : return &call_multmm_16_16_K_known<20,1,T>;
                        break;
              case 21 : return &call_multmm_16_16_K_known<21,1,T>;
                        break;
              case 22 : return &call_multmm_16_16_K_known<22,1,T>;
                        break;
              case 23 : return &call_multmm_16_16_K_known<23,1,T>;
                        break;
              case 24 : return &call_multmm_16_16_K_known<24,1,T>;
                        break;
              case 25 : return &call_multmm_16_16_K_known<25,1,T>;
                        break;
              case 26 : return &call_multmm_16_16_K_known<26,1,T>;
                        break;
              case 27 : return &call_multmm_16_16_K_known<27,1,T>;
                        break;
              case 28 : return &call_multmm_16_16_K_known<28,1,T>;
                        break;
              case 29 : return &call_multmm_16_16_K_known<29,1,T>;
                        break;
              case 30 : return &call_multmm_16_16_K_known<30,1,T>;
                        break;
              case 31 : return &call_multmm_16_16_K_known<31,1,T>;
                        break;
              case 32 : return &call_multmm_16_16_K_known<32,1,T>;
                        break;
              case 33 : return &call_multmm_16_16_K_known<33,1,T>;
                        break;
              case 34 : return &call_multmm_16_16_K_known<34,1,T>;
                        break;
              case 35 : return &call_multmm_16_16_K_known<35,1,T>;
                        break;
              case 36 : return &call_multmm_16_16_K_known<36,1,T>;
                        break;
              case 37 : return &call_multmm_16_16_K_known<37,1,T>;
                        break;
              case 38 : return &call_multmm_16_16_K_known<38,1,T>;
                        break;
              case 39 : return &call_multmm_16_16_K_known<39,1,T>;
                        break;
              case 40 : return &call_multmm_16_16_K_known<40,1,T>;
                        break;
              case 41 : return &call_multmm_16_16_K_known<41,1,T>;
                        break;
              case 42 : return &call_multmm_16_16_K_known<42,1,T>;
                        break;
              case 43 : return &call_multmm_16_16_K_known<43,1,T>;
                        break;
              case 44 : return &call_multmm_16_16_K_known<44,1,T>;
                        break;
              case 45 : return &call_multmm_16_16_K_known<45,1,T>;
                        break;
              case 46 : return &call_multmm_16_16_K_known<46,1,T>;
                        break;
              case 47 : return &call_multmm_16_16_K_known<47,1,T>;
                        break;
              case 48 : return &call_multmm_16_16_K_known<48,1,T>;
                        break;
              case 49 : return &call_multmm_16_16_K_known<49,1,T>;
                        break;
              case 50 : return &call_multmm_16_16_K_known<50,1,T>;
                        break;
              case 51 : return &call_multmm_16_16_K_known<51,1,T>;
                        break;
              case 52 : return &call_multmm_16_16_K_known<52,1,T>;
                        break;
              case 53 : return &call_multmm_16_16_K_known<53,1,T>;
                        break;
              case 54 : return &call_multmm_16_16_K_known<54,1,T>;
                        break;
              case 55 : return &call_multmm_16_16_K_known<55,1,T>;
                        break;
              case 56 : return &call_multmm_16_16_K_known<56,1,T>;
                        break;
              case 57 : return &call_multmm_16_16_K_known<57,1,T>;
                        break;
              case 58 : return &call_multmm_16_16_K_known<58,1,T>;
                        break;
              case 59 : return &call_multmm_16_16_K_known<59,1,T>;
                        break;
              case 60 : return &call_multmm_16_16_K_known<60,1,T>;
                        break;
              case 61 : return &call_multmm_16_16_K_known<61,1,T>;
                        break;
              case 62 : return &call_multmm_16_16_K_known<62,1,T>;
                        break;
              case 63 : return &call_multmm_16_16_K_known<63,1,T>;
                        break;
            }
            // Never gets here.
            return 0;
        }
    };
    template <class T>
    struct get_Kcleanup<UNKNOWN,32,T>
    {
        typedef void Kcleanup(
            const int M, const int N, const int K,
            const Scaling<1,T>& x, const T* A, const T* B, T* C);
        static inline Kcleanup* call(const int K)
        {
            TMVAssert(K > 0);
            TMVAssert(K < 32);
            switch(K) 
            {
              case 1 : return &call_multmm_16_16_K_known<1,1,T>; 
                       break;
              case 2 : return &call_multmm_16_16_K_known<2,1,T>;
                       break;
              case 3 : return &call_multmm_16_16_K_known<3,1,T>;
                       break;
              case 4 : return &call_multmm_16_16_K_known<4,1,T>;
                       break;
              case 5 : return &call_multmm_16_16_K_known<5,1,T>;
                       break;
              case 6 : return &call_multmm_16_16_K_known<6,1,T>;
                       break;
              case 7 : return &call_multmm_16_16_K_known<7,1,T>;
                       break;
              case 8 : return &call_multmm_16_16_K_known<8,1,T>;
                       break;
              case 9 : return &call_multmm_16_16_K_known<9,1,T>;
                       break;
              case 10 : return &call_multmm_16_16_K_known<10,1,T>;
                        break;
              case 11 : return &call_multmm_16_16_K_known<11,1,T>;
                        break;
              case 12 : return &call_multmm_16_16_K_known<12,1,T>;
                        break;
              case 13 : return &call_multmm_16_16_K_known<13,1,T>;
                        break;
              case 14 : return &call_multmm_16_16_K_known<14,1,T>;
                        break;
              case 15 : return &call_multmm_16_16_K_known<15,1,T>;
                        break;
              case 16 : return &call_multmm_16_16_K_known<16,1,T>;
                        break;
              case 17 : return &call_multmm_16_16_K_known<17,1,T>;
                        break;
              case 18 : return &call_multmm_16_16_K_known<18,1,T>;
                        break;
              case 19 : return &call_multmm_16_16_K_known<19,1,T>;
                        break;
              case 20 : return &call_multmm_16_16_K_known<20,1,T>;
                        break;
              case 21 : return &call_multmm_16_16_K_known<21,1,T>;
                        break;
              case 22 : return &call_multmm_16_16_K_known<22,1,T>;
                        break;
              case 23 : return &call_multmm_16_16_K_known<23,1,T>;
                        break;
              case 24 : return &call_multmm_16_16_K_known<24,1,T>;
                        break;
              case 25 : return &call_multmm_16_16_K_known<25,1,T>;
                        break;
              case 26 : return &call_multmm_16_16_K_known<26,1,T>;
                        break;
              case 27 : return &call_multmm_16_16_K_known<27,1,T>;
                        break;
              case 28 : return &call_multmm_16_16_K_known<28,1,T>;
                        break;
              case 29 : return &call_multmm_16_16_K_known<29,1,T>;
                        break;
              case 30 : return &call_multmm_16_16_K_known<30,1,T>;
                        break;
              case 31 : return &call_multmm_16_16_K_known<31,1,T>;
                        break;
            }
            // Never gets here.
            return 0;
        }
    };
    template <class T>
    struct get_Kcleanup<UNKNOWN,16,T>
    {
        typedef void Kcleanup(
            const int M, const int N, const int K,
            const Scaling<1,T>& x, const T* A, const T* B, T* C);
        static inline Kcleanup* call(const int K)
        {
            TMVAssert(K > 0);
            TMVAssert(K < 16);
            switch(K) 
            {
              case 1 : return &call_multmm_16_16_K_known<1,1,T>; 
                       break;
              case 2 : return &call_multmm_16_16_K_known<2,1,T>;
                       break;
              case 3 : return &call_multmm_16_16_K_known<3,1,T>;
                       break;
              case 4 : return &call_multmm_16_16_K_known<4,1,T>;
                       break;
              case 5 : return &call_multmm_16_16_K_known<5,1,T>;
                       break;
              case 6 : return &call_multmm_16_16_K_known<6,1,T>;
                       break;
              case 7 : return &call_multmm_16_16_K_known<7,1,T>;
                       break;
              case 8 : return &call_multmm_16_16_K_known<8,1,T>;
                       break;
              case 9 : return &call_multmm_16_16_K_known<9,1,T>;
                       break;
              case 10 : return &call_multmm_16_16_K_known<10,1,T>;
                        break;
              case 11 : return &call_multmm_16_16_K_known<11,1,T>;
                        break;
              case 12 : return &call_multmm_16_16_K_known<12,1,T>;
                        break;
              case 13 : return &call_multmm_16_16_K_known<13,1,T>;
                        break;
              case 14 : return &call_multmm_16_16_K_known<14,1,T>;
                        break;
              case 15 : return &call_multmm_16_16_K_known<15,1,T>;
                        break;
            }
            // Never gets here.
            return 0;
        }
    };
#endif

    // For these algorithms we use a little trick for complex matrices to 
    // enable the use of the optimized kernels, which are only defined for
    // real-valued matrices.
    // The SSE commands don't work out quite as nicely for complex operations,
    // so we instead recast a complex temporary matrix as two successive real 
    // matrices, one for the real part and one for the imaginary part.
    // The main algorithm is written in a manner that pretends this isn't 
    // necessary.
    // Then all the operations on the matrices are done through these
    // next few templates, which implement the funny storage scheme.

    template <bool iscomplex, int cs, int rs, class M1>
    struct MyCopy // real
    {
        typedef typename M1::real_type RT;
        typedef MatrixView<RT,1> M2;

        static void call(const ConstMatrixView<RT>& m1x,
                         RT* m2p, int M, int N, int si, int sj)
        {
            TMVAssert(si == 1);
            M1 m1(m1x);
            M2 m2(m2p,M,N,si,sj);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"MyCopy real"<<std::endl;
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"M,N = "<<M<<','<<N<<std::endl;
            std::cout<<"si,sj = "<<si<<','<<sj<<std::endl;
            std::cout<<"m2p = "<<m2p<<"..."<<
                (m2p+(M-1)*si+(N-1)*sj+1)<<std::endl;
            std::cout<<"m1 = "<<m1<<std::endl;
#endif
            CopyM_Helper<-4,cs,rs,M1,M2>::call(m1,m2); 
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"m2 = "<<m2<<std::endl;
            std::cout<<"after copy\n";
#endif
        }
    };
    template <int cs, int rs, class M1>
    struct MyCopy<true,cs,rs,M1> // complex
    {
        typedef typename M1::real_type RT;
        typedef typename M1::const_realpart_type M1r;
        typedef MatrixView<RT,1> M2r;

        static void call(
            const ConstMatrixView<RT>& m1x,
            RT* m2p, int M, int N, int si, int sj)
        {
            TMVAssert(si == 1);
            const int size = N*sj;
            M1r m1r(m1x);
            M1r m1i(m1x.cptr()+1,m1r.colsize(),m1r.rowsize(),
                    m1r.stepi(),m1r.stepj());
            M2r m2r(m2p,M,N,si,sj);
            M2r m2i(m2p+size,M,N,si,sj);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"MyCopy complex"<<std::endl;
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"m2p = "<<m2p<<"..."<<
                (m2p+(M-1)*si+(N-1)*sj+1)<<std::endl;
            //std::cout<<"m1r = "<<m1r<<std::endl;
            //std::cout<<"m1i = "<<m1i<<std::endl;
            //std::cout<<"m2r = "<<m2r<<std::endl;
            //std::cout<<"m2i = "<<m2i<<std::endl;
#endif
            const int ix1 = M1::_conj ? -1 : 1; // 1 or -1 as required
            const Scaling<ix1,RT> one;
            CopyM_Helper<-4,cs,rs,M1r,M2r>::call(m1r,m2r); 
            MultXM_Helper<-4,cs,rs,false,ix1,RT,M1r,M2r>::call(one,m1i,m2i); 
#ifdef PRINTALGO_MM_BLOCK
            //std::cout<<"m2r => "<<m2r<<std::endl;
            //std::cout<<"m2i => "<<m2i<<std::endl;
#endif
        }
    };

    template <bool iscomplex, int cs, int rs, class M1>
    struct get_copy
    {
        typedef typename M1::real_type RT;
        typedef typename M1::const_rmview_type M1rm;
        typedef typename M1::const_cmview_type M1cm;

        typedef void F(const ConstMatrixView<RT>& m1x,
                       RT* m2p, int M, int N, int si, int sj);

        static F* call(const M1& m)
        {
            if (m.isrm()) 
                return &MyCopy<iscomplex,cs,rs,M1rm>::call;
            else if (m.iscm())
                return &MyCopy<iscomplex,cs,rs,M1cm>::call;
            else
                return &MyCopy<iscomplex,cs,rs,M1>::call;
        }
    };

    template <int MB, int NB, int KB, int ix, class T> 
    struct select_multmm
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_M_N_K(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<16,16,16,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_16_16(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<16,16,32,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_16_32(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<16,16,64,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_16_64(M,N,K,x,A,B,C); }
    };
    template <int MB, int ix, class T>
    struct select_multmm<MB,16,16,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_M_16_16_known<MB>(M,N,K,x,A,B,C); }
    };
    template <int MB, int ix, class T>
    struct select_multmm<MB,16,32,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_M_16_32_known<MB>(M,N,K,x,A,B,C); }
    };
    template <int MB, int ix, class T>
    struct select_multmm<MB,16,64,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_M_16_64_known<MB>(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<UNKNOWN,16,16,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_M_16_16(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<UNKNOWN,16,32,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_M_16_32(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<UNKNOWN,16,64,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_M_16_64(M,N,K,x,A,B,C); }
    };
    template <int NB, int ix, class T>
    struct select_multmm<16,NB,16,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_N_16_known<NB>(M,N,K,x,A,B,C); }
    };
    template <int NB, int ix, class T>
    struct select_multmm<16,NB,32,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_N_32_known<NB>(M,N,K,x,A,B,C); }
    };
    template <int NB, int ix, class T>
    struct select_multmm<16,NB,64,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_N_64_known<NB>(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<16,UNKNOWN,16,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_N_16(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<16,UNKNOWN,32,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_N_32(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<16,UNKNOWN,64,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_N_64(M,N,K,x,A,B,C); }
    };
    template <int KB, int ix, class T>
    struct select_multmm<16,16,KB,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_16_K_known<KB>(M,N,K,x,A,B,C); }
    };
    template <int ix, class T>
    struct select_multmm<16,16,UNKNOWN,ix,T>
    {
        static inline void call(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const T* A, const T* B, T* C)
        { multmm_16_16_K(M,N,K,x,A,B,C); }
    };

    template <bool iscomplex1, bool iscomplex2, 
              int cs, int rs, int xs, class RT>
    struct MyProd // both real
    {
        static void call(
            const int M, const int N, const int K,
            RT* m1p, const int si1, const int sj1,
            RT* m2p, const int si2, const int sj2,
            RT* m3p, const int si3, const int sj3)
        {
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<RT> m1(m1p,M,K,si1,sj1);
            MatrixView<RT> m2(m2p,K,N,si2,sj2);
            MatrixView<RT> m3(m3p,M,N,si3,sj3);
            std::cout<<"MyProd both real"<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
#endif
            const Scaling<1,RT> one;
            select_multmm<cs,rs,xs,1,RT>::call(M,N,K,one,m1p,m2p,m3p);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"After call to multmm\n";
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };
    template <int cs, int rs, int xs, class RT>
    struct MyProd<true,false,cs,rs,xs,RT> // m1 is complex
    {
        typedef std::complex<RT> CT;

        static void call(
            const int M, const int N, const int K,
            RT* m1p, const int si1, const int sj1,
            RT* m2p, const int si2, const int sj2,
            RT* m3p, const int si3, const int sj3)
        {
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<CT> m1((CT*)m1p,M,K,si1,sj1);
            MatrixView<RT> m2(m2p,K,N,si2,sj2);
            MatrixView<CT> m3((CT*)m3p,M,N,si3,sj3);
            std::cout<<"MyProd (m1 complex)"<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"m1 = "<<m1<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            const int xn1 = M*si1;
            const int xn3 = N*sj3;
            const Scaling<1,RT> one;
            select_multmm<cs,rs,xs,1,RT>::call(M,N,K,one,m1p,m2p,m3p);
            select_multmm<cs,rs,xs,1,RT>::call(M,N,K,one,m1p+xn1,m2p,m3p+xn3);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"After call to multmm\n";
            //std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };
    template <int cs, int rs, int xs, class RT>
    struct MyProd<false,true,cs,rs,xs,RT> // m2 is complex
    {
        typedef std::complex<RT> CT;

        static void call(
            const int M, const int N, const int K,
            RT* m1p, const int si1, const int sj1,
            RT* m2p, const int si2, const int sj2,
            RT* m3p, const int si3, const int sj3)
        {
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<RT> m1(m1p,M,K,si1,sj1);
            MatrixView<CT> m2((CT*)m2p,K,N,si2,sj2);
            MatrixView<CT> m3((CT*)m3p,M,N,si3,sj3);
            std::cout<<"MyProd (m2 complex)"<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"m1 = "<<m1<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            const int xn2 = N*sj2;
            const int xn3 = N*sj3;
            const Scaling<1,RT> one;
            select_multmm<cs,rs,xs,1,RT>::call(M,N,K,one,m1p,m2p,m3p);
            select_multmm<cs,rs,xs,1,RT>::call(M,N,K,one,m1p,m2p+xn2,m3p+xn3);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"After call to multmm\n";
            //std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };
    template <int cs, int rs, int xs, class RT>
    struct MyProd<true,true,cs,rs,xs,RT> // both complex
    {
        typedef std::complex<RT> CT;

        static void call(
            const int M, const int N, const int K,
            RT* m1p, const int si1, const int sj1,
            RT* m2p, const int si2, const int sj2,
            RT* m3p, const int si3, const int sj3)
        {
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<CT> m1((CT*)m1p,M,K,si1,sj1);
            MatrixView<CT> m2((CT*)m2p,K,N,si2,sj2);
            MatrixView<CT> m3((CT*)m3p,M,N,si3,sj3);
            std::cout<<"MyProd (both complex)"<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3<<std::endl;
#endif
            const int xn1 = M*si1;
            const int xn2 = N*sj2;
            const int xn3 = N*sj3;
            const Scaling<1,RT> one;
            const Scaling<-1,RT> mone;
            select_multmm<cs,rs,xs,1,RT>::call(M,N,K,one,m1p,m2p,m3p);
            select_multmm<cs,rs,xs,-1,RT>::call(
                M,N,K,mone,m1p+xn1,m2p+xn2,m3p);
            select_multmm<cs,rs,xs,1,RT>::call(M,N,K,one,m1p+xn1,m2p,m3p+xn3);
            select_multmm<cs,rs,xs,1,RT>::call(M,N,K,one,m1p,m2p+xn2,m3p+xn3);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"After call to multmm\n";
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };

    template <bool iscomplex1, bool iscomplex2, class RT>
    struct MyCleanup // both real
    {
        typedef void Kcleanup(
            const int M, const int N, const int K,
            const Scaling<1,RT>& x, const RT* A, const RT* B, RT* C);

        static void call(
            const int M, const int N, const int K,
            RT* m1p, const int si1, const int sj1,
            RT* m2p, const int si2, const int sj2,
            RT* m3p, const int si3, const int sj3,
            Kcleanup* cleanup)
        {
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<RT> m1(m1p,M,K,si1,sj1);
            MatrixView<RT> m2(m2p,K,N,si2,sj2);
            MatrixView<RT> m3(m3p,M,N,si3,sj3);
            std::cout<<"MyCleanup (both real) \n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3<<std::endl;
#endif
            TMVAssert(cleanup);
            const Scaling<1,RT> one;
            (*cleanup)(M,N,K,one,m1p,m2p,m3p); 
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"After cleanup call\n";
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };
    template <class RT>
    struct MyCleanup<true,false,RT> // m1 is complex
    {
        typedef std::complex<RT> CT;

        typedef void Kcleanup(
            const int M, const int N, const int K,
            const Scaling<1,RT>& x, const RT* A, const RT* B, RT* C);

        static void call(
            const int M, const int N, const int K,
            RT* m1p, const int si1, const int sj1,
            RT* m2p, const int si2, const int sj2,
            RT* m3p, const int si3, const int sj3,
            Kcleanup* cleanup)
        {
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<CT> m1((CT*)m1p,M,K,si1,sj1);
            MatrixView<RT> m2(m2p,K,N,si2,sj2);
            MatrixView<CT> m3((CT*)m3p,M,N,si3,sj3);
            std::cout<<"MyCleanup (m1 complex) \n";
            //std::cout<<"m1 = "<<m1<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            TMVAssert(cleanup);
            const Scaling<1,RT> one;
            const int xn1 = M*si1;
            const int xn3 = N*sj3;
            (*cleanup)(M,N,K,one,m1p,m2p,m3p);
            (*cleanup)(M,N,K,one,m1p+xn1,m2p,m3p+xn3);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"After cleanup call\n";
            //std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };
    template <class RT>
    struct MyCleanup<false,true,RT> // m2 is complex
    {
        typedef std::complex<RT> CT;
        typedef void Kcleanup(
            const int M, const int N, const int K,
            const Scaling<1,RT>& x, const RT* A, const RT* B, RT* C);

        static void call(
            const int M, const int N, const int K,
            RT* m1p, const int si1, const int sj1,
            RT* m2p, const int si2, const int sj2,
            RT* m3p, const int si3, const int sj3,
            Kcleanup* cleanup)
        {
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<RT> m1(m1p,M,K,si1,sj1);
            MatrixView<CT> m2((CT*)m2p,K,N,si2,sj2);
            MatrixView<CT> m3((CT*)m3p,M,N,si3,sj3);
            std::cout<<"MyCleanup (m2 complex) \n";
            //std::cout<<"m1 = "<<m1<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            TMVAssert(cleanup);
            const Scaling<1,RT> one;
            const int size2 = si2 == 1 ? N*sj2 : K*si2;
            const int size3 = si3 == 1 ? N*sj3 : M*si3;
            (*cleanup)(M,N,K,one,m1p,m2p,m3p);
            (*cleanup)(M,N,K,one,m1p,m2p+size2,m3p+size3);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"After cleanup call\n";
            //std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };
    template <class RT>
    struct MyCleanup<true,true,RT> // both complex
    {
        typedef std::complex<RT> CT;
        typedef void Kcleanup(
            const int M, const int N, const int K,
            const Scaling<1,RT>& x, const RT* A, const RT* B, RT* C);

        static void call(
            const int M, const int N, const int K,
            RT* m1p, const int si1, const int sj1,
            RT* m2p, const int si2, const int sj2,
            RT* m3p, const int si3, const int sj3,
            Kcleanup* cleanup)
        {
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<CT> m1((CT*)m1p,M,K,si1,sj1);
            MatrixView<CT> m2((CT*)m2p,K,N,si2,sj2);
            MatrixView<CT> m3((CT*)m3p,M,N,si3,sj3);
            std::cout<<"MyCleanup (both complex) \n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3<<std::endl;
#endif
            TMVAssert(cleanup);
            const Scaling<1,RT> one;
            const int xn1 = M*si1;
            const int xn2 = N*sj2;
            const int xn3 = N*sj3;

            (*cleanup)(M,N,K,one,m1p,m2p,m3p);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"MyCleanup (both complex) \n";
            std::cout<<"si1, si2, si3 = "<<si1<<" "<<si2<<" "<<si3<<std::endl;
            std::cout<<"sj1, sj2, sj3 = "<<sj1<<" "<<sj2<<" "<<sj3<<std::endl;
            std::cout<<"xn1, xn2, xn3 = "<<xn1<<" "<<xn2<<" "<<xn3<<std::endl;
            MatrixView<RT> m1r(m1p,M,K,si1,sj1);
            MatrixView<RT> m2r(m2p,K,N,si2,sj2);
            MatrixView<RT> m3r(m3p,M,N,si3,sj3);
            std::cout<<"m1r = "<<TMV_Text(m1r)<<"  "<<m1r<<std::endl;
            std::cout<<"m2r = "<<TMV_Text(m2r)<<"  "<<m2r<<std::endl;
            std::cout<<"m3r => "<<TMV_Text(m3r)<<"  "<<m3r<<std::endl;
#endif

            // This next one is tricky, since we have the wrong function to 
            // call with -1 as the first argument.
            // The way I do it here works, but it loses any known xs value
            // if that's not the value in M1 or M2.
            // I couldn't figure out any way to get it with template magic,
            // since  the cleanup object doesn't have the value we need,
            // and its type is Kcleanup which is a generic call function,
            // so we can't even  do something like a Traits specialization.
            // Anyway, they right way to do this is surely to edit the calling 
            // programs to provide two functions or something like that, and
            // that's too much work for something that really doesn't have a 
            // huge impact on the speed.
            const Scaling<-1,RT> mone;
            multmm_M_N_K(M,N,K,mone,m1p+xn1,m2p+xn2,m3p);
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<RT> m1i(m1p+xn1,M,K,si1,sj1);
            MatrixView<RT> m2i(m2p+xn2,K,N,si2,sj2);
            std::cout<<"m1i = "<<TMV_Text(m1i)<<"  "<<m1i<<std::endl;
            std::cout<<"m2i = "<<TMV_Text(m2i)<<"  "<<m2i<<std::endl;
            std::cout<<"m3r => "<<TMV_Text(m3r)<<"  "<<m3r<<std::endl;
#endif

            (*cleanup)(M,N,K,one,m1p+xn1,m2p,m3p+xn3);
            (*cleanup)(M,N,K,one,m1p,m2p+xn2,m3p+xn3);
#ifdef PRINTALGO_MM_BLOCK
            MatrixView<RT> m3i(m3p+xn3,M,N,si3,sj3);
            std::cout<<"m3i => "<<TMV_Text(m3i)<<"  "<<m3i<<std::endl;
            std::cout<<"After cleanup call\n";
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };

    template <bool iscomplex_x, bool iscomplex2, bool add, 
              int ix, class T, int cs, int rs, class M2>
    struct MyAssign // both real
    {
        typedef typename M2::real_type RT;
        typedef typename M2::value_type T2; 
        // Note: T2 might be complex, but m1 isn't.
        typedef ConstMatrixView<RT,1> M1r;

        static void call(
            const RT* xp, 
            RT* m1p, const int M, const int N, const int si1, const int sj1,
            MatrixView<RT> m2x)
        {
            TMVAssert(si1 == 1);
            const Scaling<ix,T> x(*((T*)(xp)));
            M1r m1(m1p,M,N,si1,sj1);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"MyAssign real"<<std::endl;
            std::cout<<"x = "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"cs,rs,M,N,si1,sj1 = "<<cs<<','<<rs<<','<<M<<','<<N<<
                ','<<si1<<','<<sj1<<std::endl;
            std::cout<<"m2x = "<<TMV_Text(m2x)<<std::endl;
#endif
            M2 m2((T2*)m2x.ptr(),m2x.colsize(),m2x.rowsize(),
                  m2x.stepi(),m2x.stepj());
            MultXM_Helper<-4,cs,rs,add,ix,T,M1r,M2>::call(x,m1,m2); 
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"m2p = "<<m2x.ptr()<<"..."<<
                (m2x.ptr()+(M-1)*m2x.stepi()+(N-1)*m2x.stepj()+1)<<std::endl;
            std::cout<<"m2 => "<<m2<<std::endl;
#endif
        }
    };
    template <bool add, int ix, class T, int cs, int rs, class M2>
    struct MyAssign<true,false,add,ix,T,cs,rs,M2> // x is complex
    {
        typedef typename M2::real_type RT;
        typedef typename M2::value_type T2; 
        // Note: T2 might be complex, but m1 isn't.
        typedef ConstMatrixView<RT,1> M1r;

        static void call(
            const RT* xp, 
            RT* m1p, const int M, const int N, const int si1, const int sj1,
            MatrixView<RT> m2x)
        {
            TMVAssert(si1 == 1);
            const Scaling<ix,T>& x(*((T*)(xp)));
            M1r m1(m1p,M,N,si1,sj1);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"MyAssign x complex"<<std::endl;
            std::cout<<"x = "<<T(x)<<std::endl;
            //std::cout<<"m1 = "<<m1<<std::endl;
#endif
            M2 m2((T2*)m2x.ptr(),m2x.colsize(),m2x.rowsize(),
                  m2x.stepi(),m2x.stepj());
            MultXM_Helper<-4,cs,rs,add,ix,T,M1r,M2>::call(x,m1,m2); 
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"m2p = "<<m2x.ptr()<<"..."<<
                (m2x.ptr()+(M-1)*m2x.stepi()+(N-1)*m2x.stepj()+1)<<std::endl;
            //std::cout<<"m2 => "<<m2<<std::endl;
#endif
        }
    };
    template <bool add, int ix, class T, int cs, int rs, class M2>
    struct MyAssign<false,true,add,ix,T,cs,rs,M2> // m2 is complex
    {
        typedef typename M2::real_type RT;
        typedef ConstMatrixView<RT,1> M1r;
        typedef typename M2::realpart_type M2r;
        typedef typename M2::imagpart_type M2i;

        static void call(
            const RT* xp,
            RT* m1p, const int M, const int N,
            const int si1, const int sj1,
            MatrixView<RT> m2x)
        {
            TMVAssert(si1 == 1);
            const Scaling<ix,T>& x(*((T*)(xp)));
            M2r m2r(m2x);
            M2i m2i(m2x.ptr()+1,m2r.colsize(),m2r.rowsize(),
                    m2r.stepi(),m2r.stepj());
            const int size = N*sj1;
            M1r m1r(m1p,M,N,si1,sj1);
            M1r m1i(m1p+size,M,N,si1,sj1);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"MyAssign m2 complex"<<std::endl;
            std::cout<<"si1,sj1,size = "<<si1<<"  "<<sj1
                <<"  "<<size<<std::endl;
            std::cout<<"x = "<<T(x)<<std::endl;
            std::cout<<"m1r = "<<m1r<<std::endl;
            std::cout<<"m1i = "<<m1i<<std::endl;
#endif
            MultXM_Helper<-4,cs,rs,add,ix,T,M1r,M2r>::call(x,m1r,m2r); 
            MultXM_Helper<-4,cs,rs,add,ix,T,M1r,M2i>::call(x,m1i,m2i); 
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"m2p = "<<m2x.ptr()<<"..."<<
                (m2x.ptr()+(M-1)*m2x.stepi()+(N-1)*m2x.stepj()+1)<<std::endl;
            std::cout<<"m1r = "<<m1r<<std::endl;
            std::cout<<"m2r => "<<m2r<<std::endl;
#endif
        }
    };
    template <bool add, int ix, class T, int cs, int rs, class M2>
    struct MyAssign<true,true,add,ix,T,cs,rs,M2> // both complex
    {
        typedef typename M2::real_type RT;
        typedef typename M2::realpart_type M2r;
        typedef typename M2::imagpart_type M2i;
        typedef ConstMatrixView<RT,1> M1r;

        static void call(
            const RT* xp,
            RT* m1p, const int M, const int N,
            const int si1, const int sj1,
            MatrixView<RT> m2x)
        {
            TMVAssert(si1 == 1);
            const T x(*((T*)(xp)));
            M2r m2r(m2x);
            M2i m2i(m2x.ptr()+1,m2r.colsize(),m2r.rowsize(),
                    m2r.stepi(),m2r.stepj());
            const int size = N*sj1;
            M1r m1r(m1p,M,N,si1,sj1);
            M1r m1i(m1p+size,M,N,si1,sj1);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"MyAssign both complex"<<std::endl;
            std::cout<<"x = "<<T(x)<<std::endl;
            std::cout<<"m1r = "<<m1r<<std::endl;
            std::cout<<"m1i = "<<m1i<<std::endl;
#endif
            //m2r (+)= TMV_REAL(x) * m1r;
            //m2r -=   TMV_IMAG(x) * m1i;
            //m2i (+)= TMV_REAL(x) * m1i;
            //m2i +=   TMV_IMAG(x) * m1r;
            MultXM_Helper<-4,cs,rs,add,0,RT,M1r,M2r>::call(
                Scaling<0,RT>(TMV_REAL(x)),m1r,m2r); 
            MultXM_Helper<-4,cs,rs,true,0,RT,M1r,M2i>::call(
                Scaling<0,RT>(-TMV_IMAG(x)),m1i,m2r); 
            MultXM_Helper<-4,cs,rs,add,0,RT,M1r,M2r>::call(
                Scaling<0,RT>(TMV_REAL(x)),m1i,m2i); 
            MultXM_Helper<-4,cs,rs,true,0,RT,M1r,M2i>::call(
                Scaling<0,RT>(TMV_IMAG(x)),m1r,m2i); 
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"m2p = "<<m2x.ptr()<<"..."<<
                (m2x.ptr()+(M-1)*m2x.stepi()+(N-1)*m2x.stepj()+1)<<std::endl;
            std::cout<<"m1r = "<<m1r<<std::endl;
            std::cout<<"m2r => "<<m2r<<std::endl;
#endif
        }
    };

    template <bool zx, bool z2, class T, int cs, int rs, class M2>
    struct get_assign
    {
        typedef typename M2::real_type RT;
        typedef typename M2::rmview_type M2rm;
        typedef typename M2::cmview_type M2cm;

        typedef void F(
            const RT* xp,
            RT* m1p, const int M, const int N, const int si1, const int sj1,
            MatrixView<RT> m2x);

        static F* call(const bool add, const T& x, const M2& m)
        {
            if (add)
                if (TMV_IMAG(x) == RT(0))
                    if (TMV_REAL(x) == RT(1))
                        if (m.isrm()) 
                            return &MyAssign<zx,z2,true,1,RT,cs,rs,M2rm>::call;
                        else if (m.iscm())
                            return &MyAssign<zx,z2,true,1,RT,cs,rs,M2cm>::call;
                        else
                            return &MyAssign<zx,z2,true,1,RT,cs,rs,M2>::call;
                    else if (TMV_REAL(x) == RT(-1))
                        if (m.isrm()) 
                            return 
                                &MyAssign<zx,z2,true,-1,RT,cs,rs,M2rm>::call;
                        else if (m.iscm())
                            return 
                                &MyAssign<zx,z2,true,-1,RT,cs,rs,M2cm>::call;
                        else
                            return &MyAssign<zx,z2,true,-1,RT,cs,rs,M2>::call;
                    else 
                        if (m.isrm()) 
                            return &MyAssign<zx,z2,true,0,RT,cs,rs,M2rm>::call;
                        else if (m.iscm())
                            return &MyAssign<zx,z2,true,0,RT,cs,rs,M2cm>::call;
                        else
                            return &MyAssign<zx,z2,true,0,RT,cs,rs,M2>::call;
                else
                    if (m.isrm()) 
                        return &MyAssign<zx,z2,true,0,T,cs,rs,M2rm>::call;
                    else if (m.iscm())
                        return &MyAssign<zx,z2,true,0,T,cs,rs,M2cm>::call;
                    else
                        return &MyAssign<zx,z2,true,0,T,cs,rs,M2>::call;
            else
                if (TMV_IMAG(x) == RT(0))
                    if (TMV_REAL(x) == RT(1))
                        if (m.isrm()) 
                            return 
                                &MyAssign<zx,z2,false,1,RT,cs,rs,M2rm>::call;
                        else if (m.iscm())
                            return 
                                &MyAssign<zx,z2,false,1,RT,cs,rs,M2cm>::call;
                        else
                            return &MyAssign<zx,z2,false,1,RT,cs,rs,M2>::call;
                    else if (TMV_REAL(x) == RT(-1))
                        if (m.isrm()) 
                            return 
                                &MyAssign<zx,z2,false,-1,RT,cs,rs,M2rm>::call;
                        else if (m.iscm())
                            return 
                                &MyAssign<zx,z2,false,-1,RT,cs,rs,M2cm>::call;
                        else
                            return &MyAssign<zx,z2,false,-1,RT,cs,rs,M2>::call;
                    else 
                        if (m.isrm()) 
                            return 
                                &MyAssign<zx,z2,false,0,RT,cs,rs,M2rm>::call;
                        else if (m.iscm())
                            return 
                                &MyAssign<zx,z2,false,0,RT,cs,rs,M2cm>::call;
                        else
                            return &MyAssign<zx,z2,false,0,RT,cs,rs,M2>::call;
                else
                    if (m.isrm()) 
                        return &MyAssign<zx,z2,false,0,T,cs,rs,M2rm>::call;
                    else if (m.iscm())
                        return &MyAssign<zx,z2,false,0,T,cs,rs,M2cm>::call;
                    else
                        return &MyAssign<zx,z2,false,0,T,cs,rs,M2>::call;
        }
    };


    // A helper function to round a column length up to the next
    // multiple of 2 or 4 if required for the SSE commands for that type.
    template <class T>
    static inline int RoundUp(const int x)
    { return x; }
    template <>
    inline int RoundUp<double>(const int x)
    { 
        return (
#ifdef __SSE2__
            x == 0 ? 0 : (((x-1)>>1)+1)<<1
#else
            x
#endif
        );
    }
    template <>
    inline int RoundUp<float>(const int x)
    { 
        return (
#ifdef __SSE__
            x == 0 ? 0 : (((x-1)>>2)+1)<<2
#else
            x
#endif
        );
    }

    //
    // Algo 63: RecursiveBlockMultMM
    //

    // DoRecursiveBlockMultMM3 is the main workhorse routine for the 
    // RecursiveBlockMultMM function described in more detail below.  
    // It splits a section up in to (up to) 8 sub-problems and 
    // recurses back to itself to solve them.
    //
    // x,m1,m2,m3 are the original matrices (and scale) being worked on, 
    //   with the slight complication that we do this in terms of the real
    //   part if anything is complex.  The complex-ness is dealth with in
    //   the various copy and assign functions.  This complication helps
    //   to greatly reduce the compiled size.
    // MB,NB,KB are the block sizes being used on this block if it is the 
    //   final call in the recursion..
    // MB0,NB0,KB0 are the real full block sizes.
    // lgMB, etc are the base-2 log of the block sizes (used for << and >>)
    // i1,j1,k1,i2,j2,k2 are the current index ranges being worked on.
    // Mb,Nb,Kb are the number of blocks to do in each direction.  There is 
    //   a slight difference in the meaning of Mb,Nb as opposed to Kb.
    //   Mb,Nb include the last partial block in that direction if any.
    //   Kb does not -- it counts only the number of full KB blocks.
    // m1p,m2p,m3p are pointers to temporary memory used for copying
    //   blocks into packed storage which makes the block products much faster.
    // two1,two2,two3 help deal with the possibility of complex m1,m2,m3.
    //   They are 2 if the matrix is complex, and 1 if real.
    // firstm1, firstm2, firstm3, indicate whether this is the first time
    //   we are dealing with this section of m1,m2,m3.  If so, we need to
    //   copy the data into the temporary storage (or zero m3 in that case).
    // lastm3 indicates whether this is the last time we are dealing with 
    //   this section of m3.  If se, we need to assign the result back to
    //   the actual m3 matrix.
    // copy1, etc. are pointers to the functions that do the actual copying,
    //   assigning, and products.  These have template parameters that 
    //   govern whether the matrix is complex or real, what the step sizes
    //   are, and so forth.  So they are efficient at performing their tasks.
    //   But passing them as pointers here means that the compiler doesn't
    //   have to know these template parameters for this function, thus
    //   reducing the compiled size.
    //   The ones that start with "my" are the ones to use if this is the
    //   final 1x1x1 block.  The first one in each set is the one to use
    //   when we have a full block.  The others are used for partial blocks.
    //   The non "my" functions are bumped up to the "my" spot during the 
    //   recursion where appropriate, so when we get to the final recursion,
    //   the correct function is always in the "my" spot.
    //   
    template <class CopyF, class AssignF, class ProdF, class CleanupF,
              class KF, class RT>
    static void DoRecursiveBlockMultMM3(
        const RT* x, const ConstMatrixView<RT>& m1,
        const ConstMatrixView<RT>& m2, MatrixView<RT> m3,
        const int MB, const int NB, const int KB,
        const int MB0, const int NB0, 
        const int lgMB, const int lgNB, const int lgKB,
        const int i1, const int j1, int k1,
        const int i2, const int j2, const int k2,
        const int Mb, const int Nb, const int Kb,
        RT* m1p, RT* m2p, RT* m3p,
        const int two1, const int two2, const int two3,
        bool firstm1, bool firstm2, bool firstm3, bool lastm3,
        CopyF* copy1, CopyF* copy1_k, CopyF* mycopy1, CopyF* mycopy1_k,
        CopyF* copy2, CopyF* copy2_k, CopyF* mycopy2, CopyF* mycopy2_k,
        AssignF* assign, AssignF* assigna, AssignF* assignb, AssignF* myassign, 
        ProdF* prod, ProdF* proda, ProdF* prodb, ProdF* myprod,
        CleanupF* mycleanup, KF* kf, KF* kfa, KF* kfb, KF* mykf)
    {
#ifdef PRINTALGO_MM_BLOCK
        std::cout<<"Start recursion:\n";
        std::cout<<"i1,j1,k1 = "<<i1<<','<<j1<<','<<k1<<std::endl;
        std::cout<<"i2,j2,k2 = "<<i2<<','<<j2<<','<<k2<<std::endl;
        std::cout<<"Mb,Nb,Kb = "<<Mb<<','<<Nb<<','<<Kb<<std::endl;
        std::cout<<"MB,NB,KB = "<<MB<<','<<NB<<','<<KB<<std::endl;
        std::cout<<"firstm1 = "<<firstm1<<", firstm2 = "<<firstm2<<std::endl;
        std::cout<<"firstm3 = "<<firstm3<<", lastm3 = "<<lastm3<<std::endl;
        std::cout<<"m1p,m2p,m3p = "<<m1p<<','<<m2p<<','<<m3p<<std::endl;
        std::cout<<"two1,2,3 = "<<two1<<','<<two2<<','<<two3<<std::endl;
        std::cout<<"First values = "<<*m1p<<','<<*m2p<<','<<*m3p<<std::endl;
        std::cout<<"kf = "<<(size_t)kf<<"  "<<(size_t)kfa<<"  "<<(size_t)kfb
            <<"  "<<(size_t)mykf<<std::endl;
#endif

        if (Mb > 1 && Nb > 1 && Kb > 1)
        {
            const int Mx = Mb>>1; // = Mb/2
            const int Nx = Nb>>1; // = Nb/2
            const int Kx = Kb>>1; // = Kb/2
            const int Mz = Mx<<lgMB;
            const int Nz = Nx<<lgNB;
            const int Kz = Kx<<lgKB;
            const int Mw = (i2-i1-Mz);
            const int Nw = (j2-j1-Nz);
            const int Kw = RoundUp<RT>(k2-k1-Kz);

            const int im = i1 + Mz;
            const int jm = j1 + Nz;
            const int km = k1 + Kz;

            const int n1a = two1*Mz*Kz;
            const int n1b = n1a + two1*Mz*Kw;
            const int n1c = n1b + two1*Mw*Kz;

            const int n2a = two2*Nz*Kz;
            const int n2b = n2a + two2*Nz*Kw;
            const int n2c = n2b + two2*Nw*Kz;

            const int n3a = two3*Mz*Nz;
            const int n3b = n3a + two3*Mz*Nw;
            const int n3c = n3b + two3*Mw*Nz;
#ifdef PRINTALGO_MM_BLOCK
            const int My = Mb-Mx; 
            const int Ny = Nb-Nx;
            const int Ky = Kb-Kx;
            std::cout<<"Regular split (all > 1)\n";
            std::cout<<"Mxyzw = "<<Mx<<','<<My<<','<<Mz<<','<<Mw<<std::endl;
            std::cout<<"Nxyzw = "<<Nx<<','<<Ny<<','<<Nz<<','<<Nw<<std::endl;
            std::cout<<"Kxyzw = "<<Kx<<','<<Ky<<','<<Kz<<','<<Kw<<std::endl;
            std::cout<<"n1abc = "<<n1a<<','<<n1b<<','<<n1c<<std::endl;
            std::cout<<"n2abc = "<<n2a<<','<<n2b<<','<<n2c<<std::endl;
            std::cout<<"n3abc = "<<n3a<<','<<n3b<<','<<n3c<<std::endl;
#endif

            // Split matrix product up into:
            // [ A B ] * [ E F ] = [ I J ]
            // [ C D ]   [ G H ]   [ K L ]

            // AE = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB0,NB0,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,k1,im,jm,km,Mx,Nx,Kx,
                m1p,m2p,m3p, two1,two2,two3,
                firstm1,firstm2,firstm3,false,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assign,assign,
                prod,prod,prod,prod,
                mycleanup,kf,kf,kf,kf);
            // BG = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB0,NB0,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,km,im,jm,k2,Mx,Nx,Kb-Kx,
                m1p+n1a,m2p+n2a,m3p, two1,two2,two3,
                firstm1,firstm2,false,lastm3,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assign,assign,
                prod,prod,prod,prod,
                mycleanup,kf,kf,kf,kf);
            // AF = J
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB0,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,jm,k1,im,j2,km,Mx,Nb-Nx,Kx,
                m1p,m2p+n2b,m3p+n3a, two1,two2,two3,
                false,firstm2,firstm3,false,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assign,assigna,
                prod,proda,prod,proda,
                mycleanup,kf,kfa,kf,kfa);
            // BH = J
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB0,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,jm,km,im,j2,k2,Mx,Nb-Nx,Kb-Kx,
                m1p+n1a,m2p+n2c,m3p+n3a, two1,two2,two3,
                false,firstm2,false,lastm3,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assign,assigna,
                prod,proda,prod,proda,
                mycleanup,kf,kfa,kf,kfa);
            // CE = K
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB0,KB,MB0,NB0,lgMB,lgNB,lgKB,
                im,j1,k1,i2,jm,km,Mb-Mx,Nx,Kx,
                m1p+n1b,m2p,m3p+n3b, two1,two2,two3,
                firstm1,false,firstm3,false,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assignb,assignb,
                prod,prod,prodb,prodb,
                mycleanup,kf,kf,kfb,kfb);
            // DG = K
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB0,KB,MB0,NB0,lgMB,lgNB,lgKB,
                im,j1,km,i2,jm,k2,Mb-Mx,Nx,Kb-Kx,
                m1p+n1c,m2p+n2a,m3p+n3b, two1,two2,two3,
                firstm1,false,false,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assignb,assignb,
                prod,prod,prodb,prodb,
                mycleanup,kf,kf,kfb,kfb);
            // CF = L
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                im,jm,k1,i2,j2,km,Mb-Mx,Nb-Nx,Kx,
                m1p+n1b,m2p+n2b,m3p+n3c, two1,two2,two3,
                false,false,firstm3,false,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
            // DH = L
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                im,jm,km,i2,j2,k2,Mb-Mx,Nb-Nx,Kb-Kx,
                m1p+n1c,m2p+n2c,m3p+n3c, two1,two2,two3,
                false,false,false,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
        }
        else if (Mb == 1 && Nb == 1 && Kb == 1) 
        {
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"No split (all == 1)"<<std::endl;
            // For real:
            //std::cout<<"actual m1 = "<<
            //    m1.cSubMatrix(i1,i2,k1,k1+KB)<<std::endl;
            //std::cout<<"actual m2 = "<<
            //    m2.cSubMatrix(k1,k1+KB,j1,j2)<<std::endl;
            // For complex:
            std::cout<<"actual m1 = "<<
                ConstMatrixView<std::complex<RT> >(
                    (std::complex<RT>*)m1.cptr(),m1.colsize(),m1.rowsize(),
                    m1.stepi()/2,m1.stepj()/2).cSubMatrix(i1,i2,k1,k1+KB)
                <<std::endl;
            std::cout<<"actual m2 = "<<
                ConstMatrixView<std::complex<RT> >(
                    (std::complex<RT>*)m2.cptr(),m2.colsize(),m2.rowsize(),
                    m2.stepi()/2,m2.stepj()/2).cSubMatrix(k1,k1+KB,j1,j2)
                <<std::endl;
#endif
#ifdef TEST_POINTERS
            TMVAssert(m1p < glob_m1p_end);
            TMVAssert(m2p < glob_m2p_end);
            TMVAssert(m3p < glob_m3p_end);
            TMVAssert(m1p+(two1*MB*KB) <= glob_m1p_end);
            TMVAssert(m2p+(two2*KB*NB) <= glob_m2p_end);
            TMVAssert(m3p+(two3*MB*MB) <= glob_m3p_end);
            std::cout<<"end of m2 = "<<
                ConstVectorView<long,1>((const long*)glob_m2p_end,10)
                <<std::endl;
#endif
            if (firstm1) 
                (*mycopy1) (m1.cSubMatrix(i1,i2,k1,k1+KB).transpose(),
                            m1p,KB,MB,1,KB);
            if (firstm2) 
                (*mycopy2) (m2.cSubMatrix(k1,k1+KB,j1,j2),m2p,KB,NB,1,KB);
            if (firstm3)
            {
                const int size3 = two3*MB*NB;
                VectorView<RT,1> m3x(m3p,size3,1);
                m3x.setZero();
            }
            (*myprod) (MB,NB,KB, m1p,KB,1, m2p,1,KB, m3p,1,MB);

            if (lastm3) 
            {
                if (mykf) 
                {
                    const int size1 = (two1*MB)<<lgKB;
                    const int size2 = (two2*NB)<<lgKB;
                    m1p += size1;
                    m2p += size2;
                    k1 += KB;
                    const int Kc = k2-k1;
                    const int Kd = RoundUp<RT>(Kc);
#ifdef TEST_POINTERS
                    TMVAssert(m1p < glob_m1p_end);
                    TMVAssert(m2p < glob_m2p_end);
                    TMVAssert(m3p < glob_m3p_end);
                    TMVAssert(m1p+(two1*MB*Kd) <= glob_m1p_end);
                    TMVAssert(m2p+(two2*Kd*NB) <= glob_m2p_end);
                    TMVAssert(m3p+(two3*MB*MB) <= glob_m3p_end);
                    std::cout<<"end of m2 = "<<
                        ConstVectorView<long,1>((const long*)glob_m2p_end,10)
                        <<std::endl;
#endif
                    if (firstm1) 
                        (*mycopy1_k) (m1.cSubMatrix(i1,i2,k1,k2).transpose(),
                                      m1p,Kc,MB,1,Kd);
                    if (firstm2) 
                        (*mycopy2_k) (m2.cSubMatrix(k1,k2,j1,j2),
                                      m2p,Kc,NB,1,Kd);
                    (*mycleanup) (MB,NB,Kc,m1p,Kd,1,m2p,1,Kd,m3p,1,MB,mykf);
                }
                (*myassign) (x,m3p,MB,NB,1,MB,m3.cSubMatrix(i1,i2,j1,j2));
            }
        }
        else if (Mb > 1 && Nb > 1)
        {
            TMVAssert(Mb > 1 && Nb > 1 && Kb == 1);
            const int Mx = Mb>>1; // = Mb/2
            const int Nx = Nb>>1; // = Nb/2
            const int Mz = Mx<<lgMB;
            const int Nz = Nx<<lgNB;
            const int Mw = (i2-i1-Mz);
            const int Nw = (j2-j1-Nz);
            const int Kw = RoundUp<RT>(k2-k1);

            const int im = i1 + Mz;
            const int jm = j1 + Nz;

            const int n1b = two1*Mz*Kw;
            const int n2b = two2*Nz*Kw;
            const int n3a = two3*Mz*Nz;
            const int n3b = n3a + two3*Mz*Nw;
            const int n3c = n3b + two3*Mw*Nz;
#ifdef PRINTALGO_MM_BLOCK
            const int My = Mb-Mx; 
            const int Ny = Nb-Nx;
            std::cout<<"M,N > 1 split\n";
            std::cout<<"Mxyzw = "<<Mx<<','<<My<<','<<Mz<<','<<Mw<<std::endl;
            std::cout<<"Nxyzw = "<<Nx<<','<<Ny<<','<<Nz<<','<<Nw<<std::endl;
            std::cout<<"Kw = "<<Kw<<std::endl;
            std::cout<<"n1b = "<<n1b<<std::endl;
            std::cout<<"n2b = "<<n2b<<std::endl;
            std::cout<<"n3abc = "<<n3a<<','<<n3b<<','<<n3c<<std::endl;
#endif
            // Split matrix product up into:
            // [ A ] * [ E F ] = [ I J ]
            // [ C ]             [ K L ]

            // AE = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB0,NB0,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,k1,im,jm,k2,Mx,Nx,Kb,
                m1p,m2p,m3p, two1,two2,two3,
                firstm1,firstm2,firstm3,lastm3,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assign,assign,
                prod,prod,prod,prod,
                mycleanup,kf,kf,kf,kf);
            // AF = J
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB0,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,jm,k1,im,j2,k2,Mx,Nb-Nx,Kb,
                m1p,m2p+n2b,m3p+n3a, two1,two2,two3,
                false,firstm2,firstm3,lastm3,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assign,assigna,
                prod,proda,prod,proda,
                mycleanup,kf,kfa,kf,kfa);
            // CE = K
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB0,KB,MB0,NB0,lgMB,lgNB,lgKB,
                im,j1,k1,i2,jm,k2,Mb-Mx,Nx,Kb,
                m1p+n1b,m2p,m3p+n3b, two1,two2,two3,
                firstm1,false,firstm3,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assignb,assignb,
                prod,prod,prodb,prodb,
                mycleanup,kf,kf,kfb,kfb);
            // CF = L
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                im,jm,k1,i2,j2,k2,Mb-Mx,Nb-Nx,Kb,
                m1p+n1b,m2p+n2b,m3p+n3c, two1,two2,two3,
                false,false,firstm3,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
        }
        else if (Mb > 1 && Kb > 1)
        {
            TMVAssert(Mb > 1 && Kb > 1 && Nb == 1);
            const int Mx = Mb>>1; // = Mb/2
            const int Kx = Kb>>1; // = Kb/2
            const int Mz = Mx<<lgMB;
            const int Kz = Kx<<lgKB;
            const int Mw = (i2-i1-Mz);
            const int Nw = (j2-j1);
            const int Kw = RoundUp<RT>(k2-k1-Kz);

            const int im = i1 + Mz;
            const int km = k1 + Kz;

            const int n1a = two1*Mz*Kz;
            const int n1b = n1a + two1*Mz*Kw;
            const int n1c = n1b + two1*Mw*Kz;
            const int n2a = two2*Nw*Kz;
            const int n3b = two3*Mz*Nw;

#ifdef PRINTALGO_MM_BLOCK
            const int My = Mb-Mx; 
            const int Ky = Kb-Kx;
            std::cout<<"M,K > 1 split\n";
            std::cout<<"Mxyzw = "<<Mx<<','<<My<<','<<Mz<<','<<Mw<<std::endl;
            std::cout<<"Nw = "<<Nw<<std::endl;
            std::cout<<"Kxyzw = "<<Kx<<','<<Ky<<','<<Kz<<','<<Kw<<std::endl;
            std::cout<<"n1abc = "<<n1a<<','<<n1b<<','<<n1c<<std::endl;
            std::cout<<"n2a = "<<n2a<<std::endl;
            std::cout<<"n3b = "<<n3b<<std::endl;
#endif

            // Split matrix product up into:
            // [ A B ] * [ E ] = [ I ]
            // [ C D ]   [ G ]   [ K ]

            // AE = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB0,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,k1,im,j2,km,Mx,Nb,Kx,
                m1p,m2p,m3p, two1,two2,two3,
                firstm1,firstm2,firstm3,false,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assign,assigna,
                prod,proda,prod,proda,
                mycleanup,kf,kfa,kf,kfa);
            // BG = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB0,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,km,im,j2,k2,Mx,Nb,Kb-Kx,
                m1p+n1a,m2p+n2a,m3p, two1,two2,two3,
                firstm1,firstm2,false,lastm3,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assign,assigna,
                prod,proda,prod,proda,
                mycleanup,kf,kfa,kf,kfa);
            // CE = K
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                im,j1,k1,i2,j2,km,Mb-Mx,Nb,Kx,
                m1p+n1b,m2p,m3p+n3b, two1,two2,two3,
                firstm1,false,firstm3,false,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
            // DG = K
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                im,j1,km,i2,j2,k2,Mb-Mx,Nb,Kb-Kx,
                m1p+n1c,m2p+n2a,m3p+n3b, two1,two2,two3,
                firstm1,false,false,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
        }
        else if (Nb > 1 && Kb > 1)
        {
            TMVAssert(Nb > 1 && Kb > 1 && Mb == 1);
            const int Nx = Nb>>1; // = Nb/2
            const int Kx = Kb>>1; // = Kb/2
            const int Nz = Nx<<lgNB;
            const int Kz = Kx<<lgKB;
            const int Mw = (i2-i1);
            const int Nw = (j2-j1-Nz);
            const int Kw = RoundUp<RT>(k2-k1-Kz);

            const int jm = j1 + Nz;
            const int km = k1 + Kz;

            const int n1a = two1*Mw*Kz;
            const int n2a = two2*Nz*Kz;
            const int n2b = n2a + two2*Nz*Kw;
            const int n2c = n2b + two2*Nw*Kz;
            const int n3a = two3*Mw*Nz;

#ifdef PRINTALGO_MM_BLOCK
            const int Ny = Nb-Nx;
            const int Ky = Kb-Kx;
            std::cout<<"N,K > 1 split\n";
            std::cout<<"Mw = "<<Mw<<std::endl;
            std::cout<<"Nxyzw = "<<Nx<<','<<Ny<<','<<Nz<<','<<Nw<<std::endl;
            std::cout<<"Kxyzw = "<<Kx<<','<<Ky<<','<<Kz<<','<<Kw<<std::endl;
            std::cout<<"n1a = "<<n1a<<std::endl;
            std::cout<<"n2abc = "<<n2a<<','<<n2b<<','<<n2c<<std::endl;
            std::cout<<"n3a = "<<n3a<<std::endl;
#endif

            // Split matrix product up into:
            // [ A B ] * [ E F ] = [ I J ]
            //           [ G H ]          

            // AE = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB0,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,k1,i2,jm,km,Mb,Nx,Kx,
                m1p,m2p,m3p, two1,two2,two3,
                firstm1,firstm2,firstm3,false,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assignb,assignb,
                prod,prod,prodb,prodb,
                mycleanup,kf,kf,kfb,kfb);
            // BG = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB0,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,km,i2,jm,k2,Mb,Nx,Kb-Kx,
                m1p+n1a,m2p+n2a,m3p, two1,two2,two3,
                firstm1,firstm2,false,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assignb,assignb,
                prod,prod,prodb,prodb,
                mycleanup,kf,kf,kfb,kfb);
            // AF = J
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,jm,k1,i2,j2,km,Mb,Nb-Nx,Kx,
                m1p,m2p+n2b,m3p+n3a, two1,two2,two3,
                false,firstm2,firstm3,false,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
            // BH = J
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,jm,km,i2,j2,k2,Mb,Nb-Nx,Kb-Kx,
                m1p+n1a,m2p+n2c,m3p+n3a, two1,two2,two3,
                false,firstm2,false,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
        }
        else if (Mb > 1)
        {
            TMVAssert(Mb > 1 && Nb == 1 && Kb == 1);
            const int Mx = Mb>>1; // = Mb/2
            const int Mz = Mx<<lgMB;
            const int Nw = (j2-j1);
            const int Kw = RoundUp<RT>(k2-k1);

            const int im = i1 + Mz;

            const int n1b = two1*Mz*Kw;
            const int n3b = two3*Mz*Nw;

#ifdef PRINTALGO_MM_BLOCK
            const int My = Mb-Mx; 
            const int Mw = (i2-i1-Mz);
            std::cout<<"M > 1 split\n";
            std::cout<<"Mxyzw = "<<Mx<<','<<My<<','<<Mz<<','<<Mw<<std::endl;
            std::cout<<"Nw = "<<Nw<<std::endl;
            std::cout<<"Kw = "<<Kw<<std::endl;
            std::cout<<"n1b = "<<n1b<<std::endl;
            std::cout<<"n3b = "<<n3b<<std::endl;
#endif

            // Split matrix product up into:
            // [ A ] * [ E ] = [ I ]
            // [ C ]           [ K ]

            // AE = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB0,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,k1,im,j2,k2,Mx,Nb,Kb,
                m1p,m2p,m3p, two1,two2,two3,
                firstm1,firstm2,firstm3,lastm3,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assign,assigna,
                prod,proda,prod,proda,
                mycleanup,kf,kfa,kf,kfa);
            // CE = K
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                im,j1,k1,i2,j2,k2,Mb-Mx,Nb,Kb,
                m1p+n1b,m2p,m3p+n3b, two1,two2,two3,
                firstm1,false,firstm3,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
        }
        else if (Nb > 1)
        {
            TMVAssert(Nb > 1 && Kb == 1 && Mb == 1);
            const int Nx = Nb>>1; // = Nb/2
            const int Nz = Nx<<lgNB;
            const int Mw = (i2-i1);
            const int Kw = RoundUp<RT>(k2-k1);

            const int jm = j1 + Nz;

            const int n2b = two2*Nz*Kw;
            const int n3a = two3*Mw*Nz;

#ifdef PRINTALGO_MM_BLOCK
            const int Ny = Nb-Nx;
            const int Nw = (j2-j1-Nz);
            std::cout<<"N > 1 split\n";
            std::cout<<"Mw = "<<Mw<<std::endl;
            std::cout<<"Nxyzw = "<<Nx<<','<<Ny<<','<<Nz<<','<<Nw<<std::endl;
            std::cout<<"Kw = "<<Kw<<std::endl;
            std::cout<<"n2b = "<<n2b<<std::endl;
            std::cout<<"n3a = "<<n3a<<std::endl;
#endif

            // Split matrix product up into:
            // [ A ] * [ E F ] = [ I J ]
            //         

            // AE = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB0,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,k1,i2,jm,k2,Mb,Nx,Kb,
                m1p,m2p,m3p, two1,two2,two3,
                firstm1,firstm2,firstm3,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assignb,assignb,
                prod,prod,prodb,prodb,
                mycleanup,kf,kf,kfb,kfb);
            // AF = J
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,jm,k1,i2,j2,k2,Mb,Nb-Nx,Kb,
                m1p,m2p+n2b,m3p+n3a, two1,two2,two3,
                false,firstm2,firstm3,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
        }
        else // Kb > 1
        {
            TMVAssert(Kb > 1 && Mb == 1 && Nb == 1);
            const int Kx = Kb>>1; // = Kb/2
            const int Kz = Kx<<lgKB;
            const int Mw = (i2-i1);
            const int Nw = (j2-j1);

            const int km = k1 + Kz;

            const int n1a = two1*Mw*Kz;
            const int n2a = two2*Nw*Kz;

#ifdef PRINTALGO_MM_BLOCK
            const int Ky = Kb-Kx;
            const int Kw = RoundUp<RT>(k2-k1-Kz);
            std::cout<<"N,K > 1 split\n";
            std::cout<<"Mw = "<<Mw<<std::endl;
            std::cout<<"Nw = "<<Nw<<std::endl;
            std::cout<<"Kxyzw = "<<Kx<<','<<Ky<<','<<Kz<<','<<Kw<<std::endl;
            std::cout<<"n1a = "<<n1a<<std::endl;
            std::cout<<"n2a = "<<n2a<<std::endl;
#endif

            // Split matrix product up into:
            // [ A B ] * [ E ] = [ I ]
            //           [ G ]          

            // AE = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,k1,i2,j2,km,Mb,Nb,Kx,
                m1p,m2p,m3p, two1,two2,two3,
                firstm1,firstm2,firstm3,false,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
            // BG = I
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB0,NB0,lgMB,lgNB,lgKB,
                i1,j1,km,i2,j2,k2,Mb,Nb,Kb-Kx,
                m1p+n1a,m2p+n2a,m3p, two1,two2,two3,
                firstm1,firstm2,false,lastm3,
                copy1,copy1_k,mycopy1,mycopy1_k,
                copy2,copy2_k,mycopy2,mycopy2_k,
                assign,assigna,assignb,myassign,
                prod,proda,prodb,myprod,
                mycleanup,kf,kfa,kfb,mykf);
        }
    }

    // Let's look at the order of the block operations in the algo 61
    // recursion schedule.  Let's say for now that all matrices are 
    // 4 blocks by 4 blocks.
    // For m1, we use capital letters, and for m2 lowercase:
    //
    // m1 = [ A  B  C  D ]     m2 = [ a  b  c  d ]
    //      [ E  F  G  H ]          [ e  f  g  h ]
    //      [ I  J  K  L ]          [ i  j  k  l ]
    //      [ M  N  O  P ]          [ m  n  o  p ]
    //
    // If we follow through the order of the splits, then we find that
    // the operations are (the rhs indicates the m3 location in i,j notation)
    //
    // Aa Be Ab Bf Ea Fe Eb Ff    00 00 01 01 10 10 11 11
    // Ci Dm Cj Dn Gi Hm Gj Hn    00 00 01 01 10 10 11 11
    // Ac Bg Ad Bh Ec Fg Ed Fh    02 02 03 03 12 12 13 13
    // Ck Do Cl Dp Gk Ho Gl Hp    02 02 03 03 12 12 13 13
    // Ia Je Ib Jf Ma Ne Mb Nf    20 20 21 21 30 30 31 31
    // Ki Lm Kj Ln Oi Pm Oj Pn    20 20 21 21 30 30 31 31
    // Ic Jg Id Jh Mc Ng Md Nh    22 22 23 23 32 32 33 33
    // Kk Lo Kl Lp Ok Po Ol Pp    22 22 23 23 32 32 33 33
    //
    // So we can see here that we can have better cache properties if 
    // the matrices are stored in the following order:
    //
    // m1 = [ A B E F C D G H I J M N K L O P ]
    // m2 = [ a e b f i m j n c g d h k o l p ]
    //
    // m3 should be stored with the same pattern as m1.
    //
    // Obviously with more than 4x4 blocks, this pattern should be continued
    // down into the components that were called a single letter here.
    // And non powers of two are a straightforward extension, especially
    // when we define the order recursively, so we can just use the 
    // same Mx = M/2 kind of expressions as we used above in the recursive
    // algorithm.
    //
    // We do not extend the idea to the edges where we have less than full
    // blocks.  Instead, we simply call the regular non-copying recursive
    // algorithm to clean up that portion of the calculation.  
    //
    // Finally, since we are copying to new storage anyway, we now have 
    // the flexibility to use whatever storage we want.  Since the 
    // RCC algorithm is the fastest, we copy m1 into rowmajor blocks,
    // and m2 and m3 into colmajor blocks.

    template <class CopyF, class AssignF, class ProdF, class CleanupF,
              class KF, class RT>
    static void DoRecursiveBlockMultMM2(
        const RT* x, const ConstMatrixView<RT>& m1,
        const ConstMatrixView<RT>& m2, MatrixView<RT> m3,
        const int M, const int N, const int K,
        const int MB, const int NB, const int KB,
        const int lgMB, const int lgNB, const int lgKB,
        const bool m1_z, const bool m2_z, const bool m3_z,
        CopyF* copy1, CopyF* copy1_k, CopyF* copy1b, CopyF* copy1b_k,
        CopyF* copy2, CopyF* copy2_k, CopyF* copy2a, CopyF* copy2a_k,
        AssignF* assign, AssignF* assigna, AssignF* assignb, AssignF* assignc, 
        ProdF* prod, ProdF* proda, ProdF* prodb, ProdF* prodc,
        CleanupF* cleanup, KF* kf, KF* kfa, KF* kfb, KF* kfc)
    {
        const int Mb = (M>>lgMB); // = M/MB
        const int Nb = (N>>lgNB); // = N/NB
        const int Kb = (K>>lgKB); // = K/KB
        const int Ma = (Mb<<lgMB); // = M/MB*MB
        const int Na = (Nb<<lgNB); // = N/NB*NB
        const int Ka = (Kb<<lgKB); // = K/KB*KB
        const int Mc = M-Ma; // = M%MB
        const int Nc = N-Na; // = N%NB
        const int Kc = K-Ka; // = K%KB
        // *bx and *cx are the values to pass to the DoRecursiveMultMM3
        // function.  They take accound of whether there is a partial block
        // that needs to be dealt with.  If there is we add one to the number
        // of block and use the real Mc or Nc.  But if not, we don't augment
        // the block count, and we use the full block size for Mc or Nc, which
        // indicate the size of the final block in the recursion.
        const int Mbx = Mc == 0 ? Mb : Mb+1;
        const int Nbx = Nc == 0 ? Nb : Nb+1;
        const int Mcx = Mc == 0 ? MB : Mc;
        const int Ncx = Nc == 0 ? NB : Nc;
        // Similiarly the following "x" functiosn pass the right function
        // according to what Mcx or Ncx is.
        CopyF* copy1bx = Mc==0 ? copy1 : copy1b;
        CopyF* copy1bx_k = Mc==0 ? copy1_k : copy1b_k;
        CopyF* copy2ax = Nc==0 ? copy2 : copy2a;
        CopyF* copy2ax_k = Nc==0 ? copy2_k : copy2a_k;
        AssignF* assignax = Nc==0 ? assign : assigna;
        AssignF* assignbx = Mc==0 ? assign : assignb;
        AssignF* assigncx = Mc==0 ? assignax : Nc==0 ? assignb : assignc;
        ProdF* prodax = Nc==0 ? prod : proda;
        ProdF* prodbx = Mc==0 ? prod : prodb;
        ProdF* prodcx = Mc==0 ? prodax : Nc==0 ? prodb : prodc;
        KF* kfax = Nc==0 ? kf : kfa;
        KF* kfbx = Mc==0 ? kf : kfb;
        KF* kfcx = Mc==0 ? kfax : Nc==0 ? kfb : kfc;

        const int Kd = RoundUp<RT>(Kc);
        // = Kc rounded up to multiple of 2 or 4 as required for SSE commands
        const int Ktot_d = (Kb<<lgKB) + Kd;
        TMVAssert(Ktot_d == RoundUp<RT>(K));

        const int two1 = m1_z ? 2 : 1;
        const int two2 = m2_z ? 2 : 1;
        const int two3 = m3_z ? 2 : 1;

#ifdef PRINTALGO_MM_BLOCK
        std::cout<<"RecursiveBlockMultMM2:\n";
        std::cout<<"m1,2,3_z = "<<m1_z<<','<<m2_z<<','<<m3_z<<std::endl;
        std::cout<<"two1,2,3 = "<<two1<<','<<two2<<','<<two3<<std::endl;
        std::cout<<"Ma,Na,Ka = "<<Ma<<','<<Na<<','<<Ka<<std::endl;
        std::cout<<"Mb,Nb,Kb = "<<Mb<<','<<Nb<<','<<Kb<<std::endl;
        std::cout<<"Mc,Nc,Kc,Kd = "<<Mc<<','<<Nc<<','<<Kc<<','<<Kd<<std::endl;
        std::cout<<"Mbx,Mcx,Nbx,Ncx = "<<Mbx<<','<<Mcx<<
            ','<<Nbx<<','<<Ncx<<std::endl;
        std::cout<<"Ktot_d = "<<Ktot_d<<std::endl;
#endif

        TMVAssert(Mb >= 2);
        TMVAssert(Nb >= 2);
        TMVAssert(Kb >= 2);

        if (Mb >= 3*std::max(Nb,Kb)/2) // M is largest
        {
            TMVAssert(Mb >= 3);
            TMVAssert(Nb >= 2 || Kb >= 2);
            int nsplit = Mb / (std::max(Nb,Kb)/2);
            const int Mx = (Mb-1)/nsplit+1;
            int My = Mb - (nsplit-1) * Mx;
            while (My <= 0) {
                // This doesn't happen often, but it is possible.
                // e.g. Mb = 19, Nb/2 = 3 -> nsplit = 6, Mx = 4, My = -1
                // Dropping nsplit by 1 solves the problem:
                // nsplit = 5, Mx = 4, My = 3
                --nsplit;
                My += Mx;
            }
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"nsplit = "<<nsplit<<std::endl;
            std::cout<<"Mx,My = "<<Mx<<','<<My<<std::endl;
#endif
            int i1m = 0;
            TMVAssert(nsplit >= 2);
            int Mz = Mx<<lgMB;

            int size1 = two1*Mz*Ktot_d;
            int size2 = two2*Ktot_d*N;
            int size3 = two3*Mz*N;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"m1 size = "<<size1<<std::endl;
            std::cout<<"m2 size = "<<size2<<std::endl;
            std::cout<<"m3 size = "<<size3<<std::endl;
#endif
            AlignedArray<RT> m1_temp(size1);
            AlignedArray<RT> m2_temp(size2);
            AlignedArray<RT> m3_temp(size3);
            RT*const m1p = m1_temp;
            RT*const m2p = m2_temp;
            RT*const m3p = m3_temp;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"M >= 3*max(N,K)/2\n";
            std::cout<<"m1p = "<<m1p<<" end = "<<(m1p+size1)<<std::endl;
            std::cout<<"m2p = "<<m2p<<" end = "<<(m2p+size2)<<std::endl;
            std::cout<<"m3p = "<<m3p<<" end = "<<(m3p+size3)<<std::endl;
#endif
#ifdef TEST_POINTERS
            glob_m1p_end = m1p+size1;
            glob_m2p_end = m2p+size2;
            glob_m3p_end = m3p+size3;
#endif

            for(int n=0;n<nsplit-1;++n) {
                const int i2m = i1m + Mz;
                DoRecursiveBlockMultMM3(
                    x,m1,m2,m3,
                    MB,Ncx,KB,MB,NB,lgMB,lgNB,lgKB,
                    i1m,0,0,i2m,N,K,Mx,Nbx,Kb,
                    m1p,m2p,m3p, two1,two2,two3,
                    true,n==0,true,true,
                    copy1,copy1_k,copy1,copy1_k,
                    copy2,copy2_k,copy2ax,copy2ax_k,
                    assign,assignax,assign,assignax,
                    prod,prodax,prod,prodax,
                    cleanup,kf,kfax,kf,kfax);
                i1m = i2m;
            }
            const int i2m = i1m + (My<<lgMB);
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,Ncx,KB,MB,NB,lgMB,lgNB,lgKB,
                i1m,0,0,i2m,N,K,My,Nbx,Kb,
                m1p,m2p,m3p, two1,two2,two3,
                true,false,true,true,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,copy2ax,copy2ax_k,
                assign,assignax,assign,assignax,
                prod,prodax,prod,prodax,
                cleanup,kf,kfax,kf,kfax);
            i1m = i2m;
            if (Mc)
                DoRecursiveBlockMultMM3(
                    x,m1,m2,m3,
                    Mc,Ncx,KB,MB,NB,lgMB,lgNB,lgKB,
                    i1m,0,0,M,N,K,1,Nbx,Kb,
                    m1p,m2p,m3p, two1,two2,two3,
                    true,false,true,true,
                    copy1,copy1_k,copy1b,copy1b_k,
                    copy2,copy2_k,copy2ax,copy2ax_k,
                    assign,assignax,assignbx,assigncx,
                    prod,prodax,prodbx,prodcx,
                    cleanup,kf,kfax,kfbx,kfcx);
        }
        else if (Nb >= 3*std::max(Mb,Kb)/2) // N is largest
        {
            TMVAssert(Nb >= 3);
            TMVAssert(Mb >= 2 || Kb >= 2);
            int nsplit = Nb / (std::max(Mb,Kb)/2);
            const int Nx = (Nb-1)/nsplit+1; 
            int Ny = Nb - (nsplit-1) * Nx;
            while (Ny <= 0) { --nsplit; Ny += Nx; }
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"nsplit = "<<nsplit<<std::endl;
            std::cout<<"Nx,Ny = "<<Nx<<','<<Ny<<std::endl;
#endif
            int j1m = 0;
            TMVAssert(nsplit >= 2);
            int Nz = Nx<<lgNB;

            int size1 = two1*M*Ktot_d;
            int size2 = two2*Ktot_d*Nz;
            int size3 = two3*M*Nz;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"m1 size = "<<size1<<std::endl;
            std::cout<<"m2 size = "<<size2<<std::endl;
            std::cout<<"m3 size = "<<size3<<std::endl;
#endif
            AlignedArray<RT> m1_temp(size1);
            AlignedArray<RT> m2_temp(size2);
            AlignedArray<RT> m3_temp(size3);
            RT*const m1p = m1_temp;
            RT*const m2p = m2_temp;
            RT*const m3p = m3_temp;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"N >= 3*max(M,K)/2\n";
            std::cout<<"m1p = "<<m1p<<" end = "<<(m1p+size1)<<std::endl;
            std::cout<<"m2p = "<<m2p<<" end = "<<(m2p+size2)<<std::endl;
            std::cout<<"m3p = "<<m3p<<" end = "<<(m3p+size3)<<std::endl;
#endif
#ifdef TEST_POINTERS
            glob_m1p_end = m1p+size1;
            glob_m2p_end = m2p+size2;
            glob_m3p_end = m3p+size3;
#endif

            for(int n=0;n<nsplit-1;++n) {
                const int j2m = j1m + Nz;
                DoRecursiveBlockMultMM3(
                    x,m1,m2,m3,
                    Mcx,NB,KB,MB,NB,lgMB,lgNB,lgKB,
                    0,j1m,0,M,j2m,K,Mbx,Nx,Kb,
                    m1p,m2p,m3p, two1,two2,two3,
                    n==0,true,true,true,
                    copy1,copy1_k,copy1bx,copy1bx_k,
                    copy2,copy2_k,copy2,copy2_k,
                    assign,assign,assignbx,assignbx,
                    prod,prod,prodbx,prodbx,
                    cleanup,kf,kf,kfbx,kfbx);
                j1m = j2m;
            }
            const int j2m = j1m + (Ny<<lgNB);
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                Mcx,NB,KB,MB,NB,lgMB,lgNB,lgKB,
                0,j1m,0,M,j2m,K,Mbx,Ny,Kb,
                m1p,m2p,m3p, two1,two2,two3,
                false,true,true,true,
                copy1,copy1_k,copy1bx,copy1bx_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assignbx,assignbx,
                prod,prod,prodbx,prodbx,
                cleanup,kf,kf,kfbx,kfbx);
            j1m = j2m;
            if (Nc)
                DoRecursiveBlockMultMM3(
                    x,m1,m2,m3,
                    Mcx,Nc,KB,MB,NB,lgMB,lgNB,lgKB,
                    0,j1m,0,M,N,K,Mbx,1,Kb,
                    m1p,m2p,m3p, two1,two2,two3,
                    false,true,true,true,
                    copy1,copy1_k,copy1bx,copy1bx_k,
                    copy2,copy2_k,copy2ax,copy2ax_k,
                    assign,assignax,assignbx,assigncx,
                    prod,prodax,prodbx,prodcx,
                    cleanup,kf,kfax,kfbx,kfcx);
        }
        else if (Kb >= 2*std::max(Mb,Nb)) // K is largest
        {
            TMVAssert(Kb >= 4);
            TMVAssert(Mb >= 2 || Nb >= 2);
            int nsplit = Kb / std::max(Mb,Nb);
            const int Kx = (Kb-1)/nsplit+1;
            int Ky = Kb - (nsplit-1) * Kx;
            while (Ky <= 0) { --nsplit; Ky += Kx; }
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"nsplit = "<<nsplit<<std::endl;
            std::cout<<"Kx,Ky = "<<Kx<<','<<Ky<<std::endl;
#endif
            int k1m = 0;
            TMVAssert(nsplit >= 2);
            int Kz = Kx<<lgKB;

            // Make room for last partial block if necessary
            int Kzz = (Ky == Kx) ? Kz + Kd : Kz;

            int size1 = two1*M*Kzz;
            int size2 = two2*Kzz*N;
            int size3 = two3*M*N;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"m1 size = "<<size1<<std::endl;
            std::cout<<"m2 size = "<<size2<<std::endl;
            std::cout<<"m3 size = "<<size3<<std::endl;
#endif
            AlignedArray<RT> m1_temp(size1);
            AlignedArray<RT> m2_temp(size2);
            AlignedArray<RT> m3_temp(size3);
            RT*const m1p = m1_temp;
            RT*const m2p = m2_temp;
            RT*const m3p = m3_temp;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"K >= 2*max(M,N)\n";
            std::cout<<"m1p = "<<m1p<<" end = "<<(m1p+size1)<<std::endl;
            std::cout<<"m2p = "<<m2p<<" end = "<<(m2p+size2)<<std::endl;
            std::cout<<"m3p = "<<m3p<<" end = "<<(m3p+size3)<<std::endl;
#endif
#ifdef TEST_POINTERS
            glob_m1p_end = m1p+size1;
            glob_m2p_end = m2p+size2;
            glob_m3p_end = m3p+size3;
#endif

            for(int n=0;n<nsplit-1;++n) {
                const int k2m = k1m + Kz;
                DoRecursiveBlockMultMM3(
                    x,m1,m2,m3,
                    Mcx,Ncx,KB,MB,NB,lgMB,lgNB,lgKB,
                    0,0,k1m,M,N,k2m,Mbx,Nbx,Kx,
                    m1p,m2p,m3p, two1,two2,two3,
                    true,true,n==0,false,
                    copy1,copy1_k,copy1bx,copy1bx_k,
                    copy2,copy2_k,copy2ax,copy2ax_k,
                    assign,assignax,assignbx,assigncx,
                    prod,prodax,prodbx,prodcx,
                    cleanup,kf,kfax,kfbx,kfcx);
                k1m = k2m;
            }
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                Mcx,Ncx,KB,MB,NB,lgMB,lgNB,lgKB,
                0,0,k1m,M,N,K,Mbx,Nbx,Ky,
                m1p,m2p,m3p, two1,two2,two3,
                true,true,false,true,
                copy1,copy1_k,copy1bx,copy1bx_k,
                copy2,copy2_k,copy2ax,copy2ax_k,
                assign,assignax,assignbx,assigncx,
                prod,prodax,prodbx,prodcx,
                cleanup,kf,kfax,kfbx,kfcx);
        }
        else // None is much larger than the others (or at least one other)
        {
            TMVAssert(Mb >= 2);
            TMVAssert(Nb >= 2);
            TMVAssert(Kb >= 2);
            const int Mx = Mb-(Mb>>1); // = M/2, rounding up if odd
            const int Nx = Nb-(Nb>>1); // = N/2, rounding up if odd
            const int Kx = Kb-(Kb>>1); // = K/2, rounding up if odd
            const int im = Mx<<lgMB;
            const int jm = Nx<<lgNB;
            const int km = Kx<<lgKB;
            const int Mz = Mx<<lgMB;
            const int Nz = Nx<<lgNB;
            const int Kz = Kx<<lgKB;
            int Mw = (M-Mz);
            int Nw = (N-Nz);
            int Kw = RoundUp<RT>(K-Kz);
            int size1 = two1*std::max(Mz,Mw)*Ktot_d;
            int size2 = two2*Ktot_d*N;
            int size3 = two3*std::max(Mz,Mw)*std::max(Nz,Nw);
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"Mx,Nx,Kx = "<<Mx<<','<<Nx<<','<<Kx<<std::endl;
            std::cout<<"Mz,Nz,Kz = "<<Mz<<','<<Nz<<','<<Kz<<std::endl;
            std::cout<<"Mw,Nw,Kw = "<<Mw<<','<<Nw<<','<<Kw<<std::endl;
            std::cout<<"im,jm,km = "<<im<<','<<jm<<','<<km<<std::endl;
            std::cout<<"m1 size = "<<size1<<std::endl;
            std::cout<<"m2 size = "<<size2<<std::endl;
            std::cout<<"m3 size = "<<size3<<std::endl;
#endif
            AlignedArray<RT> m1_temp(size1);
            AlignedArray<RT> m2_temp(size2);
            AlignedArray<RT> m3_temp(size3);
            RT*const m1p = m1_temp;
            RT*const m2p = m2_temp;
            RT*const m3p = m3_temp;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"None much larger than the others\n";
            std::cout<<"m1p = "<<m1p<<" end = "<<(m1p+size1)<<std::endl;
            std::cout<<"m2p = "<<m2p<<" end = "<<(m2p+size2)<<std::endl;
            std::cout<<"m3p = "<<m3p<<" end = "<<(m3p+size3)<<std::endl;
            std::cout<<"two1,2,3 = "<<two1<<','<<two2<<','<<two3<<std::endl;
#endif
#ifdef TEST_POINTERS
            glob_m1p_end = m1p+size1;
            glob_m2p_end = m2p+size2;
            glob_m3p_end = m3p+size3;
            std::cout<<"end of m2 = "<<
                ConstVectorView<long,1>((const long*)glob_m2p_end,10)
                <<std::endl;
#endif

            const int n1a = two1*Mz*Kz;
            const int n1c = two1*Mw*Kz;

            const int n2a = two2*Nz*Kz;
            const int n2b = n2a + two2*Nz*Kw;
            const int n2c = n2b + two2*Nw*Kz;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"n1ac = "<<n1a<<','<<n1c<<std::endl;
            std::cout<<"n2abc = "<<n2a<<','<<n2b<<','<<n2c<<std::endl;
#endif

            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB,NB,lgMB,lgNB,lgKB,
                0,0,0,im,jm,km,Mx,Nx,Kx,
                m1p,m2p,m3p, two1,two2,two3,
                true,true,true,false,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assign,assign,
                prod,prod,prod,prod,
                cleanup,kf,kf,kf,kf);
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,NB,KB,MB,NB,lgMB,lgNB,lgKB,
                0,0,km,im,jm,K,Mx,Nx,Kb-Kx,
                m1p+n1a,m2p+n2a,m3p, two1,two2,two3,
                true,true,false,true,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assign,assign,
                prod,prod,prod,prod,
                cleanup,kf,kf,kf,kf);
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,Ncx,KB,MB,NB,lgMB,lgNB,lgKB,
                0,jm,0,im,N,km,Mx,Nbx-Nx,Kx,
                m1p,m2p+n2b,m3p, two1,two2,two3,
                false,true,true,false,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,copy2ax,copy2ax_k,
                assign,assignax,assign,assignax,
                prod,prodax,prod,prodax,
                cleanup,kf,kfax,kf,kfax);
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                MB,Ncx,KB,MB,NB,lgMB,lgNB,lgKB,
                0,jm,km,im,N,K,Mx,Nbx-Nx,Kb-Kx,
                m1p+n1a,m2p+n2c,m3p, two1,two2,two3,
                false,true,false,true,
                copy1,copy1_k,copy1,copy1_k,
                copy2,copy2_k,copy2ax,copy2ax_k,
                assign,assignax,assign,assignax,
                prod,prodax,prod,prodax,
                cleanup,kf,kfax,kf,kfax);
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                Mcx,NB,KB,MB,NB,lgMB,lgNB,lgKB,
                im,0,0,M,jm,km,Mbx-Mx,Nx,Kx,
                m1p,m2p,m3p, two1,two2,two3,
                true,false,true,false,
                copy1,copy1_k,copy1bx,copy1bx_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assignbx,assignbx,
                prod,prod,prodbx,prodbx,
                cleanup,kf,kf,kfbx,kfbx);
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                Mcx,NB,KB,MB,NB,lgMB,lgNB,lgKB,
                im,0,km,M,jm,K,Mbx-Mx,Nx,Kb-Kx,
                m1p+n1c,m2p+n2a,m3p, two1,two2,two3,
                true,false,false,true,
                copy1,copy1_k,copy1bx,copy1bx_k,
                copy2,copy2_k,copy2,copy2_k,
                assign,assign,assignbx,assignbx,
                prod,prod,prodbx,prodbx,
                cleanup,kf,kf,kfbx,kfbx);
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                Mcx,Ncx,KB,MB,NB,lgMB,lgNB,lgKB,
                im,jm,0,M,N,km,Mbx-Mx,Nbx-Nx,Kx,
                m1p,m2p+n2b,m3p, two1,two2,two3,
                false,false,true,false,
                copy1,copy1_k,copy1bx,copy1bx_k,
                copy2,copy2_k,copy2ax,copy2ax_k,
                assign,assignax,assignbx,assigncx,
                prod,prodax,prodbx,prodcx,
                cleanup,kf,kfax,kfbx,kfcx);
            DoRecursiveBlockMultMM3(
                x,m1,m2,m3,
                Mcx,Ncx,KB,MB,NB,lgMB,lgNB,lgKB,
                im,jm,km,M,N,K,Mbx-Mx,Nbx-Nx,Kb-Kx,
                m1p+n1c,m2p+n2c,m3p, two1,two2,two3,
                false,false,false,true,
                copy1,copy1_k,copy1bx,copy1bx_k,
                copy2,copy2_k,copy2ax,copy2ax_k,
                assign,assignax,assignbx,assigncx,
                prod,prodax,prodbx,prodcx,
                cleanup,kf,kfax,kfbx,kfcx);
        }
    }

    template <int algo, int cs, int rs, int xs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_RecursiveBlock_Helper;

    // algo 1: The normal recursive block algorithm.
    template <int cs, int rs, int xs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_RecursiveBlock_Helper<1,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Traits2<T1,T2>::type PT3;
            typedef typename M3::real_type RT;
            const bool x_z = Traits<T>::iscomplex;
            const bool m1_z = Traits<T1>::iscomplex;
            const bool m2_z = Traits<T2>::iscomplex;
            const bool m3_z = Traits<PT3>::iscomplex;
            const int MB    = 16;
            const int NB    = 16;
#if 0
            const int KB    = 16;
#else
            const int KB    = 
#ifdef __SSE__
                Traits2<RT,float>::sametype ? (
                    ( xs == UNKNOWN || xs >= 128 ) ? 64 :
                    ( xs >= 96 ) ? 32 :
                    ( xs >= 64 ) ? 64 :
                    ( xs >= 48 ) ? 16 :
                    ( xs >= 32 ) ? 32 :
                    16 ) :
#endif
#ifdef __SSE2__
                Traits2<RT,double>::sametype ? (
                    ( xs == UNKNOWN || xs >= 64 ) ? 32 :
                    ( xs >= 48 ) ? 16 :
                    ( xs >= 32 ) ? 32 :
                    16 ) :
#endif
                16;
#endif
            const int lgMB  = IntTraits<MB>::log;
            const int lgNB  = IntTraits<NB>::log;
            const int lgKB  = IntTraits<KB>::log;
            const int csx = cs == UNKNOWN ? UNKNOWN : (cs-((cs>>lgMB)<<lgMB));
            const int rsx = rs == UNKNOWN ? UNKNOWN : (rs-((rs>>lgNB)<<lgNB));
            const int xsx = xs == UNKNOWN ? UNKNOWN : (xs-((xs>>lgKB)<<lgKB));
            typedef typename M1::const_submatrix_type M1sub;
            typedef typename M1sub::const_transpose_type M1sub_t;
            typedef typename M2::const_submatrix_type M2sub;
            typedef typename M3::submatrix_type M3sub;

            // We will determine the value of K once and call a known K
            // algorithm for the K cleanup calls.
            typedef void KF(
                const int M, const int N, const int K,
                const Scaling<1,RT>& x, const RT* A, const RT* B, RT* C);

            typedef void CopyF(const ConstMatrixView<RT>& m1x,
                               RT* m2p, int M, int N, int si, int sj);
            typedef void AssignF(
                const RT* x, 
                RT* m1p, const int M, const int N, const int si1, const int sj1,
                MatrixView<RT> m2x);
            typedef void ProdF(
                const int M, const int N, const int K,
                RT* m1p, const int si1, const int sj1,
                RT* m2p, const int si2, const int sj2,
                RT* m3p, const int si3, const int sj3);
            typedef void CleanupF(
                const int M, const int N, const int K,
                RT* m1p, const int si1, const int sj1,
                RT* m2p, const int si2, const int sj2,
                RT* m3p, const int si3, const int sj3,
                KF* kf);

            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;

            const int Mb = (M>>lgMB); // = M/MB
            const int Nb = (N>>lgNB); // = N/NB
            const int Kb = (K>>lgKB); // = K/KB

            // The following requirement should be checked by algo 71.
            // It should only call this function when these relations are true.
            TMVAssert(Mb>=2);
            TMVAssert(Nb>=2);
            TMVAssert(Kb>=2);
            // TODO: Make sure this is ok.

            const int Mc = M-(Mb<<lgMB); // = M%MB
            const int Nc = N-(Nb<<lgNB); // = N%NB
            const int Kc = K-(Kb<<lgKB); // = K%KB

            KF* kf = 0;
            KF* kfa = 0;
            KF* kfb = 0;
            KF* kfc = 0;
            if (Kc) {
                kf = get_Kcleanup<xsx,KB,RT>::call(Kc);
                if (Nc)
                    kfa = &call_multmm_M_N_K<1,RT>;
                if (Mc)
                    kfb = &call_multmm_M_N_K<1,RT>;
                if (Mc && Nc)
                    kfc = &call_multmm_M_N_K<1,RT>;
            }

            CopyF* copy1 = &MyCopy<m1_z,KB,MB,M1sub_t>::call;
            CopyF* copy1_k = &MyCopy<m1_z,xsx,MB,M1sub_t>::call;
            CopyF* copy1b = &MyCopy<m1_z,KB,csx,M1sub_t>::call;
            CopyF* copy1b_k = &MyCopy<m1_z,xsx,csx,M1sub_t>::call;
            CopyF* copy2 = &MyCopy<m2_z,KB,NB,M2sub>::call;
            CopyF* copy2_k = &MyCopy<m2_z,xsx,NB,M2sub>::call;
            CopyF* copy2a = &MyCopy<m2_z,KB,rsx,M2sub>::call;
            CopyF* copy2a_k = &MyCopy<m2_z,xsx,rsx,M2sub>::call;

            AssignF* assign = &MyAssign<x_z,m3_z,add,ix,T,MB,NB,M3sub>::call;
            AssignF* assigna = &MyAssign<x_z,m3_z,add,ix,T,MB,rsx,M3sub>::call;
            AssignF* assignb = &MyAssign<x_z,m3_z,add,ix,T,csx,NB,M3sub>::call;
            AssignF* assignc = &MyAssign<x_z,m3_z,add,ix,T,csx,rsx,M3sub>::call;
            T xx = T(x);

            ProdF* prod = &MyProd<m1_z,m2_z,MB,NB,KB,RT>::call;
            ProdF* proda = &MyProd<m1_z,m2_z,MB,rsx,KB,RT>::call;
            ProdF* prodb = &MyProd<m1_z,m2_z,csx,NB,KB,RT>::call;
            ProdF* prodc = &MyProd<m1_z,m2_z,csx,rsx,KB,RT>::call;

            CleanupF* cleanup = &MyCleanup<m1_z,m2_z,RT>::call;

            DoRecursiveBlockMultMM2(
                (RT*)(&xx),
                m1.realPart().xView(),m2.realPart().xView(),
                m3.realPart().xView(),
                M,N,K, MB,NB,KB, lgMB,lgNB,lgKB, m1_z,m2_z,m3_z,
                copy1,copy1_k,copy1b,copy1b_k,
                copy2,copy2_k,copy2a,copy2a_k,
                assign,assigna,assignb,assignc,
                prod,proda,prodb,prodc,
                cleanup,kf,kfa,kfb,kfc);
        }
    };

    // algo 2: m1, m2, or m3 is nomajor (but might really be rm or cm).
    // Use get_copy and get_assign to figure out which copy or assign
    // function to use.
    template <int cs, int rs, int xs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_RecursiveBlock_Helper<2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Traits2<T1,T2>::type PT3;
            typedef typename M3::real_type RT;

            const bool x_z = Traits<T>::iscomplex;
            const bool m1_z = Traits<T1>::iscomplex;
            const bool m2_z = Traits<T2>::iscomplex;
            const bool m3_z = Traits<PT3>::iscomplex;
            const int MB    = 16;
            const int NB    = 16;
#if 0
            const int KB    = 16;
#else
            const int KB    = 
#ifdef __SSE__
                Traits2<RT,float>::sametype ? 64 :
#endif
#ifdef __SSE2__
                Traits2<RT,double>::sametype ? 32 :
#endif
                16;
#endif
            const int lgMB  = IntTraits<MB>::log;
            const int lgNB  = IntTraits<NB>::log;
            const int lgKB  = IntTraits<KB>::log;

            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_view_type M2v;
            typedef typename M3::view_type M3v;

            typedef void KF(
                const int M, const int N, const int K,
                const Scaling<1,RT>& x, const RT* A, const RT* B, RT* C);

            typedef void CopyF(const ConstMatrixView<RT>& m1x,
                               RT* m2p, int M, int N, int si, int sj);
            typedef void AssignF(
                const RT* x, 
                RT* m1p, const int M, const int N, const int si1, const int sj1,
                MatrixView<RT> m2x);
            typedef void ProdF(
                const int M, const int N, const int K,
                RT* m1p, const int si1, const int sj1,
                RT* m2p, const int si2, const int sj2,
                RT* m3p, const int si3, const int sj3);
            typedef void CleanupF(
                const int M, const int N, const int K,
                RT* m1p, const int si1, const int sj1,
                RT* m2p, const int si2, const int sj2,
                RT* m3p, const int si3, const int sj3,
                KF* kf);

            const int M = m3.colsize();
            const int N = m3.rowsize();
            const int K = m1.rowsize();

            const int XX = UNKNOWN;

            const int Mb = (M>>lgMB); // = M/MB
            const int Nb = (N>>lgNB); // = N/NB
            const int Kb = (K>>lgKB); // = K/KB

            TMVAssert(Mb>=2);
            TMVAssert(Nb>=2);
            TMVAssert(Kb>=2);

            const int Mc = M-(Mb<<lgMB); // = M%MB
            const int Nc = N-(Nb<<lgNB); // = N%NB
            const int Kc = K-(Kb<<lgKB); // = K%KB

            KF* kf = 0;
            KF* kfa = 0;
            KF* kfb = 0;
            KF* kfc = 0;
            if (Kc) {
                kf = get_Kcleanup<XX,KB,RT>::call(Kc);
                if (Nc)
                    kfa = &call_multmm_M_N_K<1,RT>;
                if (Mc)
                    kfb = &call_multmm_M_N_K<1,RT>;
                if (Mc && Nc)
                    kfc = &call_multmm_M_N_K<1,RT>;
            }

            CopyF* copy1 = get_copy<m1_z,KB,MB,M1t>::call(m1.transpose());
            CopyF* copy1_k = get_copy<m1_z,XX,MB,M1t>::call(m1.transpose());
            CopyF* copy1b = get_copy<m1_z,KB,XX,M1t>::call(m1.transpose());
            CopyF* copy1b_k = get_copy<m1_z,XX,XX,M1t>::call(m1.transpose());
            CopyF* copy2 = get_copy<m2_z,KB,NB,M2v>::call(m2);
            CopyF* copy2_k = get_copy<m2_z,XX,NB,M2v>::call(m2);
            CopyF* copy2a = get_copy<m2_z,KB,XX,M2v>::call(m2);
            CopyF* copy2a_k = get_copy<m2_z,XX,XX,M2v>::call(m2);

            AssignF* assign = get_assign<x_z,m3_z,T,MB,NB,M3v>::call(add,x,m3);
            AssignF* assigna = get_assign<x_z,m3_z,T,MB,XX,M3v>::call(add,x,m3);
            AssignF* assignb = get_assign<x_z,m3_z,T,XX,NB,M3v>::call(add,x,m3);
            AssignF* assignc = get_assign<x_z,m3_z,T,XX,XX,M3v>::call(add,x,m3);

            T xx = T(x);

            ProdF* prod = &MyProd<m1_z,m2_z,MB,NB,KB,RT>::call;
            ProdF* proda = &MyProd<m1_z,m2_z,MB,XX,KB,RT>::call;
            ProdF* prodb = &MyProd<m1_z,m2_z,XX,NB,KB,RT>::call;
            ProdF* prodc = &MyProd<m1_z,m2_z,XX,XX,KB,RT>::call;

            CleanupF* cleanup = &MyCleanup<m1_z,m2_z,RT>::call;

            DoRecursiveBlockMultMM2(
                (RT*)(&xx),
                m1.realPart().xView(),m2.realPart().xView(),
                m3.realPart().xView(),
                M,N,K, MB,NB,KB, lgMB,lgNB,lgKB, m1_z,m2_z,m3_z,
                copy1,copy1_k,copy1b,copy1b_k,
                copy2,copy2_k,copy2a,copy2a_k,
                assign,assigna,assignb,assignc,
                prod,proda,prodb,prodc,
                cleanup,kf,kfa,kfb,kfc);
        }
    };

    // algo -2: Call correct version depending on whether majority is known
    template <int cs, int rs, int xs, bool add,
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_RecursiveBlock_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo =
                ( (M1::_rowmajor || M1::_colmajor) &&
                  (M2::_rowmajor || M2::_colmajor) &&
                  (M3::_rowmajor || M3::_colmajor) ) ? 1 : 2;
            MultMM_RecursiveBlock_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
        }
    };

    // algo 98: Call inst
    template <int cs, int rs, int xs,
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_RecursiveBlock_Helper<98,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstMultMM_RecursiveBlock(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int xs,
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_RecursiveBlock_Helper<98,cs,rs,xs,true,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstAddMultMM_RecursiveBlock(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo -1: Check for inst
    template <int cs, int rs, int xs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_RecursiveBlock_Helper<-1,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst =
                M1::unknownsizes &&
                M2::unknownsizes &&
                M3::unknownsizes &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                inst ? 98 :
                -2;
            MultMM_RecursiveBlock_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_RecursiveBlock(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        MultMM_RecursiveBlock_Helper<-1,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(
            x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM_RecursiveBlock(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        MultMM_RecursiveBlock_Helper<-2,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(
            x,m1v,m2v,m3v);
    }



    //
    // Algo 64: BlockMultMM
    //

    // Algo 64 is similar to 63, but not recursive divide and conquer.
    // Instead we loop over M then N (or N then M) then K in blocks.
    //
    // This order of operations also reduces the memory required for
    // copying into blocks.  We still copy all of one of the matrices
    // (either m1 or m3 depending on the values of M,N), but then
    // we only need one block row of m1 or one block column of m2 that
    // we reuse over and over.  Also, we only need a single temporary
    // block for M3 that we reuse each time as well.
    //
    // Finally, since this ordering is so simple, it is not too hard to 
    // incorporate the edge calculations into the same model.  So rather
    // than call some other algorithm for the edges, we do the edges here.


    // DoBlockMultMM3 is the routine that does the inner K loop of the below
    // BlockMultMM calculation.
    //
    // It calculates a single MBxNB block of m3 using a temporary
    // compact matrix, and then stores it into the real m3 when the
    // calculation is done.
    //
    // x,m1,m2,m3 are the real matrices being worked on (and the scaling).
    // i1, j1, i2, j2 indicate which block of m3 we are calculating.
    // MB,NB are the block sizes.
    // Kb is the number of blocks in the K direction to loop over.
    // m1p, m2p are temporary storage to use for m1, m2 blocks
    // firstm1, firstm2 say whether this is the first time we are 
    //   accessing these blocks of m1, m2, and therefore whether the
    //   data need to be loaded into the temporary matrices.
    // m3s is the temporary compact matrix to store the result.
    // cleanup is the function to call for the the final partial block
    //   after the end of the full-block loop.
    // k1x, k2x are the k indices for the final partial block.
    // Kc is the k2x-k1x, and Kd is the stepj to use for the partial block.
    //
    // I originally wrote this for MB, NB, etc., but the edge cleanup
    // works exactly the same way, so the way we chose to use the same 
    // code for the edges is to pass in the functions to call for each
    // aspect of the calculation.  This has the advantage of avoiding a 
    // lot of different template instantiations of what amount to really
    // the same algorithm.  So it greatly decreases the compiled size
    // compared to what we would get by doing this with templates, and
    // it doesn't appreciably reduce the speed.
    //
    // Also, to reduce the compiled library size, I write all this 
    // in terms of just RT, rather than T.  So m1,m2,m3 are the realPart()
    // portions of the full m1,m2,m3.  In the functions calls that 
    // do each part of the calculation, the variables are recast as complex
    // if appropriate.
    template <class CopyF, class AssignF, class ProdF, class CleanupF,
              class KF, class RT>
    static void DoBlockMultMM3(
        const RT* x, const ConstMatrixView<RT>& m1,
        const ConstMatrixView<RT>& m2, MatrixView<RT> m3,
        const int i1, const int j1, const int i2, const int j2,
        const int MB, const int NB, const int KB, const int Kb,
        const int size1, const int size2, const int size3,
        RT* m1p, RT* m2p, RT* m3p, bool firstm1, bool firstm2, 
        CopyF* mycopy1, CopyF* mycopy2, CopyF* mycopy1_k, CopyF* mycopy2_k,
        AssignF* myassign, ProdF* myprod,
        CleanupF* mycleanup, KF* kf,
        const int Ka, const int K, const int Kc, const int Kd)
    {
#ifdef PRINTALGO_MM_BLOCK
        std::cout<<"DoBlockMultMM3\n";
        std::cout<<"i1,j1 = "<<i1<<','<<j1<<std::endl;
        std::cout<<"i2,j2 = "<<i2<<','<<j2<<std::endl;
        std::cout<<"MB,NB,Kb = "<<MB<<','<<NB<<','<<Kb<<std::endl;
        std::cout<<"firstm1 = "<<firstm1<<", firstm2 = "<<firstm2<<std::endl;
        std::cout<<"Ka,K,Kc,Kd = "<<Ka<<','<<K<<','<<Kc<<','<<Kd<<std::endl;
        std::cout<<"m1p = "<<m1p<<"  "<<*m1p<<std::endl;
        std::cout<<"m2p = "<<m2p<<"  "<<*m2p<<std::endl;
        std::cout<<"m3p = "<<m3p<<"  "<<*m3p<<std::endl;
        std::cout<<"sizes = "<<size1<<"  "<<size2<<"  "<<size3<<std::endl;
#endif

        VectorView<RT,1> m3x(m3p,size3,1);
        m3x.setZero();
        for (int k=0;k<Kb;++k,m1p+=size1,m2p+=size2)
        {
            if (firstm1) 
                (*mycopy1) (m1.cSubMatrix(i1,i2,k*KB,(k+1)*KB).transpose(),
                            m1p,KB,MB,1,KB);
            if (firstm2) 
                (*mycopy2) (m2.cSubMatrix(k*KB,(k+1)*KB,j1,j2),m2p,KB,NB,1,KB);
            (*myprod) (MB,NB,KB, m1p,KB,1, m2p,1,KB, m3p,1,MB);
        }

        if (kf) {
            if (firstm1) 
                (*mycopy1_k) (m1.cSubMatrix(i1,i2,Ka,K).transpose(),
                              m1p,Kc,MB,1,Kd);
            if (firstm2) 
                (*mycopy2_k) (m2.cSubMatrix(Ka,K,j1,j2),m2p,Kc,NB,1,Kd);
            (*mycleanup) (MB,NB,Kc, m1p,Kd,1, m2p,1,Kd, m3p,1,MB, kf);
        }
        (*myassign) (x,m3p,MB,NB,1,MB,m3.cSubMatrix(i1,i2,j1,j2));
    }

    // This is the meat of the calculation for algo 64.
    // We loop over M and N in blocks, but do the full length of K
    // for that block all at once.  This tends to be better than the 
    // above recursive block algorith for medium-sized matrices.
    template <class CopyF, class AssignF, class ProdF, class CleanupF,
              class KF, class RT>
    static void DoBlockMultMM2(
        const RT* x, const ConstMatrixView<RT>& m1,
        const ConstMatrixView<RT>& m2, MatrixView<RT> m3,
        const int M, const int N, const int K,
        const int MB, const int NB, const int KB,
        const int lgMB, const int lgNB, const int lgKB,
        const bool m1_z, const bool m2_z, const bool m3_z,
        CopyF* copy1, CopyF* copy1_k, CopyF* copy1b, CopyF* copy1b_k,
        CopyF* copy2, CopyF* copy2_k, CopyF* copy2a, CopyF* copy2a_k,
        AssignF* assign, AssignF* assigna, AssignF* assignb, AssignF* assignc, 
        ProdF* prod, ProdF* proda, ProdF* prodb, ProdF* prodc,
        CleanupF* cleanup, KF* kf, KF* kfa, KF* kfb, KF* kfc)
    {
        const int Mb = (M>>lgMB); // = M/MB
        const int Nb = (N>>lgNB); // = N/NB
        const int Kb = (K>>lgKB); // = K/KB
        const int Ma = (Mb<<lgMB); // = M/MB*MB
        const int Na = (Nb<<lgNB); // = N/NB*NB
        const int Ka = (Kb<<lgKB); // = K/KB*KB
        const int Mc = M-(Mb<<lgMB); // = M%MB
        const int Nc = N-(Nb<<lgNB); // = N%NB
        const int Kc = K-(Kb<<lgKB); // = K%KB
        const int Kd = RoundUp<RT>(Kc);
        // = Kc rounded up to multiple of 2 or 4 as required for SSE commands
        const int Ktot_d = (Kb<<lgKB) + Kd;
        TMVAssert(Ktot_d == RoundUp<RT>(K));
#ifdef PRINTALGO_MM_BLOCK
        std::cout<<"Mb,Nb,Kb = "<<Mb<<','<<Nb<<','<<Kb<<std::endl;
        std::cout<<"Mc,Nc,Kc,Kd = "<<Mc<<','<<Nc<<','<<Kc<<','<<Kd<<std::endl;
#endif

        const int two1 = m1_z ? 2 : 1;
        const int two2 = m2_z ? 2 : 1;
        const int two3 = m3_z ? 2 : 1;

        const int size1 = two1*MB*KB;
        const int size1y = two1*Mc*KB;
        const int size2 = two2*NB*KB;
        const int size2y = two2*Nc*KB;
        const int size3 = two3*MB*NB;

        AlignedArray<RT> m3_temp(size3);
        RT*const m3p = m3_temp;
        if (N >= M) {
            // Then we loop over the columns of m2 (in blocks).
            // We store a full copy of m1 in block format.
            // Each block column of m2 is copied one at a time into 
            // temporary storage.
            AlignedArray<RT> m1_temp(two1*M*Ktot_d);
            AlignedArray<RT> m2_temp(two2*NB*Ktot_d);
            RT*const m1p0 = m1_temp;
            RT*const m2p0 = m2_temp;
            const int fullsize1 = two1*MB*Ktot_d;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"N >= M\n";
            std::cout<<"m1p0 = "<<m1p0<<" end = "<<(m1p0+size1)<<std::endl;
            std::cout<<"m2p0 = "<<m2p0<<" end = "<<(m2p0+size2)<<std::endl;
            std::cout<<"m3p = "<<m3p<<" end = "<<(m3p+size3)<<std::endl;
#endif
#ifdef TEST_POINTERS
            glob_m1p_end = m1p0+size1;
            glob_m2p_end = m2p0+size2;
            glob_m3p_end = m3p+size3;
#endif

            for(int j=0;j<Nb;++j)
            {
                RT* m1p = m1p0;
                for(int i=0;i<Mb;++i)
                {
                    DoBlockMultMM3(
                        x,m1,m2,m3,
                        i*MB,j*NB,(i+1)*MB,(j+1)*NB, MB,NB,KB,Kb,
                        size1,size2,size3,
                        m1p,m2p0,m3p, j==0,i==0, 
                        copy1,copy2,copy1_k,copy2_k,
                        assign, prod, cleanup, kf,
                        Ka,K,Kc,Kd);
                    m1p += fullsize1;
                }
                if (Mc) {
                    DoBlockMultMM3(
                        x,m1,m2,m3,
                        Ma,j*NB,M,(j+1)*NB, Mc,NB,KB,Kb,
                        size1y,size2,size3,
                        m1p,m2p0,m3p, j==0,Mb==0, 
                        copy1b,copy2,copy1b_k,copy2_k,
                        assignb, prodb, cleanup, kfb,
                        Ka,K,Kc,Kd);
                }
            }
            if (Nc) {
                RT* m1p = m1p0;
                for(int i=0;i<Mb;++i)
                {
                    DoBlockMultMM3(
                        x,m1,m2,m3,
                        i*MB,Na,(i+1)*MB,N, MB,Nc,KB,Kb,
                        size1,size2y,size3,
                        m1p,m2p0,m3p, Nb==0,i==0,
                        copy1,copy2a,copy1_k,copy2a_k,
                        assigna, proda, cleanup, kfa,
                        Ka,K,Kc,Kd);
                    m1p += fullsize1;
                }
                if (Mc) {
                    DoBlockMultMM3(
                        x,m1,m2,m3,
                        Ma,Na,M,N, Mc,Nc,KB,Kb,
                        size1y,size2y,size3,
                        m1p,m2p0,m3p, Nb==0,Mb==0,
                        copy1b,copy2a,copy1b_k,copy2a_k,
                        assignc, prodc, cleanup, kfc,
                        Ka,K,Kc,Kd);
                }
            }
        } else {
            // Then we loop over the rows of m1 and store a full copy of m2.
            AlignedArray<RT> m1_temp(two1*MB*Ktot_d);
            AlignedArray<RT> m2_temp(two2*N*Ktot_d);
            RT*const m1p0 = m1_temp;
            RT*const m2p0 = m2_temp;
            const int fullsize2 = two2*NB*Ktot_d;
#ifdef PRINTALGO_MM_BLOCK
            std::cout<<"N <= M\n";
            std::cout<<"m1p0 = "<<m1p0<<" end = "<<(m1p0+size1)<<std::endl;
            std::cout<<"m2p0 = "<<m2p0<<" end = "<<(m2p0+size2)<<std::endl;
            std::cout<<"m3p = "<<m3p<<" end = "<<(m3p+size3)<<std::endl;
#endif
#ifdef TEST_POINTERS
            glob_m1p_end = m1p0+size1;
            glob_m2p_end = m2p0+size2;
            glob_m3p_end = m3p+size3;
#endif

            for(int i=0;i<Mb;++i)
            {
                RT* m2p = m2p0;
                for(int j=0;j<Nb;++j)
                {
                    DoBlockMultMM3(
                        x,m1,m2,m3,
                        i*MB,j*NB,(i+1)*MB,(j+1)*NB, MB,NB,KB,Kb,
                        size1,size2,size3,
                        m1p0,m2p,m3p, j==0,i==0,
                        copy1,copy2,copy1_k,copy2_k,
                        assign, prod, cleanup, kf,
                        Ka,K,Kc,Kd);
                    m2p += fullsize2;
                }
                if (Nc) {
                    DoBlockMultMM3(
                        x,m1,m2,m3,
                        i*MB,Na,(i+1)*MB,N, MB,Nc,KB,Kb,
                        size1,size2y,size3,
                        m1p0,m2p,m3p, Nb==0,i==0, 
                        copy1,copy2a,copy1_k,copy2a_k,
                        assigna, proda, cleanup, kfa,
                        Ka,K,Kc,Kd);
                }
            }
            if (Mc) {
                RT* m2p = m2p0;
                for(int j=0;j<Nb;++j)
                {
                    DoBlockMultMM3(
                        x,m1,m2,m3,
                        Ma,j*NB,M,(j+1)*NB, Mc,NB,KB,Kb,
                        size1y,size2,size3,
                        m1p0,m2p,m3p, j==0,Mb==0, 
                        copy1b,copy2,copy1b_k,copy2_k,
                        assignb, prodb, cleanup, kfb,
                        Ka,K,Kc,Kd);
                    m2p += fullsize2;
                }
                if (Nc) {
                    DoBlockMultMM3(
                        x,m1,m2,m3,
                        Ma,Na,M,N, Mc,Nc,KB,Kb,
                        size1y,size2y,size3,
                        m1p0,m2p,m3p, Nb==0,Mb==0, 
                        copy1b,copy2a,copy1b_k,copy2a_k,
                        assignc, prodc, cleanup, kfc,
                        Ka,K,Kc,Kd);
                }
            }
        }
    }

    template <int algo, int cs, int rs, int xs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_Block_Helper;

    // algo 1: The normal block algorithm.
    template <int cs, int rs, int xs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_Block_Helper<1,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Traits2<T1,T2>::type PT3;
            typedef typename M3::real_type RT;
            const bool x_z = Traits<T>::iscomplex;
            const bool m1_z = Traits<T1>::iscomplex;
            const bool m2_z = Traits<T2>::iscomplex;
            const bool m3_z = Traits<PT3>::iscomplex;
            const int MB    = 16;
            const int NB    = 16;
#if 0
            const int KB    = 16;
#else
            const int KB    = 
#ifdef __SSE__
                Traits2<RT,float>::sametype ? (
                    ( xs == UNKNOWN || xs >= 128 ) ? 64 :
                    ( xs >= 96 ) ? 32 :
                    ( xs >= 64 ) ? 64 :
                    ( xs >= 48 ) ? 16 :
                    ( xs >= 32 ) ? 32 :
                    16 ) :
#endif
#ifdef __SSE2__
                Traits2<RT,double>::sametype ? (
                    ( xs == UNKNOWN || xs >= 64 ) ? 32 :
                    ( xs >= 48 ) ? 16 :
                    ( xs >= 32 ) ? 32 :
                    16 ) :
#endif
                16;
#endif
            const int lgMB  = IntTraits<MB>::log;
            const int lgNB  = IntTraits<NB>::log;
            const int lgKB  = IntTraits<KB>::log;
            const int csx = cs == UNKNOWN ? UNKNOWN : (cs-((cs>>lgMB)<<lgMB));
            const int rsx = rs == UNKNOWN ? UNKNOWN : (rs-((rs>>lgNB)<<lgNB));
            const int xsx = xs == UNKNOWN ? UNKNOWN : (xs-((xs>>lgKB)<<lgKB));
            typedef typename M1::const_submatrix_type M1sub;
            typedef typename M1sub::const_transpose_type M1sub_t;
            typedef typename M2::const_submatrix_type M2sub;
            typedef typename M3::submatrix_type M3sub;

            // We will determine the value of K once and call a known K
            // algorithm for the K cleanup calls.
            typedef void KF(
                const int M, const int N, const int K,
                const Scaling<1,RT>& x, const RT* A, const RT* B, RT* C);

            typedef void CopyF(const ConstMatrixView<RT>& m1x,
                               RT* m2p, int M, int N, int si, int sj);
            typedef void AssignF(
                const RT* x, 
                RT* m1p, const int M, const int N, const int si1, const int sj1,
                MatrixView<RT> m2x);
            typedef void ProdF(
                const int M, const int N, const int K,
                RT* m1p, const int si1, const int sj1,
                RT* m2p, const int si2, const int sj2,
                RT* m3p, const int si3, const int sj3);
            typedef void CleanupF(
                const int M, const int N, const int K,
                RT* m1p, const int si1, const int sj1,
                RT* m2p, const int si2, const int sj2,
                RT* m3p, const int si3, const int sj3,
                KF* kf);

            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;

            const int Mb = (M>>lgMB); // = M/MB
            const int Nb = (N>>lgNB); // = N/NB
            const int Kb = (K>>lgKB); // = K/KB
            const int Mc = M-(Mb<<lgMB); // = M%MB
            const int Nc = N-(Nb<<lgNB); // = N%NB
            const int Kc = K-(Kb<<lgKB); // = K%KB

            KF* kf = 0;
            KF* kfa = 0;
            KF* kfb = 0;
            KF* kfc = 0;
            if (Kc) {
                kf = get_Kcleanup<xsx,KB,RT>::call(Kc);
                if (Nc)
                    kfa = &call_multmm_M_N_K<1,RT>;
                if (Mc)
                    kfb = &call_multmm_M_N_K<1,RT>;
                if (Mc && Nc)
                    kfc = &call_multmm_M_N_K<1,RT>;
            }

            CopyF* copy1 = &MyCopy<m1_z,KB,MB,M1sub_t>::call;
            CopyF* copy1_k = &MyCopy<m1_z,xsx,MB,M1sub_t>::call;
            CopyF* copy1b = &MyCopy<m1_z,KB,csx,M1sub_t>::call;
            CopyF* copy1b_k = &MyCopy<m1_z,xsx,csx,M1sub_t>::call;
            CopyF* copy2 = &MyCopy<m2_z,KB,NB,M2sub>::call;
            CopyF* copy2_k = &MyCopy<m2_z,xsx,NB,M2sub>::call;
            CopyF* copy2a = &MyCopy<m2_z,KB,rsx,M2sub>::call;
            CopyF* copy2a_k = &MyCopy<m2_z,xsx,rsx,M2sub>::call;

            AssignF* assign = &MyAssign<x_z,m3_z,add,ix,T,MB,NB,M3sub>::call;
            AssignF* assigna = &MyAssign<x_z,m3_z,add,ix,T,MB,rsx,M3sub>::call;
            AssignF* assignb = &MyAssign<x_z,m3_z,add,ix,T,csx,NB,M3sub>::call;
            AssignF* assignc = &MyAssign<x_z,m3_z,add,ix,T,csx,rsx,M3sub>::call;

            T xx = T(x);

            ProdF* prod = &MyProd<m1_z,m2_z,MB,NB,KB,RT>::call;
            ProdF* proda = &MyProd<m1_z,m2_z,MB,rsx,KB,RT>::call;
            ProdF* prodb = &MyProd<m1_z,m2_z,csx,NB,KB,RT>::call;
            ProdF* prodc = &MyProd<m1_z,m2_z,csx,rsx,KB,RT>::call;

            CleanupF* cleanup = &MyCleanup<m1_z,m2_z,RT>::call;

            DoBlockMultMM2(
                (RT*)(&xx),
                m1.realPart().xView(),m2.realPart().xView(),
                m3.realPart().xView(),
                M,N,K, MB,NB,KB, lgMB,lgNB,lgKB, m1_z,m2_z,m3_z,
                copy1,copy1_k,copy1b,copy1b_k,
                copy2,copy2_k,copy2a,copy2a_k,
                assign,assigna,assignb,assignc,
                prod,proda,prodb,prodc,
                cleanup,kf,kfa,kfb,kfc);
        }
    };

    // algo 2: m1, m2, or m3 is nomajor (but might really be rm or cm).
    // Use get_copy and get_assign to figure out which copy or assign
    // function to use.
    template <int cs, int rs, int xs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_Block_Helper<2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Traits2<T1,T2>::type PT3;
            typedef typename M3::real_type RT;
            const bool x_z = Traits<T>::iscomplex;
            const bool m1_z = Traits<T1>::iscomplex;
            const bool m2_z = Traits<T2>::iscomplex;
            const bool m3_z = Traits<PT3>::iscomplex;
            const int MB    = 16;
            const int NB    = 16;
#if 0
            const int KB    = 16;
#else
            const int KB    = 
#ifdef __SSE__
                Traits2<RT,float>::sametype ? 64 :
#endif
#ifdef __SSE2__
                Traits2<RT,double>::sametype ? 32 :
#endif
                16;
#endif
            const int lgMB  = IntTraits<MB>::log;
            const int lgNB  = IntTraits<NB>::log;
            const int lgKB  = IntTraits<KB>::log;

            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_view_type M2v;
            typedef typename M3::view_type M3v;

            typedef void KF(
                const int M, const int N, const int K,
                const Scaling<1,RT>& x, const RT* A, const RT* B, RT* C);

            typedef void CopyF(const ConstMatrixView<RT>& m1x,
                               RT* m2p, int M, int N, int si, int sj);
            typedef void AssignF(
                const RT* x, 
                RT* m1p, const int M, const int N, const int si1, const int sj1,
                MatrixView<RT> m2x);
            typedef void ProdF(
                const int M, const int N, const int K,
                RT* m1p, const int si1, const int sj1,
                RT* m2p, const int si2, const int sj2,
                RT* m3p, const int si3, const int sj3);
            typedef void CleanupF(
                const int M, const int N, const int K,
                RT* m1p, const int si1, const int sj1,
                RT* m2p, const int si2, const int sj2,
                RT* m3p, const int si3, const int sj3,
                KF* kf);

            const int M = m3.colsize();
            const int N = m3.rowsize();
            const int K = m1.rowsize();

            const int XX = UNKNOWN;

            const int Mb = (M>>lgMB); // = M/MB
            const int Nb = (N>>lgNB); // = N/NB
            const int Kb = (K>>lgKB); // = K/KB
            const int Mc = M-(Mb<<lgMB); // = M%MB
            const int Nc = N-(Nb<<lgNB); // = N%NB
            const int Kc = K-(Kb<<lgKB); // = K%KB

            KF* kf = 0;
            KF* kfa = 0;
            KF* kfb = 0;
            KF* kfc = 0;
            if (Kc) {
                kf = get_Kcleanup<XX,KB,RT>::call(Kc);
                if (Nc)
                    kfa = &call_multmm_M_N_K<1,RT>;
                if (Mc)
                    kfb = &call_multmm_M_N_K<1,RT>;
                if (Mc && Nc)
                    kfc = &call_multmm_M_N_K<1,RT>;
            }

            CopyF* copy1 = get_copy<m1_z,KB,MB,M1t>::call(m1.transpose());
            CopyF* copy1_k = get_copy<m1_z,XX,MB,M1t>::call(m1.transpose());
            CopyF* copy1b = get_copy<m1_z,KB,XX,M1t>::call(m1.transpose());
            CopyF* copy1b_k = get_copy<m1_z,XX,XX,M1t>::call(m1.transpose());
            CopyF* copy2 = get_copy<m2_z,KB,NB,M2v>::call(m2);
            CopyF* copy2_k = get_copy<m2_z,XX,NB,M2v>::call(m2);
            CopyF* copy2a = get_copy<m2_z,KB,XX,M2v>::call(m2);
            CopyF* copy2a_k = get_copy<m2_z,XX,XX,M2v>::call(m2);

            AssignF* assign = get_assign<x_z,m3_z,T,MB,NB,M3v>::call(add,x,m3);
            AssignF* assigna = get_assign<x_z,m3_z,T,MB,XX,M3v>::call(add,x,m3);
            AssignF* assignb = get_assign<x_z,m3_z,T,XX,NB,M3v>::call(add,x,m3);
            AssignF* assignc = get_assign<x_z,m3_z,T,XX,XX,M3v>::call(add,x,m3);

            T xx = T(x);

            ProdF* prod = &MyProd<m1_z,m2_z,MB,NB,KB,RT>::call;
            ProdF* proda = &MyProd<m1_z,m2_z,MB,XX,KB,RT>::call;
            ProdF* prodb = &MyProd<m1_z,m2_z,XX,NB,KB,RT>::call;
            ProdF* prodc = &MyProd<m1_z,m2_z,XX,XX,KB,RT>::call;

            CleanupF* cleanup = &MyCleanup<m1_z,m2_z,RT>::call;

            DoBlockMultMM2(
                (RT*)(&xx),
                m1.realPart().xView(),m2.realPart().xView(),
                m3.realPart().xView(),
                M,N,K, MB,NB,KB, lgMB,lgNB,lgKB, m1_z,m2_z,m3_z,
                copy1,copy1_k,copy1b,copy1b_k,
                copy2,copy2_k,copy2a,copy2a_k,
                assign,assigna,assignb,assignc,
                prod,proda,prodb,prodc,
                cleanup,kf,kfa,kfb,kfc);
        }
    };

    // algo -2: Call correct version depending on whether majority is known
    template <int cs, int rs, int xs, bool add,
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_Block_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo =
                ( (M1::_rowmajor || M1::_colmajor) &&
                  (M2::_rowmajor || M2::_colmajor) &&
                  (M3::_rowmajor || M3::_colmajor) ) ? 1 : 2;
            MultMM_Block_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
        }
    };

    // algo 98: Call inst
    template <int cs, int rs, int xs, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_Block_Helper<98,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstMultMM_Block(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int xs, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_Block_Helper<98,cs,rs,xs,true,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstAddMultMM_Block(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo -1: Check for inst
    template <int cs, int rs, int xs, bool add,
              int ix, class T, class M1, class M2, class M3>
    struct MultMM_Block_Helper<-1,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst =
                M1::unknownsizes &&
                M2::unknownsizes &&
                M3::unknownsizes &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                inst ? 98 :
                -2;
            MultMM_Block_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_Block(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        MultMM_Block_Helper<-1,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(
            x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM_Block(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        MultMM_Block_Helper<-2,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(
            x,m1v,m2v,m3v);
    }



} // namespace tmv

#endif 
