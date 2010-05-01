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


#ifndef TMV_MultMM_Winograd_H
#define TMV_MultMM_Winograd_H

#include "TMV_Matrix.h"

// Q4 is the minimum size to use a recursive Winograd algorithm
#ifdef TMV_USE_RECURSIVE_BLOCK
#define TMV_Q4 2048
#else
#define TMV_Q4 1024
#endif

namespace tmv {

    // Defined in TMV_MultMM_Winograd.cpp
    template <class T1, bool CC1, class T2, bool CC2, class T3>
    void InstMultMM_Winograd(
        const T3 x, 
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,CC1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,CC2>& m2,
        MatrixView<T3> m3);
    template <class T1, bool CC1, class T2, bool CC2, class T3>
    void InstAddMultMM_Winograd(
        const T3 x, 
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,CC1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,CC2>& m2,
        MatrixView<T3> m3);

    // Defined below:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_Winograd(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM_Winograd(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);

    // Need this for final pass.
    template <int algo, int cs, int rs, int xs, bool add, 
              int ix, class T, class M1, class M2, class M3> 
    struct MultMM_Helper;

    template <int algo, bool add, 
              int ix, class T, class M1, class M2, class M3> 
    struct MultMM_Winograd_Helper;

    // algo 0: Final pass: call another algorithm
    template <bool add, int ix, class T, class M1, class M2, class M3> 
    struct MultMM_Winograd_Helper<0,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            const int xx = UNKNOWN;
            MultMM_Helper<73,xx,xx,xx,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
        }
    };

    // algo 10: Normal recursion
    template <bool add, int ix, class T, class M1, class M2, class M3> 
    struct MultMM_Winograd_Helper<1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3 m3)
        {
            // This is a recursive algorithm that avoids some of the 
            // matrix multiplication steps by using some clever algebra.
            // This particular implementation of the Winograd algorithm
            // is applicable to non-square matrices and transitions to 
            // a regular algorithm at a particular recursion threshold.
            //
            // In this code, each matrix is divided into 4 parts:
            // m1 == A = ( A0  A1 )
            //           ( A2  A3 )
            // m2 == B = ( B0  B1 )
            //           ( B2  B3 )
            // m3 == C = ( C0  C1 )
            //           ( C2  C3 )
            // 
            // The basic algorithm in pseudo-code, taken directly from the
            // Alberto and Nicolau paper (Algorithm 1) is technically
            // only applicable if M,N,K are all even.  But odd values are
            // allowed if the definition of + for matrices is extended to 
            // apply to different sized matrices, where the extra row or 
            // column is just copied to the destination matrix.
            //
            // Note: There are some typos in both their listing in Algorithm 1
            // and the version found in the appendix, although the 
            // non-psuedo-code listing in the appendix is accurate.
            // The listing below has (I think) corrected all the typos.
            //
            // Notation: (+=) means = if add=false, or += if add=true
            //
            // if (min(m,n,k) < TMV_Q4) { Call direct algorithm }
            // else {
            //   X = A2 + A3
            //   Y = B1 - B0
            //   Z = X * Y
            //   C3 (+=) Z
            //   C1 (+=) Z
            // 
            //   Z = A0 * B0
            //   C0 (+=) Z
            //
            //   C0 += A1 * B2
            //   X = X - A0
            //   Y = B3 - Y
            //   Z += X * Y
            //   C1 += Z
            //   
            //   X = A1 - X
            //   C1 += X * B3
            // 
            //   Y = B2 - Y
            //   C2 (+=) A3 * Y
            //
            //   X = A0 - A2
            //   Y = B3 - B1
            //   Z += X * Y
            //   C3 += Z
            //   C2 += Z
            // }
            //
            // X = A0-A2
            // Y = B3-B1
            // Z = A0*B0 + (A2+A3-A0)*(B0+B3-B1) + (A0-A2)*(B3-B1)
            // If you track through the calculation, the values of C are:
            // C0 (+=) A0*B0 + A1*B2
            // C1 (+=) (A2+A3)*(B1-B0) + A0*B0 + (A2+A3-A0)*(B0+B3-B1)
            //            + (A0+A1-A2-A3)*B3
            // C2 (+=) A3*(B1+B2-B0-B3) + A0*B0 + (A2+A3-A0)*(B0+B3-B1) 
            //            + (A0-A2)*(B3-B1)
            // C3 (+=) (A2+A3)*(B1-B0) + A0*B0 + (A2+A3-A0)*(B0+B3-B1) 
            //            + (A0-A2)*(B3-B1)
            //
            // It is straightforward to check that these equations simplify
            // to the right answers:
            // C0 (+=) A0*B0 + A1*B2
            // C1 (+=) (A2+A3)*(B1-B0) + A0*B0 + (A2+A3-A0)*(B0+B3-B1)
            //            + (A0+A1-A2-A3)*B3
            //       = A2*B1 + A3*B1 - A2*B0 - A3*B0 + A0*B0 + A2*B0 + A2*B3 
            //         - A2*B1 + A3*B0 + A3*B3 - A3*B1 - A0*B0 - A0*B3 + A0*B1
            //         + A0*B3 + A1*B3 - A2*B3 - A3*B3
            //       = A0*B1 + A1*B3
            // C2 (+=) A3*(B1+B2-B0-B3) + A0*B0 + (A2+A3-A0)*(B0+B3-B1) 
            //            + (A0-A2)*(B3-B1)
            //       = A3*B1 + A3*B2 - A3*B0 - A3*B3 + A0*B0 + A2*B0 + A2*B3 
            //         - A2*B1 + A3*B0 + A3*B3 - A3*B1 - A0*B0 - A0*B3 + A0*B1
            //         + A0*B3 - A0*B1 - A2*B3 + A2*B1
            //       = A2*B0 + A3*B2
            // C3 (+=) (A2+A3)*(B1-B0) + A0*B0 + (A2+A3-A0)*(B0+B3-B1) 
            //            + (A0-A2)*(B3-B1)
            //       = A2*B1 + A3*B1 - A2*B0 - A3*B0 + A0*B0 + A2*B0 + A2*B3 
            //         - A2*B1 + A3*B0 + A3*B3 - A3*B1 - A0*B0 - A0*B3 + A0*B1
            //         + A0*B3 - A0*B1 - A2*B3 + A2*B1
            //       = A2*B1 + A3*B3
            //

            TMVAssert(m1.colsize() == m3.colsize());
            TMVAssert(m1.rowsize() == m2.colsize());
            TMVAssert(m2.rowsize() == m3.rowsize());

            const int M = m3.colsize();
            const int N = m3.rowsize();
            const int K = m1.rowsize();
            TMVAssert(int(m1.colsize()) == M);
            TMVAssert(int(m2.rowsize()) == N);
            TMVAssert(int(m2.colsize()) == K);

            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Traits<T3>::real_type RT;
            Scaling<1,RT> one;
            Scaling<-1,RT> mone;

            if (M < TMV_Q4 && N < TMV_Q4 && K < TMV_Q4) {
                return MultMM_Winograd_Helper<0,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
            }

            const int Mb = M>>1;
            const int Ma = M-Mb;
            const int Nb = N>>1;
            const int Na = N-Nb;
            const int Kb = K>>1;
            const int Ka = K-Kb;

#ifdef PRINTALGO_MM_WIN
            std::cout<<"Winograd recursive branch: \n";
            std::cout<<"M = "<<Ma<<" + "<<Mb<<std::endl;
            std::cout<<"N = "<<Na<<" + "<<Nb<<std::endl;
            std::cout<<"K = "<<Ka<<" + "<<Kb<<std::endl;
#endif

            // X,Y,Z are temporaries.
            // We use the subscripts 1,2,3 to match the size of the 
            // corrsponding A, B or C submatrix.  0 matches the full size.
            Matrix<T1,RowMajor> X(Ma,Ka,T1(0)); 
            Matrix<T2,ColMajor> Y(Ka,Na,T2(0));
            Matrix<T3,ColMajor> Z(Ma,Na,T3(0));

            typedef MatrixView<T1,UNKNOWN,1> M1r;
            typedef MatrixView<T2,1> M2c;
            typedef MatrixView<T3,1> M3c;

            M1r X0 = X.view();
            M1r X1 = X.colRange(0,Kb).xView();
            M1r X2 = X.rowRange(0,Mb);
            M1r X3 = X.subMatrix(0,Mb,0,Kb);

            M2c Y0 = Y.view();
            M2c Y1 = Y.colRange(0,Nb);
            M2c Y2 = Y.rowRange(0,Kb);
            M2c Y3 = Y.subMatrix(0,Kb,0,Nb);

            M3c Z0 = Z.view();
            M3c Z1 = Z.colRange(0,Nb);
            M3c Z2 = Z.rowRange(0,Mb);
            M3c Z3 = Z.subMatrix(0,Mb,0,Nb);

            M1 A0 = m1.subMatrix(0,Ma,0,Ka);
            M1 A1 = m1.subMatrix(0,Ma,Ka,K);
            M1 A2 = m1.subMatrix(Ma,M,0,Ka);
            M1 A3 = m1.subMatrix(Ma,M,Ka,K);
            M2 B0 = m2.subMatrix(0,Ka,0,Na);
            M2 B1 = m2.subMatrix(0,Ka,Na,N);
            M2 B2 = m2.subMatrix(Ka,K,0,Na);
            M2 B3 = m2.subMatrix(Ka,K,Na,N);
            M3 C0 = m3.subMatrix(0,Ma,0,Na);
            M3 C1 = m3.subMatrix(0,Ma,Na,N);
            M3 C2 = m3.subMatrix(Ma,M,0,Na);
            M3 C3 = m3.subMatrix(Ma,M,Na,N);

            // X = A2 + A3
            NoAliasCopy(A2,X2);
            NoAliasMultXM<true>(one,A3,X3);

            // Y = B1 - B0
            NoAliasCopy(B1,Y1);
            NoAliasMultXM<true>(mone,B0,Y);

            // Z = X * Y
            InlineMultMM_Winograd<false>(one,X2,Y0,Z2);

            // C3 (+=) Z
            NoAliasMultXM<add>(x,X3,C3);

            // C1 (+=) Z
            NoAliasMultXM<add>(x,Z1,C1);

            // Z = A0 * B0
            InlineMultMM_Winograd<false>(one,A0,B0,Z0);

            // C0 (+=) Z
            NoAliasMultXM<add>(x,Z,C0);

            // C0 += A1 * B2
            InlineMultMM_Winograd<true>(x,A1,B2,C0);

            // X = X - A0
            NoAliasMultXM<true>(mone,A0,X);

            // Y = B3 - Y
            Scale(mone,Y);
            NoAliasMultXM<true>(one,B3,Y3);

            // Z += X * Y
            InlineMultMM_Winograd<true>(one,X0,Y0,Z0);

            // C1 += Z
            NoAliasMultXM<true>(x,Z1,C1);

            // X = A1 - X
            Scale(mone,X);
            NoAliasMultXM<true>(one,A1,X1);

            // C1 += X * B3
            InlineMultMM_Winograd<true>(x,X1,B3,C1);

            // Y = B2 - Y
            Scale(mone,Y);
            NoAliasMultXM<true>(one,B2,Y2);

            // C2 (+=) A3 * Y
            InlineMultMM_Winograd<add>(x,A3,Y2,C2);

            // X = A0 - A2
            NoAliasCopy(A0,X);
            NoAliasMultXM<true>(mone,A2,X2);

            // Y = B3 - B1
            NoAliasMultXM<false>(mone,B1,Y1);
            NoAliasMultXM<true>(one,B3,Y3);

            // Z += X * Y
            InlineMultMM_Winograd<true>(one,X0,Y1,Z1);

            // C3 += Z
            NoAliasMultXM<true>(x,Z3,C3);

            // C2 += Z
            NoAliasMultXM<true>(x,Z2,C2);
        }
    };

    // algo 98: Call inst
    template <int ix, class T, class M1, class M2, class M3> 
    struct MultMM_Winograd_Helper<98,false,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstMultMM_Winograd(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int ix, class T, class M1, class M2, class M3> 
    struct MultMM_Winograd_Helper<98,true,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstAddMultMM_Winograd(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo -1: Check for inst
    template <bool add, int ix, class T, class M1, class M2, class M3> 
    struct MultMM_Winograd_Helper<-1,add,ix,T,M1,M2,M3>
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
                1;
            MultMM_Winograd_Helper<algo,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_Winograd(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        typedef typename M1::const_xview_type M1v;
        typedef typename M2::const_xview_type M2v;
        typedef typename M3::xview_type M3v;
        M1v m1v = m1.xView();
        M2v m2v = m2.xView();
        M3v m3v = m3.xView();
        MultMM_Winograd_Helper<-1,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM_Winograd(
        const Scaling<ix,T>& x,
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        typedef typename M1::const_xview_type M1v;
        typedef typename M2::const_xview_type M2v;
        typedef typename M3::xview_type M3v;
        M1v m1v = m1.xView();
        M2v m2v = m2.xView();
        M3v m3v = m3.xView();
        MultMM_Winograd_Helper<1,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }





} // namespace tmv

#endif 
