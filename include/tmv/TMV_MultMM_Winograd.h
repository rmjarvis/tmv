

#ifndef TMV_MultMM_Winograd_H
#define TMV_MultMM_Winograd_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_MultXM_Funcs.h"

#ifndef TMV_MM_MIN_WINOGRAD
#ifdef TMV_USE_RECURSIVE_BLOCK
#define TMV_MM_MIN_WINOGRAD 2048
#else
#define TMV_MM_MIN_WINOGRAD 1024
#endif
#endif

namespace tmv {

    // Defined in TMV_MultMM_Winograd.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM_Winograd(
        const T3 x, 
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM_Winograd(
        const T3 x, 
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_Winograd(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);


    // Need this for final pass.
    template <int algo, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper;

    template <int algo, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Winograd_Helper;

    // algo 11: Normal recursion
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Winograd_Helper<11,add,ix,T,M1,M2,M3>
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
            // if (min(m,n,k) < TMV_MM_MIN_WINOGRAD) {
            //     Call direct algorithm 
            // } else {
            //     X = A2 + A3
            //     Y = B1 - B0
            //     Z = X * Y
            //     C3 (+=) Z
            //     C1 (+=) Z
            // 
            //     Z = A0 * B0
            //     C0 (+=) Z
            //
            //     C0 += A1 * B2
            //     X = X - A0
            //     Y = B3 - Y
            //     Z += X * Y
            //     C1 += Z
            //     
            //     X = A1 - X
            //     C1 += X * B3
            // 
            //     Y = B2 - Y
            //     C2 (+=) A3 * Y
            //
            //     X = A0 - A2
            //     Y = B3 - B1
            //     Z += X * Y
            //     C3 += Z
            //     C2 += Z
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

            const ptrdiff_t M = m3.colsize();
            const ptrdiff_t N = m3.rowsize();
            const ptrdiff_t K = m1.rowsize();
            TMVAssert(m1.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            TMVAssert(m2.colsize() == K);

            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Traits<T3>::real_type RT;

            if (M < TMV_MM_MIN_WINOGRAD && N < TMV_MM_MIN_WINOGRAD && 
                K < TMV_MM_MIN_WINOGRAD) {
                return MultMM_Winograd_Helper<12,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
            }

            const ptrdiff_t Mb = M>>1;
            const ptrdiff_t Ma = M-Mb;
            const ptrdiff_t Nb = N>>1;
            const ptrdiff_t Na = N-Nb;
            const ptrdiff_t Kb = K>>1;
            const ptrdiff_t Ka = K-Kb;

#ifdef PRINTALGO_MM_WIN
            std::cout<<"Winograd recursive branch: \n";
            std::cout<<"M = "<<Ma<<" + "<<Mb<<std::endl;
            std::cout<<"N = "<<Na<<" + "<<Nb<<std::endl;
            std::cout<<"K = "<<Ka<<" + "<<Kb<<std::endl;
#endif

            // X,Y,Z are temporaries.
            // We use the subscripts 1,2,3 to match the size of the 
            // corrsponding A, B or C submatrix.  0 matches the full size.
            Matrix<T1,RowMajor | NoDivider | NoAlias> X(Ma,Ka,T1(0)); 
            Matrix<T2,ColMajor | NoDivider | NoAlias> Y(Ka,Na,T2(0));
            Matrix<T3,ColMajor | NoDivider | NoAlias> Z(Ma,Na,T3(0));

            // Only use the versions that maintain knowledge of the
            // majority of X,Y,Z submatrices if the input m1,m2, and m3
            // all have such knowledge.
            const bool majority_known = 
                (M1::_colmajor || M1::_rowmajor) &&
                (M2::_colmajor || M2::_rowmajor) &&
                (M3::_colmajor || M3::_rowmajor);
            const int cmA = (majority_known ? ColMajor : NonMajor) | NoAlias;
            const int rmA = (majority_known ? RowMajor : NonMajor) | NoAlias;
            const int maybe_one = majority_known ? 1 : 0;
            const int maybe_mone = majority_known ? -1 : 0;
            typedef MatrixView<T1,rmA> M1r;
            typedef MatrixView<T2,cmA> M2c;
            typedef MatrixView<T3,cmA> M3c;
            Scaling<maybe_one,RT> one(RT(1));
            Scaling<maybe_mone,RT> mone(RT(-1));

            M1r X0 = X.view();
            M1r X1 = X.colRange(0,Kb);
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
            Copy(A2,X2);
            MultXM<true>(one,A3,X3);

            // Y = B1 - B0
            Copy(B1,Y1);
            MultXM<true>(mone,B0,Y0);

            // Z = X * Y
            // Z is already 0, so can use add=true.
            MultMM_Winograd<true>(one,X2,Y0,Z2);

            // C3 (+=) Z
            MultXM<add>(x,Z3,C3);

            // C1 (+=) Z
            MultXM<add>(x,Z1,C1);

            // Z = A0 * B0
            Z.setZero();
            MultMM_Winograd<true>(one,A0,B0,Z0);

            // C0 (+=) Z
            MultXM<add>(x,Z0,C0);

            // C0 += A1 * B2
            MultMM_Winograd<true>(x,A1,B2,C0);

            // X = X - A0
            MultXM<true>(mone,A0,X0);

            // Y = B3 - Y
            Scale(mone,Y0);
            MultXM<true>(one,B3,Y3);

            // Z += X * Y
            MultMM_Winograd<true>(one,X0,Y0,Z0);

            // C1 += Z
            MultXM<true>(x,Z1,C1);

            // X = A1 - X
            Scale(mone,X0);
            MultXM<true>(one,A1,X1);

            // C1 += X * B3
            MultMM_Winograd<true>(x,X1,B3,C1);

            // Y = B2 - Y
            Scale(mone,Y0);
            MultXM<true>(one,B2,Y2);

            // C2 (+=) A3 * Y
            MultMM_Winograd<add>(x,A3,Y2,C2);

            // X = A0 - A2
            Copy(A0,X0);
            MultXM<true>(mone,A2,X2);

            // Y = B3 - B1
            MultXM<false>(mone,B1,Y1);
            MultXM<true>(one,B3,Y3);

            // Z += X * Y
            MultMM_Winograd<true>(one,X0,Y1,Z1);

            // C3 += Z
            MultXM<true>(x,Z3,C3);

            // C2 += Z
            MultXM<true>(x,Z2,C2);
        }
    };

    // algo 12: Final pass: call another algorithm
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Winograd_Helper<12,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t xx = Unknown;
            MultMM_Helper<73,xx,xx,xx,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
        }
    };

    // algo 90: Call inst
    template <int ix, class T, class M1, class M2, class M3>
    struct MultMM_Winograd_Helper<90,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMM_Winograd(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int ix, class T, class M1, class M2, class M3>
    struct MultMM_Winograd_Helper<90,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMM_Winograd(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo -3: Only one algorithm here, so do it.
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Winograd_Helper<-3,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            MultMM_Winograd_Helper<11,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -2: Check for inst
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Winograd_Helper<-2,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst =
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                inst ? 90 :
                -3;
            MultMM_Winograd_Helper<algo,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Winograd_Helper<-1,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        { MultMM_Winograd_Helper<-2,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3); }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_Winograd(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        typedef typename M1::const_xview_type M1v;
        typedef typename M2::const_xview_type M2v;
        typedef typename M3::xview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.xView();
        TMV_MAYBE_CREF(M1,M2v) m2v = m2.xView();
        TMV_MAYBE_REF(M1,M3v) m3v = m3.xView();
        MultMM_Winograd_Helper<-2,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM_Winograd(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        typedef typename M1::const_xview_type M1v;
        typedef typename M2::const_xview_type M2v;
        typedef typename M3::xview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.xView();
        TMV_MAYBE_CREF(M1,M2v) m2v = m2.xView();
        TMV_MAYBE_REF(M1,M3v) m3v = m3.xView();
        MultMM_Winograd_Helper<-3,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }





} // namespace tmv

#endif 
