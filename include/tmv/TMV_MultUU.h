

#ifndef TMV_MultUU_H
#define TMV_MultUU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_MultXV.h"
#include "TMV_MultUV.h"
#include "TMV_Rank1VVM.h"
#include "TMV_MultUM.h"
#include "TMV_MultXU.h"

#ifdef PRINTALGO_UU
#include <iostream>
#endif

// UNROLL is the maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_UU_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_UU_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_UU_UNROLL 9
#else
#define TMV_UU_UNROLL 0
#endif

// Inline the MV (MultUV, Rank1VVM) calls.
#if TMV_OPT >= 1
#define TMV_UU_INLINE_MV
#endif

// The size to stop recursing.
#if TMV_OPT >= 3
#define TMV_UU_RECURSE 8
#else
#define TMV_UU_RECURSE 1
#endif

namespace tmv {

    // Defined in TMV_MultUU.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        const ConstUpperTriMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        const ConstUpperTriMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m2,
        const ConstLowerTriMatrixView<T2,C2>& m1, LowerTriMatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        const ConstUpperTriMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        const ConstUpperTriMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m2,
        const ConstLowerTriMatrixView<T2,C2>& m1, LowerTriMatrixView<T3> m3);

    template <int algo, ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper;

    // algo 0: Trivial, nothing to do.
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<0,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {}
    };

    // algo 1: s == 1, so reduces to scalar product
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<1,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 1: N,s,x = "<<1<<','<<1<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(!M3::_unit);
            Maybe<add>::add(m3.ref(0,0), x*m1.cref(0,0)*m2.cref(0,0));
        }
    };

    // algo 11: UpperTri loop over n
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<11,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 11: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            const bool u2 = M2::_unit;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_subtrimatrix_type M1s;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int ix2 = u2 ? ix : 0;
            const ptrdiff_t xx = Unknown;
#ifdef TMV_UU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(ptrdiff_t j=N-1;j>=0;--j) {
                // m3.col(j,0,j+1) = m1.subTriMatrix(0,j+1) * m2.col(j,0,j+1)
                // ==>
                // m3.col(j,0,j) = m1.col(j,0,j) * m2(j,j) +
                //                 m1.subTriMatrix(0,j) * m2.col(j,0,j)
                // m3(j,j) = m1(j,j) * m2(j,j)
                const Scaling<ix2,PT2> xd(Maybe<!u2>::prod(m2.cref(j,j),x));
                M1s m1s = m1.cSubTriMatrix(0,j);
                M1c m1c = m1.get_col(j,0,j);
                M2c m2c = m2.get_col(j,0,j);
                M3c m3c = m3.get_col(j,0,j);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd,m1c,m3c);
                MultUV_Helper<algo2,xx,true,ix,T,M1s,M2c,M3c>::call(
                    x,m1s,m2c,m3c);
                Maybe2<!M3::_unit,add>::add(
                    m3.ref(j,j), Maybe<!u1>::prod(m1.cref(j,j),xd));
            }
        }
    };

    // algo 12: UpperTri loop over m
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<12,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 12: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            const bool u2 = M2::_unit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M2::const_subtrimatrix_type M2s;
            typedef typename M2s::const_transpose_type M2t;
            typedef typename M3::row_sub_type M3r;
            const int ix1 = u1 ? ix : 0;
            const ptrdiff_t xx = Unknown;
#ifdef TMV_UU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(ptrdiff_t i=0;i<N;++i) {
                // m3.row(i,i,N) = m1.row(i,i,N) * m2.subTriMatrix(i,N)
                // ==>
                // m3.row(i,i+1,N) = m1.row(i,i+1,N) * m2.subTriMatrix(i+1,N) +
                //                   m1(i,i) * m2.row(i,i+1,N)
                // m3(i,i) = m1(i,i) * m2(i,i)
                const Scaling<ix1,PT1> xd(Maybe<!u1>::prod(m1.cref(i,i),x));
                M1r m1r = m1.get_row(i,i+1,N);
                M2t m2t = m2.cSubTriMatrix(i+1,N).transpose();
                M2r m2r = m2.get_row(i,i+1,N);
                M3r m3r = m3.get_row(i,i+1,N);
                MultUV_Helper<algo2,xx,add,ix,T,M2t,M1r,M3r>::call(
                    x,m2t,m1r,m3r);
                MultXV_Helper<-4,xx,true,ix1,PT1,M2r,M3r>::call(xd,m2r,m3r);
                Maybe2<!M3::_unit,add>::add(
                    m3.ref(i,i), Maybe<!u2>::prod(m2.cref(i,i),xd));
            }
        }
    };

    // algo 13: UpperTri loop over k
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<13,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 13: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            const bool u2 = M2::_unit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::row_sub_type M3r;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3::submatrix_type M3s;
            const int ix1 = u1 ? ix : 0;
            const int ix2 = u2 ? ix : 0;
            const ptrdiff_t xx = Unknown;
#ifdef TMV_UU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(ptrdiff_t k=N-1;k>=0;--k) {
                // m3.subMatrix(0,k+1,k,N) = m1.col(k,0,k+1) ^ m2.row(k,k,N)
                // ==>
                // m3.subMatrix(0,k,k+1,N) = m1.col(k,0,k) ^ m2.row(k,k+1,N)
                // m3.col(k,0,k) = m1.col(k,0,k) * m2(k,k)
                // m3.row(k,k+1,N) = m1(k,k) * m2.row(k,k+1,N)
                // m3(k,k) = m1(k,k) * m2(k,k)
                const Scaling<ix1,PT1> xd1(Maybe<!u1>::prod(m1.cref(k,k),x));
                const Scaling<ix2,PT2> xd2(Maybe<!u2>::prod(m2.cref(k,k),x));
                M1c m1c = m1.get_col(k,0,k);
                M2r m2r = m2.get_row(k,k+1,N);
                M3c m3c = m3.get_col(k,0,k);
                M3r m3r = m3.get_row(k,k+1,N);
                M3s m3s = m3.cSubMatrix(0,k,k+1,N);
                Rank1VVM_Helper<algo2,xx,xx,true,ix,T,M1c,M2r,M3s>::call(
                    x,m1c,m2r,m3s);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd2,m1c,m3c);
                MultXV_Helper<-4,xx,true,ix1,PT1,M2r,M3r>::call(xd1,m2r,m3r);
                Maybe2<!M3::_unit,add>::add(
                    m3.ref(k,k), Maybe<!u1>::prod(m1.cref(k,k),xd2));
            }
        }
    };

    // algo 16: UpperTri: Unroll small case
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<16,s,add,ix,T,M1,M2,M3>
    {
        template <ptrdiff_t I, ptrdiff_t N>
        struct Unroller
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                // [ C00 C01 ] = [ A00 A01 ] [ B00 B01 ]
                // [  0  C11 ]   [  0  A11 ] [  0  B11 ]

                const ptrdiff_t Nx = N/2;
                const ptrdiff_t Ny = N-Nx;
                const ptrdiff_t I1 = I+Nx;
                const ptrdiff_t I2 = I+N;
                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1st A00 = m1.cSubTriMatrix(I,I1);
                M1smt A01t = m1.cSubMatrix(I,I1,I1,I2).transpose();
                M2sm B01 = m2.cSubMatrix(I,I1,I1,I2);
                M2stt B11t = m2.cSubTriMatrix(I1,I2).transpose();
                M3sm C01 = m3.cSubMatrix(I,I1,I1,I2);
                M3smt C01t = C01.transpose();

                // C01 (+)= A01 B11
                MultUM_Helper<-4,Ny,Nx,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B11t,A01t,C01t);

                // C01 += A00 B01
                MultUM_Helper<-4,Nx,Ny,true,ix,T,M1st,M2sm,M3sm>::call(
                    x,A00,B01,C01);

                // C00 (+)= A00 B00
                Unroller<I,Nx>::unroll(x,m1,m2,m3);

                // C11 (+)= A11 B11
                Unroller<I1,Ny>::unroll(x,m1,m2,m3);
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                const bool u1 = M1::_unit;
                const bool u2 = M2::_unit;
                Maybe2<!M3::_unit,add>::add(
                    m3.ref(I,I), Maybe<!u1>::prod(
                            m1.cref(I,I),Maybe<!u2>::prod(m2.cref(I,I),x)));
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& , const M1& , const M2& , M3& ) {}
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UU
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
            std::cout<<"UU algo 16: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,m2,m3); 
        }
    };
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<16,Unknown,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = m3.size();
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 16: N,s,x = "<<N<<','<<Unknown<<
                ','<<T(x)<<std::endl;
#endif
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const int algo2 = rxr ? 12 : crx ? 13 : xcc ? 11 : 13;
            const int algo1 = M3::_unit ? 0 : 1;
            switch (N) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   MultUU_Helper<algo1,1,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 2 :
                   MultUU_Helper<16,2,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 3 :
                   MultUU_Helper<16,3,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              default :
                   MultUU_Helper<algo2,Unknown,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
            }
        }
    };

    // algo 17: Split the UpperTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<17,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 17: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

            const ptrdiff_t sp1 = IntTraits<s>::Sp1;
            const ptrdiff_t sp2 = IntTraits<sp1>::Sp1;
            // nops = 1/6 n(n+1)(n+2)
            const ptrdiff_t nops = 
                IntTraits2<IntTraits2<s,sp1>::safeprod,sp2>::safeprod / 6;
            const bool unroll = 
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_UU_UNROLL;
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? ( M3::_unit ? 0 : 1 ) :
                unroll ? 16 : 
                // For known s, always recurse down to unroll size
                s != Unknown ? 0 :
                (TMV_UU_RECURSE == 1) ? ( M3::_unit ? 0 : 1) :
                rxr ? 12 : crx ? 13 : xcc ? 11 : 13;
            const int algo3 =  // The algorithm for N > UU_RECURSE
                unroll || s == 1 ? 0 : 17;
            const int algo4 =  // The algorithm for MultUM
                unroll || s == 1 ? 0 : -2;
#ifdef PRINTALGO_UU
            std::cout<<"algo2,3,4 = "<<algo2<<"  "<<algo3<<"  "<<algo4<<std::endl;
#endif

            if (s==Unknown ? (N > TMV_UU_RECURSE) : (s > 1 && !unroll)) {
                // [ C00 C01 ] = [ A00 A01 ] [ B00 B01 ]
                // [  0  C11 ]   [  0  A11 ] [  0  B11 ]

                const ptrdiff_t Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round N/2 up to a multiple of 16.)
                const ptrdiff_t sx = IntTraits<s>::half_roundup;
                const ptrdiff_t sy = IntTraits2<s,sx>::diff;

                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::subtrimatrix_type M3st;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1st A00 = m1.cSubTriMatrix(0,Nx);
                M1smt A01t = m1.cSubMatrix(0,Nx,Nx,N).transpose();
                M1st A11 = m1.cSubTriMatrix(Nx,N);
                M2st B00 = m2.cSubTriMatrix(0,Nx);
                M2sm B01 = m2.cSubMatrix(0,Nx,Nx,N);
                M2st B11 = m2.cSubTriMatrix(Nx,N);
                M2stt B11t = B11.transpose();
                M3st C00 = m3.cSubTriMatrix(0,Nx);
                M3sm C01 = m3.cSubMatrix(0,Nx,Nx,N);
                M3smt C01t = C01.transpose();
                M3st C11 = m3.cSubTriMatrix(Nx,N);

                // C01 (+)= A01 B11
                MultUM_Helper<algo4,sy,sx,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B11t,A01t,C01t);

                // C01 += A00 B01
                MultUM_Helper<algo4,sx,sy,true,ix,T,M1st,M2sm,M3sm>::call(
                    x,A00,B01,C01);

                // C00 (+)= x A00 B00
                MultUU_Helper<algo3,sx,add,ix,T,M1st,M2st,M3st>::call(
                    x,A00,B00,C00);

                // C11 (+)= x A11 B11
                MultUU_Helper<algo3,sy,add,ix,T,M1st,M2st,M3st>::call(
                    x,A11,B11,C11);
            } else {
                MultUU_Helper<algo2,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 21: LowerTri loop over n
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<21,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 21: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            const bool u2 = M2::_unit;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_subtrimatrix_type M1s;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int ix2 = u2 ? ix : 0;
            const ptrdiff_t xx = Unknown;
#ifdef TMV_UU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(ptrdiff_t j=0;j<N;++j) {
                // m3.col(j,j,N) = m1.subTriMatrix(j,N) * m2.col(j,j,N)
                // ==>
                // m3.col(j,j+1,N) = m1.col(j,j+1,N) * m2(j,j) +
                //                   m1.subTriMatrix(j+1,N) * m2.col(j,j+1,N)
                // m3(j,j) = m1(j,j) * m2(j,j)
                const Scaling<ix2,PT2> xd(Maybe<!u2>::prod(m2.cref(j,j),x));
                M1s m1s = m1.cSubTriMatrix(j+1,N);
                M1c m1c = m1.get_col(j,j+1,N);
                M2c m2c = m2.get_col(j,j+1,N);
                M3c m3c = m3.get_col(j,j+1,N);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd,m1c,m3c);
                MultUV_Helper<algo2,xx,true,ix,T,M1s,M2c,M3c>::call(
                    x,m1s,m2c,m3c);
                Maybe2<!M3::_unit,add>::add(
                    m3.ref(j,j), Maybe<!u1>::prod(m1.cref(j,j),xd));
            }
        }
    };

    // algo 22: LowerTri loop over m
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<22,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 22: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            const bool u2 = M2::_unit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M2::const_subtrimatrix_type M2s;
            typedef typename M2s::const_transpose_type M2t;
            typedef typename M3::row_sub_type M3r;
            const int ix1 = u1 ? ix : 0;
            const ptrdiff_t xx = Unknown;
#ifdef TMV_UU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(ptrdiff_t i=N;i--;) {
                // m3.row(i,0,i+1) = m1.row(i,0,i+1) * m2.subTriMatrix(0,i+1)
                // ==>
                // m3.row(i,0,i) = m1.row(i,0,i) * m2.subTriMatrix(0,i) +
                //                 m1(i,i) * m2.row(i,0,i)
                // m3(i,i) = m1(i,i) * m2(i,i)
                const Scaling<ix1,PT1> xd(Maybe<!u1>::prod(m1.cref(i,i),x));
                M1r m1r = m1.get_row(i,0,i);
                M2t m2t = m2.cSubTriMatrix(0,i).transpose();
                M2r m2r = m2.get_row(i,0,i);
                M3r m3r = m3.get_row(i,0,i);
                MultUV_Helper<algo2,xx,add,ix,T,M2t,M1r,M3r>::call(
                    x,m2t,m1r,m3r);
                MultXV_Helper<-4,xx,true,ix1,PT1,M2r,M3r>::call(xd,m2r,m3r);
                Maybe2<!M3::_unit,add>::add(
                    m3.ref(i,i), Maybe<!u2>::prod(m2.cref(i,i),xd));
            }
        }
    };

    // algo 23: LowerTri loop over k
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<23,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 23: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            const bool u2 = M2::_unit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::row_sub_type M3r;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3::submatrix_type M3s;
            const int ix1 = u1 ? ix : 0;
            const int ix2 = u2 ? ix : 0;
            const ptrdiff_t xx = Unknown;
#ifdef TMV_UU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(ptrdiff_t k=0;k<N;++k) {
                // m3.subMatrix(k,N,0,k+1) = m1.col(k,k,N) ^ m2.row(k,0,k+1)
                // ==>
                // m3.subMatrix(k+1,N,0,k) = m1.col(k,k+1,N) ^ m2.row(k,0,k)
                // m3.col(k,k+1,N) = m1.col(k,k+1,N) * m2(k,k)
                // m3.row(k,0,k) = m1(k,k) * m2.row(k,0,k)
                // m3(k,k) = m1(k,k) * m2(k,k)
                const Scaling<ix1,PT1> xd1(Maybe<!u1>::prod(m1.cref(k,k),x));
                const Scaling<ix2,PT2> xd2(Maybe<!u2>::prod(m2.cref(k,k),x));
                M1c m1c = m1.get_col(k,k+1,N);
                M2r m2r = m2.get_row(k,0,k);
                M3c m3c = m3.get_col(k,k+1,N);
                M3r m3r = m3.get_row(k,0,k);
                M3s m3s = m3.cSubMatrix(k+1,N,0,k);
                Rank1VVM_Helper<algo2,xx,xx,true,ix,T,M1c,M2r,M3s>::call(
                    x,m1c,m2r,m3s);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd2,m1c,m3c);
                MultXV_Helper<-4,xx,true,ix1,PT1,M2r,M3r>::call(xd1,m2r,m3r);
                Maybe2<!M3::_unit,add>::add(
                    m3.ref(k,k), Maybe<!u1>::prod(m1.cref(k,k),xd2));
            }
        }
    };

    // algo 26: LowerTri: Unroll small case
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<26,s,add,ix,T,M1,M2,M3>
    {
        template <ptrdiff_t I, ptrdiff_t N>
        struct Unroller
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                // [ C00  0  ] = [ A00  0  ] [ B00  0  ]
                // [ C10 C11 ]   [ A10 A11 ] [ B10 B11 ]

                const ptrdiff_t Nx = N/2;
                const ptrdiff_t Ny = N-Nx;
                const ptrdiff_t I1 = I+Nx;
                const ptrdiff_t I2 = I+N;
                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1smt A10t = m1.cSubMatrix(I1,I2,I,I1).transpose();
                M1st A11 = m1.cSubTriMatrix(I1,I2);
                M2sm B10 = m2.cSubMatrix(I1,I2,I,I1);
                M2stt B00t = m2.cSubTriMatrix(I,I1).transpose();
                M3sm C10 = m3.cSubMatrix(I1,I2,I,I1);
                M3smt C10t = C10.transpose();

                // C10 (+)= A10 B00
                MultUM_Helper<-4,Nx,Ny,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B00t,A10t,C10t);

                // C10 += A11 B10
                MultUM_Helper<-4,Ny,Nx,true,ix,T,M1st,M2sm,M3sm>::call(
                    x,A11,B10,C10);

                // C00 (+)= A00 B00
                Unroller<I,Nx>::unroll(x,m1,m2,m3);

                // C11 (+)= A11 B11
                Unroller<I1,Ny>::unroll(x,m1,m2,m3);
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                const bool u1 = M1::_unit;
                const bool u2 = M2::_unit;
                Maybe2<!M3::_unit,add>::add(
                    m3.ref(I,I), Maybe<!u1>::prod(
                            m1.cref(I,I),Maybe<!u2>::prod(m2.cref(I,I),x)));
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& , const M1& , const M2& , M3& ) {}
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 26: N,s,x = "<<2<<','<<2<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,m2,m3); 
        }
    };
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<26,Unknown,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = m3.size();
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 26: N,s,x = "<<N<<','<<Unknown<<
                ','<<T(x)<<std::endl;
#endif
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const int algo2 = rxr ? 22 : crx ? 23 : xcc ? 21 : 23;
            const int algo1 = M3::_unit ? 0 : 1;
            switch (N) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   MultUU_Helper<algo1,1,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 2 :
                   MultUU_Helper<26,2,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 3 :
                   MultUU_Helper<26,3,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              default :
                   MultUU_Helper<algo2,Unknown,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
            }
        }
    };

    // algo 27: Split the LowerTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<27,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 27: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

            const ptrdiff_t sp1 = IntTraits<s>::Sp1;
            const ptrdiff_t sp2 = IntTraits<sp1>::Sp1;
            // nops = 1/6 n(n+1)(n+2)
            const ptrdiff_t nops = 
                IntTraits2<IntTraits2<s,sp1>::safeprod,sp2>::safeprod / 6;
            const bool unroll = 
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_UU_UNROLL;
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? ( M3::_unit ? 0 : 1 ) :
                unroll ? 26 : 
                // For known s, always recurse down to unroll size
                s != Unknown ? 0 :
                (TMV_UU_RECURSE == 1) ? ( M3::_unit ? 0 : 1) :
                rxr ? 22 : crx ? 23 : xcc ? 21 : 23;
            const int algo3 =  // The algorithm for N > UU_RECURSE
                unroll || s == 1 ? 0 : 27;
            const int algo4 =  // The algorithm for MultUM
                unroll || s == 1 ? 0 : -2;
#ifdef PRINTALGO_UU
            std::cout<<"algo2,3,4 = "<<algo2<<"  "<<algo3<<"  "<<algo4<<std::endl;
#endif

            if (s==Unknown ? (N > TMV_UU_RECURSE) : (s > 1 && !unroll)) {
                // [ C00  0  ] = [ A00  0  ] [ B00  0  ]
                // [ C10 C11 ]   [ A10 A11 ] [ B10 B11 ]

                const ptrdiff_t Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round N/2 up to a multiple of 16.)
                const ptrdiff_t sx = IntTraits<s>::half_roundup;
                const ptrdiff_t sy = IntTraits2<s,sx>::diff;

                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::subtrimatrix_type M3st;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1st A00 = m1.cSubTriMatrix(0,Nx);
                M1smt A10t = m1.cSubMatrix(Nx,N,0,Nx).transpose();
                M1st A11 = m1.cSubTriMatrix(Nx,N);
                M2st B00 = m2.cSubTriMatrix(0,Nx);
                M2stt B00t = B00.transpose();
                M2sm B10 = m2.cSubMatrix(Nx,N,0,Nx);
                M2st B11 = m2.cSubTriMatrix(Nx,N);
                M3st C00 = m3.cSubTriMatrix(0,Nx);
                M3sm C10 = m3.cSubMatrix(Nx,N,0,Nx);
                M3smt C10t = C10.transpose();
                M3st C11 = m3.cSubTriMatrix(Nx,N);

                // C10 (+)= A10 B00
                MultUM_Helper<algo4,sx,sy,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B00t,A10t,C10t);

                // C10 += A11 B10
                MultUM_Helper<algo4,sy,sx,true,ix,T,M1st,M2sm,M3sm>::call(
                    x,A11,B10,C10);

                // C00 (+)= x A00 B00
                MultUU_Helper<algo3,sx,add,ix,T,M1st,M2st,M3st>::call(
                    x,A00,B00,C00);

                // C11 (+)= x A11 B11
                MultUU_Helper<algo3,sy,add,ix,T,M1st,M2st,M3st>::call(
                    x,A11,B11,C11);
            } else {
                MultUU_Helper<algo2,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 83: Use temporary for m1*m2
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<83,s,add,ix,T,M1,M2,M3>
    {
        // This algo is slightly complicated by the fact that we allow
        // UnknownDiag.  If m3 is UnknownDiag, then its _shape dictates
        // NonUnitDiag.  But that means at m3c is not copyable back to m3 
        // directly (can't copy NonUnitDiag -> UnitDiag).
        // So we have another layer of indirection at the end to make sure 
        // that an UnknownDiag m3 is copied correctly.
        template <bool unknowndiag, int dummy>
        struct copyBack
        { // unknowndiag = false
            template <class M3c>
            static TMV_INLINE void call(
                const Scaling<ix,T>& x, const M3c& m3c, M3& m3)
            { MultXU_Helper<-2,s,add,ix,T,M3c,M3>::call(x,m3c,m3); }
        };
        template <int dummy>
        struct copyBack<true,dummy>
        {
            template <class M3c>
            static void call(
                const Scaling<ix,T>& x, const M3c& m3c, M3& m3)
            {
                TMVStaticAssert(!add);
                if (m3.isunit()) {
                    TMVAssert(T(x) == T(1));
                    typedef typename M3c::const_unitdiag_type M3cu;
                    M3cu m3cu = m3c.viewAsUnitDiag();
                    CopyU_Helper<-2,s,M3cu,M3>::call(m3cu,m3);
                } else {
                    MultXU_Helper<-2,s,add,ix,T,M3c,M3>::call(x,m3c,m3);
                }
            }
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UU
            const ptrdiff_t N = s == Unknown ? m3.size() : s;
            std::cout<<"UU algo 83: N,s,x = "<<N<< ','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type PT3;
            const int A = M1::_rowmajor && M2::_rowmajor ? RowMajor : ColMajor;
            const ptrdiff_t s3 = M3::_shape;
            typedef typename MCopyHelper<PT3,s3,s,s,A>::type M3c;
            // can't be unitdiag unless x == 1 and !add and m1,m2 are unit
            const bool unknowndiag = 
                M3::_unknowndiag && ix != -1 && !add &&
                (M1::_unit || M1::_unknowndiag) && 
                (M2::_unit || M2::_unknowndiag);
            copyBack<unknowndiag,1>::call(x,M3c(m1*m2),m3);
        }
    };

    // algo 90: call inst
    template <ptrdiff_t s, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<90,s,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <ptrdiff_t s, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<90,s,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 91: call inst alias
    template <ptrdiff_t s, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<91,s,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <ptrdiff_t s, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<91,s,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 96: Transpose
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<96,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            MultUU_Helper<-2,s,add,ix,T,M2t,M1t,M3t>::call(x,m2t,m1t,m3t);
        }
    };

    // algo 196: Transpose
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<196,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            MultUU_Helper<99,s,add,ix,T,M2t,M1t,M3t>::call(x,m2t,m1t,m3t);
        }
    };

    // algo 97: Conjugate
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<97,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultUU_Helper<-2,s,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<197,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultUU_Helper<99,s,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: Inline check for aliases
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<98,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // We set up the algorithm so that m1 can be in the 
            // same (exact or opposite) storage as m3.
            // Also, m2 can be in the opposite storage..
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_UU
            std::cout<<"UU Check aliases:\n";
            std::cout<<"s1, s2 = "<<s1<<"  "<<s2<<std::endl;
            std::cout<<"Exact13 = "<<(s1 && ExactSameStorage(m1,m3))<<std::endl;
            std::cout<<"Opp13 = "<<(s1 && OppositeStorage(m1,m3))<<std::endl;
            std::cout<<"Exact23 = "<<(s2 && ExactSameStorage(m2,m3))<<std::endl;
            std::cout<<"Opp23 = "<<(s2 && OppositeStorage(m2,m3))<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
#endif
            if ( (!s1 || ExactSameStorage(m1,m3) || OppositeStorage(m1,m3)) &&
                 (!s2 || OppositeStorage(m2,m3)) ) {
                // No aliasing (or no clobbering)
#ifdef PRINTALGO_UU
                std::cout<<"No clobber\n";
#endif
                MultUU_Helper<-2,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if ( 
                (!s1 || OppositeStorage(m1,m3)) &&
                (!s2 || ExactSameStorage(m2,m3) || OppositeStorage(m2,m3)) ) {
                // Can transpose to get no clobbering storage
#ifdef PRINTALGO_UU
                std::cout<<"Transpose to get no clobber.\n";
#endif
                MultUU_Helper<96,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (!s2) {
                // Copy m1 to m3, becomes MultEq op
#ifdef PRINTALGO_UU
                std::cout<<"Copy m1 to m3\n";
#endif
                Copy(m1,m3);
                typedef typename M3::const_view_type M3c;
                M3c m3c = m3;
                MultUU_Helper<-2,s,add,ix,T,M3c,M2,M3>::call(x,m3c,m2,m3);
            } else if (!s1) {
                // Copy m2 to m3, becomes MultEq op
#ifdef PRINTALGO_UU
                std::cout<<"Copy m2 to m3\n";
#endif
                Copy(m2,m3);
                typedef typename M1::const_transpose_type M1t;
                typedef typename M3::const_transpose_type M3ct;
                typedef typename M3::transpose_type M3t;
                M1t m1t = m1.transpose();
                M3t m3t = m3.transpose();
                M3ct m3ct = m3t;
                MultUU_Helper<-2,s,add,ix,T,M3ct,M1t,M3t>::call(x,m3ct,m1t,m3t);
            } else {
                // Use temporary for m1*m2
#ifdef PRINTALGO_UU
                std::cout<<"Use temporary\n";
#endif
                MultUU_Helper<83,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<99,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? ( M3::_unit ? 0 : 1 ) :
                M3::_conj ? 197 :
                inst ? 91 : 
                98;
            MultUU_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -4: No branches or copies.
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<-4,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool upper = M1::_upper;
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const ptrdiff_t s2 = s > 20 ? Unknown : s;
            const ptrdiff_t s2p1 = IntTraits<s2>::Sp1;
            const ptrdiff_t s2p2 = IntTraits<s2p1>::Sp1;
            // nops = 1/6 n(n+1)(n+2)
            const ptrdiff_t nops = 
                IntTraits2<IntTraits2<s2,s2p1>::safeprod,s2p2>::safeprod / 6;
            const bool unroll = 
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_UU_UNROLL;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? ( M3::_unit ? 0 : 1 ) :
                unroll ? ( upper ? 16 : 26 ) :

                upper ? 
                TMV_OPT >= 1 ? 17 :
                rxr ? 12 : crx ? 13 : xcc ? 11 : 13 :
                
                // lower
                TMV_OPT >= 1 ? 27 :
                rxr ? 22 : crx ? 23 : xcc ? 21 : 23;
            MultUU_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<-3,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // Possible algorithms are:
            //
            // Trivial and special for small TriMatrix sizes.
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar product function
            //
            // UpperTri:
            // 11 = loop over n: MultUV
            // 12 = loop over m: MultUV 
            // 13 = loop over k: Rank1
            // 16 = unroll small case
            // 17 = split each trimatrix into 3 submatrices and recurse
            // 
            // LowerTri:
            // 21 = loop over n: MultUV
            // 22 = loop over m: MultUV 
            // 23 = loop over k: Rank1
            // 26 = unroll small case
            // 27 = split each trimatrix into 3 submatrices and recurse
            //
            // Copy matrices to new storage
            // 83 = temp m1*m2
            
            const int algo = 
                // TODO: Add checks for bad majority on m1,m2,m3.
                -4;
#ifdef PRINTALGO_UU
            const ptrdiff_t N = s==Unknown ? m3.size() : s;
            std::cout<<"InlineMultUU: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"N = "<<N<<std::endl;
            std::cout<<"s = "<<s<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_UU
            typedef typename M3::real_type RT;
            typedef typename M3::value_type T3;
            Matrix<T3> m1c = m1;
            Matrix<T3> m2c = m2;
            Matrix<T3> m3i = m3;
            Matrix<T3> m3c = m3;
            MultMM<add>(x,m1c,m2c,m3c);
#endif
            MultUU_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#ifdef XDEBUG_UU
            if (Norm(m3-m3c) > 1.e-3*(Norm(m1c)*Norm(m2c)+(add?Norm(m3i):RT(0)))) {
                std::cout<<"m1 = "<<m1c<<std::endl;
                std::cout<<"m2 = "<<m2c<<std::endl;
                std::cout<<"m3 = "<<m3i<<std::endl;
                std::cout<<"m3 => "<<m3<<std::endl;
                std::cout<<"Correct m3 = "<<m3c<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<-2,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? ( M3::_unit ? 0 : 1 ) :
                M3::_conj ? 97 :
                inst ? 90 : 
                -3;
            MultUU_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<-1,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? ( M3::_unit ? 0 : 1 ) :
                M3::_checkalias ? 99 : 
                -2;
            MultUU_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M1::_size,M3::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVStaticAssert(!M3::_unit || ix == 1);
        TMVAssert(m1.size() == m3.size());
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m3.isunit() || T(x) == T(1));

        const ptrdiff_t s = Sizes<Sizes<M1::_size,M2::_size>::size,M3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultUU_Helper<-1,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M1::_size,M3::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVStaticAssert(!M3::_unit || ix == 1);
        TMVAssert(m1.size() == m3.size());
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m3.isunit() || T(x) == T(1));

        const ptrdiff_t s = Sizes<Sizes<M1::_size,M2::_size>::size,M3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultUU_Helper<-3,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M1::_size,M3::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVStaticAssert(!M3::_unit || ix == 1);
        TMVAssert(m1.size() == m3.size());
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m3.isunit() || T(x) == T(1));

        const ptrdiff_t s = Sizes<Sizes<M1::_size,M2::_size>::size,M3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultUU_Helper<98,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <class M1, int ix, class T, class M2>
    inline void MultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(!M1::_unit || ix == 1);
        TMVAssert(!m1.isunit() || T(x) == T(1));
        MultMM<false>(x,m1.mat(),m2.mat(),m1.mat());
    }

} // namespace tmv

#endif 
