

#ifndef TMV_MultBM_H
#define TMV_MultBM_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_MultBV.h"
#include "TMV_MultMV.h"
#include "TMV_Rank1VVM.h"
#include "TMV_MultMD.h"
#include "TMV_SmallDiagMatrix.h"
#include "TMV_MultMM_Funcs.h"
#include "TMV_MultXM_Funcs.h"

#ifdef PRINTALGO_BM
#include <iostream>
#endif

#ifdef XDEBUG_BM
#include <iostream>
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_BandMatrixIO.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_SumMM.h"
#endif

// ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
#define TMV_BM_ZeroIX (ix==0)
//#define TMV_BM_ZeroIX (ix!=1)

namespace tmv {

    // Defined in TMV_MultMM.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);

    //
    // BandMatrix * Matrix
    //

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper;

    // algo 0: Trivial, nothing to do.
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<0,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {}
    };

    // algo 1: Trivial, just m3.setZero();
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<1,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        { m3.setZero(); }
    };

    // algo 2: cs == 1, so reduces to MultBV
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
            std::cout<<"MB algo 2: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<
                ','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_colrange_type M2cr;
            typedef typename M2cr::const_transpose_type M2crt;
            typedef typename M3::row_type M3r;
            const ptrdiff_t xx = Unknown;
            M1r m1r = m1.get_row(0,0,m1.nhi()+1);
            M2crt m2t = m2.cColRange(0,m1.nhi()+1).transpose();
            M3r m3r = m3.get_row(0);
            MultMV_Helper<-4,xx,xs,add,ix,T,M2crt,M1r,M3r>::call(x,m2t,m1r,m3r);
        }
    };

    // algo 102: same as 2, but use -1 algo
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<102,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
            std::cout<<"MB algo 102: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<
                ','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_colrange_type M2cr;
            typedef typename M2cr::const_transpose_type M2crt;
            typedef typename M3::row_type M3r;
            const ptrdiff_t xx = Unknown;
            M1r m1r = m1.get_row(0,0,m1.nhi()+1);
            M2crt m2t = m2.cColRange(0,m1.nhi()+1).transpose();
            M3r m3r = m3.get_row(0);
            MultMV_Helper<-1,xx,xs,add,ix,T,M2crt,M1r,M3r>::call(x,m2t,m1r,m3r);
        }
    };

    // algo 202: same as 2, but use -2 algo
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<202,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
            std::cout<<"MB algo 202: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<
                ','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_colrange_type M2cr;
            typedef typename M2cr::const_transpose_type M2crt;
            typedef typename M3::row_type M3r;
            const ptrdiff_t xx = Unknown;
            M1r m1r = m1.get_row(0,0,m1.nhi()+1);
            M2crt m2t = m2.cColRange(0,m1.nhi()+1).transpose();
            M3r m3r = m3.get_row(0);
            MultMV_Helper<-2,xx,xs,add,ix,T,M2crt,M1r,M3r>::call(x,m2t,m1r,m3r);
        }
    };

    // algo 3: rs == 1, so reduces to MultBV
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<3,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
            std::cout<<"MB algo 3: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<
                ','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            M2c m2c = m2.get_col(0);
            M3c m3c = m3.get_col(0);
            MultBV_Helper<-4,cs,xs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
        }
    };

    // algo 103: same as 3, but use -1 algo
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<103,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
            std::cout<<"MB algo 103: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<
                ','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            M2c m2c = m2.get_col(0);
            M3c m3c = m3.get_col(0);
            MultBV_Helper<-1,cs,xs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
        }
    };

    // algo 203: same as 3, but use -2 algo
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<203,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
            std::cout<<"MB algo 203: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<
                ','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            M2c m2c = m2.get_col(0);
            M3c m3c = m3.get_col(0);
            MultBV_Helper<-2,cs,xs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
        }
    };

    // algo 4: xs == 1, so reduces to Rank1Update
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<4,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            std::cout<<"MM algo 4: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<
                ','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::rowrange_type M3rr;
            const ptrdiff_t xx = Unknown;
            M1c m1c = m1.get_col(0,0,m1.nlo()+1);
            M2r m2r = m2.get_row(0);
            M3rr m3rr = m3.cRowRange(0,m1.nlo()+1);
            Rank1VVM_Helper<-4,xx,rs,add,ix,T,M1c,M2r,M3rr>::call(
                x,m1c,m2r,m3rr);
        }
    };

    // algo 104: same as 4, but use -1 algo
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<104,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            std::cout<<"MM algo 104: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<
                ','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::rowrange_type M3rr;
            const ptrdiff_t xx = Unknown;
            M1c m1c = m1.get_col(0,0,m1.nlo()+1);
            M2r m2r = m2.get_row(0);
            M3rr m3rr = m3.cRowRange(0,m1.nlo()+1);
            Rank1VVM_Helper<-1,xx,rs,add,ix,T,M1c,M2r,M3rr>::call(
                x,m1c,m2r,m3rr);
        }
    };

    // algo 204: same as 4, but use -2 algo
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<204,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            std::cout<<"MM algo 204: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<
                ','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::rowrange_type M3rr;
            const ptrdiff_t xx = Unknown;
            M1c m1c = m1.get_col(0,0,m1.nlo()+1);
            M2r m2r = m2.get_row(0);
            M3rr m3rr = m3.cRowRange(0,m1.nlo()+1);
            Rank1VVM_Helper<-2,xx,rs,add,ix,T,M1c,M2r,M3rr>::call(
                x,m1c,m2r,m3rr);
        }
    };


    // algo 11: ?CC -- Loop over N
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<11,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
            std::cout<<"MB algo 11: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            for(ptrdiff_t j=0;j<N;++j) {
                // m3.get_col(j) += x * m1 * m2.get_col(j)
                M2c m2j = m2.get_col(j);
                M3c m3j = m3.get_col(j);
                MultBV_Helper<-4,cs,xs,add,ix,T,M1,M2c,M3c>::call(
                    x,m1,m2j,m3j);
            }
        }
    };

    // algo 21: R?R -- Loop over M
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<21,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
#ifdef PRINTALGO_BM
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            std::cout<<"MB algo 21: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_rowrange_type M2rr;
            typedef typename M2rr::const_transpose_type M2rrt;
            typedef typename M3::row_type M3r;
            const ptrdiff_t xx = Unknown;

            const ptrdiff_t lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const ptrdiff_t i1 = m1.nlo();
            const ptrdiff_t i2 = TMV_MIN(M,K-m1.nhi());
            const ptrdiff_t i3 = TMV_MIN(M,K+m1.nlo());
            ptrdiff_t j1 = 0;
            ptrdiff_t j2 = m1.nhi()+1;

            ptrdiff_t i=0;
            for(;i<i1;++i) {
                // m3.get_row(j) += x * m1.get_row(j) * m2
                M1r m1i = m1.get_row(i,j1,j2);
                M2rrt m2it = m2.cRowRange(j1,j2).transpose();
                M3r m3i = m3.get_row(i);
                MultMV_Helper<-4,rs,xx,add,ix,T,M2rrt,M1r,M3r>::call(
                    x,m2it,m1i,m3i);
                if (j2 < K) ++j2;
            }
            if (i1 < i2) TMVAssert(j2 == m1.nlo()+m1.nhi()+1);
            for(;i<i2;++i) {
                M1r m1i = m1.get_row(i,j1,j2);
                M2rrt m2it = m2.cRowRange(j1,j2).transpose();
                M3r m3i = m3.get_row(i);
                MultMV_Helper<-4,rs,lh,add,ix,T,M2rrt,M1r,M3r>::call(
                    x,m2it,m1i,m3i);
                ++j1; ++j2;
            }
            for(;i<i3;++i) {
                M1r m1i = m1.get_row(i,j1,K);
                M2rrt m2it = m2.cRowRange(j1,K).transpose();
                M3r m3i = m3.get_row(i);
                MultMV_Helper<-4,rs,xx,add,ix,T,M2rrt,M1r,M3r>::call(
                    x,m2it,m1i,m3i);
                ++j1; 
            }
            if (i3 < M) Maybe<!add>::zero2(m3.cRowRange(i3,M));
        }
    };

    // algo 31: CR? -- Loop over K
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<31,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
#ifdef PRINTALGO_BM
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            std::cout<<"MB algo 31: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::rowrange_type M3rr;
            const ptrdiff_t xx = Unknown;

            const ptrdiff_t lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const ptrdiff_t j1 = m1.nhi();
            const ptrdiff_t j2 = TMV_MIN(K,M-m1.nlo());
            const ptrdiff_t j3 = TMV_MIN(K,M+m1.nhi());
            ptrdiff_t i1 = 0;
            ptrdiff_t i2 = m1.nlo()+1;

            Maybe<!add>::zero(m3);
            ptrdiff_t j=0;
            for(;j<j1;++j) {
                M1c m1j = m1.get_col(j,i1,i2);
                M2r m2j = m2.get_row(j);
                M3rr m3j = m3.cRowRange(i1,i2);
                Rank1VVM_Helper<-4,xx,rs,true,ix,T,M1c,M2r,M3rr>::call(
                    x,m1j,m2j,m3j);
                if (i2 < M) ++i2;
            }
            if (j1 < j2) TMVAssert(i2-i1 == m1.nlo()+m1.nhi()+1);
            for(;j<j2;++j) {
                M1c m1j = m1.get_col(j,i1,i2);
                M2r m2j = m2.get_row(j);
                M3rr m3j = m3.cRowRange(i1,i2);
                Rank1VVM_Helper<-4,lh,rs,true,ix,T,M1c,M2r,M3rr>::call(
                    x,m1j,m2j,m3j);
                ++i1; ++i2;
            }
            for(;j<j3;++j) {
                M1c m1j = m1.get_col(j,i1,M);
                M2r m2j = m2.get_row(j);
                M3rr m3j = m3.cRowRange(i1,M);
                Rank1VVM_Helper<-4,xx,rs,true,ix,T,M1c,M2r,M3rr>::call(
                    x,m1j,m2j,m3j);
                ++i1; 
            }
        }
    };

    // algo 41: DRR, DCC -- Loop over diagonals
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<41,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
#ifdef PRINTALGO_BM
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            std::cout<<"MB algo 41: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_diag_sub_type M1d;
            typedef typename DMVO<M1d>::cv M1dm;
            typedef typename M2::const_rowrange_type M2rr;
            typedef typename M2rr::const_transpose_type M2rrt;
            typedef typename M3::rowrange_type M3rr;
            typedef typename M3rr::transpose_type M3rrt;

            const ptrdiff_t xx = Unknown;
            ptrdiff_t len = TMV_MIN(M-m1.nlo(),K);
            const ptrdiff_t ds = IntTraits2<cs,xs>::min;

            Maybe<!add>::zero(m3);
            for(ptrdiff_t k=m1.nlo();k;--k) {
                M1dm m1k = DiagMatrixViewOf(m1.diag(-k,0,len));
                M2rrt m2k = m2.cRowRange(0,len).transpose();
                M3rrt m3k = m3.cRowRange(k,k+len).transpose();
                MultMD_Helper<-4,rs,xx,true,ix,T,M2rrt,M1dm,M3rrt>::call(
                    x,m2k,m1k,m3k);
                if (len < K) ++len;
            }
            TMVAssert(len == TMV_MIN(M,K));
            M1dm m10 = DiagMatrixViewOf(m1.diag());
            M2rrt m20 = m2.cRowRange(0,len).transpose();
            M3rrt m30 = m3.cRowRange(0,len).transpose();
            MultMD_Helper<-4,rs,ds,true,ix,T,M2rrt,M1dm,M3rrt>::call(
                x,m20,m10,m30);
            for(ptrdiff_t k=1;k<=m1.nhi();++k) {
                if (k+len > K) --len;
                M1dm m1k = DiagMatrixViewOf(m1.diag(k,0,len));
                M2rrt m2k = m2.cRowRange(k,k+len).transpose();
                M3rrt m3k = m3.cRowRange(0,len).transpose();
                MultMD_Helper<-4,rs,xx,true,ix,T,M2rrt,M1dm,M3rrt>::call(
                    x,m2k,m1k,m3k);
            }
        }
    };

    template <int ix, class T, class M> class ProdXM;

    // algo 81: copy x*m1
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<81,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs == Unknown ? m3.colsize() : cs;
            const ptrdiff_t K = xs == Unknown ? m1.rowsize() : xs;
            const ptrdiff_t N = rs == Unknown ? m3.rowsize() : rs;
            std::cout<<"MB algo 81: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M3::real_type RT;
            const Scaling<1,RT> one;
            typedef typename Traits2<T,T1>::type PT1;
            const int A = M2::_colmajor || M3::_rowmajor ? RowMajor : ColMajor;
            typedef typename MCopyHelper<PT1,Band,cs,xs,A>::type M1c;
            typename M3::noalias_type m3na = m3.noAlias();
            MultMM<add>(one,M1c(ProdXM<ix,T,M1>(x,m1)),m2,m3na);
        }
    };

    // algo 82: copy x*m2
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<82,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t N = rs == Unknown ? m3.rowsize() : rs;
            const ptrdiff_t K = xs == Unknown ? m1.rowsize() : rs;
            const ptrdiff_t M = cs == Unknown ? m3.colsize() : cs;
            std::cout<<"MB algo 82: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::value_type T2;
            typedef typename M3::real_type RT;
            const Scaling<1,RT> one;
            typedef typename Traits2<T,T2>::type PT2;
            const int A = M1::_colmajor && M3::_rowmajor ? RowMajor : ColMajor;
            typedef typename MCopyHelper<PT2,Rec,xs,rs,A>::type M2c;
            typename M3::noalias_type m3na = m3.noAlias();
            MultMM<add>(one,m1,M2c(ProdXM<ix,T,M2>(x,m2)),m3na);
        }
    };

    // algo 83: Use temporary for m1*m2
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<83,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs == Unknown ? m3.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m3.rowsize() : rs;
            const ptrdiff_t K = xs == Unknown ? m1.rowsize() : rs;
            std::cout<<"MB algo 83: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type PT3;
            const int A = M1::_rowmajor && M2::_rowmajor ? RowMajor : ColMajor;
            typedef typename MCopyHelper<PT3,Rec,cs,rs,A>::type M3c;
            typename M3::noalias_type m3na = m3.noAlias();
            MultXM<add>(x,M3c(m1*m2),m3na);
        }
    };

    // algo 90: call inst
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<90,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<90,cs,rs,xs,true,ix,T,M1,M2,M3>
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
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<91,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<91,cs,rs,xs,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 97: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<97,cs,rs,xs,add,ix,T,M1,M2,M3>
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
            MultBM_Helper<-2,cs,rs,xs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<197,cs,rs,xs,add,ix,T,M1,M2,M3>
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
            MultBM_Helper<99,cs,rs,xs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: Inline check for aliases
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<98,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_BM
            std::cout<<"MB algo 98: \n";
            std::cout<<"s1 = "<<s1<<std::endl;
            std::cout<<"s2 = "<<s2<<std::endl;
#endif
            if (!s1 && !s2) {
                // No aliasing
                MultBM_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (s1 && !s2) {
                // Use temporary for m1
                MultBM_Helper<81,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (!s1 && s2) {
                // Use temporary for m2
                MultBM_Helper<82,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else {
                // Use temporary for m1*m2
                MultBM_Helper<83,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<99,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                (xs == Unknown || xs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                M3::_conj ? 197 :
                inst ? 91 : 
                98;
            MultBM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -4: No branches or copies
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<-4,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xrr = M2::_rowmajor && M3::_rowmajor;
            const bool dcc = M1::_diagmajor && M2::_colmajor && M3::_colmajor;

            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 2 :
                rs == 1 ? 3 :
                dcc ? 41 : xcc ? 11 : rxr ? 21 : crx ? 31 : xrr ? 41 : 
                M3::_rowmajor ? 21 :
                M2::_rowmajor ? 31 :
                11;
            MultBM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<-3,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;

            // Possible algorithms are:
            //
            // Trivial:
            //  0 = cs or rs == 0, or xs == 0 && add, so nothing to do
            //  1 = xs == 0 && !add, so m3.setZero();
            //  2 = cs == 1: reduces to trivial MultMV function
            //  3 = rs == 1: reduces to trivial MultMV function
            //
            // TODO: Add possibility that m1 is upper or lower banded.
            // 11 = ?CC, loop over n
            // 21 = R?R, loop over m
            // 31 = CR?, loop over k
            // 41 = DRR, DCC, loop over diagonals
            //
            // 81 = copy x*m1
            // 82 = copy x*m2
            // 83 = temp m1*m2

            const bool xcc = M2::_colmajor && M3::_colmajor;
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xrr = M2::_rowmajor && M3::_rowmajor;
            const bool dcc = M1::_diagmajor && M2::_colmajor && M3::_colmajor;

            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 202 :
                rs == 1 ? 203 :
                !(M1::_rowmajor || M1::_colmajor || M1::_diagmajor) ? 81 :
                !(M2::_rowmajor || M2::_colmajor) ? 82 :
                !(M3::_rowmajor || M3::_colmajor) ? 83 :
                dcc ? 41 : xcc ? 11 : rxr ? 21 : crx ? 31 : xrr ? 41 : 
                M3::_rowmajor ? 21 :
                M2::_rowmajor ? 31 :
                // The above should have caught all possibilities.
                // So selection should never get here.
                -999;
                TMVStaticAssert(algo != -999);
#ifdef PRINTALGO_BM
            const ptrdiff_t M = cs==Unknown ? m3.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? m3.rowsize() : rs;
            const ptrdiff_t K = xs==Unknown ? m1.rowsize() : xs;
            std::cout<<"InlineMultBM: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"M = "<<M<<"  N = "<<N<<"  K = "<<K<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<"  xs = "<<xs<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_BM
            typedef typename M3::real_type RT;
            typedef typename M3::value_type T3;
            Matrix<T3> m1c = m1;
            Matrix<T3> m2c = m2;
            Matrix<T3> m3i = m3;
            Matrix<T3> m3c = m3;
            MultMM<add>(x,m1c,m2c,m3c);
#endif
            MultBM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
#ifdef XDEBUG_BM
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

    // algo -402: Same as algo -2, but use -4 if no Inst version
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<-402,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                (xs == Unknown || xs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 2 :
                rs == 1 ? 3 :
                xs == 1 ? 4 :
                inst ? 90 : 
                -4;
            MultBM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                (xs == Unknown || xs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 202 :
                rs == 1 ? 203 :
                xs == 1 ? 204 :
                M3::_conj ? 97 :
                inst ? 90 : 
                -3;
            MultBM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBM_Helper<-1,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 102 :
                rs == 1 ? 103 :
                xs == 1 ? 104 :
                M3::_checkalias ? 99 : 
                -2;
            MultBM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        // cs = m3.colsize
        // rs = m3.rowsize
        // xs = "extra size" = m1.rowsize, m2.colsize
        const ptrdiff_t cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const ptrdiff_t rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const ptrdiff_t xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultBM_Helper<-1,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const ptrdiff_t cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const ptrdiff_t rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const ptrdiff_t xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultBM_Helper<-3,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const ptrdiff_t cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const ptrdiff_t rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const ptrdiff_t xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultBM_Helper<98,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();;
        MultMM<add>(x,m2.transpose(),m1.transpose(),m3t);
    }


} // namespace tmv

#endif 
