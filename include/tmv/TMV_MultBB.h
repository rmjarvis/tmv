

#ifndef TMV_MultBB_H
#define TMV_MultBB_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_MultBV.h"
#include "TMV_MultMV.h"
#include "TMV_MultBD.h"
#include "TMV_Rank1VVM.h"
#include "TMV_MultMM_Funcs.h"
#include "TMV_MultXM_Funcs.h"
#include "TMV_SmallDiagMatrix.h"
#include "TMV_SmallBandMatrix.h"

#ifdef PRINTALGO_BB
#include <iostream>
#endif

#ifdef XDEBUG_BB
#include <iostream>
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_BandMatrixIO.h"
#include "tmv/TMV_NormB.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_AddBB.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_MultMM.h"
#endif

// ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
#define TMV_BB_ZeroIX (ix==0)
//#define TMV_BB_ZeroIX (ix!=1)

namespace tmv {

    // Defined in TMV_MultBB.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstBandMatrixView<T2,C2>& m2, BandMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstBandMatrixView<T2,C2>& m2, BandMatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstBandMatrixView<T2,C2>& m2, BandMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstBandMatrixView<T2,C2>& m2, BandMatrixView<T3> m3);

    //
    // BandMatrix * BandMatrix
    //

    template <int algo, int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper;

    // algo 0: Trivial, nothing to do.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<0,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {}
    };

    // algo 1: Trivial, just m3.setZero();
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<1,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        { m3.setZero(); }
    };

    // algo 2: cs == 1, so reduces to MultBV
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            std::cout<<"BB algo 2: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<
                ','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_colrange_type M2cr;
            typedef typename M2cr::const_transpose_type M2crt;
            typedef typename M3::row_sub_type M3r;
            const int xx = TMV_UNKNOWN;
            M1r m1r = m1.get_row(0,0,m1.nhi()+1);
            M2crt m2t = m2.cColRange(0,m1.nhi()+1).transpose();
            M3r m3r = m3.get_row(0,0,m2t.colsize());
            MultBV_Helper<-4,xx,xx,add,ix,T,M2crt,M1r,M3r>::call(x,m2t,m1r,m3r);
        }
    };

    // algo 102: same as 2, but use -1 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<102,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            std::cout<<"BB algo 102: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<
                ','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_colrange_type M2cr;
            typedef typename M2cr::const_transpose_type M2crt;
            typedef typename M3::row_sub_type M3r;
            const int xx = TMV_UNKNOWN;
            M1r m1r = m1.get_row(0,0,m1.nhi()+1);
            M2crt m2t = m2.cColRange(0,m1.nhi()+1).transpose();
            M3r m3r = m3.get_row(0,0,m2t.colsize());
            MultBV_Helper<-1,xx,xx,add,ix,T,M2crt,M1r,M3r>::call(x,m2t,m1r,m3r);
        }
    };

    // algo 202: same as 2, but use -2 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<202,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            std::cout<<"BB algo 202: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<
                ','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_colrange_type M2cr;
            typedef typename M2cr::const_transpose_type M2crt;
            typedef typename M3::row_sub_type M3r;
            const int xx = TMV_UNKNOWN;
            M1r m1r = m1.get_row(0,0,m1.nhi()+1);
            M2crt m2t = m2.cRowRange(0,m1.nhi()+1).transpose();
            M3r m3r = m3.get_row(0,0,m2t.colsize());
            MultBV_Helper<-2,xx,xx,add,ix,T,M2crt,M1r,M3r>::call(x,m2t,m1r,m3r);
        }
    };

    // algo 3: rs == 1, so reduces to MultBV
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<3,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            std::cout<<"BB algo 3: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<
                ','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_colrange_type M1cr;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int xx = TMV_UNKNOWN;
            M2c m2c = m2.get_col(0,0,m2.nlo()+1);
            M1cr m1cr = m1.cColRange(0,0,m2.nlo()+1);
            M3c m3c = m3.get_col(0,0,m1cr.colsize());
            MultBV_Helper<-4,xx,xx,add,ix,T,M1cr,M2c,M3c>::call(x,m1cr,m2c,m3c);
        }
    };

    // algo 103: same as 3, but use -1 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<103,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            std::cout<<"BB algo 103: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<
                ','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_colrange_type M1cr;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int xx = TMV_UNKNOWN;
            M2c m2c = m2.get_col(0,0,m2.nlo()+1);
            M1cr m1cr = m1.cColRange(0,0,m2.nlo()+1);
            M3c m3c = m3.get_col(0,0,m1cr.colsize());
            MultBV_Helper<-1,xx,xx,add,ix,T,M1cr,M2c,M3c>::call(x,m1cr,m2c,m3c);
        }
    };

    // algo 203: same as 3, but use -2 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<203,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            std::cout<<"BB algo 203: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<
                ','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_colrange_type M1cr;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int xx = TMV_UNKNOWN;
            M2c m2c = m2.get_col(0,0,m2.nlo()+1);
            M1cr m1cr = m1.cColRange(0,0,m2.nlo()+1);
            M3c m3c = m3.get_col(0,0,m1cr.colsize());
            MultBV_Helper<-2,xx,xx,add,ix,T,M1cr,M2c,M3c>::call(x,m1cr,m2c,m3c);
        }
    };

    // algo 4: xs == 1, so reduces to Rank1Update
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<4,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"MM algo 4: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<
                ','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::submatrix_type M3s;
            const int xx = TMV_UNKNOWN;
            M1c m1c = m1.get_col(0,0,m1.nlo()+1);
            M2r m2r = m2.get_row(0,0,m2.nhi()+1);
            M3s m3s = m3.cSubMatrix(0,m1.nlo()+1,0,m2.nhi()+1);
            Rank1VVM_Helper<-4,xx,rs,add,ix,T,M1c,M2r,M3s>::call(x,m1c,m2r,m3s);
        }
    };

    // algo 104: same as 4, but use -1 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<104,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"MM algo 104: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<
                ','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::submatrix_type M3s;
            const int xx = TMV_UNKNOWN;
            M1c m1c = m1.get_col(0,0,m1.nlo()+1);
            M2r m2r = m2.get_row(0,0,m2.nhi()+1);
            M3s m3s = m3.cSubMatrix(0,m1.nlo()+1,0,m2.nhi()+1);
            Rank1VVM_Helper<-1,xx,rs,add,ix,T,M1c,M2r,M3s>::call(x,m1c,m2r,m3s);
        }
    };

    // algo 204: same as 4, but use -2 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<204,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"MM algo 204: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<
                ','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::submatrix_type M3s;
            const int xx = TMV_UNKNOWN;
            M1c m1c = m1.get_col(0,0,m1.nlo()+1);
            M2r m2r = m2.get_row(0,0,m2.nhi()+1);
            M3s m3s = m3.cSubMatrix(0,m1.nlo()+1,0,m2.nhi()+1);
            Rank1VVM_Helper<-2,xx,rs,add,ix,T,M1c,M2r,M3s>::call(x,m1c,m2r,m3s);
        }
    };

    // algo 5: ??R -- transpose
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<5,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            std::cout<<"BB algo 5: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            MultBB_Helper<-3,rs,cs,xs,add,ix,T,M2t,M1t,M3t>::call(
                x,m2t,m1t,m3t);
        }
    };

    // algo 405: same as 5 but go back to -4, rather than -3
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<405,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            std::cout<<"BB algo 405: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            MultBB_Helper<-4,rs,cs,xs,add,ix,T,M2t,M1t,M3t>::call(
                x,m2t,m1t,m3t);
        }
    };

    // algo 11: ?CC -- Loop over N
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<11,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
#ifdef PRINTALGO_BB
            std::cout<<"BB algo 11: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_subbandmatrix_type M1s;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int xx = TMV_UNKNOWN;

            const int hi3 = TMV_MIN(N-1 , m1.nhi()+m2.nhi());
            const int lo3 = TMV_MIN(M-1 , m1.nlo()+m2.nlo());
            TMVAssert(hi3 <= m3.nhi());
            TMVAssert(lo3 <= m3.nlo());

            const int lh2 = IntTraits<IntTraits2<M2::_nlo,M2::_nhi>::sum>::Sp1;
            const int xhi3 = IntTraits2<M1::_nhi,M2::_nhi>::sum;
            const int xlo3 = IntTraits2<M1::_nlo,M2::_nlo>::sum;
            const int lh3 = IntTraits<IntTraits2<xlo3,xhi3>::sum>::Sp1;
            const int j1 = hi3;
            int j2 = TMV_MIN(TMV_MIN(N,M-lo3),K-m2.nlo());
            const int j3 = TMV_MIN(TMV_MIN(N,M+hi3),K+m2.nhi());
            if ((hi3 == N-1 || lo3 == M-1) && j1<j2) j2 = j1;
            int i1 = 0;
            int i2 = lo3+1;
            int k1 = 0;
            int k2 = m2.nlo()+1;
            int xlo = m1.nlo();
            int xhi = TMV_MIN(m1.nhi(),m2.nlo());

            int j=0;
            for(;j<j1;++j) {
                M1s m1j = m1.cSubBandMatrix(i1,i2,k1,k2,xlo,xhi);
                M2c m2j = m2.get_col(j,k1,k2);
                M3c m3j = m3.get_col(j,i1,i2);
                MultBV_Helper<-4,xx,xx,add,ix,T,M1s,M2c,M3c>::call(
                    x,m1j,m2j,m3j);
                if (i2 < M) ++i2;
                if (k2 < K) ++k2;
                if (j + m2.nlo() < m1.nhi()) ++xhi;
                if (j >= m2.nhi()) { ++k1; --xhi; if (xlo < M-1) ++xlo; }
            }
            if (j1 < j2) { 
                TMVAssert(i2 == lo3+hi3+1); 
                TMVAssert(hi3 == m1.nhi() + m2.nhi());
                TMVAssert(lo3 == m1.nlo() + m2.nlo());
                TMVAssert(k2-k1 == m2.nlo()+m2.nhi()+1);
            }
            for(;j<j2;++j) {
                M1s m1j = m1.cSubBandMatrix(i1,i2,k1,k2,xlo,xhi);
                M2c m2j = m2.get_col(j,k1,k2);
                M3c m3j = m3.get_col(j,i1,i2);
                MultBV_Helper<-4,lh3,lh2,add,ix,T,M1s,M2c,M3c>::call(
                    x,m1j,m2j,m3j);
                ++i1; ++i2; ++k1; ++k2;
            }
            if (k2 > K) --k2;
            if (M-i1 == xlo) --xlo;
            TMVAssert(k2 <= K);
            for(;j<j3;++j) {
                M1s m1j = m1.cSubBandMatrix(i1,M,k1,k2,xlo,xhi);
                M2c m2j = m2.get_col(j,k1,k2);
                M3c m3j = m3.get_col(j,i1,M);
                MultBV_Helper<-4,xx,xx,add,ix,T,M1s,M2c,M3c>::call(
                    x,m1j,m2j,m3j);
                ++i1; ++k1;
                if (k2 < K) ++k2;
                else if (k2-k1 == xhi) --xhi;
                if (M-i1 == xlo) --xlo;
                if (k1 == K) { 
                    if (++j < N) 
                        Maybe<!add>::zero2(
                            m3.cSubBandMatrix(i1,M,j,N,M-i1-1,0));
                    break;
                }
            }
        }
    };

    // algo 31: CR? -- Loop over K
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<31,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_BB
            std::cout<<"BB algo 31: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::submatrix_type M3s;
            const int xx = TMV_UNKNOWN;

            const int lh1 = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int lh2 = IntTraits<IntTraits2<M2::_nlo,M2::_nhi>::sum>::Sp1;
            const int k1 = TMV_MAX(m1.nhi(),m2.nlo());
            const int k2 = TMV_MIN(K,TMV_MIN(M-m1.nlo(),N-m2.nhi()));
            const int k3 = TMV_MIN(K,TMV_MIN(M+m1.nhi(),N+m2.nlo()));
            int i1 = 0;
            int i2 = m1.nlo()+1;
            int j1 = 0;
            int j2 = m2.nhi()+1;

            Maybe<!add>::zero(m3);
            int k=0;
            for(;k<k1;++k) {
                M1c m1k = m1.get_col(k,i1,i2);
                M2r m2k = m2.get_row(k,j1,j2);
                M3s m3k = m3.cSubMatrix(i1,i2,j1,j2);
                Rank1VVM_Helper<-4,xx,xx,true,ix,T,M1c,M2r,M3s>::call(
                    x,m1k,m2k,m3k);
                if (i2 < M) ++i2;
                if (j2 < N) ++j2;
                if (k >= m1.nhi()) ++i1; 
                if (k >= m2.nlo()) ++j1; 
            }
            if (k1 < k2) {
                TMVAssert(i2-i1 == m1.nlo()+m1.nhi()+1);
                TMVAssert(j2-j1 == m2.nlo()+m2.nhi()+1);
            }
            for(;k<k2;++k) {
                M1c m1k = m1.get_col(k,i1,i2);
                M2r m2k = m2.get_row(k,j1,j2);
                M3s m3k = m3.cSubMatrix(i1,i2,j1,j2);
                Rank1VVM_Helper<-4,lh1,lh2,true,ix,T,M1c,M2r,M3s>::call(
                    x,m1k,m2k,m3k);
                ++i1; ++i2; ++j1; ++j2;
            }
            if (i2 > M) --i2;
            if (j2 > N) --j2;
            TMVAssert(i2 <= M);
            TMVAssert(j2 <= N);
            for(;k<k3;++k) {
                M1c m1k = m1.get_col(k,i1,i2);
                M2r m2k = m2.get_row(k,j1,j2);
                M3s m3k = m3.cSubMatrix(i1,i2,j1,j2);
                Rank1VVM_Helper<-4,xx,xx,true,ix,T,M1c,M2r,M3s>::call(
                    x,m1k,m2k,m3k);
                if (++i1 == M) break;
                if (++j1 == N) break;
                if (i2 < M) ++i2;
                if (j2 < N) ++j2;
            }
        }
    };

    // algo 41: CDC -- Loop over diagonals
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<41,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
#ifdef PRINTALGO_BB
            std::cout<<"BB algo 41: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_subbandmatrix_type M1s;
            typedef typename M2::const_diag_sub_type M2d;
            typedef typename DMVO<M2d>::cv M2dm;
            typedef typename M3::subbandmatrix_type M3s;

            const int xx = TMV_UNKNOWN;
            int len = TMV_MIN(K-m2.nlo(),N);
            const int ds = IntTraits2<xs,rs>::min;
            int i1 = TMV_MAX(m2.nlo()-m1.nhi(),0);

            Maybe<!add>::zero(m3);
            for(int k=m2.nlo();k;--k) {
                // TODO: It's probably more efficient to keep track of
                // xlo,xhi values as we go rather than let them be 
                // calculated by cSubBandMatrix.
                M1s m1k = m1.cSubBandMatrix(i1,M,k,k+len);
                M2dm m2k = DiagMatrixViewOf(m2.diag(-k,0,len));
                M3s m3k = m3.cSubBandMatrix(i1,M,0,len);
                MultBD_Helper<-4,xx,xx,true,ix,T,M1s,M2dm,M3s>::call(
                    x,m1k,m2k,m3k);
                if (i1 > 0) --i1;
                if (len < N) ++len;
            }
            TMVAssert(len == TMV_MIN(K,N));
            M1s m10 = m1.cSubBandMatrix(0,M,0,len);
            M2dm m20 = DiagMatrixViewOf(m2.diag());
            M3s m30 = m3.cSubBandMatrix(0,M,0,len);
            MultBD_Helper<-4,xx,ds,true,ix,T,M1s,M2dm,M3s>::call(
                x,m10,m20,m30);
            for(int k=1;k<=m2.nhi();++k) {
                if (k+len > N) --len;
                M1s m1k = m1.cSubBandMatrix(0,M,0,len);
                M2dm m2k = DiagMatrixViewOf(m2.diag(k,0,len));
                M3s m3k = m3.cSubBandMatrix(0,M,k,k+len);
                MultBD_Helper<-4,xx,xx,true,ix,T,M1s,M2dm,M3s>::call(
                    x,m1k,m2k,m3k);
            }
        }
    };

    // algo 81: copy x*m1
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<81,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int M = cs == TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int K = xs == TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            const int N = rs == TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"BB algo 81: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M3::real_type RT;
            const Scaling<1,RT> one;
            typedef typename Traits2<T,T1>::type PT1;
            const int A = M2::_colmajor || M3::_rowmajor ? RowMajor : ColMajor;
            typedef typename MCopyHelper<PT1,Band,cs,xs,A>::type M1c;
            typename M3::noalias_type m3na = m3.noAlias();
            MultMM<add>(one,M1c(x*m1),m2,m3na);
        }
    };

    // algo 82: copy x*m2
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<82,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            const int N = rs == TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            const int M = cs == TMV_UNKNOWN ? int(m3.colsize()) : cs;
            std::cout<<"BB algo 82: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::value_type T2;
            typedef typename M3::real_type RT;
            const Scaling<1,RT> one;
            typedef typename Traits2<T,T2>::type PT2;
            // TODO: Once algo 41,42 are done, check for diagmajor options.
            const int A = M1::_colmajor && M3::_rowmajor ? RowMajor : ColMajor;
            typedef typename MCopyHelper<PT2,Band,xs,rs,A>::type M2c;
            typename M3::noalias_type m3na = m3.noAlias();
            MultMM<add>(one,m1,M2c(x*m2),m3na);
        }
    };

    // algo 83: Use temporary for m1*m2
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<83,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_BB
            std::cout<<"BB algo 83: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type PT3;
            const int A = M1::_rowmajor && M2::_rowmajor ? RowMajor : ColMajor;
            typedef typename MCopyHelper<PT3,Band,cs,rs,A>::type M3c;
            typename M3::noalias_type m3na = m3.noAlias();
            MultXM<add>(x,M3c(m1*m2),m3na);
        }
    };

    // algo 90: call inst
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<90,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo 90\n";
#endif
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<90,cs,rs,xs,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo 90\n";
#endif
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 91: call inst alias
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<91,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo 91\n";
#endif
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<91,cs,rs,xs,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo 91\n";
#endif
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<97,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo 97\n";
#endif
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultBB_Helper<-2,cs,rs,xs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<197,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo 197\n";
#endif
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultBB_Helper<99,cs,rs,xs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<98,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_BB
            std::cout<<"BB algo 98: \n";
            std::cout<<"s1 = "<<s1<<std::endl;
            std::cout<<"s2 = "<<s2<<std::endl;
#endif
            if (!s1 && !s2) {
                // No aliasing
                MultBB_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (s1 && !s2) {
                // Use temporary for m1
                MultBB_Helper<81,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (!s1 && s2) {
                // Use temporary for m2
                MultBB_Helper<82,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else {
                // Use temporary for m1*m2
                MultBB_Helper<83,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<99,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                (xs == TMV_UNKNOWN || xs > 16) &&
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
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo 99\n";
#endif
            MultBB_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<-4,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor && 
                (M3::_colmajor || M3::_rowmajor);
            const bool xdx = M2::_diagmajor && (
                (M1::_colmajor && M3::_colmajor) ||
                (M1::_rowmajor && M3::_rowmajor) ||
                (M1::_diagmajor && M3::_diagmajor) );
            const bool dxx = M1::_diagmajor && (
                (M2::_colmajor && M3::_colmajor) ||
                (M2::_rowmajor && M3::_rowmajor) ||
                (M2::_diagmajor && M3::_diagmajor) );

            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 2 :
                rs == 1 ? 3 :
                M3::_rowmajor && !M3::_colmajor ? 405 :
                xcc ? 11 : crx ? 31 : xdx ? 41 : dxx ? 405 :
                M2::_diagmajor ? 41 :
                M1::_diagmajor ? 405 :
                41;
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo -4\n";
#endif
            MultBB_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<-3,cs,rs,xs,add,ix,T,M1,M2,M3>
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
            //  5 = transpose to m3t += x * m2t * m1t
            //
            // TODO: Add in possibility that m1,m2 are upper or lower banded.
            // 11 = ?CC, loop over n
            // 31 = CR?, loop over k
            // 41 = ?D?, loop over diagonals of m2
            //
            // 81 = copy x*m1
            // 82 = copy x*m2
            // 83 = temp m1*m2

            const bool xcc = M2::_colmajor && M3::_colmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor && 
                (M3::_colmajor || M3::_rowmajor);
            const bool xdx = M2::_diagmajor && (
                (M1::_colmajor && M3::_colmajor) ||
                (M1::_rowmajor && M3::_rowmajor) ||
                (M1::_diagmajor && M3::_diagmajor) );
            const bool dxx = M1::_diagmajor && (
                (M2::_colmajor && M3::_colmajor) ||
                (M2::_rowmajor && M3::_rowmajor) ||
                (M2::_diagmajor && M3::_diagmajor) );

            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 202 :
                rs == 1 ? 203 :
                !(M1::_rowmajor || M1::_colmajor || M1::_diagmajor) ? 81 :
                !(M2::_rowmajor || M2::_colmajor || M2::_diagmajor) ? 82 :
                !(M3::_rowmajor || M3::_colmajor || M3::_diagmajor) ? 83 :
                M3::_rowmajor && !M3::_colmajor ? 5 :
                xcc ? 11 : crx ? 31 : xdx ? 41 : dxx ? 5 :
                M3::_colmajor ? 82 :
                // m3 must be diagmajor at this point
                M2::_diagmajor ? 81 :
                M1::_diagmajor ? 82 : 
                83;
#ifdef PRINTALGO_BB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            const int K = xs==TMV_UNKNOWN ? int(m1.rowsize()) : xs;
            std::cout<<"InlineMultBB: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"M = "<<M<<"  N = "<<N<<"  K = "<<K<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<"  xs = "<<xs<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_BB
            typedef typename M3::real_type RT;
            typedef typename M3::value_type T3;
            Matrix<T3> m1c = m1;
            Matrix<T3> m2c = m2;
            Matrix<T3> m3i = m3;
            Matrix<T3> m3c = m3;
            MultMM<add>(x,m1c,m2c,m3c);
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
#endif
            MultBB_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
#ifdef XDEBUG_BB
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
            const int lo3 = TMV_MIN(int(m3.colsize())-1 , m1.nlo()+m2.nlo());
            const int hi3 = TMV_MIN(int(m3.rowsize())-1 , m1.nhi()+m2.nhi());
            if (m3.nlo() > lo3)
                Maybe<!add>::zero2(m3.diagRange(-m3.nlo(),-lo3));
            if (m3.nhi() > hi3)
                Maybe<!add>::zero2(m3.diagRange(m3.nhi()+1,hi3+1));
#ifdef XDEBUG_BB
            std::cout<<"lo3, hi3 = "<<lo3<<','<<hi3<<std::endl;
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
#ifdef XDEBUG_BB
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
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<-402,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                (xs == TMV_UNKNOWN || xs > 16) &&
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
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo -402\n";
#endif
            MultBB_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                (xs == TMV_UNKNOWN || xs > 16) &&
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
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo -2\n";
#endif
            MultBB_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBB_Helper<-1,cs,rs,xs,add,ix,T,M1,M2,M3>
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
#ifdef PRINTALGO_BB
            std::cout<<"MultBB algo -1\n";
#endif
            MultBB_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Band_Mutable<M3>& m3)
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
        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultBB_Helper<-1,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Band_Mutable<M3>& m3)
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
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultBB_Helper<-3,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Band_Mutable<M3>& m3)
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
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultBB_Helper<98,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    //
    // M = B * B
    //

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        const int lo = TMV_MIN(int(m3.colsize())-1 , m1.nlo()+m2.nlo());
        const int hi = TMV_MIN(int(m3.rowsize())-1 , m1.nhi()+m2.nhi());
        typename BMVO<M3>::b b3 = BandMatrixViewOf(m3,lo,hi);
        MultMM<add>(x,m1,m2,b3);
        if (lo < int(m3.colsize())-1) {
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m3.cRowRange(lo+1,m3.colsize()),m3.colsize()-lo-2,0));
        }
        if (hi < int(m3.rowsize())-1) {
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m3.cColRange(hi+1,m3.rowsize()),0,m3.rowsize()-hi-2));
        }
    }

    //
    // U = B * B
    //

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert(!M3::_unit);
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        const int lo = TMV_MIN(int(m3.colsize())-1 , m1.nlo()+m2.nlo());
        const int hi = TMV_MIN(int(m3.rowsize())-1 , m1.nhi()+m2.nhi());
        const int lo3 = Maybe<M3::_upper>::select(lo,hi);
        const int hi3 = Maybe<M3::_upper>::select(hi,lo);
        TMVAssert(lo3 == 0);
        typename BMVOTri<M3>::b b3 = BandMatrixViewOf(m3,hi3);
        MultMM<add>(x,m1,m2,b3);
        if (hi3 < int(m3.size())-1)
            Maybe<!add>::zero2(m3.offDiag(hi3+1));
    }

    //
    // D = B * B
    //

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        TMVAssert(m1.nlo() == 0);
        TMVAssert(m1.nhi() == 0);
        TMVAssert(m2.nlo() == 0);
        TMVAssert(m2.nhi() == 0);
        MultMM<add>(
            x,DiagMatrixViewOf(m1.diag()),DiagMatrixViewOf(m2.diag()),m3);
    }

    //
    // B * U, U * B
    //

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
#ifdef PRINTALGO_BB
        std::cout<<"MultMM U*B:\n";
        std::cout<<"add = "<<add<<std::endl;
        std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
        std::cout<<"m1 = "<<m1.mat()<<std::endl;
        std::cout<<"m2 = "<<m2.mat()<<std::endl;
        std::cout<<"m3 = "<<m3.mat()<<std::endl;
#endif
        if (m2.isunit()) {
            // TODO: Worth writing a version that doesn't just copy?
            const int s = ShapeTraits<M2::_shape>::nonunit_shape;
            typedef typename M2::value_type T2;
            typename MCopyHelper<T2,s,M2::_size,M2::_size>::type m2c = m2;
            MultMM<add>(x,m1.mat(),BandMatrixViewOf(m2c),m3.mat());
        } else {
            MultMM<add>(
                x,m1.mat(),BandMatrixViewOf(m2.viewAsNonUnitDiag()),m3.mat()); 
        }
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
#ifdef PRINTALGO_BB
        std::cout<<"MultMM U*B:\n";
        std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
        std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
        std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
        std::cout<<"add = "<<add<<std::endl;
        std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
        std::cout<<"m1 = "<<m1.mat()<<std::endl;
        std::cout<<"m2 = "<<m2.mat()<<std::endl;
        std::cout<<"m3 = "<<m3.mat()<<std::endl;
#endif
        if (m1.isunit()) {
            const int s = ShapeTraits<M1::_shape>::nonunit_shape;
            typedef typename M1::value_type T1;
            typename MCopyHelper<T1,s,M1::_size,M1::_size>::type m1c = m1;
            MultMM<add>(x,BandMatrixViewOf(m1c),m2.mat(),m3.mat());
        } else {
            MultMM<add>(
                x,BandMatrixViewOf(m1.viewAsNonUnitDiag()),m2.mat(),m3.mat()); 
        }
    }

} // namespace tmv

#endif 
