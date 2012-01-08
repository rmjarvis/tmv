
//---------------------------------------------------------------------------
//
// This file defines the PackedQ class.
//
// A PackedQ object represents a unitary matrix that is the 
// result of a QR (or similar) decomposition.  
// The QR decomposition is able to store almost all the information about
// the Q matrix in the lower trapezoid part of the matrix, with R being 
// stored in the upper triangle part.  The only extra information needed
// is a vector of beta values.
//
// Basically each column of the matrix and the corresponding beta value
// define a Householder reflection.  The full Q matrix is the product 
// of these reflections.  The main use of a PackedQ object is to multiply
// Q by some other vector or matrix.  This is implemented by multiplying
// each Householder matrix in turn.  This is more efficient than directly
// constructing the full Q matrix and then multiplying by that.
//
// The constructor takes as arguments a reference to the matrix
// where the columns are stored and a reference to the beta vector.
// The PackedQ object stores views of these.  So if they go out of scope,
// then the PackedQ object will be invalid.  
//
// The PackedQ class is a template with two template parameters.
// These parameters are the types of the Q and beta objects which
// described above.  These templates are called M and V below.
//
// Constructors:
//
//    PackedQ(const M& Q, const V& beta)
//        Create a PackedQ matrix using the values in Q and beta.
//
// Access Functions:
//
//     int colsize() const
//     int rowsize() const
//         Return the size of the matrix.
//     
//     const M& getQ() const
//     const V& getBeta() const
//         Return the component elements
//
//     int operator()(int i, int j) const
//         This is required for BaseMatrix, but it's extremely inefficient
//         for PackedQ.  (You basically need to construct the full Q
//         matrix up to the j column.)
//         So this function is private to make it a compiler error if anyone
//         tries to use it.
//
// Functions:
//
//     int det() const
//         Returns the determinant of Q.  This is always either 1 or -1.
//     int logDet(int* sign=0) const
//         Returns 0.  If requested, the sign is either 1 or -1.
//
//
// Operators:
//
//     q*v
//     v*q
//     v/q
//     v%q
//     v *= q
//     v /= q
//     v %= q
//
//     q*m
//     m*q
//     m/q
//     m%q
//     m *= q
//     m /= q
//     m %= q
//     
//     These are the main reason for having a PackedQ class.
//     The multiplication and division operations can be performed 
//     without having to construct the actual Q matrix.
//     This is done by successively applying the Householder reflections
//     to the vector or matrix.
//

#ifndef TMV_PackedQ_H
#define TMV_PackedQ_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseVector.h"
#include "TMV_UnpackQ.h"
#include "TMV_Householder.h"
#include "TMV_SmallMatrix.h"

//#define PRINTALGO_QR

#ifdef PRINTALGO_QR
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#endif

// BLOCKSIZE is the block size to use in algo 21
#define TMV_QR_BLOCKSIZE 48

namespace tmv {

    // 
    // First the basic functionality of doing the PackedQ multiplication
    // or division without having to unpack the full Q matrix.
    //

    // Defined in TMV_PackedQ.cpp
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_MultEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        MatrixView<T2> m2);
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_LDivEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        MatrixView<T2> m2);
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_MultEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        VectorView<T2> v2);
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_LDivEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        VectorView<T2> v2);

    // cs,rs refer to M1.  xs is the other dimension of M2
    template <int algo, bool div, int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0 or K == 0)
    template <bool div, int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<0,div,cs,rs,xs,M1,V1,M2>
    { static TMV_INLINE void call(const M1& , const V1& , M2& ) {} };

    // algo 11: Normal case
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<11,false,cs,rs,xs,M1,V1,M2>
    {
        static void call(const M1& Q, const V1& beta, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 11: div,cs,rs,xs = "<<
                false<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"Q = "<<Q<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //typename M1::copy_type QQ = Q;
            //UnpackQ(QQ,beta);
            //std::cout<<"QQ = "<<QQ<<std::endl;
            //std::cout<<"QQ * m2 = "<<QQ*m2<<std::endl;
#endif
            typedef typename M1::real_type RT;

            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? Q.rowsize() : rs;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::row_type M2r;
            typedef typename M2::rowrange_type M2rr;

            typedef typename M2::value_type T2;
            typedef typename VCopyHelper<T2,xs>::type V3;
            V3 temp = VectorSizer<T2>(m2.rowsize());

            for(int j=N-1;j>=0;--j) if (beta(j) != RT(0)) {
                M1c u = Q.col(j,j+1,M);
                M2r m2a = m2.row(j);
                M2rr m2b = m2.rowRange(j+1,M);
#ifdef PRINTALGO_QR
                //std::cout<<"u = "<<u<<std::endl;
                //std::cout<<"beta = "<<beta(j)<<std::endl;
                //std::cout<<"m2a = "<<m2a<<std::endl;
                //std::cout<<"m2b = "<<m2b<<std::endl;
#endif
                HouseholderMultEq(u,beta(j),m2a,m2b,temp);
#ifdef PRINTALGO_QR
                //std::cout<<"m2 => "<<m2<<std::endl;
#endif
            }
        }
    };
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<11,true,cs,rs,xs,M1,V1,M2>
    {
        static void call(const M1& Q, const V1& beta, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 11: div,cs,rs,xs = "<<
                true<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"Q = "<<Q<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //typename M1::copy_type QQ = Q;
            //UnpackQ(QQ,beta);
            //std::cout<<"QQ = "<<QQ<<std::endl;
            //std::cout<<"QQt * m2 = "<<QQ.transpose()*m2<<std::endl;
#endif
            typedef typename M1::real_type RT;

            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? Q.rowsize() : rs;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::row_type M2r;
            typedef typename M2::rowrange_type M2rr;

            typedef typename M2::value_type T2;
            typedef typename VCopyHelper<T2,xs>::type V3;
            V3 temp = VectorSizer<T2>(m2.rowsize());

            for(int j=0;j<N;++j) if (beta(j) != RT(0)) {
                M1c u = Q.col(j,j+1,M);
                M2r m2a = m2.row(j);
                M2rr m2b = m2.rowRange(j+1,M);
#ifdef PRINTALGO_QR
                //std::cout<<"u = "<<u<<std::endl;
                //std::cout<<"beta = "<<beta(j)<<std::endl;
                //std::cout<<"m2a = "<<m2a<<std::endl;
                //std::cout<<"m2b = "<<m2b<<std::endl;
#endif
                HouseholderMultEq(u,beta(j),m2a,m2b,temp);
#ifdef PRINTALGO_QR
                //std::cout<<"m2 => "<<m2<<std::endl;
#endif
            }
        }
    };

    // algo 13: Construct Q directly and multiply.
    template <bool div, int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<13,div,cs,rs,xs,M1,V1,M2>
    {
        static void call(const M1& Q, const V1& beta, M2& m2)
        {
            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? Q.rowsize() : rs;
#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 13: div,cs,rs,xs = "<<
                true<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"Q = "<<Q<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            typedef typename M1::real_type RT;
            typedef typename M1::value_type T1;

            typedef typename MCopyHelper<T1,Rec,cs,cs>::type M1f;
            M1f Qfull = MatrixSizer<T1>(M,M);
            Qfull.colRange(0,N) = Q;
            Qfull.colRange(N,M).setZero();

            UnpackQ(Qfull,beta);

            typename M2::copy_type m2c = m2;
            m2 = Maybe<div>::transpose(Qfull) * m2c;
        }
    };

    // algo 21: Block algorithm
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<21,false,cs,rs,xs,M1,V1,M2>
    {
        static void call(const M1& Q, const V1& beta, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 21: div,cs,rs,xs = "<<
                false<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"Q = "<<Q<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //typedef typename M2::copy_type M2c;
            //M2c m2c = m2;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;

            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? Q.rowsize() : rs;
            const int Nb = TMV_QR_BLOCKSIZE;
            const int s1 = IntTraits2<Nb,rs>::min;
            const int N1 = TMV_MIN(Nb,N);
            typedef typename MCopyHelper<T1,UpperTri,s1,s1>::type Ztype;
            Ztype BaseZ(MatrixSizer<T1>(N1,N1));

            typedef typename MCopyHelper<T2,Rec,s1,M2::_rowsize>::type M3;
            typedef typename M3::rowrange_type M3r;
            M3 tempBase = MatrixSizer<T2>(N1,m2.rowsize());

            typedef typename Ztype::subtrimatrix_type Zs;
            typedef typename M1::const_submatrix_type M1s;
            typedef typename V1::const_subvector_type V1s;
            typedef typename M2::rowrange_type M2r;

            for(int j2=N;j2>0;) {
                int j1 = j2 > Nb ? j2-Nb : 0;
                M1s Y = Q.subMatrix(j1,M,j1,j2);
                Zs Z = BaseZ.subTriMatrix(0,j2-j1);
                V1s b1 = beta.subVector(j1,j2);
                M2r m2r = m2.rowRange(j1,M);
                M3r temp = tempBase.rowRange(0,j2-j1);

                BlockHouseholderMakeZ(Y,Z,b1);
                BlockHouseholderLMult(Y,Z,m2r,temp);
                j2 = j1;
            }
#ifdef PRINTALGO_QR
            //std::cout<<"m2 => "<<m2<<std::endl;
            //PackedQ_MultEq_Helper<11,false,cs,rs,xs,M1,V1,M2c>::call(
                //Q,beta,m2c);
            //std::cout<<"Correct m2 = "<<m2c<<std::endl;
#endif
        }
    };
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<21,true,cs,rs,xs,M1,V1,M2>
    {
        static void call(const M1& Q, const V1& beta, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 21: div,cs,rs,xs = "<<
                true<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"Q = "<<Q<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;

            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? Q.rowsize() : rs;
            const int Nb = TMV_QR_BLOCKSIZE;
            const int s1 = IntTraits2<Nb,rs>::min;
            const int N1 = TMV_MIN(Nb,N);
            typedef typename MCopyHelper<T1,UpperTri,s1,s1>::type Ztype;
            Ztype BaseZ(MatrixSizer<T1>(N1,N1));

            typedef typename MCopyHelper<T2,Rec,s1,M2::_rowsize>::type M3;
            typedef typename M3::rowrange_type M3r;
            M3 tempBase = MatrixSizer<T2>(N1,m2.rowsize());

            typedef typename Ztype::subtrimatrix_type Zs;
            typedef typename M1::const_submatrix_type M1s;
            typedef typename V1::const_subvector_type V1s;
            typedef typename M2::rowrange_type M2r;

            for(int j1=0;j1<N;) {
                int j2 = TMV_MIN(j1+Nb,N);
                M1s Y = Q.subMatrix(j1,M,j1,j2);
                Zs Z = BaseZ.subTriMatrix(0,j2-j1);
                V1s b1 = beta.subVector(j1,j2);
                M2r m2r = m2.rowRange(j1,M);
                M3r temp = tempBase.rowRange(0,j2-j1);

                BlockHouseholderMakeZ(Y,Z,b1);
                BlockHouseholderLDiv(Y,Z,m2r,temp);
                j1 = j2;
            }
        }
    };

    // algo 27: Block algorithm -- single block
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<27,false,cs,rs,xs,M1,V1,M2>
    {
        static void call(const M1& Q, const V1& beta, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 27: div,cs,rs,xs = "<<
                false<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"Q = "<<Q<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;

            const int N = rs==Unknown ? Q.rowsize() : rs;

            typedef typename MCopyHelper<T1,UpperTri,rs,rs>::type Ztype;
            Ztype Z(MatrixSizer<T1>(N,N));

            typedef typename MCopyHelper<T2,Rec,rs,M2::_rowsize>::type M3;
            M3 temp = MatrixSizer<T2>(N,m2.rowsize());

            BlockHouseholderMakeZ(Q,Z,beta);
            BlockHouseholderLMult(Q,Z,m2,temp);
        }
    };
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<27,true,cs,rs,xs,M1,V1,M2>
    {
        static void call(const M1& Q, const V1& beta, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 27: div,cs,rs,xs = "<<
                true<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"Q = "<<Q<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;

            const int N = rs==Unknown ? Q.rowsize() : rs;

            typedef typename MCopyHelper<T1,UpperTri,rs,rs>::type Ztype;
            Ztype Z(MatrixSizer<T1>(N,N));

            typedef typename MCopyHelper<T2,Rec,rs,xs>::type M3;
            M3 temp = MatrixSizer<T2>(N,m2.rowsize());

            BlockHouseholderMakeZ(Q,Z,beta);
            BlockHouseholderLDiv(Q,Z,m2,temp);
        }
    };

    // algo 31: Decide which algorithm to use from runtime size
    template <bool div, int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<31,div,cs,rs,xs,M1,V1,M2>
    {
        static void call(const M1& Q, const V1& beta, M2& m2)
        {
            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? Q.rowsize() : rs;
            const int K = xs==Unknown ? m2.rowsize() : xs;
            typedef typename M1::value_type T;
            const int l2cache = TMV_L2_CACHE*1024/sizeof(T);
            const int csrs = IntTraits2<cs,rs>::prod;
            const int algo27 =
                csrs != Unknown && csrs <= l2cache ? 0 :
                (rs == Unknown || rs <= 128) ? 27 : 0;
            const int algo21 =
                csrs != Unknown && csrs <= l2cache ? 0 :
                (rs == Unknown || rs > 128) ? 21 : 0;
#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 31: div,cs,rs,xs = "<<
                div<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"l2cache = "<<l2cache<<std::endl;
            std::cout<<"csrs = "<<csrs<<std::endl;
            std::cout<<"algo27 = "<<algo27<<std::endl;
            std::cout<<"algo21 = "<<algo21<<std::endl;
            std::cout<<"M*N <= l2cache ? "<<(M*N <= l2cache)<<std::endl;
            std::cout<<"use 11 ? "<<(
                M < 16 ||
                2*M*N*(K-M) < N*N*(K-2*M+2*N/3) + (
                    M < 32 ? M*M*K*3/4 : M*M*K/2 ) ) <<std::endl;
#endif

            // algo 13 uses more floating point operations than algo 11.
            // However, algo 11 is all Level 2 calculations, while algo 13
            // ends with a direct matrix-matrix product, so can use Level 3
            // operations.  This can make it faster for some matrix sizes.
            // We calculate the total ops for algo 11 and 13, but we discount
            // the level 3 ops by 25% when M is < 32 and 50% for higher M's.
            // If M < 16, level 3 isn't really any faster, so just use algo 11.
            //
            // Algo 11: nops = (2M-N)NK
            //
            // Also 13: nops = 2(M-N)MN + 2N^3/3 + x M^2K
            // x is our discount factor.  0.75 or 0.5.
            //
            // nops_11 <? nops_13
            // 2MNK - N^2K <? 2M^2N - 2MN^2 + 2N^3/3 + xM^2K
            //
            // 2MN(K-M) <? N^2*(K-2M+2N/3) + xM^2K
            
            if (M*N <= l2cache)
                if (M < 16 ||
                    2*M*N*(K-M) < N*N*(K-2*M+2*N/3) + (
                        M < 32 ? M*M*K*3/4 : M*M*K/2 ) ) 
                    PackedQ_MultEq_Helper<11,div,cs,rs,xs,M1,V1,M2>::call(
                        Q,beta,m2);
                else
                    PackedQ_MultEq_Helper<13,div,cs,rs,xs,M1,V1,M2>::call(
                        Q,beta,m2);
            else if (N <= 128)
                PackedQ_MultEq_Helper<algo27,div,cs,rs,xs,M1,V1,M2>::call(
                    Q,beta,m2);
            else
                PackedQ_MultEq_Helper<algo21,div,cs,rs,xs,M1,V1,M2>::call(
                    Q,beta,m2);
        }
    };

    // algo 32: Decide which algorithm to use from runtime size,
    // all sizes known, csrs <= l2cache, cs >= 16
    template <bool div, int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<32,div,cs,rs,xs,M1,V1,M2>
    {
        static void call(const M1& Q, const V1& beta, M2& m2)
        {
            TMVStaticAssert(cs != Unknown);
            TMVStaticAssert(rs != Unknown);
            TMVStaticAssert(xs != Unknown);
            TMVStaticAssert(cs >= 16);

#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 32: div,cs,rs,xs = "<<
                div<<','<<cs<<','<<rs<<','<<xs<<std::endl;
#endif
            typedef typename M1::value_type T;
            const int algo1 = 
                2*cs*rs*(xs-cs) < rs*rs*(xs-2*cs+2*rs/3) + (
                    cs < 32 ? cs*cs*xs*3/4 : cs*cs*xs/2 ) ? 11 : 13;
            PackedQ_MultEq_Helper<algo1,div,cs,rs,xs,M1,V1,M2>::call(
                Q,beta,m2);
        }
    };

    // algo 90: call InstPackedQ_MultEq
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<90,false,cs,rs,xs,M1,V1,M2>
    {
        static TMV_INLINE void call(const M1& Q, const V1& beta, M2& m2)
        { InstPackedQ_MultEq(Q.xView(),beta.xView(),m2.xView()); }
    };
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<90,true,cs,rs,xs,M1,V1,M2>
    {
        static TMV_INLINE void call(const M1& Q, const V1& beta, M2& m2)
        { InstPackedQ_LDivEq(Q.xView(),beta.xView(),m2.xView()); }
    };

    // algo 91: call InstPackedQ_MultEq -- M2 is a vector
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<91,false,cs,rs,xs,M1,V1,M2>
    {
        static TMV_INLINE void call(const M1& Q, const V1& beta, M2& m2)
        { InstPackedQ_MultEq(Q.xView(),beta.xView(),m2.col(0).xView()); }
    };
    template <int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<91,true,cs,rs,xs,M1,V1,M2>
    {
        static TMV_INLINE void call(const M1& Q, const V1& beta, M2& m2)
        { InstPackedQ_LDivEq(Q.xView(),beta.xView(),m2.col(0).xView()); }
    };

    // algo 97: Conjugate
    template <bool div, int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<97,div,cs,rs,xs,M1,V1,M2>
    {
        static TMV_INLINE void call(const M1& Q, const V1& beta, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c Qc = Q.conjugate();
            M2c m2c = m2.conjugate();
            PackedQ_MultEq_Helper<-2,div,cs,rs,xs,M1c,V1,M2c>::call(Qc,beta,m2c);
        }
    };

    // algo -3: Determine which algorithm to use
    template <bool div, int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<-3,div,cs,rs,xs,M1,V1,M2>
    {
        static TMV_INLINE void call(const M1& Q, const V1& beta, M2& m2)
        {
            typedef typename M2::value_type T;
            const int l2cache = TMV_L2_CACHE*1024/sizeof(T);
            const int csrs = IntTraits2<cs,rs>::prod;
            const int algo =
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                TMV_OPT == 0 ? 11 :
                rs == Unknown ? 31 :
                cs == Unknown ? 31 :
                csrs <= l2cache ? (
                    cs < 16 ? 11 : xs == Unknown ? 31 : 
                    32 ) :
                rs <= 128 ? 27 : 21;
#ifdef PRINTALGO_QR
            std::cout<<"Inline PackedQ_MultEq\n";
            std::cout<<"div,cs,rs,xs = "<<div<<','<<cs<<','<<rs<<std::endl;
            std::cout<<"Q = "<<TMV_Text(Q)<<std::endl;
            std::cout<<"beta = "<<TMV_Text(beta)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            std::cout<<"Q = "<<Q<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
#endif
            PackedQ_MultEq_Helper<algo,div,cs,rs,xs,M1,V1,M2>::call(Q,beta,m2);
#ifdef PRINTALGO_QR
            std::cout<<"m2 => "<<m2<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <bool div, int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<-2,div,cs,rs,xs,M1,V1,M2>
    {
        static TMV_INLINE void call(const M1& Q, const V1& beta, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                (xs == Unknown || xs > 16 || xs == 1) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                cs == 0 || rs == 0 || xs == 0 ? 0 : 
                M2::_conj ? 97 :
                inst ? (xs == 1 ? 91 : 90) :
                -3;
#ifdef PRINTALGO_QR
            std::cout<<"PackedQ_MultEq algo 32: div,cs,rs,xs = "<<
                div<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"inst = "<<inst<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            PackedQ_MultEq_Helper<algo,div,cs,rs,xs,M1,V1,M2>::call(Q,beta,m2);
        }
    };

    template <bool div, int cs, int rs, int xs, class M1, class V1, class M2>
    struct PackedQ_MultEq_Helper<-1,div,cs,rs,xs,M1,V1,M2>
    {
        static TMV_INLINE void call(const M1& Q, const V1& beta, M2& m2)
        { PackedQ_MultEq_Helper<-2,div,cs,rs,xs,M1,V1,M2>::call(Q,beta,m2); }
    };

    template <class M1, class V1, class M2>
    inline void InlinePackedQ_MultEq(
        const BaseMatrix_Rec<M1>& Q, const BaseVector_Calc<V1>& beta,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,V1::_size>::size;
        const int xs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) Qv = Q.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        PackedQ_MultEq_Helper<-3,false,cs,rs,xs,M1v,V1,M2v>::call(Qv,betav,m2v);
    }

    template <class M1, class V1, class M2>
    inline void PackedQ_MultEq(
        const BaseMatrix_Rec<M1>& Q, const BaseVector_Calc<V1>& beta,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,V1::_size>::size;
        const int xs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) Qv = Q.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        PackedQ_MultEq_Helper<-2,false,cs,rs,xs,M1v,V1,M2v>::call(Qv,betav,m2v);
    }

    template <class M1, class V1, class V2>
    inline void InlinePackedQ_MultEq(
        const BaseMatrix_Rec<M1>& Q, const BaseVector_Calc<V1>& beta,
        BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V2::_size,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == v2.size());
        const int cs = Sizes<V2::_size,M1::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,V1::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::value_type T2;
        const int xs = 1;
        const int s2 = V2::_step;
        const int xx = Unknown;
        const int c = V2::_conj ? Conj : NonConj;
        typedef typename MViewHelper<T2,Rec,cs,xs,s2,xx,c>::type V2v;
        TMV_MAYBE_CREF(M1,M1v) Qv = Q.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        V2v v2v = ColVectorViewOf(v2);
        PackedQ_MultEq_Helper<-3,false,cs,rs,xs,M1v,V1,V2v>::call(Qv,betav,v2v);
    }

    template <class M1, class V1, class V2>
    inline void PackedQ_MultEq(
        const BaseMatrix_Rec<M1>& Q, const BaseVector_Calc<V1>& beta,
        BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V2::_size,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == v2.size());
        const int cs = Sizes<V2::_size,M1::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,V1::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::value_type T2;
        const int xs = 1;
        const int s2 = V2::_step;
        const int xx = Unknown;
        const int c = V2::_conj ? Conj : NonConj;
        typedef typename MViewHelper<T2,Rec,cs,xs,s2,xx,c>::type V2v;
        TMV_MAYBE_CREF(M1,M1v) Qv = Q.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        V2v v2v = ColVectorViewOf(v2);
        PackedQ_MultEq_Helper<-2,false,cs,rs,xs,M1v,V1,V2v>::call(Qv,betav,v2v);
    }

    template <class M1, class V1, class M2>
    inline void InlinePackedQ_LDivEq(
        const BaseMatrix_Rec<M1>& Q, const BaseVector_Calc<V1>& beta,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,V1::_size>::size;
        const int xs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) Qv = Q.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        PackedQ_MultEq_Helper<-3,true,cs,rs,xs,M1v,V1,M2v>::call(Qv,betav,m2v);
    }

    template <class M1, class V1, class M2>
    inline void PackedQ_LDivEq(
        const BaseMatrix_Rec<M1>& Q, const BaseVector_Calc<V1>& beta,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,V1::_size>::size;
        const int xs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) Qv = Q.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        PackedQ_MultEq_Helper<-2,true,cs,rs,xs,M1v,V1,M2v>::call(Qv,betav,m2v);
    }

    template <class M1, class V1, class V2>
    inline void InlinePackedQ_LDivEq(
        const BaseMatrix_Rec<M1>& Q, const BaseVector_Calc<V1>& beta,
        BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V2::_size,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == v2.size());
        const int cs = Sizes<V2::_size,M1::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,V1::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::value_type T2;
        const int xs = 1;
        const int s2 = V2::_step;
        const int xx = Unknown;
        const int c = V2::_conj ? Conj : NonConj;
        typedef typename MViewHelper<T2,Rec,cs,xs,s2,xx,c>::type V2v;
        TMV_MAYBE_CREF(M1,M1v) Qv = Q.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        V2v v2v = ColVectorViewOf(v2);
        PackedQ_MultEq_Helper<-3,true,cs,rs,xs,M1v,V1,V2v>::call(Qv,betav,v2v);
    }

    template <class M1, class V1, class V2>
    inline void PackedQ_LDivEq(
        const BaseMatrix_Rec<M1>& Q, const BaseVector_Calc<V1>& beta,
        BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V2::_size,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == v2.size());
        const int cs = Sizes<V2::_size,M1::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,V1::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::value_type T2;
        const int xs = 1;
        const int s2 = V2::_step;
        const int xx = Unknown;
        const int c = V2::_conj ? Conj : NonConj;
        typedef typename MViewHelper<T2,Rec,cs,xs,s2,xx,c>::type V2v;
        TMV_MAYBE_CREF(M1,M1v) Qv = Q.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        V2v v2v = ColVectorViewOf(v2);
        PackedQ_MultEq_Helper<-2,true,cs,rs,xs,M1v,V1,V2v>::call(Qv,betav,v2v);
    }


    // 
    // PackedQ class 
    //

    template <class M, class V> 
    class PackedQ;

    template <class M, class V> 
    struct Traits<PackedQ<M,V> >
    {
        typedef typename M::value_type value_type;
        typedef typename M::real_type real_type;
        typedef typename M::complex_type complex_type;

        enum { isreal = Traits<M>::isreal };
        enum { iscomplex = Traits<M>::iscomplex };

        typedef PackedQ<M,V> type;
        typedef Matrix<value_type> copy_type;
        typedef copy_type calc_type;
        typedef copy_type eval_type;
        typedef InvalidType inverse_type;

        enum { _colsize = M::_colsize };
        enum { _rowsize = M::_rowsize };
        enum { _nlo = IntTraits2<IntTraits<_colsize>::Sm1,0>::max };
        enum { _nhi = IntTraits2<IntTraits<_rowsize>::Sm1,0>::max };
        enum { _shape = Rec };
        enum { _fort = M::_fort };
        enum { _calc = false };
    };

    template <class M, class V> 
    class PackedQ : 
        public BaseMatrix<PackedQ<M,V> >
    {
    public :

        typedef typename M::value_type T;
        typedef typename M::real_type RT;

        PackedQ(const BaseMatrix_Rec<M>& _Q, const BaseVector<V>& _beta) :
            Q(_Q.mat()), beta(_beta.vec()) 
        {
            TMVStaticAssert((Traits2<typename V::value_type,RT>::sametype));
            TMVAssert(beta.size() == Q.rowsize()); 
        }

        PackedQ(const PackedQ<M,V>& rhs) : Q(rhs.Q), beta(rhs.beta) {}

        ~PackedQ() {}

        //
        // Accesss
        //

        TMV_INLINE const M& getQ() const { return Q; }
        TMV_INLINE const V& getBeta() const { return beta; }

        int det() 
        { 
            int d=1;
            const int n = beta.size();
            for(int i=0; i<n; ++i) if (beta!=RT(0)) d = -d;
            return d;
        }

        int logDet(int* sign=0) const
        { if (sign) *sign = det(); return 0; }


        //
        // Create matrix version
        //

        template <class M2>
        void assignTo(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            Q.assignTo(m2);
            UnpackQ(m2,beta);
        }

        template <class M2>
        void makeFullQ(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            TMVAssert(m2.colsize() == Q.colsize());
            TMVAssert(m2.rowsize() == Q.colsize());

            m2.colRange(0,Q.rowsize()) = Q;
            UnpackQ(m2,beta);
        }


        // 
        // Auxilliary functions
        //

        TMV_INLINE int colsize() const { return Q.colsize(); }
        TMV_INLINE int rowsize() const { return Q.rowsize(); }
        TMV_INLINE int nlo() const { return TMV_MAX(colsize()-1,0); }
        TMV_INLINE int nhi() const { return TMV_MAX(rowsize()-1,0); }

    private : 

        const M& Q;
        const V& beta;

        // Don't allow op=
        void operator=(const PackedQ<M,V>& rhs);
    };

    template <class M, class V>
    TMV_INLINE int Det(const PackedQ<M,V>& Q)
    { return Q.det(); }
    template <class M, class V>
    TMV_INLINE int LogDet(const PackedQ<M,V>& Q)
    { return Q.logDet(); }


    //
    // v3 = Q * v2
    //

    template <bool add, int ix, class T, class M1, class V1, class V2, class V3>
    inline void MultMV(
        const Scaling<ix,T>& x, const PackedQ<M1,V1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        if (add) {
            v3 += (x*m1*v2).calc();
        } else if (m1.isSquare()) {
            v3 = v2;
            PackedQ_MultEq(m1.getQ(),m1.getBeta(),v3.vec());
            Scale(x,v3);
        } else {
            TMVAssert(v3.size() > v2.size());
            v3.subVector(0,v2.size()) = v2;
            v3.subVector(v2.size(),v3.size()).setZero();
            PackedQ_MultEq(m1.getQ(),m1.getBeta(),v3.vec());
            Scale(x,v3);
        }
    }


    //
    // v3 = v1 * Q
    // v3 = QT * v1
    // v3* = Qt * v1*
    //

    template <bool add, int ix, class T, class V1, class M2, class V2, class V3>
    inline void MultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    {
        typename V3::conjugate_type v3c = v3.conjugate();
        if (add) {
            v3 += (x*v1*m2).calc();
        } else if (m2.isSquare()) {
            v3 = v1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v3c);
            Scale(x,v3);
        } else {
            TMVAssert(v1.size() > v3.size());
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            typename VCopyHelper<T12,V1::_size>::type v1c =
                v1.conjugate();
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v1c);
            v3 = x * v1c.subVector(0,v3.size());
        }
    }


    //
    // v *= Q
    //

    template <int ix, class T, class V1, class M2, class V2>
    inline void MultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const PackedQ<M2,V2>& m2)
    {
        typename V1::conjugate_type v1c = v1.conjugate();
        PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v1c);
        m2.divEq(v1c);
        Scale(x,v1);
    }


    //
    // v / Q
    //

    template <int ix, class T, class V1, class M2, class V2, class V3>
    inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    {
        if (m2.isSquare()) {
            v3 = v1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v3);
            Scale(x,v3);
        } else {
            TMVAssert(v1.size() > v3.size());
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            typename VCopyHelper<T12,V1::_size>::type v1c = v1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v1c);
            v3 = x * v1c.subVector(0,v3.size());
        }
    }


    //
    // v /= Q
    //

    template <class V1, class M2, class V2>
    inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const PackedQ<M2,V2>& m2)
    { PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v1); }


    // 
    // v % Q
    // v3 = v1 * Qt
    // v3 = Q* v1
    // v3* = Q v1*
    //

    template <int ix, class T, class V1, class M2, class V2, class V3>
    inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    {
        typename V3::conjugate_type v3c = v3.conjugate();
        if (m2.isSquare()) {
            v3 = v1;
            PackedQ_MultEq(m2.getQ(),m2.getBeta(),v3c);
            Scale(x,v3);
        } else {
            TMVAssert(v3.size() > v1.size());
            v3.subVector(0,v1.size()) = v1;
            v3.subVector(v1.size(),v3.size()).setZero();
            PackedQ_MultEq(m2.getQ(),m2.getBeta(),v3c);
            Scale(x,v3);
        }
    }


    // 
    // v %= Q
    //

    template <class V1, class M2, class V2>
    inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const PackedQ<M2,V2>& m2)
    {
        typename V1::conjugate_type v1c = v1.conjugate();
        PackedQ_MultEq(m2.getQ(),m2.getBeta(),v1c);
    }



    //
    // m3 = Q * m2
    //

    template <bool add, int ix, class T, class M1, class V1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const PackedQ<M1,V1>& m1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else if (m1.isSquare()) {
            m3 = m2;
            PackedQ_MultEq(m1.getQ(),m1.getBeta(),m3.mat());
            Scale(x,m3);
        } else {
            TMVAssert(m3.size() > m2.size());
            m3.rowRange(0,m2.colsize()) = m2;
            m3.rowRange(m2.colsize(),m3.colsize()).setZero();
            PackedQ_MultEq(m1.getQ(),m1.getBeta(),m3.mat());
            Scale(x,m3);
        }
    }


    //
    // m3 = m1 * Q
    // m3t = Qt * m1t
    //

    template <bool add, int ix, class T, class M1, class M2, class V2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1.mat()*m2.mat()).calc();
        } else if (m2.isSquare()) {
            m3 = m1;
            typename M3::adjoint_type m3a = m3.adjoint();
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m3a);
            Scale(x,m3);
        } else {
            TMVAssert(m1.rowsize() > m3.rowsize());
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            const int cs = M1::_colsize;
            const int rs = M1::_rowsize;
            const int A = M1::_rowmajor ? ColMajor : RowMajor;
            typename MCopyHelper<T12,Rec,rs,cs,A>::type m1a = m1.adjoint();
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m1a);
            m3.adjoint() = x * m1a.rowRange(0,m3.colsize());
        }
    }

    //
    // m *= Q
    //

    template <int ix, class T, class M1, class M2, class V2>
    inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1, const Scaling<ix,T>& x,
        const PackedQ<M2,V2>& m2)
    {
        typename M1::adjoint_type m1a = m1.adjoint();
        PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m1a);
        m2.divEq(m1a);
        Scale(x,m1);
    }


    //
    // m / Q
    //

    template <int ix, class T, class M1, class M2, class V2, class M3>
    inline void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (m2.isSquare()) {
            m3 = m1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m3);
            Scale(x,m3);
        } else {
            TMVAssert(m1.colsize() > m3.colsize());
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            const int cs = M1::_colsize;
            const int rs = M1::_rowsize;
            const int A = M1::_rowmajor ? RowMajor : ColMajor;
            typename MCopyHelper<T12,Rec,cs,rs,A>::type m1c = m1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m1c);
            m3 = x * m1c.rowRange(0,m3.colsize());
        }
    }


    //
    // m /= Q
    //

    template <class M1, class M2, class V2>
    inline void LDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const PackedQ<M2,V2>& m2)
    { PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m1); }


    // 
    // m % Q
    // m3 = m1 Qt
    // m3t = Q m1t
    //

    template <int ix, class T, class M1, class M2, class V2, class M3>
    inline void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::adjoint_type m3a = m3.adjoint();
        if (m2.isSquare()) {
            m3 = m1;
            PackedQ_MultEq(m2.getQ(),m2.getBeta(),m3a);
            Scale(x,m3);
        } else {
            TMVAssert(m3.size() > m1.size());
            m3.rowRange(0,m1.colsize()) = m1;
            m3.rowRange(m1.colsize(),m3.colsize()).setZero();
            PackedQ_MultEq(m2.getQ(),m2.getBeta(),m3a);
            Scale(x,m3);
        }
    }


    // 
    // m %= Q
    //

    template <class M1, class M2, class V2>
    inline void RDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const PackedQ<M2,V2>& m2)
    {
        typename M1::adjoint_type m1a = m1.adjoint();
        PackedQ_MultEq(m2.getQ(),m2.getBeta(),m1a);
    }



    //
    // TMV_Text
    //

    template <class M, class V>
    inline std::string TMV_Text(const PackedQ<M,V>& Q)
    {
        std::ostringstream s;
        s << "PackedQ< "<<TMV_Text(Q.getQ())<< " , ";
        s << TMV_Text(Q.getBeta())<<" >";
        return s.str();
    }


} // namespace mv

#undef TMV_QR_BLOCKSIZE


#endif
