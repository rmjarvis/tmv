

#ifndef TMV_UnpackQ_H
#define TMV_UnpackQ_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_MultMV.h"
#include "TMV_MultMM.h"
#include "TMV_DivVU.h"
#include "TMV_DivMU.h"
#include "TMV_Householder.h"

#ifdef PRINTALGO_QR
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_TriMatrixIO.h"
#include "TMV_ProdMM.h"
#include "TMV_MultUL.h"
#endif

// BLOCKSIZE is the block size to use in algo 21
#define TMV_QR_BLOCKSIZE 48

namespace tmv {

    // Defined in TMV_UnpackQ.cpp
    template <class T, class RT>
    void InstUnpackQ(MatrixView<T> Q, const ConstVectorView<RT>& beta);

    template <int algo, int cs, int rs, class M, class V>
    struct UnpackQ_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or 1, or N == 0)
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<0,cs,rs,M,V>
    { static TMV_INLINE void call(M& A, const V& beta) {} };

    // algo 11: Non-block algorithm, loop over n
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<11,cs,rs,M1,V>
    {
        static void call(M1& Q, const V& beta)
        {
            // This is essentially the reverse of the QR Decomposition
            // algorithm where R is taken to be an identity matrix.
            
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? beta.size() : rs;
            const int K = Q.rowsize();
#ifdef PRINTALGO_QR
            std::cout<<"UnpackQ algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename V::const_reverse_iterator IT;
            typedef typename M1::reference Mref;
            typedef typename M1::col_sub_type Mc;
            typedef typename M1::row_sub_type Mr;
            typedef typename M1::submatrix_type Ms;

            typedef Vector<T,NoAlias> V2;
            typedef typename V2::subvector_type V2s;
            V2 tempBase(K);

            if (N > 1) Q.colRange(0,N).upperTri().offDiag().setZero();

            IT bj = beta.rbegin(); // rbegin, so iterate from end.
            for (int j=N-1;j>=0;--j,++bj) {
                // Work on the lower part of column j
                Mref Qjj = Q.ref(j,j);
                Mc u = Q.col(j,j+1,M);
                Mr Q1a = Q.row(j,j+1,K);
                Ms Q1b = Q.subMatrix(j+1,M,j+1,K);
                V2s temp = tempBase.subVector(0,K-j-1);
                // Reflect the rest of the matrix to the right of this column.
                HouseholderMultEq(u,*bj,Q1a,Q1b,temp);
                // Unpack this column into H * ej
                HouseholderUnpack(u,*bj,Qjj);
            }
        }
    };

    // algo 21: Block algorithm
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<21,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        static void call(M1& Q, const V& beta)
        {
            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? beta.size() : rs;
            const int K = Q.rowsize();
#ifdef PRINTALGO_QR
            std::cout<<"UnpackQ algo 21: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const int Nx = TMV_QR_BLOCKSIZE;
            const int s1 = IntTraits2<Nx,rs>::min;
            const int N1 = TMV_MIN(Nx,N);
            typedef typename MCopyHelper<T,UpperTri,s1,s1>::type Ztype;
            Ztype BaseZ = MatrixSizer<T>(N1,N1);

            typedef typename Ztype::subtrimatrix_type Zs;
            typedef typename M1::submatrix_type M1s;
            typedef typename V::const_subvector_type Vs;

            typedef Matrix<T,NoDivider|NoAlias> M2;
            typedef typename M2::submatrix_type M2s;
            typedef typename M2s::uppertri_type M2t;
            M2 tempBase = MatrixSizer<T>(N1,K);

            Q.colRange(0,N).upperTri().setZero();
            for(int j2=N;j2>0;) {
                int j1 = j2 > Nx ? j2-Nx : 0;
                M1s Y = Q.subMatrix(j1,M,j1,j2);
                Zs Z = BaseZ.subTriMatrix(0,j2-j1);
                Vs b1 = beta.subVector(j1,j2);
                M1s Q2 = Q.subMatrix(j1,M,j2,K);
                M2s temp1 = tempBase.subMatrix(0,j2-j1,0,K-j2);
                M2t temp2 = tempBase.subMatrix(0,j2-j1,0,j2-j1).upperTri();

                BlockHouseholderMakeZ(Y,Z,b1);
                BlockHouseholderLMult(Y,Z,Q2,temp1);
                BlockHouseholderUnpack(Y,Z,temp2);
                j2 = j1;
            }
        }
    };

    // algo 27: Block algorithm -- single block
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<27,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        static void call(M1& Q, const V& beta)
        {
            const int N = rs==Unknown ? beta.size() : rs;
            const int K = Q.rowsize();
#ifdef PRINTALGO_QR
            const int M = cs==Unknown ? Q.colsize() : cs;
            std::cout<<"UnpackQ algo 27: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename MCopyHelper<T,UpperTri,rs,rs>::type Ztype;
            Ztype Z = MatrixSizer<T>(N,N);

            typedef typename M1::colrange_type M1c;

            const int xx = Unknown;
            typedef typename MCopyHelper<T,Rec,rs,xx>::type M2;
            typedef typename M2::colrange_type M2c;
            typedef typename M2c::uppertri_type M2t;
            M2 tempBase(N,TMV_MAX(K-N,N));
            M2t temp2 = tempBase.rowRange(0,N).upperTri();

            M1c Q1 = Q.colRange(0,N);

            Q1.upperTri().setZero();
            BlockHouseholderMakeZ(Q1,Z,beta);
            if (K > N) {
                M1c Q2 = Q.colRange(N,K);
                M2c temp1 = tempBase.colRange(0,K-N);
                BlockHouseholderLMult(Q1,Z,Q2,temp1);
            }
            BlockHouseholderUnpack(Q1,Z,temp2);
        }
    };

    // algo 31: Decide which algorithm to use from runtime size
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<31,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        static void call(M1& Q, const V& beta)
        {
            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? beta.size() : rs;
#ifdef PRINTALGO_QR
            std::cout<<"UnpackQ algo 31: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const int algo27 =
                (rs == Unknown || rs <= 128) ? 27 : 0;
            const int algo21 =
                (rs == Unknown || rs > 128) ? 21 : 0;
            const int l2cache = TMV_L2_CACHE*1024/sizeof(T);

            if (M*N <= l2cache)
                UnpackQ_Helper<11,cs,rs,M1,V>::call(Q,beta);
            else if (N <= 128)
                UnpackQ_Helper<algo27,cs,rs,M1,V>::call(Q,beta);
            else
                UnpackQ_Helper<algo21,cs,rs,M1,V>::call(Q,beta);
        }
    };

    // algo 81: Copy to colmajor
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<81,cs,rs,M,V>
    {
        static inline void call(M& Q, const V& beta)
        {
#ifdef PRINTALGO_QR
            std::cout<<"UnpackQ algo 81: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M::value_type T;
            typedef typename MCopyHelper<T,Rec,cs,rs>::type Mcm;
            Mcm Qcm = Q;
            UnpackQ_Helper<-2,cs,rs,Mcm,V>::call(Qcm,beta);
            typename M::noalias_type Qna = Q.noAlias();
            Copy(Qcm,Qna);
        }
    };

    // algo 90: call InstUnpackQ
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<90,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& Q, const V& beta)
        { InstUnpackQ(Q.xView(),beta.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<97,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& Q, const V& beta)
        {
            typedef typename M::conjugate_type Mc;
            Mc Qc = Q.conjugate();
            UnpackQ_Helper<-2,cs,rs,Mc,V>::call(Qc,beta);
        }
    };

    // algo -4: No copies or branches
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<-4,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        static TMV_INLINE void call(M1& Q, const V& beta)
        {
            typedef typename M1::value_type T;
            const int csrs = IntTraits2<cs,rs>::prod;
            const int l2cache = TMV_L2_CACHE*1024/sizeof(T);
            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                TMV_OPT == 0 ? 11 :
                rs == Unknown ? 31 :
                cs == Unknown ? 31 :
                csrs <= l2cache ? 11 :
                rs <= 128 ? 27 : 21;
#ifdef PRINTALGO_QR
            std::cout<<"Inline UnpackQ: \n";
            std::cout<<"Q = "<<TMV_Text(Q)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<Q.colsize()<<"  "<<beta.size()<<std::endl;
            std::cout<<"csrs = "<<csrs<<std::endl;
            std::cout<<"l2cache = "<<l2cache<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            UnpackQ_Helper<algo,cs,rs,M1,V>::call(Q,beta);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<-3,cs,rs,M1,V>
    {
        static TMV_INLINE void call(M1& Q, const V& beta)
        {
            const int algo = (
                ( cs != Unknown && rs != Unknown &&
                  cs <= 16 && rs <= 16 ) ? -4 :
                ( TMV_OPT >= 2 && !M1::_colmajor ) ? 81 :
                -4 );
#ifdef PRINTALGO_QR
            const int M = cs==Unknown ? Q.colsize() : cs;
            const int N = rs==Unknown ? beta.size() : rs;
            std::cout<<"UnpackQ algo -3: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            UnpackQ_Helper<algo,cs,rs,M1,V>::call(Q,beta);
#ifdef PRINTALGO_QR
            std::cout<<"Done UnpackQ\n";
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<-2,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& Q, const V& beta)
        {
            typedef typename M::value_type T;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                M::_conj ? 97 :
                inst ? 90 :
                -3;
            UnpackQ_Helper<algo,cs,rs,M,V>::call(Q,beta);
        }
    };

    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<-1,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& Q, const V& beta)
        { UnpackQ_Helper<-2,cs,rs,M,V>::call(Q,beta); }
    };

    template <class M, class V>
    inline void InlineUnpackQ(
        BaseMatrix_Rec_Mutable<M>& Q, const BaseVector_Calc<V>& beta)
    {
        TMVStaticAssert(V::isreal);
        TMVStaticAssert((Traits2<
                         typename M::value_type,
                         typename V::value_type>::samebase));
        TMVAssert(Q.colsize() >= beta.size());
        TMVAssert(Q.rowsize() == beta.size() ||
                  Q.rowsize() == Q.colsize());

        const int cs = M::_colsize;
        const int rs = V::_size;
        typedef typename M::cview_type Mv;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_REF(M,Mv) Qv = Q.cView();
        TMV_MAYBE_CREF(V,Vv) betav = beta.cView();
        UnpackQ_Helper<-3,cs,rs,Mv,Vv>::call(Qv,betav);
    }

    template <class M, class V>
    inline void UnpackQ(
        BaseMatrix_Rec_Mutable<M>& Q, const BaseVector_Calc<V>& beta)
    {
        TMVStaticAssert(V::isreal);
        TMVStaticAssert((Traits2<
                         typename M::value_type,
                         typename V::value_type>::samebase));
        TMVAssert(Q.colsize() >= beta.size());
        TMVAssert(Q.rowsize() == beta.size() ||
                  Q.rowsize() == Q.colsize());

        const int cs = M::_colsize;
        const int rs = V::_size;
        typedef typename M::cview_type Mv;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_REF(M,Mv) Qv = Q.cView();
        TMV_MAYBE_CREF(V,Vv) betav = beta.cView();
        UnpackQ_Helper<-2,cs,rs,Mv,Vv>::call(Qv,betav);
    }

} // namespace tmv

#undef TMV_QR_BLOCKSIZE

#endif

