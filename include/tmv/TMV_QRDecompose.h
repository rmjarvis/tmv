

#ifndef TMV_QRDecompose_H
#define TMV_QRDecompose_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Householder.h"
//#include "TMV_MultMM.h"
//#include "TMV_DivVU.h"
//#include "TMV_DivMU.h"

#ifdef PRINTALGO_QR
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_TriMatrixIO.h"
#include "TMV_ProdMM.h"
#include "TMV_MultUL.h"
#include "TMV_PermuteM.h"
#endif

// BLOCKSIZE is the block size to use in algo 21, etc.
#define TMV_QR_BLOCKSIZE 48


namespace tmv {

    // Defined in TMV_QRDecompose.cpp
    template <class T, class RT>
    void InstQR_Decompose(MatrixView<T> m, VectorView<RT> beta);

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, class M, class V>
    struct QRDecompose_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or 1, or N == 0)
    template <ptrdiff_t cs, ptrdiff_t rs, class M, class V>
    struct QRDecompose_Helper<0,cs,rs,M,V>
    { static TMV_INLINE void call(M& A, V& beta) {} };

    // algo 11: Non-block algorithm, loop over n
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<11,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta)
        {
            // We want to decompose the matrix (input as A) into Q * R
            // where Q is unitary and R is upper triangular.
            // Q and R are stored in the same matrix (output of A),
            // with the beta values for the Householder matrices returned
            // in beta.
            
            typedef typename M1::value_type T;

            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename V::iterator IT;
            typedef typename M1::col_sub_type M1c;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;
            typedef typename VCopyHelper<T,rs>::type V2;
            typedef typename V2::subvector_type V2s;
            V2 tempBase = VectorSizer<T>(N);

            IT bj = beta.begin();
            for (ptrdiff_t j=0;j<N;++j,++bj) {
                // Work on the lower part of column j
                M1c u = A.col(j,j+1,M);
                // Compute the Householder reflection of this column
                // (and perform the reflection on this column).
                HouseholderReflect(A.ref(j,j),u,*bj);
                // Reflect the rest of the matrix to the right of this column.
                M1r A1a = A.row(j,j+1,N);
                M1s A1b = A.subMatrix(j+1,M,j+1,N);
                V2s temp = tempBase.subVector(0,N-j-1);
                HouseholderMultEq(u,*bj,A1a,A1b,temp);
            }
        }
    };

    // algo 16: Unrolled 
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<16,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        typedef typename M1::real_type RT;
        enum { AA = (
                (M1::_conj ? Conj : NonConj) |
                (M1::_colmajor ? ColMajor : M1::_rowmajor ? RowMajor : NonMajor)
        ) };
        enum { AS1 = M1::_stepi };
        enum { AS2 = M1::_stepj };
        enum { C = M1::_conj };

        // This prevents large matrices from unrolling a lot of 
        // intermediate MultMM, etc. calls.
        enum { csx = cs == Unknown || cs > 64 ? Unknown : cs };


        template <ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t N1>
        struct Unroll_Helper // N1 > 1
        {
            template <class Vt>
            static inline void unroll(M1& A, V& beta, Vt& tempBase) 
            {
                const ptrdiff_t jmid = IntTraits<IntTraits2<j1,j2>::sum>::halfS;
                Unroll_Helper<j1,jmid,jmid-j1>::unroll(A,beta,tempBase);
                Unroll_Helper<jmid,j2,j2-jmid>::unroll(A,beta,tempBase);
            }
        };

        template <ptrdiff_t j1, ptrdiff_t j2>
        struct Unroll_Helper<j1,j2,1> // N == 1
        {
            template <class Vt>
            static inline void unroll(M1& A, V& beta, Vt& tempBase) 
            {
                TMVStaticAssert(j2 == j1+1);
                const ptrdiff_t Mx = IntTraits2<csx,j2>::diff;
                const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
                typename VViewHelper<T,Mx,AS1,C>::type u = A.col(j1,j2,M);
                RT b0;
                HouseholderReflect(A.ref(j1,j1),u,b0);
                beta.ref(j1) = b0;

                typename VViewHelper<T,rs-j2,AS2,C>::type A1a = A.row(j1,j2,rs);
                typename MViewHelper<T,Rec,Mx,rs-j2,AS1,AS2,C>::type A1b =
                    A.subMatrix(j2,M,j2,rs);
                typename VViewHelper<T,rs-j2,1>::type temp =
                    tempBase.subVector(0,rs-j1-1);
                HouseholderMultEq(u,b0,A1a,A1b,temp);
            }
        };

        static inline void call(M1& A, V& beta)
        {
            TMVStaticAssert(rs != Unknown);
#ifdef PRINTALGO_QR
            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
            std::cout<<"QRDecompose algo 16: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            //typedef typename M1::copy_type M1x;
            //typedef typename V::copy_type Vx;
            //M1x Ac = A;
            //Vx bc = beta;
#endif
            SmallVector<T,rs> tempBase;
            Unroll_Helper<0,rs,rs>::unroll(A,beta,tempBase);
#ifdef PRINTALGO_QR
            //std::cout<<"A -> "<<A<<std::endl;
            //std::cout<<"beta -> "<<beta<<std::endl;
            //QRDecompose_Helper<11,cs,rs,M1x,Vx>::call(Ac,bc);
            //std::cout<<"correct A = "<<Ac<<std::endl;
            //std::cout<<"beta = "<<bc<<std::endl;
#endif
        }
    };


    // algo 21: Block algorithm
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<21,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta)
        {
            typedef typename M1::value_type T;

            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 21: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            //typedef typename M1::copy_type M1x;
            //typedef typename V::copy_type Vx;
            //M1x Ac = A;
            //Vx bc = beta;
#endif

            typedef typename M1::col_sub_type M1c;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;

            const ptrdiff_t NB = 4;
            typedef typename MCopyHelper<T,UpperTri,NB,NB>::type Ztype;
            typedef typename Ztype::subtrimatrix_type Zs;
            Ztype BaseZ = MatrixSizer<T>(NB,NB);
            typename Ztype::view_type Z = BaseZ.view();

            typedef typename V::iterator IT;

            typedef typename MCopyHelper<T,Rec,NB,Unknown>::type M3;
            typedef typename M3::col_sub_type M3c;
            const ptrdiff_t Si = M3::_stepi;
            const ptrdiff_t Sj = M3::_stepj;
            typedef typename MViewHelper<T,Rec,NB,Unknown,Si,Sj>::type M3cr;

            M3 tempBase = MatrixSizer<T>(NB,TMV_MAX(ptrdiff_t(1),N-NB));

            IT bj = beta.begin();
            ptrdiff_t j1=0;
            for(ptrdiff_t j2=j1+NB; j2<N; j1=j2,j2+=NB) {
                M1s A1 = A.subMatrix(j1,M,j1,j2);
                for(ptrdiff_t j=j1;j<j2;++j,++bj) {
                    M1c u = A.col(j,j+1,M);
                    HouseholderReflect(A.ref(j,j),u,*bj);
                    M1r A2a = A.row(j,j+1,j2);
                    M1s A2b = A.subMatrix(j+1,M,j+1,j2);
                    M3c temp = tempBase.col(0,0,j2-j-1);
                    HouseholderMultEq(u,*bj,A2a,A2b,temp);
                    M1s A3 = A.subMatrix(j1,M,j1,j+1);
                    Zs Z1 = Z.subTriMatrix(0,j-j1+1);
                    BlockHouseholderAugment(A3,Z1,*bj);
                }
                M1s A4 = A.subMatrix(j1,M,j2,N);
                M3cr temp = tempBase.colRange(0,N-j2);
                BlockHouseholderLDiv(A1,Z,A4,temp);
            }
            M1s A1 = A.subMatrix(j1,M,j1,N);
            for(ptrdiff_t j=j1;j<N;++j,++bj) {
                M1c u = A.col(j,j+1,M);
                HouseholderReflect(A.ref(j,j),u,*bj);
                M1r A2a = A.row(j,j+1,N);
                M1s A2b = A.subMatrix(j+1,M,j+1,N);
                M3c temp = tempBase.col(0,0,N-j-1);
                HouseholderMultEq(u,*bj,A2a,A2b,temp);
            }
#ifdef PRINTALGO_QR
            //std::cout<<"A -> "<<A<<std::endl;
            //std::cout<<"beta -> "<<beta<<std::endl;
            //QRDecompose_Helper<11,cs,rs,M1x,Vx>::call(Ac,bc);
            //std::cout<<"correct A = "<<Ac<<std::endl;
            //std::cout<<"beta = "<<bc<<std::endl;
            //std::cout<<"diff = "<<(A-Ac).copy().clip(1.e-5)<<std::endl;
            //std::cout<<(beta-bc).copy().clip(1.e-5)<<std::endl;
#endif
        }
    };

    // Used by both algo 22 and 27.
    template <class M1, class M2>
    inline void RecursiveQRDecompose(M1& A, M2& Z, bool makeZ)
    {
        // This is very similar to the BlockHouseholder_MakeZ function
        // in Householder.cpp.  The difference is the addition of the 
        // Householder_Reflects.
        // The makeZ parameter should be set to true if you want the Z
        // matrix to be correct on output.  If you don't need the Z matrix
        // after this call, setting makeZ to false will speed it up slightly.
        // (In either case, the diagonal of Z is correctly set to beta.)

        const ptrdiff_t cs = M1::_colsize;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_size>::size;
        const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
        const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;

        typedef typename M1::value_type T;
        typedef typename M1::real_type RT;
        typedef typename M1::col_sub_type M1c;
        typedef typename M1::row_sub_type M1r;
        typedef typename M1::submatrix_type M1s;
        typedef typename M2::subtrimatrix_type M2s;
        typedef typename M2::submatrix_type M2sm;
        typedef typename M2::col_sub_type M2c;

        if (N > 2) {
            ptrdiff_t j1 = N/2;
            M1s A1 = A.colRange(0,j1);
            M2s Z1 = Z.subTriMatrix(0,j1);
            RecursiveQRDecompose(A1,Z1,true);

            M1s A2x = A.colRange(j1,N);
            M2sm Z3 = Z.subMatrix(0,j1,j1,N);
            // Use Z3 as the temporary, since it happens to be the right
            // shape and we don't need to write to it yet.
            BlockHouseholderLDiv(A1,Z1,A2x,Z3);

            M1s A2 = A.subMatrix(j1,M,j1,N);
            M2s Z2 = Z.subTriMatrix(j1,N);
            RecursiveQRDecompose(A2,Z2,makeZ);

            if (makeZ) {
                Z3 = A1.rowRange(j1,N).adjoint() *
                    A.subMatrix(j1,N,j1,N).unitLowerTri();
                Z3 += A1.rowRange(N,M).adjoint() * A.subMatrix(N,M,j1,N);
                Z3 = -Z1*Z3;
                Z3 *= Z2;
            }
        } else if (N==2) {
            M1c u0 = A.col(0,1,M);
            RT b0;
            HouseholderReflect(A.ref(0,0),u0,b0);
            Z.ref(0,0) = b0;
            M1r A1a = A.row(0,1,N);
            M1s A1b = A.subMatrix(1,M,1,N);
            M2c temp = Z.col(1,0,1);
            HouseholderMultEq(u0,b0,A1a,A1b,temp);
            M1c u1 = A.col(1,2,M);
            RT b1;
            HouseholderReflect(A.ref(1,1),u1,b1);
            Z.ref(1,1) = b1;
            if (makeZ) {
                T temp = A.col(0,2,M).conjugate()*u1;
                temp += TMV_CONJ(A.cref(1,0));
                Z.ref(0,1) = -b0*b1*temp;
            }
        } else { // N == 1
            M1c u = A.col(0,1,M);
            RT b0;
            HouseholderReflect(A.ref(0,0),u,b0);
            Z.ref(0,0) = b0;
        }
    }

    // algo 22: Block algorithm, using recursive within each block
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<22,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta)
        {
            typedef typename M1::value_type T;

            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 22: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M1::col_sub_type M1c;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;
            const ptrdiff_t Si1 = M1::_stepi;
            const ptrdiff_t Sj1 = M1::_stepj;
            const int C = M1::_conj ? Conj : NonConj;

            const ptrdiff_t NB = TMV_QR_BLOCKSIZE;

            const ptrdiff_t s1 = IntTraits2<NB,rs>::min;
            const ptrdiff_t N1 = TMV_MIN(NB,N);
            typedef typename MCopyHelper<T,UpperTri,s1,s1>::type Ztype;
            Ztype BaseZ = MatrixSizer<T>(N1,N1);
            typename Ztype::view_type Z = BaseZ.view();

            const ptrdiff_t s2 = IntTraits2<rs,NB>::diff;
            const ptrdiff_t s3 = IntTraits2<s1,s2>::max;
            const ptrdiff_t N3 = TMV_MAX(N1,N-NB);
            typedef typename MCopyHelper<T,Rec,s1,s3>::type M3;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3::colrange_type M3cr;
            M3 tempBase = MatrixSizer<T>(N1,N3);

            typedef typename MViewHelper<T1,Rec,Unknown,NB,Si1,Sj1,C>::type M1sa;
            const ptrdiff_t s4 = (
                rs == Unknown ? Unknown : 
                rs - NB*(rs/NB) );
            const ptrdiff_t s5 = (
                (cs == Unknown || rs == Unknown) ? Unknown :
                cs - NB*(rs/NB) );
            typedef typename MViewHelper<T1,Rec,s5,s4,Si1,Sj1,C>::type M1sb;

            typedef typename V::iterator IT;

            ptrdiff_t j1=0;
            for(ptrdiff_t j2=j1+NB; j2<N; j1=j2,j2+=NB) {
                M1sa A1 = A.subMatrix(j1,M,j1,j2);

                RecursiveQRDecompose(A1,Z,true);

                M1s A2 = A.subMatrix(j1,M,j2,N);
                M3cr temp = tempBase.colRange(0,N-j2);
                BlockHouseholderLDiv(A1,Z,A2,temp);
                beta.subVector(j1,j2) = Z.diag().realPart();
            }

            M1sb A1 = A.subMatrix(j1,M,j1,N);
            IT bj = beta.subVector(j1,N).begin();
            for(ptrdiff_t j=j1;j<N;++j,++bj) {
                M1c u = A.col(j,j+1,M);
                HouseholderReflect(A.ref(j,j),u,*bj);
                M1r A2a = A.row(j,j+1,N);
                M1s A2b = A.subMatrix(j+1,M,j+1,N);
                M3c temp = tempBase.col(0,0,N-j-1);
                HouseholderMultEq(u,*bj,A2a,A2b,temp);
            }
        }
    };

    // algo 26: Unrolled recursive
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<26,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        typedef typename M1::real_type RT;
        typedef typename MCopyHelper<T,UpperTri,rs,rs>::type M2;
        enum { ZS1 = 1 };
        enum { ZS2 = rs };
        enum { ZA = ColMajor };
        enum { AA = (
                (M1::_conj ? Conj : NonConj) |
                (M1::_colmajor ? ColMajor : M1::_rowmajor ? RowMajor : NonMajor)
        ) };
        enum { AS1 = M1::_stepi };
        enum { AS2 = M1::_stepj };
        enum { C = M1::_conj };
        typedef typename MCopyHelper<T,Rec,rs/2,(rs+1)/2>::type M3;
        enum { TS1 = M3::_stepi };
        enum { TS2 = M3::_stepj };

        // This prevents large matrices from unrolling a lot of 
        // intermediate MultMM, etc. calls.
        enum { csx = cs == Unknown || cs > 64 ? Unknown : cs };

        template <bool makeZ, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t N1>
        struct MakeZ_Helper  // makeZ = false, any N1
        { static TMV_INLINE void call(const M1& , M2& ) {} };

        template <ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t N1>
        struct MakeZ_Helper<true,j1,j2,N1> // N1 > 2
        {
            static inline void call(const M1& A, M2& Z)
            {
                TMVStaticAssert(j2 == j1+N1);
                //std::cout<<"MakeZ:\n";
                const ptrdiff_t jmid = IntTraits<IntTraits2<j1,j2>::sum>::halfS;
                const ptrdiff_t Na = jmid-j1;
                const ptrdiff_t Nb = j2-jmid;
                const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
                const ptrdiff_t Mx = IntTraits2<csx,j2>::diff;

                typename MViewHelper<T,UpperTri,Na,Na,ZS1,ZS2>::type Z1 = 
                    Z.subTriMatrix(j1,jmid);
                //std::cout<<"Z1 = "<<Z1<<std::endl;
                typename MViewHelper<T,UpperTri,Nb,Nb,ZS1,ZS2>::type Z2 = 
                    Z.subTriMatrix(jmid,j2);
                //std::cout<<"Z2 = "<<Z2<<std::endl;
                typename MViewHelper<T,Rec,Na,Nb,ZS1,ZS2>::type Z3 = 
                    Z.subMatrix(j1,jmid,jmid,j2);
                //std::cout<<"Z3 = "<<Z3<<std::endl;

                typename MViewHelper<T,Rec,Nb,Na,AS1,AS2,C>::ctype A1a = 
                    A.subMatrix(jmid,j2,j1,jmid);
                //std::cout<<"A1a = "<<A1a<<std::endl;
                typename MViewHelper<T,Rec,Mx,Na,AS1,AS2,C>::ctype A1b = 
                    A.subMatrix(j2,M,j1,jmid);
                //std::cout<<"A1b = "<<A1b<<std::endl;
                typename MViewHelper<T,Rec,Nb,Nb,AS1,AS2,C>::ctype A2a = 
                    A.subMatrix(jmid,j2,jmid,j2);
                //std::cout<<"A2a = "<<A2a<<std::endl;
                typename MViewHelper<T,Rec,Mx,Nb,AS1,AS2,C>::ctype A2b = 
                    A.subMatrix(j2,M,jmid,j2);
                //std::cout<<"A2b = "<<A2b<<std::endl;

                //std::cout<<"Z = "<<Z.subTriMatrix(j1,j2)<<std::endl;
                Z3 = A1a.adjoint() * A2a.unitLowerTri();
                //std::cout<<"Z3 = "<<Z3<<std::endl;
                Z3 += A1b.adjoint() * A2b;
                //std::cout<<"Z3 => "<<Z3<<std::endl;
                Z3 = -Z1*Z3;
                //std::cout<<"Z3 => "<<Z3<<std::endl;
                //std::cout<<"Z2 = "<<Z2<<std::endl;
                Z3 *= Z2;
                //std::cout<<"Z3 => "<<Z3<<std::endl;
            }
        };
        template <ptrdiff_t j1, ptrdiff_t j2>
        struct MakeZ_Helper<true,j1,j2,2> // N1 == 2
        {
            static inline void call(const M1& A, M2& Z)
            {
                TMVStaticAssert(j2 == j1+2);
                //std::cout<<"N==2 MakeZ:\n";
                const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
                const ptrdiff_t Mx = IntTraits2<csx,j2>::diff;
                typename VViewHelper<T,Mx,AS1,C>::ctype A1b = A.col(j1,j2,M);
                //std::cout<<"A1b = "<<A1b<<std::endl;
                typename VViewHelper<T,Mx,AS1,C>::ctype A2b = A.col(j1+1,j2,M);
                //std::cout<<"A2b = "<<A2b<<std::endl;
                T temp = A1b.conjugate()*A2b;
                //std::cout<<"temp = "<<temp<<std::endl;
                temp += TMV_CONJ(A.cref(j1+1,j1));
                //std::cout<<"temp => "<<temp<<std::endl;
                Z.ref(j1,j1+1) = -Z.cref(j1,j1)*Z.cref(j1+1,j1+1)*temp;
                //std::cout<<"Z3 -> "<<Z.cref(j1,j1+1)<<std::endl;
            }
        };
        template <ptrdiff_t j1, ptrdiff_t j2>
        struct MakeZ_Helper<true,j1,j2,1> // N1 == 1
        { static TMV_INLINE void call(const M1& , M2& ) {} };


        template <bool makeZ, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t N1>
        struct Unroll_Helper // N1 > 2
        {
            static inline void unroll(M1& A, M2& Z, M3& tempBase) 
            {
                TMVStaticAssert(j2 == j1+N1);
                const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
                //std::cout<<"Start unroll: N = "<<N1<<"  "<<j1<<','<<j2<<std::endl;
                const ptrdiff_t jmid = IntTraits<IntTraits2<j1,j2>::sum>::halfS;
                const ptrdiff_t Na = jmid-j1;
                const ptrdiff_t Nb = j2-jmid;
                const ptrdiff_t Mx1 = IntTraits2<csx,j1>::diff;

                Unroll_Helper<true,j1,jmid,jmid-j1>::unroll(A,Z,tempBase);
                //std::cout<<"After first recurse: "<<N1<<std::endl;

                typename MViewHelper<T,Rec,Mx1,Na,AS1,AS2,C>::type A1 = 
                    A.subMatrix(j1,M,j1,jmid);
                //std::cout<<"A1 = "<<A1<<std::endl;
                typename MViewHelper<T,UpperTri,Na,Na,ZS1,ZS2>::type Z1 = 
                    Z.subTriMatrix(j1,jmid);
                //std::cout<<"Z1 = "<<Z1<<std::endl;
                typename MViewHelper<T,Rec,Mx1,Nb,AS1,AS2,C>::type A2x = 
                    A.subMatrix(j1,M,jmid,j2);
                //std::cout<<"A2x = "<<A2x<<std::endl;
                typename MViewHelper<T,Rec,Na,Nb,TS1,TS2>::type temp = 
                    tempBase.subMatrix(0,jmid-j1,0,j2-jmid);
                //std::cout<<"temp = "<<temp<<std::endl;
                BlockHouseholderLDiv(A1,Z1,A2x,temp);
                //std::cout<<"After LDiv: "<<N1<<std::endl;

                Unroll_Helper<makeZ,jmid,j2,j2-jmid>::unroll(A,Z,tempBase);
                //std::cout<<"After second recurse: "<<N1<<std::endl;

                MakeZ_Helper<makeZ,j1,j2,N1>::call(A,Z);
                //std::cout<<"After MakeZ: "<<N1<<std::endl;
            }
        };

        template <bool makeZ, ptrdiff_t j1, ptrdiff_t j2>
        struct Unroll_Helper<makeZ,j1,j2,2> // N==2
        {
            static inline void unroll(M1& A, M2& Z, M3& tempBase) 
            {
                TMVStaticAssert(j2 == j1+2);
                const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
                //std::cout<<"Start unroll: N = "<<2<<"  "<<j1<<','<<j2<<std::endl;
                const ptrdiff_t Mx1 = IntTraits2<csx,j1+1>::diff;
                const ptrdiff_t Mx2 = IntTraits2<csx,j2>::diff;

                typename VViewHelper<T,Mx1,AS1,C>::type u0 = A.col(j1,j1+1,M);
                //std::cout<<"u0 = "<<u0<<std::endl;
                RT b0;
                HouseholderReflect(A.ref(j1,j1),u0,b0);
                Z.ref(j1,j1) = b0;
                //std::cout<<"After first reflect: "<<2<<std::endl;

                typename VViewHelper<T,1,AS2,C>::type A1a = A.row(j1,j1+1,j2);
                //std::cout<<"A1a = "<<A1a<<std::endl;
                typename MViewHelper<T,Rec,Mx1,1,AS1,AS2,C>::type A1b =
                    A.subMatrix(j1+1,M,j1+1,j2);
                //std::cout<<"A1b = "<<A1b<<std::endl;
                SmallVectorView<T,1,1> temp = tempBase.col(0,0,1);
                HouseholderMultEq(u0,b0,A1a,A1b,temp);
                //std::cout<<"After MultEq: "<<2<<std::endl;

                RT b1;
                typename VViewHelper<T,Mx2,AS1,C>::type u1 = A.col(j1+1,j2,M);
                //std::cout<<"u1 = "<<u1<<std::endl;
                HouseholderReflect(A.ref(j1+1,j1+1),u1,b1);
                Z.ref(j1+1,j1+1) = b1;
                //std::cout<<"After second reflect: "<<2<<std::endl;

                MakeZ_Helper<makeZ,j1,j2,2>::call(A,Z);
                //std::cout<<"After MakeZ: "<<2<<std::endl;
            }
        };

        template <bool makeZ, ptrdiff_t j1, ptrdiff_t j2>
        struct Unroll_Helper<makeZ,j1,j2,1> // N == 1
        {
            static inline void unroll(M1& A, M2& Z, M3& tempBase) 
            {
                TMVStaticAssert(j2 == j1+1);
                const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
                //std::cout<<"Start unroll: N = "<<1<<"  "<<j1<<','<<j2<<std::endl;
                const ptrdiff_t Mx = IntTraits2<csx,j2>::diff;
                typename VViewHelper<T,Mx,AS1,C>::type u0 = A.col(j1,j1+1,M);
                //std::cout<<"u0 = "<<u0<<std::endl;
                RT b0;
                HouseholderReflect(A.ref(j1,j1),u0,b0);
                Z.ref(j1,j1) = b0;
                //std::cout<<"Done N = "<<1<<std::endl;
            }
        };
        static void call(M1& A, V& beta)
        {
            TMVStaticAssert(rs != Unknown);
#ifdef PRINTALGO_QR
            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
            std::cout<<"QRDecompose algo 26: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            //typedef typename M1::copy_type M1x;
            //typedef typename V::copy_type Vx;
            //M1x Ac = A;
            //Vx bc = beta;
#endif
            M2 Z;
            M3 tempBase;
            Unroll_Helper<false,0,rs,rs>::unroll(A,Z,tempBase);
            beta = Z.diag().realPart();
#ifdef PRINTALGO_QR
            //std::cout<<"A -> "<<A<<std::endl;
            //std::cout<<"beta -> "<<beta<<std::endl;
            //QRDecompose_Helper<27,cs,rs,M1x,Vx>::call(Ac,bc);
            //std::cout<<"correct A = "<<Ac<<std::endl;
            //std::cout<<"beta = "<<bc<<std::endl;
#endif
        }
    };

    // algo 27: Recursive algorithm
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<27,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta)
        {
            typedef typename M1::value_type T;

            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
#ifdef PRINTALGO_QR
            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            std::cout<<"QRDecompose algo 27: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename MCopyHelper<T,UpperTri,rs,rs>::type Ztype;
            Ztype Z = MatrixSizer<T>(N,N);
            typename Ztype::view_type Zv = Z.view();

            RecursiveQRDecompose(A,Zv,false);
            beta = Z.diag().realPart();
        }
    };

    // algo 31: Decide which algorithm to use from runtime size
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<31,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta)
        {
            typedef typename M1::value_type T;

            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 31: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const int algo27 = 
                (rs == Unknown || rs <= 128) ? 27 : 0;
            const int algo22 = 
                (rs == Unknown || rs > 128) ? 22 : 0;
            const ptrdiff_t l2cache = TMV_L2_CACHE*1024/sizeof(T);

            // I'm assuming that this first transition is primarily 
            // memory related.  If the full matrix fits in the L2 cache,
            // then all the recursion and blocking stuff don't help much.
            // And the extra calculations for doing the Z matrix are
            // just wasted extra flops.  
            if (M*N <= l2cache)
                QRDecompose_Helper<11,cs,rs,M1,V>::call(A,beta);

            // I'm not sure if this next transition should really be 
            // in terms of one of the caches...  
            // I did my testing for T = double, so sizeof(T) = 8,
            // and this transition seemed pretty independent of M.
            // But perhaps (N*N*2 <= l2cache) ?
            else if (N <= 128)
                QRDecompose_Helper<algo27,cs,rs,M1,V>::call(A,beta);
            else
                QRDecompose_Helper<algo22,cs,rs,M1,V>::call(A,beta);
        }
    };

    // algo 32: Decide which algorithm to use from runtime column size
    // given that rs is known and <= 32
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<32,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta)
        {
            TMVStaticAssert(rs != Unknown);
            typedef typename M1::value_type T;

            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 32: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const ptrdiff_t l1cache = TMV_L1_CACHE*1024/sizeof(T);
            // For this one, the L1 cache seems to be the relevant value.  
            if (M*N <= l1cache)
                QRDecompose_Helper<16,cs,rs,M1,V>::call(A,beta);
            else
                QRDecompose_Helper<26,cs,rs,M1,V>::call(A,beta);
        }
    };

    // algo 33: Decide which algorithm to use from runtime column size
    // given that rs is known and > 32
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<33,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta)
        {
            TMVStaticAssert(rs != Unknown);

            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
#ifdef PRINTALGO_QR
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
            std::cout<<"QRDecompose algo 33: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            // Again, I don't think this transition is memory related, but
            // my timing tests were only done for double, so sizeof(T) = 8.
            if (M <= 128)
                QRDecompose_Helper<11,cs,rs,M1,V>::call(A,beta);
            else
                QRDecompose_Helper<26,cs,rs,M1,V>::call(A,beta);
        }
    };

    // algo 81: Copy to colmajor
    template <ptrdiff_t cs, ptrdiff_t rs, class M, class V>
    struct QRDecompose_Helper<81,cs,rs,M,V>
    {
        static inline void call(M& m, V& beta)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 81: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M::value_type T;
            typedef typename MCopyHelper<T,Rec,cs,rs>::type Mcm;
            Mcm mcm = m;
            QRDecompose_Helper<-2,cs,rs,Mcm,V>::call(mcm,beta);
            typename M::noalias_type mna = m.noAlias();
            Copy(mcm,mna);
        }
    };

    // algo 90: call InstQR_Decompose
    template <ptrdiff_t cs, ptrdiff_t rs, class M, class V>
    struct QRDecompose_Helper<90,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta)
        { InstQR_Decompose(m.xView(),beta.xView()); }
    };

    // algo 97: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, class M, class V>
    struct QRDecompose_Helper<97,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            QRDecompose_Helper<-2,cs,rs,Mc,V>::call(mc,beta);
        }
    };

    // algo -4: No copies or branches
    template <ptrdiff_t cs, ptrdiff_t rs, class M, class V>
    struct QRDecompose_Helper<-4,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta)
        {
            typedef typename M::value_type T;
            const ptrdiff_t maxunroll = 
                TMV_OPT==0 ? 0 : 
                TMV_OPT==1 ? 4 : 
                TMV_OPT==2 ? 16 : 
                100;
            const ptrdiff_t csrs = IntTraits2<cs,rs>::prod;
            const ptrdiff_t rsrs = IntTraits2<rs,rs>::prod;
            const ptrdiff_t l1cache = TMV_L1_CACHE*1024/sizeof(T);
            const ptrdiff_t l2cache = TMV_L2_CACHE*1024/sizeof(T);
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                rs == Unknown ? 31 :
                rs <= maxunroll ? (
                    rs <= 32 ? (
                        cs == Unknown ? 32 :
                        (csrs <= l1cache ? 16 : 26) ) :
                    cs == Unknown ? 33 : 
                    cs <= 128 ? 11 : 26) :
                rsrs > l2cache ? (rs <= 128 ? 27 : 22) :
                cs == Unknown ? 31 : 
                csrs <= l2cache ? 11 :
                rs <= 128 ? 27 : 22;
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRDecompose: \n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<m.colsize()<<"  "<<m.rowsize()<<std::endl;
            //std::cout<<"maxunroll = "<<maxunroll<<std::endl;
            //std::cout<<"csrs = "<<csrs<<std::endl;
            //std::cout<<"l1cache = "<<l1cache<<std::endl;
            //std::cout<<"l2cache = "<<l2cache<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
#endif
            QRDecompose_Helper<algo,cs,rs,M,V>::call(m,beta);
            //std::cout<<"m => "<<m<<std::endl;
            //std::cout<<"beta => "<<beta<<std::endl;
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class V>
    struct QRDecompose_Helper<-3,cs,rs,M1,V>
    {
        static TMV_INLINE void call(M1& m, V& beta)
        {
            const int algo = (
                ( cs != Unknown && rs != Unknown &&
                  cs <= 16 && rs <= 16 ) ? -4 :
                ( TMV_OPT >= 2 && !M1::_colmajor ) ? 81 :
                -4 );
#ifdef PRINTALGO_QR
            const ptrdiff_t M = cs==Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? m.rowsize() : rs;
            std::cout<<"QRDecompose algo -3: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            QRDecompose_Helper<algo,cs,rs,M1,V>::call(m,beta);
#ifdef PRINTALGO_QR
            std::cout<<"Done QRDecompose\n";
#endif
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M, class V>
    struct QRDecompose_Helper<-2,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta)
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
            QRDecompose_Helper<algo,cs,rs,M,V>::call(m,beta);
        }
    };

    template <ptrdiff_t cs, ptrdiff_t rs, class M, class V>
    struct QRDecompose_Helper<-1,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta)
        { QRDecompose_Helper<-2,cs,rs,M,V>::call(m,beta); }
    };

    template <class M, class V>
    inline void InlineQR_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta)
    {
        const ptrdiff_t cs = M::_colsize;
        const ptrdiff_t rs = Sizes<M::_rowsize,V::_size>::size;
        typedef typename M::cview_type Mv;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TMV_MAYBE_REF(V,Vv) betav = beta.cView();
        QRDecompose_Helper<-3,cs,rs,Mv,Vv>::call(mv,betav);
    }

    // This is the basic functionality
    template <class M, class V>
    inline void QR_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta)
    {
        TMVStaticAssert(V::isreal);
        TMVStaticAssert((Traits2<
                         typename M::value_type,
                         typename V::value_type>::samebase));
        TMVAssert(m.colsize() >= m.rowsize());
        TMVAssert(beta.size() == m.rowsize());
        TMVStaticAssert((Sizes<M::_rowsize,V::_size>::same));
        TMVAssert(m.rowsize() == beta.size());

        const ptrdiff_t cs = M::_colsize;
        const ptrdiff_t rs = Sizes<M::_rowsize,V::_size>::size;
        typedef typename M::cview_type Mv;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TMV_MAYBE_REF(V,Vv) betav = beta.cView();
        QRDecompose_Helper<-2,cs,rs,Mv,Vv>::call(mv,betav);
    }

    // The rest of these below are basically convenience functions
    // to allow the user to provide fewer or different arguments.
    template <class M1, class M2>
    inline void QR_Decompose(
        BaseMatrix_Rec_Mutable<M1>& Q, BaseMatrix_Tri_Mutable<M2>& R)
    {
        TMVStaticAssert((Traits2<
                         typename M1::value_type,
                         typename M2::value_type>::sametype));
        TMVStaticAssert(M2::_upper);
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(R.size() == Q.rowsize());

        typedef typename M1::real_type RT;
        Vector<RT> beta(Q.rowsize());
        QR_Decompose(Q,beta);
        Copy(Q.upperTri(),R);
        UnpackQ(Q,beta);
    }

    template <class M>
    inline void QR_Decompose(BaseMatrix_Rec_Mutable<M>& A)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        typedef typename M::real_type RT;
        Vector<RT> beta(A.rowsize());
        QR_Decompose(A,beta);
    }


    // Allow views as an argument by value (for convenience)
    template <class T, int A, int A2>
    inline void QR_Decompose(
        MatrixView<T,A> Q, UpperTriMatrixView<T,A2> R)
    {
        typedef MatrixView<T,A> M1;
        typedef UpperTriMatrixView<T,A2> M2;
        QR_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R)); 
    }

    template <class T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
    inline void QR_Decompose(
        SmallMatrixView<T,M,N,Si,Sj,A> Q,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> R)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        typedef SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> M2;
        QR_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R));
    }

    template <class T, ptrdiff_t N, int A, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
    inline void QR_Decompose(
        MatrixView<T,A> Q,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> R)
    {
        typedef MatrixView<T,A> M1;
        typedef SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> M2;
        QR_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R));
    }

    template <class T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A, int A2>
    inline void QR_Decompose(
        SmallMatrixView<T,M,N,Si,Sj,A> Q,
        UpperTriMatrixView<T,A2> R)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        typedef UpperTriMatrixView<T,A2> M2;
        QR_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R));
    }

    template <class T, int A>
    inline void QR_Decompose(MatrixView<T,A> m)
    {
        typedef MatrixView<T,A> M1;
        QR_Decompose(static_cast<BaseMatrix_Rec_Mutable<M1>&>(m));
    }

    template <class T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    inline void QR_Decompose(SmallMatrixView<T,M,N,Si,Sj,A> m)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        QR_Decompose(static_cast<BaseMatrix_Rec_Mutable<M1>&>(m));
    }

    // Don't forget the ones that mix *MatrixView with BaseMatrix_*_Mutable
    template <class T, int A, class M2>
    inline void QR_Decompose(
        MatrixView<T,A> Q, BaseMatrix_Tri_Mutable<M2>& R)
    {
        typedef MatrixView<T,A> M1;
        QR_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),R);
    }

    template <class M1, class T, int A2>
    inline void QR_Decompose(
        BaseMatrix_Rec_Mutable<M1>& Q, UpperTriMatrixView<T,A2> R)
    {
        typedef UpperTriMatrixView<T,A2> M2;
        QR_Decompose(
            Q,static_cast<BaseMatrix_Tri_Mutable<M2>&>(R)); 
    }

    template <class T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A, class M2>
    inline void QR_Decompose(
        SmallMatrixView<T,M,N,Si,Sj,A> Q, BaseMatrix_Tri_Mutable<M2>& R)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        QR_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),R);
    }

    template <class M1, class T, ptrdiff_t N, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
    inline void QR_Decompose(
        BaseMatrix_Rec_Mutable<M1>& Q,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> R)
    {
        typedef SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> M2;
        QR_Decompose(
            Q,static_cast<BaseMatrix_Tri_Mutable<M2>&>(R));
    }

} // namespace tmv

#undef TMV_QR_BLOCKSIZE

#endif

