

#ifndef TMV_QRDowndate_H
#define TMV_QRDowndate_H

//#define PRINTALGO_QR

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Householder.h"
#include "TMV_MatrixIO.h"
#include "TMV_TriMatrixIO.h"

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

#ifndef TMV_NO_THROW
    template <class T>
    class BadQRDowndate : 
        public NonPosDef
    {
    public:
        UpperTriMatrix<T> R;
        Matrix<T,NoDivider> A;

        template <class M1, class M2>
        BadQRDowndate(
            const BaseMatrix_Tri<M1>& _R, 
            const BaseMatrix_Rec<M2>& _A) throw() :
            NonPosDef("QR Downdate."), R(_R), A(_A) {}
        BadQRDowndate() throw() :
            NonPosDef("QR Downdate."), R(0), A(0,0) {}
        BadQRDowndate(const BadQRDowndate<T>& rhs) throw() :
            NonPosDef("QR Downdate."), R(rhs.R), A(rhs.A) {}
        ~BadQRDowndate() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV NonPosDef: QR Downdate found that the resulting \n";
            os<<"down-dated RtR is not positive definite. \n";
            os<<"(and hence the down date is impossible)\n";
            os<<"The partially downdated matrix is \n"<<R<<std::endl;
            os<<"The matrix attempting to be down-dated was \n"<<A<<std::endl;
        }
    };
#endif


    // Defined in TMV_QRDowndate.cpp
    template <class T, int C>
    void InstQR_Downdate(UpperTriMatrixView<T> R, MatrixView<T,C> A);

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or 1, or N == 0)
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<0,cs,rs,M1,M2>
    { static TMV_INLINE void call(M1& , M2& ) {} };

    // algo 11: Non-block algorithm, loop over n
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<11,cs,rs,M1,M2>
    {
        static void call(M1& R, M2& A)
        {
            const ptrdiff_t N = rs == Unknown ? R.size() : rs;
#ifdef PRINTALGO_QR
            const ptrdiff_t M = cs == Unknown ? A.colsize() : cs;
            std::cout<<"QRDowndate algo 11: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            // Given that [ A0 ] = Q1 R1
            //            [ A  ] 
            // Find R0 so that A0 = Q0 R0
            // Input R is R1, output is R0

            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;
            typedef typename M2::col_type M2c;
            typedef typename M1::row_sub_type M1r;
            typedef typename M2::colrange_type M2cr;
            typedef typename VCopyHelper<T,rs>::type V2;
            typedef typename V2::subvector_type V2s;
            V2 tempBase = VectorSizer<T>(N);

            RT beta;
            for(ptrdiff_t j=0;j<N;++j) {
                // Apply the Householder Reflection for this column
                M2c u = A.col(j);
                if (!HouseholderUnReflect(R.ref(j,j),u,beta)) {
#ifdef TMV_NO_THROW
                    std::cerr<<"Bad QR Downdate\n";
                    exit(1);
#else
                    throw BadQRDowndate<T>(R,A);
#endif
                }
                M1r Rj = R.row(j,j+1,N);
                M2cr Aj = A.colRange(j+1,N);

                // Rj' = Rj - beta Rj - beta ut Aj
                // Rj = (Rj' + beta ut Aj) / (1-beta)
                TMVAssert(beta != RT(1));
                V2s temp = tempBase.subVector(0,N-j-1);
                temp = beta * u.conjugate() * Aj;
                Rj += temp;
                Rj /= T(1) - beta;

                // Aj' = Aj - beta u (Rj + ut Aj)
                temp += beta * Rj;
                Aj -= u ^ temp;
            }
        }
    };
    
    // TODO: Add a cs == 1 algorithm.  I think Givens rotations are 
    // probably faster if there is only one row to remove.

    // algo 21: Block algorithm
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<21,cs,rs,M1,M2>
    {
        static void call(M1& R, M2& A)
        {
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
#ifdef PRINTALGO_QR
            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            std::cout<<"QRDowndate algo 21: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;

            typedef typename M2::col_type M2c;
            typedef typename M2::colrange_type M2cr;

            const ptrdiff_t NB = 4;
            typedef typename MCopyHelper<T,UpperTri,NB,NB>::type Ztype;
            typedef typename Ztype::subtrimatrix_type Zs;
            typedef typename Ztype::col_sub_type Zc;
            Ztype BaseZ = MatrixSizer<T>(NB,NB);
            typename Ztype::view_type Z = BaseZ.view();

            typedef typename MCopyHelper<T,Rec,NB,Unknown>::type M3;
            typedef typename M3::col_sub_type M3c;
            const ptrdiff_t Si = M3::_stepi;
            const ptrdiff_t Sj = M3::_stepj;
            typedef typename MViewHelper<T,Rec,NB,Unknown,Si,Sj>::type M3cr;

            M3 tempBase = MatrixSizer<T>(NB,TMV_MAX(ptrdiff_t(1),N-NB));

            RT beta;
            ptrdiff_t j1=0;
            for(ptrdiff_t j2=j1+NB; j2<N; j1=j2,j2+=NB) {
                M2cr A1 = A.colRange(j1,j2);
                for(ptrdiff_t j=j1;j<j2;++j) {
                    M2c u = A.col(j);
                    if (!HouseholderUnReflect(R.ref(j,j),u,beta)) {
#ifdef TMV_NO_THROW
                        std::cerr<<"Bad QR Downdate\n";
                        exit(1);
#else
                        throw BadQRDowndate<T>(R,A);
#endif
                    }
                    M1r Rj = R.row(j,j+1,j2);
                    M2cr Aj = A.colRange(j+1,j2);
                    M3c temp = tempBase.col(0,0,j2-j-1);

                    TMVAssert(beta != RT(1));
                    temp = beta * u.conjugate() * Aj;
                    Rj += temp;
                    Rj /= T(1) - beta;
                    temp += beta * Rj;
                    Aj -= u ^ temp;

                    M2cr A3 = A.colRange(j1,j+1);
                    Zs Z1 = Z.subTriMatrix(0,j-j1+1);
                    Block2HouseholderAugment(A3,Z1,beta);
                }
                M1s R4 = R.subMatrix(j1,j2,j2,N);
                M2cr A4 = A.colRange(j2,N);
                M3cr temp = tempBase.colRange(0,N-j2);

                // R4' = R4 - Zt(Yt A4 + R4)
                // R4' + Zt Yt A4 = (I-Zt) R4
                // R4 = (I-Zt)^-1 (R4' + Zt Yt A4)
                temp = A1.adjoint() * A4;
                temp = Z.adjoint() * temp;
                R4 += temp;
                Z -= T(1); // Now Z is Z-I
                R4 /= -Z.adjoint();
                Z += T(1); // Now Z is Z again.
                temp += Z.adjoint() * R4;
                A4 -= A1 * temp;
            }
            M2cr A1 = A.colRange(j1,N);
            for(ptrdiff_t j=j1;j<N;++j) {
                M2c u = A.col(j);
                if (!HouseholderUnReflect(R.ref(j,j),u,beta))
                {
#ifdef TMV_NO_THROW
                    std::cerr<<"Bad QR Downdate\n";
                    exit(1);
#else
                    throw BadQRDowndate<T>(R,A);
#endif
                }
                M1r Rj = R.row(j,j+1,N);
                M2cr A2 = A.colRange(j+1,N);
                M3c temp = tempBase.col(0,0,N-j-1);

                TMVAssert(beta != RT(1));
                temp = beta * u.conjugate() * A2;
                Rj += temp;
                Rj /= T(1) - beta;
                temp += beta * Rj;
                A2 -= u ^ temp;
            }
        }
    };

    // Used by both algo 22 and 27
    template <class M1, class M2, class M3>
    inline void RecursiveQRDowndate(M1& R, M2& A, M3& Z, bool makeZ)
    {
        //const ptrdiff_t cs = M2::_colsize;
        //const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
        const ptrdiff_t rs = Sizes<M2::_rowsize,M1::_size>::size;
        const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;

        typedef typename M1::value_type T;
        typedef typename M1::real_type RT;
        typedef typename M1::subtrimatrix_type M1s;
        typedef typename M1::submatrix_type M1sm;
        typedef typename M2::col_sub_type M2c;
        typedef typename M1::row_sub_type M1r;
        typedef typename M2::colrange_type M2cr;
        typedef typename M3::subtrimatrix_type M3s;
        typedef typename M3::submatrix_type M3sm;
        typedef typename M3::col_sub_type M3c;

        if (N > 2) {
            ptrdiff_t j1 = N/2;
            M1s R1 = R.subTriMatrix(0,j1);
            M2cr A1 = A.colRange(0,j1);
            M3s Z1 = Z.subTriMatrix(0,j1);
            RecursiveQRDowndate(R1,A1,Z1,true);

            M1sm R3 = R.subMatrix(0,j1,j1,N);
            M2cr A2 = A.colRange(j1,N);
            M3sm Z3 = Z.subMatrix(0,j1,j1,N);

            // R3' = R3 - Zt(Yt A2 + R3)
            // R3' + Zt Yt A2 = (I-Zt) R3
            // R3 = (I-Zt)^-1 (R3' + Zt Yt A2)
            // Use Z3 as the temporary, since it happens to be the right
            // shape and we don't need to write to it yet.
            Z3 = A1.adjoint() * A2;
            Z3 = Z1.adjoint() * Z3;
            R3 += Z3;
            Z1 -= T(1); // Now Z1 is Z1-I
            R3 /= -Z1.adjoint();
            Z1 += T(1); // Now Z1 is Z1 again.
            Z3 += Z1.adjoint() * R3;
            A2 -= A1 * Z3;

            M1s R2 = R.subTriMatrix(j1,N);
            M3s Z2 = Z.subTriMatrix(j1,N);
            RecursiveQRDowndate(R2,A2,Z2,makeZ);

            if (makeZ) {
                Z3 = A1.adjoint() * A2;
                Z3 = -Z1*Z3;
                Z3 *= Z2;
            }
        } else if (N==2) {
            M2c u0 = A.col(0);
            M2c u1 = A.col(1);
            RT b0,b1;
            if (!HouseholderUnReflect(R.ref(0,0),u0,b0)) {
#ifdef TMV_NO_THROW
                std::cerr<<"Bad QR Downdate\n";
                exit(1);
#else
                throw BadQRDowndate<T>();
#endif
            }
            Z.ref(0,0) = b0;
            if (b0 != RT(0)) {
                TMVAssert(b0 != RT(1));
                T temp = b0 * u0.conjugate() * u1;
                R.ref(0,1) += temp;
                R.ref(0,1) /= T(1) - b0;

                temp += b0 * R.cref(0,1);
                A.col(1) -= temp * u0;
            }
            if (!HouseholderUnReflect(R.ref(1,1),u1,b1)) {
#ifdef TMV_NO_THROW
                std::cerr<<"Bad QR Downdate\n";
                exit(1);
#else
                throw BadQRDowndate<T>();
#endif
            }
            Z.ref(1,1) = b1;
            if (makeZ) {
                T temp = u0.conjugate()*u1;
                Z.ref(0,1) = -b0*b1*temp;
            }
        } else { // N == 1
            M2c u = A.col(0);
            RT b0;
            if (!HouseholderUnReflect(R.ref(0,0),u,b0)) {
#ifdef TMV_NO_THROW
                std::cerr<<"Bad QR Downdate\n";
                exit(1);
#else
                throw BadQRDowndate<T>();
#endif
            }
            Z.ref(0,0) = b0;
        }
    }

    // algo 22: Block algorithm, using recursive within each block
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<22,cs,rs,M1,M2>
    {
        static void call(M1& R, M2& A)
        {
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
#ifdef PRINTALGO_QR
            std::cout<<"QRDowndate algo 22: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const ptrdiff_t NB = TMV_QR_BLOCKSIZE;

            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;
            const ptrdiff_t Si1 = M1::_stepi;
            const ptrdiff_t Sj1 = M1::_stepj;
            const int C1 = M1::_conj ? Conj : NonConj | NonUnitDiag;
            // Note: The C term is usually just used to set the 
            // Conj or NonConj status.  But any attribute may be put here,
            // so put the NonUnitDiag attribute to keep it from being
            // UnknownDiag (the default).
            typedef typename MViewHelper<T,UpperTri,NB,NB,Si1,Sj1,C1>::type M1a;

            typedef typename M2::col_sub_type M2c;
            typedef typename M2::colrange_type M2cr;
            const ptrdiff_t Si2 = M2::_stepi;
            const ptrdiff_t Sj2 = M2::_stepj;
            const int C2 = M2::_conj ? Conj : NonConj;
            typedef typename MViewHelper<T,Rec,cs,NB,Si2,Sj2,C2>::type M2a;
            const ptrdiff_t s4 = (
                rs == Unknown ? Unknown :
                rs - NB*(rs/NB) );
            typedef typename MViewHelper<T,Rec,cs,s4,Si2,Sj2,C2>::type M2b;

            const ptrdiff_t s1 = IntTraits2<NB,rs>::min;
            const ptrdiff_t N1 = TMV_MIN(NB,N);
            typedef typename MCopyHelper<T,UpperTri,s1,s1>::type Ztype;
            typedef typename Ztype::subtrimatrix_type Zs;
            Ztype BaseZ = MatrixSizer<T>(N1,N1);
            typename Ztype::view_type Z = BaseZ.view();

            const ptrdiff_t s2 = IntTraits2<rs,NB>::diff;
            const ptrdiff_t s3 = IntTraits2<s1,s2>::max;
            const ptrdiff_t N3 = TMV_MAX(N1,N-NB);
            typedef typename MCopyHelper<T,Rec,s1,s3>::type M3;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3::colrange_type M3cr;
            typedef typename M3::submatrix_type M3s;
            M3 tempBase = MatrixSizer<T>(N1,N3);

            RT beta;
            ptrdiff_t j1=0;
            for(ptrdiff_t j2=j1+NB; j2<N; j1=j2,j2+=NB) {
                M1a R1 = R.subTriMatrix(j1,j2);
                M2a A1 = A.colRange(j1,j2);

#ifndef TMV_NOTHROW
                try {
#endif
                    RecursiveQRDowndate(R1,A1,Z,true);
#ifndef TMV_NOTHROW
                } catch(BadQRDowndate<T>) {
                    throw BadQRDowndate<T>(R,A);
                }
#endif

                M1s R2 = R.subMatrix(j1,j2,j2,N);
                M2cr A2 = A.colRange(j2,N);
                M3cr temp = tempBase.colRange(0,N-j2);

                // R2' = R2 - Zt(Yt A2 + R2)
                // R2' + Zt Yt A2 = (I-Zt) R2
                // R2 = (I-Zt)^-1 (R2' + Zt Yt A2)
                temp = A1.adjoint() * A2;
                temp = Z.adjoint() * temp;
                R2 += temp;
                Z -= T(1); // Now Z is Z-I
                R2 /= -Z.adjoint();
                Z += T(1); // Now Z is Z again.
                temp += Z.adjoint() * R2;
                A2 -= A1 * temp;
            }

            M2b A1 = A.subMatrix(j1,M,j1,N);
            for(ptrdiff_t j=j1;j<N;++j) {
                M2c u = A.col(j);
                if (!HouseholderUnReflect(R.ref(j,j),u,beta)) {
#ifdef TMV_NO_THROW
                    std::cerr<<"Bad QR Downdate\n";
                    exit(1);
#else
                    throw BadQRDowndate<T>(R,A);
#endif
                }
                M1r Rj = R.row(j,j+1,N);
                M2cr A2 = A.colRange(j+1,N);
                M3c temp = tempBase.col(0,0,N-j-1);

                TMVAssert(beta != RT(1));
                temp = beta * u.conjugate() * A2;
                Rj += temp;
                Rj /= T(1) - beta;
                temp += beta * Rj;
                A2 -= u ^ temp;
            }
        }
    };

    // algo 27: Recursive algorithm
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<27,cs,rs,M1,M2>
    {
        static void call(M1& R, M2& A)
        {
            typedef typename M1::value_type T;

            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
#ifdef PRINTALGO_QR
            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
            std::cout<<"QRDowndate algo 27: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename MCopyHelper<T,UpperTri,rs,rs>::type Ztype;
            Ztype Z = MatrixSizer<T>(N,N);
            typename Ztype::view_type Zv = Z.view();

#ifndef TMV_NO_THROW
            try {
#endif
                RecursiveQRDowndate(R,A,Zv,false);
#ifndef TMV_NO_THROW
            } catch (BadQRDowndate<T>) {
                throw BadQRDowndate<T>(R,A);
            }
#endif
        }
    };

    // algo 31: Decide which algorithm to use from runtime size
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<31,cs,rs,M1,M2>
    {
        static void call(M1& R, M2& A)
        {
            typedef typename M1::value_type T;

            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
            const ptrdiff_t M = cs==Unknown ? A.colsize() : cs;
#ifdef PRINTALGO_QR
            std::cout<<"QRDowndate algo 31: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const int algo27 = 
                (rs == Unknown || rs <= 128) ? 27 : 0;
            const int algo22 = 
                (rs == Unknown || rs > 128) ? 22 : 0;
            const ptrdiff_t l2cache = TMV_L2_CACHE*1024/sizeof(T);

            if (N*(M+N) <= l2cache)
                QRDowndate_Helper<11,cs,rs,M1,M2>::call(R,A);
            else if (N <= 128)
                QRDowndate_Helper<algo27,cs,rs,M1,M2>::call(R,A);
            else
                QRDowndate_Helper<algo22,cs,rs,M1,M2>::call(R,A);
        }
    };

    // algo 90: call InstQR_Downdate
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<90,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& R, M2& A)
        { InstQR_Downdate(R.xView(),A.xView()); }
    };

    // algo 97: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<97,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& R, M2& A)
        {
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c Rc = R.conjugate();
            M2c Ac = A.conjugate();
            QRDowndate_Helper<-2,cs,rs,M1c,M2c>::call(Rc,Ac);
        }
    };

    // algo -4: No copies or branches
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<-4,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& R, M2& A)
        {
            typedef typename M1::value_type T;
            const ptrdiff_t csrs = IntTraits2<cs,rs>::prod;
            const ptrdiff_t rsrs = IntTraits<IntTraits2<rs,rs>::prod>::halfS;
            const ptrdiff_t totmem = IntTraits2<csrs,rsrs>::sum;
            const ptrdiff_t l2cache = TMV_L2_CACHE*1024/sizeof(T);
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                rs == Unknown ? 31 :
                rsrs > l2cache ? (rs <= 128 ? 27 : 22) :
                cs == Unknown ? 31 : 
                totmem <= l2cache ? 11 :
                rs <= 128 ? 27 : 22;
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRDowndate: \n";
            std::cout<<"R = "<<TMV_Text(R)<<std::endl;
            std::cout<<"A = "<<TMV_Text(A)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            QRDowndate_Helper<algo,cs,rs,M1,M2>::call(R,A);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<-3,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& R, M2& A)
        {
            const int algo = -4;
            QRDowndate_Helper<algo,cs,rs,M1,M2>::call(R,A);
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<-2,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& R, M2& A)
        {
            typedef typename M1::value_type T;
            const bool inst = 
                (cs == Unknown || cs > 16 || cs == 1) &&
                (rs == Unknown || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            QRDowndate_Helper<algo,cs,rs,M1,M2>::call(R,A);
        }
    };

    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct QRDowndate_Helper<-1,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& R, M2& A)
        { QRDowndate_Helper<-2,cs,rs,M1,M2>::call(R,A); }
    };

    template <class M1, class M2>
    inline void InlineQR_Downdate(
        BaseMatrix_Tri_Mutable<M1>& R, BaseMatrix_Rec_Mutable<M2>& A)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        TMVStaticAssert((Traits2<T1,T2>::sametype));
        TMVStaticAssert((Sizes<M2::_rowsize,M1::_size>::same));
        TMVStaticAssert(M1::_upper);
        TMVStaticAssert(!M1::_unit);
        TMVAssert(A.rowsize() == R.size());
        TMVAssert(!R.isunit());
        const ptrdiff_t cs = M2::_colsize;
        const ptrdiff_t rs = Sizes<M2::_rowsize,M1::_size>::size;
        typedef typename M1::nonunitdiag_type::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) Rv = R.viewAsNonUnitDiag().cView();
        TMV_MAYBE_REF(M2,M2v) Av = A.cView();
        QRDowndate_Helper<-3,cs,rs,M1v,M2v>::call(Rv,Av);
    }

    template <class M1, class M2>
    inline void QR_Downdate(
        BaseMatrix_Tri_Mutable<M1>& R, BaseMatrix_Rec_Mutable<M2>& A)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        TMVStaticAssert((Traits2<T1,T2>::sametype));
        TMVStaticAssert((Sizes<M2::_rowsize,M1::_size>::same));
        TMVStaticAssert(M1::_upper);
        TMVStaticAssert(!M1::_unit);
        TMVAssert(A.rowsize() == R.size());
        TMVAssert(!R.isunit());
        const ptrdiff_t cs = M2::_colsize;
        const ptrdiff_t rs = Sizes<M2::_rowsize,M1::_size>::size;
        typedef typename M1::nonunitdiag_type::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) Rv = R.viewAsNonUnitDiag().cView();
        TMV_MAYBE_REF(M2,M2v) Av = A.cView();
        QRDowndate_Helper<-2,cs,rs,M1v,M2v>::call(Rv,Av);
    }


    // Allow views as an argument by value (for convenience)
    template <class T, int A1, int A2>
    inline void QR_Downdate(
        UpperTriMatrixView<T,A1> R, MatrixView<T,A2> A)
    {
        typedef UpperTriMatrixView<T,A1> M1;
        typedef MatrixView<T,A2> M2;
        QR_Downdate(
            static_cast<BaseMatrix_Tri_Mutable<M1>&>(R),
            static_cast<BaseMatrix_Rec_Mutable<M2>&>(A)); 
    }

    template <class T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t Si1, ptrdiff_t Sj1, int A1, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
    inline void QR_Downdate(
        SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> R,
        SmallMatrixView<T,M,N,Si2,Sj2,A2> A)
    {
        typedef SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> M1;
        typedef SmallMatrixView<T,M,N,Si2,Sj2,A2> M2;
        QR_Downdate(
            static_cast<BaseMatrix_Tri_Mutable<M1>&>(R),
            static_cast<BaseMatrix_Rec_Mutable<M2>&>(A));
    }

    template <class T, ptrdiff_t M, ptrdiff_t N, int A1, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
    inline void QR_Downdate(
        UpperTriMatrixView<T,A1> R,
        SmallMatrixView<T,M,N,Si2,Sj2,A2> A)
    {
        typedef UpperTriMatrixView<T,A1> M1;
        typedef SmallMatrixView<T,M,N,Si2,Sj2,A2> M2;
        QR_Downdate(
            static_cast<BaseMatrix_Tri_Mutable<M1>&>(R),
            static_cast<BaseMatrix_Rec_Mutable<M2>&>(A));
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si1, ptrdiff_t Sj1, int A1, int A2>
    inline void QR_Downdate(
        SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> R,
        MatrixView<T,A2> A)
    {
        typedef SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> M1;
        typedef MatrixView<T,A2> M2;
        QR_Downdate(
            static_cast<BaseMatrix_Tri_Mutable<M1>&>(R),
            static_cast<BaseMatrix_Rec_Mutable<M2>&>(A));
    }


    // Also versions with A as a single vector:
    template <class M1, class V2>
    inline void InlineQR_Downdate(
        BaseMatrix_Tri_Mutable<M1>& R, BaseVector_Mutable<V2>& A)
    {
        typedef typename M1::value_type T1;
        typedef typename V2::value_type T2;
        TMVStaticAssert((Traits2<T1,T2>::sametype));
        TMVStaticAssert((Sizes<V2::_size,M1::_size>::same));
        TMVStaticAssert(M1::_upper);
        TMVStaticAssert(!M1::_unit);
        TMVAssert(A.size() == R.size());
        TMVAssert(!R.isunit());
        const ptrdiff_t cs = 1;
        const ptrdiff_t rs = Sizes<V2::_size,M1::_size>::size;
        typedef typename M1::nonunitdiag_type::cview_type M1v;
        typedef typename MViewHelper<T2,Rec,1,rs,Unknown,V2::_step>::type V2v;
        TMV_MAYBE_REF(M1,M1v) Rv = R.viewAsNonUnitDiag().cView();
        V2v Av = RowVectorViewOf(A);
        QRDowndate_Helper<-3,cs,rs,M1v,V2v>::call(Rv,Av);
    }

    template <class M1, class V2>
    inline void QR_Downdate(
        BaseMatrix_Tri_Mutable<M1>& R, BaseVector_Mutable<V2>& A)
    {
        typedef typename M1::value_type T1;
        typedef typename V2::value_type T2;
        TMVStaticAssert((Traits2<T1,T2>::sametype));
        TMVStaticAssert((Sizes<V2::_size,M1::_size>::same));
        TMVStaticAssert(M1::_upper);
        TMVStaticAssert(!M1::_unit);
        TMVAssert(A.size() == R.size());
        TMVAssert(!R.isunit());
        const ptrdiff_t cs = 1;
        const ptrdiff_t rs = Sizes<V2::_size,M1::_size>::size;
        typedef typename M1::nonunitdiag_type::cview_type M1v;
        typedef typename MViewHelper<T2,Rec,1,rs,Unknown,V2::_step>::type V2v;
        TMV_MAYBE_REF(M1,M1v) Rv = R.viewAsNonUnitDiag().cView();
        V2v Av = RowVectorViewOf(A);
        QRDowndate_Helper<-2,cs,rs,M1v,V2v>::call(Rv,Av);
    }

    template <class T, int A1, int A2>
    inline void QR_Downdate(
        UpperTriMatrixView<T,A1> R, VectorView<T,A2> A)
    {
        typedef UpperTriMatrixView<T,A1> M1;
        typedef VectorView<T,A2> V2;
        QR_Downdate(
            static_cast<BaseMatrix_Tri_Mutable<M1>&>(R),
            static_cast<BaseVector_Mutable<V2>&>(A)); 
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si1, ptrdiff_t Sj1, int A1, ptrdiff_t S2, int A2>
    inline void QR_Downdate(
        SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> R,
        SmallVectorView<T,N,S2,A2> A)
    {
        typedef SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> M1;
        typedef SmallVectorView<T,N,S2,A2> V2;
        QR_Downdate(
            static_cast<BaseMatrix_Tri_Mutable<M1>&>(R),
            static_cast<BaseVector_Mutable<V2>&>(A));
    }

    template <class T, ptrdiff_t N, int A1, ptrdiff_t S2, int A2>
    inline void QR_Downdate(
        UpperTriMatrixView<T,A1> R,
        SmallVectorView<T,N,S2,A2> A)
    {
        typedef UpperTriMatrixView<T,A1> M1;
        typedef SmallVectorView<T,N,S2,A2> V2;
        QR_Downdate(
            static_cast<BaseMatrix_Tri_Mutable<M1>&>(R),
            static_cast<BaseVector_Mutable<V2>&>(A));
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si1, ptrdiff_t Sj1, int A1, int A2>
    inline void QR_Downdate(
        SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> R,
        VectorView<T,A2> A)
    {
        typedef SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> M1;
        typedef VectorView<T,A2> V2;
        QR_Downdate(
            static_cast<BaseMatrix_Tri_Mutable<M1>&>(R),
            static_cast<BaseVector_Mutable<V2>&>(A));
    }




} // namespace tmv

#undef TMV_QR_BLOCKSIZE

#endif

