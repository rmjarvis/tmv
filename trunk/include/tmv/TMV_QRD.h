

//---------------------------------------------------------------------------
//
// This file contains the code for doing division using 
// QR Decomposition.
//
// The basic idea of an QR decomposition is that any 
// matrix A can be decomposed into a unitary matrix
// times an upper triangle matrix.
//
// A = Q R
//
// We do this using Householder transformations, which can
// be stored in a lower triangle matrix, thus other than the 
// diagonal, which they both need, we can store them in the 
// place of the original matrix. It is more convenient to 
// keep the diagonal of Q in place and take the diagonal of R
// separately.
//
// If R is not singular, the solution to A x = b is found from
// Q R x = b
// R x = Qt b
// which can be solved by back-substitutaion.
//
// If m > n, this does not actually give a solution to A x = b,
// since Q Qt != I  (Q is only column orthogonal if m > n.)
// But it does give the value of x which minimizes the 2-norm
// of (A x - b)
//
// The 2-norm of v is the square root of vt v, so
// |A x - b|^2 =
// (A x - b)t (A x - b) = 
// (xt At - bt) (A x - b) =
// xt At A x - bt A x - xt At b + bt b =
// xt Rt Qt Q R x - bt Q R x - xt Rt Qt b + bt b =
// |R x - Qt b|^2 + |b|^2 - |Qt b|^2
// Clearly the x which minimizes this is the solution of R x = Qt b.
//
// If R is singular, then you need QRP Decomposition (see TMV_QRPD.h).
//


#ifndef TMV_QRD_H
#define TMV_QRD_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_BaseVector.h"
#include "TMV_Divider.h"
#include "TMV_PackedQ.h"
#include "TMV_Array.h"

namespace tmv {

    // In TMV_QRDecompose.h
    template <class M, class V>
    inline void QR_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta);

    // In TMV_QRInverse.h
    template <class M1, class V1, class M2>
    inline void QR_Inverse(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& minv);
    template <class M1, class V1, class M2>
    inline void QR_InverseATA(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& ata);

    // In TMV_QRDiv.h
    template <class M1, class V1, class M2, class M3>
    inline void QR_Solve(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <class M1, class V1, class V2, class V3>
    inline void QR_Solve(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <class M1, class V1, class M2, class M3>
    inline void QR_SolveTranspose(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <class M1, class V1, class V2, class V3>
    inline void QR_SolveTranspose(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <class M1, class V1, class M2>
    inline void QR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class V1, class V2>
    inline void QR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2);
    template <class M1, class V1, class M2>
    inline void QR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class V1, class V2>
    inline void QR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2);




    // The point of the Impl class here is to implement the transfer of 
    // ownership copy semantics.
    // It also differentiates between small and non-small implementations.
    template <bool small, class M>
    struct QRD_Impl;

    template <class M>
    class QRD 
    {
    public :

        typedef typename M::value_type T;
        typedef typename M::real_type RT;
        typedef typename M::complex_type CT;
        typedef typename M::float_type FT;
        typedef typename M::zfloat_type ZFT;

        // This next bit finds the storage type to use for the qr matrix
        // regardless of what kind of matrix M is.  e.g. this should
        // work even if M is a TriMatrix or a BandMatrix, etc.
        enum { cs = IntTraits2<M::_colsize,M::_rowsize>::max };
        enum { rs = IntTraits2<M::_colsize,M::_rowsize>::min };

        enum { small = (
                M::_colsize != Unknown && M::_rowsize != Unknown
                && M::_colsize <= 32 && M::_rowsize <= 32 ) };

        typedef typename QRD_Impl<small,M>::qrx_type qrx_type;
        typedef typename QRD_Impl<small,M>::beta_type beta_type;

        typedef typename qrx_type::const_view_type getqr_type;
        typedef PackedQ<qrx_type,beta_type> getq_type;
        typedef typename qrx_type::const_uppertri_type getr_type;
        typedef const beta_type& getbeta_type;

        //
        // Constructors
        //

        // Sets up the internal storage and does the decomposition.
        template <class M2>
        QRD(const BaseMatrix<M2>& A, bool _inplace=false);

        // The copy constructor has transfer of ownership semantics.
        // This way an QRD object can be returned by value, and the 
        // copy is cheap.  I don't think there is any reason to use
        // a more sophisticated technique like shared_ptr or something
        // similar, since there are no non-const methods.
        QRD(const QRD<M>& rhs);

        // Clean up the internal storage
        ~QRD();


        //
        // Division: (not in place)
        // 

        template <class M1, class M2>
        void doSolve(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class M1, class M2>
        void solve(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const
        { doSolve(m1,m2); }
        template <class M1, class T2, int A2>
        void solve(const BaseMatrix<M1>& m1, MatrixView<T2,A2> m2) const
        { doSolve(m1,m2); }
        template <class M1, class T2, int M2, int N2, int Si2, int Sj2, int A2>
        void solve(
            const BaseMatrix<M1>& m1, SmallMatrixView<T2,M2,N2,Si2,Sj2,A2> m2) const
        { doSolve(m1,m2); }

        template <class V1, class V2>
        void doSolve(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const;
        template <class V1, class V2>
        void solve(const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const
        { doSolve(v1,v2); }
        template <class V1, class T2, int A2>
        void solve(const BaseVector<V1>& v1, VectorView<T2,A2> v2) const
        { doSolve(v1,v2); }
        template <class V1, class T2, int N2, int S2, int A2>
        void solve(const BaseVector<V1>& v1, SmallVectorView<T2,N2,S2,A2> v2) const
        { doSolve(v1,v2); }

        template <class M1, class M2>
        void doSolveTranspose(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class M1, class M2>
        void solveTranspose(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const
        { doSolveTranspose(m1,m2); }
        template <class M1, class T2, int A2>
        void solveTranspose(
            const BaseMatrix<M1>& m1, MatrixView<T2,A2> m2) const
        { doSolveTranspose(m1,m2); }
        template <class M1, class T2, int M2, int N2, int Si2, int Sj2, int A2>
        void solveTranspose(
            const BaseMatrix<M1>& m1, SmallMatrixView<T2,M2,N2,Si2,Sj2,A2> m2) const
        { doSolveTranspose(m1,m2); }

        template <class V1, class V2>
        void doSolveTranspose(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const;
        template <class V1, class V2>
        void solveTranspose(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const
        { doSolveTranspose(v1,v2); }
        template <class V1, class T2, int A2>
        void solveTranspose(
            const BaseVector<V1>& v1, VectorView<T2,A2> v2) const
        { doSolveTranspose(v1,v2); }
        template <class V1, class T2, int N2, int S2, int A2>
        void solveTranspose(
            const BaseVector<V1>& v1, SmallVectorView<T2,N2,S2,A2> v2) const
        { doSolveTranspose(v1,v2); }

        //
        // Perform the division in place
        //

        template <class M2>
        void doSolveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class M2>
        void solveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
        { doSolveInPlace(m2); }
        template <class T2, int A2>
        void solveInPlace(MatrixView<T2,A2> m2) const
        { doSolveInPlace(m2); }
        template <class T2, int M2, int N2, int Si2, int Sj2, int A2>
        void solveInPlace(SmallMatrixView<T2,M2,N2,Si2,Sj2,A2> m2) const
        { doSolveInPlace(m2); }

        template <class V2>
        void doSolveInPlace(BaseVector_Mutable<V2>& v2) const;
        template <class V2>
        void solveInPlace(BaseVector_Mutable<V2>& v2) const
        { doSolveInPlave(v2); }
        template <class T2, int A2>
        void solveInPlace(VectorView<T2,A2> v2) const
        { doSolveInPlace(v2); }
        template <class T2, int N2, int S2, int A2>
        void solveInPlace(SmallVectorView<T2,N2,S2,A2> v2) const
        { doSolveInPlace(v2); }

        template <class M2>
        void doSolveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class M2>
        void solveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
        { doSolveTransposeInPlace(m2); }
        template <class T2, int A2>
        void solveTransposeInPlace(MatrixView<T2,A2> m2) const
        { doSolveTransposeInPlace(m2); }
        template <class T2, int M2, int N2, int Si2, int Sj2, int A2>
        void solveTransposeInPlace(SmallMatrixView<T2,M2,N2,Si2,Sj2,A2> m2) const
        { doSolveTransposeInPlace(m2); }

        template <class V2>
        void doSolveTransposeInPlace(BaseVector_Mutable<V2>& v2) const;
        template <class V2>
        void solveTransposeInPlace(BaseVector_Mutable<V2>& v2) const
        { doSolveTransposeInPlave(v2); }
        template <class T2, int A2>
        void solveTransposeInPlace(VectorView<T2,A2> v2) const
        { doSolveTransposeInPlace(v2); }
        template <class T2, int N2, int S2, int A2>
        void solveTransposeInPlace(SmallVectorView<T2,N2,S2,A2> v2) const
        { doSolveTransposeInPlace(v2); }

        //
        // Determinant
        //

        T det() const;
        FT logDet(ZFT* sign) const;
        bool isSingular() const;


        //
        // Inverse
        //

        template <class M2>
        void doMakeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const;
        template <class M2>
        void makeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const
        { doMakeInverse(minv); }
        template <class T2, int A2>
        void MakeInverse(MatrixView<T2,A2> minv) const
        { doMakeInverse(minv); }
        template <class T2, int M2, int N2, int Si2, int Sj2, int A2>
        void MakeInverse(SmallMatrixView<T2,M2,N2,Si2,Sj2,A2> minv) const
        { doMakeInverse(minv); }


        //
        // InverseATA
        //

        template <class M2>
        void doMakeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const;
        template <class M2>
        void makeInverseATA(BaseMatrix_Rec_Mutable<M2>& minv) const
        { doMakeInverseATA(minv); }
        template <class T2, int A2>
        void MakeInverseATA(MatrixView<T2,A2> minv) const
        { doMakeInverseATA(minv); }
        template <class T2, int M2, int N2, int Si2, int Sj2, int A2>
        void MakeInverseATA(SmallMatrixView<T2,M2,N2,Si2,Sj2,A2> minv) const
        { doMakeInverseATA(minv); }


        // 
        // Condition (kappa_inf)
        //

        RT condition(RT normInf) const;


        //
        // Access Decomposition
        //

        bool isTrans() const;
        getq_type getQ() const;
        getr_type getR() const;
        getqr_type getQR() const;
        getbeta_type getBeta() const;

        bool preferInPlace() const { return false; }

    private :

        // mutable so the normal copy constructor with the argument
        // const QRD<M>& can release the memory.
        mutable std::auto_ptr<QRD_Impl<small,M> > pimpl;

        int colsize() const;
        int rowsize() const;

        // op= not allowed.
        QRD<M>& operator=(const QRD<M>&);
    };


    template <class T>
    class InstQRD :
        public QRD<Matrix<T,ColMajor> >,
        public Divider<T>
    {
    public :
        typedef QRD<Matrix<T,ColMajor> > base;
        typedef typename base::RT RT;
        typedef typename base::CT CT;
        typedef typename base::FT FT;
        typedef typename base::ZFT ZFT;

        // Sets up the internal storage and does the decomposition.
        template <int C>
        InstQRD(const ConstMatrixView<T,C>& A, bool _inplace=false);
        InstQRD(const InstQRD<T>& rhs);
        ~InstQRD();

        // These are the virtual functions from the Divider base class.
        void doSolveInPlace(MatrixView<RT> m2) const;
        void doSolveInPlace(MatrixView<CT> m2) const;
        void doSolveInPlace(MatrixView<CT,Conj> m2) const;
        void doSolveInPlace(VectorView<RT> v2) const;
        void doSolveInPlace(VectorView<CT> v2) const;
        void doSolveInPlace(VectorView<CT,Conj> v2) const;

        void doSolveTransposeInPlace(MatrixView<RT> m2) const;
        void doSolveTransposeInPlace(MatrixView<CT> m2) const;
        void doSolveTransposeInPlace(MatrixView<CT,Conj> m2) const;
        void doSolveTransposeInPlace(VectorView<RT> v2) const;
        void doSolveTransposeInPlace(VectorView<CT> v2) const;
        void doSolveTransposeInPlace(VectorView<CT,Conj> v2) const;

        void doSolve(
            const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const;
        void doSolve(
            const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const;
        void doSolve(
            const ConstMatrixView<RT>& m1, MatrixView<CT,Conj> m2) const;
        void doSolve(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const;
        void doSolve(
            const ConstMatrixView<CT>& m1, MatrixView<CT,Conj> m2) const;
        void doSolve(
            const ConstMatrixView<CT,Conj>& m1, MatrixView<CT> m2) const;
        void doSolve(
            const ConstMatrixView<CT,Conj>& m1, MatrixView<CT,Conj> m2) const;
        void doSolve(
            const ConstVectorView<RT>& v1, VectorView<RT> v2) const;
        void doSolve(
            const ConstVectorView<RT>& v1, VectorView<CT> v2) const;
        void doSolve(
            const ConstVectorView<RT>& v1, VectorView<CT,Conj> v2) const;
        void doSolve(
            const ConstVectorView<CT>& v1, VectorView<CT> v2) const;
        void doSolve(
            const ConstVectorView<CT>& v1, VectorView<CT,Conj> v2) const;
        void doSolve(
            const ConstVectorView<CT,Conj>& v1, VectorView<CT> v2) const;
        void doSolve(
            const ConstVectorView<CT,Conj>& v1, VectorView<CT,Conj> v2) const;

        void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<CT,Conj> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<CT>& m1, MatrixView<CT,Conj> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<CT,Conj>& m1, MatrixView<CT> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<CT,Conj>& m1, MatrixView<CT,Conj> m2) const;
        void doSolveTranspose(
            const ConstVectorView<RT>& v1, VectorView<RT> v2) const;
        void doSolveTranspose(
            const ConstVectorView<RT>& v1, VectorView<CT> v2) const;
        void doSolveTranspose(
            const ConstVectorView<RT>& v1, VectorView<CT,Conj> v2) const;
        void doSolveTranspose(
            const ConstVectorView<CT>& v1, VectorView<CT> v2) const;
        void doSolveTranspose(
            const ConstVectorView<CT>& v1, VectorView<CT,Conj> v2) const;
        void doSolveTranspose(
            const ConstVectorView<CT,Conj>& v1, VectorView<CT> v2) const;
        void doSolveTranspose(
            const ConstVectorView<CT,Conj>& v1, VectorView<CT,Conj> v2) const;

        T det() const;
        FT logDet(ZFT* sign) const;
        bool isSingular() const;

        void doMakeInverse(MatrixView<RT> minv) const;
        void doMakeInverse(MatrixView<CT> minv) const;
        void doMakeInverse(MatrixView<CT,Conj> minv) const;

        void doMakeInverseATA(MatrixView<RT> ata) const;
        void doMakeInverseATA(MatrixView<CT> ata) const;
        void doMakeInverseATA(MatrixView<CT,Conj> ata) const;
        
        RT condition(RT normInf) const;
        bool preferInPlace() const;

    private :
        // op= not allowed.
        InstQRD<T>& operator=(const InstQRD<T>&);
    };

    template <bool isvalid, bool istrans>
    struct QRHelper;

    template <>
    struct QRHelper<true,false>
    {
        template <class M1, class V1, class M2, class M3>
        static TMV_INLINE void solve(
            const M1& QRx, const V1& beta, const M2& m2, M3& m3)
        { QR_Solve(QRx,beta,0,QRx.rowsize(),m2,m3); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void solveInPlace(
            const M1& QRx, const V1& beta, M2& m2)
        { QR_SolveInPlace(QRx,beta,0,QRx.rowsize(),m2); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverse(
            const M1& QRx, const V1& beta, M2& m2)
        { QR_Inverse(QRx,beta,0,QRx.rowsize(),m2); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverseATA(
            const M1& QRx, const V1& beta, M2& m2)
        { QR_InverseATA(QRx,beta,0,QRx.rowsize(),m2); }
    };
    template <>
    struct QRHelper<true,true>
    {
        template <class M1, class V1, class M2, class M3>
        static TMV_INLINE void solve(
            const M1& QRx, const V1& beta, const M2& m2, M3& m3)
        { QR_SolveTranspose(QRx,beta,0,QRx.rowsize(),m2,m3); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void solveInPlace(
            const M1& QRx, const V1& beta, M2& m2)
        { QR_SolveTransposeInPlace(QRx,beta,0,QRx.rowsize(),m2); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverse(
            const M1& QRx, const V1& beta, M2& m2)
        { 
            typename M2::transpose_type m2t = m2.transpose();
            QR_Inverse(QRx,beta,0,QRx.rowsize(),m2t);
        }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverseATA(
            const M1& QRx, const V1& beta, M2& m2)
        { QR_InverseATA(QRx,beta,0,QRx.rowsize(),m2); }
    };
    template <bool istrans>
    struct QRHelper<false,istrans>
    {
        template <class M1, class V1, class M2, class M3>
        static TMV_INLINE void solve(const M1& , const V1& , const M2& , M3& ) 
        { TMVAssert(false && "Calling invalid QRHelper::solve\n"); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void solveInPlace(const M1& , const V1& , M2& )
        { TMVAssert(false && "Calling invalid QRHelper::solveInPlace\n"); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverse(const M1& , const V1& , M2& )
        { TMVAssert(false && "Calling invalid QRHelper::makeInverse\n"); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverseATA(const M1& , const V1& , M2& ) 
        { TMVAssert(false && "Calling invalid QRHelper::makeInverseATA\n"); }
    };

    template <class M>
    struct QRD_Impl<true,M> 
    // small = true, so cs,rs both known 
    {
        typedef typename M::real_type RT;
        typedef typename M::value_type T;
        enum { cs1 = M::_colsize };
        enum { rs1 = M::_rowsize };
        enum { istrans = int(cs1) < int(rs1) };
        enum { cs = IntTraits2<cs1,rs1>::max };
        enum { rs = IntTraits2<cs1,rs1>::min };
        enum { A = (istrans ? RowMajor : ColMajor) | NoAlias };
        typedef typename MCopyHelper<T,Rec,cs1,rs1,A>::type Mc;
        typedef typename TypeSelect< istrans ,
                typename Mc::transpose_type ,
                typename Mc::view_type >::type qrx_type;
        typedef typename VCopyHelper<RT,rs>::type beta_type;

        template <class M2>
        QRD_Impl(const BaseMatrix<M2>& A, bool ) : 
            QRx(Maybe<istrans>::transposeview(SmallQRx) )
        {
            TMVStaticAssert(M::_colsize != Unknown);
            TMVStaticAssert(M::_rowsize != Unknown);
            TMVAssert(A.colsize() == istrans ? int(rs) : int(cs));
            TMVAssert(A.rowsize() == istrans ? int(cs) : int(rs));
            //std::cout<<"QRD_Impl small\n";
            //std::cout<<"istrans = "<<istrans<<std::endl;
            //std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            //std::cout<<"A = "<<A<<std::endl;
            //std::cout<<"SmallQRx = "<<SmallQRx<<std::endl;
            //std::cout<<"QRx = "<<QRx<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            SmallQRx = A;
            //std::cout<<"SmallQRx => "<<SmallQRx<<std::endl;
            //std::cout<<"QRx => "<<QRx<<std::endl;
            QR_Decompose(QRx,beta);
            //std::cout<<"SmallQRx => "<<SmallQRx<<std::endl;
            //std::cout<<"QRx => "<<QRx<<std::endl;
            //std::cout<<"beta => "<<beta<<std::endl;
        }
        template <class M2, class M3>
        void solve(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            QRHelper<isvalid,istrans>::solve(QRx,beta,m2,m3);
        }
        template <class M2, class M3>
        void solveTranspose(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            QRHelper<isvalid,!istrans>::solve(QRx,beta,m2,m3);
        }
        template <class M2>
        void solveInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRHelper<isvalid,istrans>::solveInPlace(QRx,beta,m2);
        }
        template <class M2>
        void solveTransposeInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRHelper<isvalid,!istrans>::solveInPlace(QRx,beta,m2);
        }
        template <class M2>
        void makeInverse(M2& minv)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRHelper<isvalid,istrans>::makeInverse(QRx,beta,minv);
        }
        template <class M2>
        void makeInverseATA(M2& ata)
        { 
            const bool isvalid = M::isreal || M2::iscomplex;
            QRHelper<isvalid,istrans>::makeInverseATA(QRx,beta,ata);
        }

        Mc SmallQRx;
        qrx_type QRx;
        beta_type beta;
    };
    
    template <class M>
    struct QRD_Impl<false,M>
    {
        typedef typename M::value_type T;
        typedef typename M::real_type RT;
        enum { cs1 = M::_colsize };
        enum { rs1 = M::_rowsize };
        enum { knownsizes = cs1 != Unknown && rs1 != Unknown };
        enum { istrans1 = knownsizes && cs1 < int(rs1) };
        enum { cs = IntTraits2<cs1,rs1>::max };
        enum { rs = IntTraits2<cs1,rs1>::min };
        typedef typename MViewHelper<T,Rec,cs,rs,1,Unknown,NoAlias>::type qrx_type;
        typedef Vector<RT> beta_type;

        template <class M2>
        QRD_Impl(const BaseMatrix_Rec<M2>& A, bool _inplace) :
            // if A is short, need to transpose
            istrans(knownsizes ? istrans1 : A.colsize() < A.rowsize()),
            // inplace only if it works with a ColMajor QRx object
            inplace( _inplace && 
                     ((A.iscm() && !istrans) || (A.isrm() && istrans)) ),
            // Aptr is the pointer to new storage if any
            Aptr( inplace ? 0 : A.rowsize()*A.colsize() ),
            // QRx views this memory as the QR matrix
            QRx(
                inplace ? A.nonConst().ptr() : Aptr.get() , // ptr
                istrans ? A.rowsize() : A.colsize() ,  // colsize
                istrans ? A.colsize() : A.rowsize() ,  // rowsize
                inplace ? (istrans ? A.stepi() : A.stepj()) : 1 , // stepi
                ( inplace ? (istrans ? A.stepj() : A.stepi()) : 
                  (istrans ? A.rowsize() : A.colsize()) ) // stepj
            ),
            beta(istrans ? A.colsize() : A.rowsize())
            {
                if (!inplace) {
                    if (istrans) {
                        typename qrx_type::transpose_type QRxt = QRx.transpose();
                        Maybe<!knownsizes||istrans1>::assignTo(A,QRxt);
                    } else {
                        Maybe<!knownsizes||!istrans1>::assignTo(A,QRx);
                    }
                } else {
                    Maybe<M2::_conj>::conjself(QRx);
                }
                QR_Decompose(QRx,beta);
            }

        // If A is not a BaseMatrix_Rec, can't do it in place.
        template <class M2>
        QRD_Impl(const BaseMatrix<M2>& A, bool _inplace) :
            istrans(A.colsize() < A.rowsize()), inplace(false),
            Aptr( A.rowsize()*A.colsize() ),
            QRx(
                Aptr.get() , // ptr
                istrans ? A.rowsize() : A.colsize() ,  // colsize
                istrans ? A.colsize() : A.rowsize() ,  // rowsize
                1 , // stepi
                (istrans ? A.rowsize() : A.colsize()) // stepj
            ),
            beta(istrans ? A.colsize() : A.rowsize())
            {
                if (istrans) {
                    typename qrx_type::transpose_type QRxt = QRx.transpose();
                    Maybe<!knownsizes||istrans1>::assignTo(A,QRxt);
                } else {
                    Maybe<!knownsizes||!istrans1>::assignTo(A,QRx);
                }
                QR_Decompose(QRx,beta);
            }

        template <class M2, class M3>
        void solve(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRHelper<isvalid1,true>::solve(QRx,beta,m2,m3);
            else
                QRHelper<isvalid2,false>::solve(QRx,beta,m2,m3);
        }
        template <class M2, class M3>
        void solveTranspose(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRHelper<isvalid1,false>::solve(QRx,beta,m2,m3);
            else
                QRHelper<isvalid2,true>::solve(QRx,beta,m2,m3);
        }
        template <class M2>
        void solveInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRHelper<isvalid1,true>::solveInPlace(QRx,beta,m2);
            else
                QRHelper<isvalid2,false>::solveInPlace(QRx,beta,m2);
        }
        template <class M2>
        void solveTransposeInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRHelper<isvalid1,false>::solveInPlace(QRx,beta,m2);
            else
                QRHelper<isvalid2,true>::solveInPlace(QRx,beta,m2);
        }
        template <class M2>
        void makeInverse(M2& minv)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRHelper<isvalid1,true>::makeInverse(QRx,beta,minv);
            else
                QRHelper<isvalid2,false>::makeInverse(QRx,beta,minv);
        }
        template <class M2>
        void makeInverseATA(M2& ata)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            // For this one, it doesn't actually matter if istrans = true
            // So just use trans = false arbitrarily.
            QRHelper<isvalid,false>::makeInverseATA(QRx,beta,ata);
        }


        const bool istrans;
        const bool inplace;
        AlignedArray<typename M::value_type> Aptr;
        qrx_type QRx;
        beta_type beta;
    };

    template <class M> template <class M2>
    QRD<M>::QRD(const BaseMatrix<M2>& A, bool inplace) :
        pimpl(new QRD_Impl<small,M>(A.mat(),inplace)) {}

    template <class M>
    QRD<M>::QRD(const QRD<M>& rhs) : pimpl(rhs.pimpl.release()) {}

    template <class M>
    QRD<M>::~QRD() {}

    template <class M> template <class M2, class M3>
    void QRD<M>::doSolve(
        const BaseMatrix<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3) const
    {
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M::_rowsize>::same));
        TMVAssert(m2.rowsize() == m3.rowsize());
        TMVAssert(m2.colsize() == colsize());
        TMVAssert(m3.colsize() == rowsize());
        pimpl->solve(m2.mat(),m3.mat());
    }

    template <class M> template <class V2, class V3>
    void QRD<M>::doSolve(
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<V3::_size,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(v3.size() == rowsize());
        pimpl->solve(v2.vec(),v3.vec());
    }

    template <class M> template <class M2>
    void QRD<M>::doSolveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == colsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveInPlace(m2.mat());
    }

    template <class M> template <class V2>
    void QRD<M>::doSolveInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveInPlace(v2.vec());
    }

    template <class M> template <class M2, class M3>
    void QRD<M>::doSolveTranspose(
        const BaseMatrix<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3) const
    {
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M2::_colsize,M::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M::_colsize>::same));
        TMVAssert(m2.rowsize() == m3.rowsize());
        TMVAssert(m2.colsize() == rowsize());
        TMVAssert(m3.colsize() == colsize());
        pimpl->solveTranspose(m2.mat(),m3.mat());
    }

    template <class M> template <class V2, class V3>
    void QRD<M>::doSolveTranspose(
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_rowsize>::same));
        TMVStaticAssert((Sizes<V3::_size,M::_colsize>::same));
        TMVAssert(v2.size() == rowsize());
        TMVAssert(v3.size() == colsize());
        pimpl->solveTranspose(v2.vec(),v3.vec());
    }

    template <class M> template <class M2>
    void QRD<M>::doSolveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == rowsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveTransposeInPlace(m2.mat());
    }

    template <class M> template <class V2>
    void QRD<M>::doSolveTransposeInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveTransposeInPlace(v2.vec());
    }

    template <class M>
    typename M::value_type QRD<M>::det() const
    { return typename M::real_type(CalculateDetQ(pimpl->beta)) * getR().det(); }

    template <class M>
    typename M::float_type QRD<M>::logDet(typename M::zfloat_type* sign) const
    {
        typename M::float_type ret = getR().logDet(sign);
        if (sign) *sign *= typename M::float_type(CalculateDetQ(pimpl->beta));
        return ret;
    }                  

    template <class M>
    bool QRD<M>::isSingular() const 
    { return getR().isSingular(); }

    template <class M> template <class M2>
    void QRD<M>::doMakeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const
    {
        TMVStaticAssert((Sizes<M::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M::_rowsize,M2::_colsize>::same));
        TMVAssert(minv.colsize() == rowsize());
        TMVAssert(minv.rowsize() == colsize());
        pimpl->makeInverse(minv.mat());
    }

    template <class M> template <class M2>
    void QRD<M>::doMakeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M2::_rowsize>::same));
        TMVAssert(ata.rowsize() == rowsize());
        TMVAssert(ata.colsize() == rowsize());
        pimpl->makeInverseATA(ata.mat());
    }

    template <class M>
    bool QRD<M>::isTrans() const 
    { return pimpl->istrans; }

    template <class M>
    typename QRD<M>::getq_type QRD<M>::getQ() const 
    { return typename QRD<M>::getq_type(pimpl->QRx,pimpl->beta); }

    template <class M>
    typename QRD<M>::getr_type QRD<M>::getR() const 
    { return pimpl->QRx.upperTri(); }

    template <class M>
    typename QRD<M>::getqr_type QRD<M>::getQR() const 
    { return pimpl->QRx; }

    template <class M>
    typename QRD<M>::getbeta_type QRD<M>::getBeta() const 
    { return pimpl->beta; }

    template <class M>
    typename M::real_type QRD<M>::condition(RT normInf) const 
    {
        // FIXME: This is a placeholder until I write the real function.
        // Make sure to do this before releasing the code!
        // See page 129 of Golub and van Loan.
        //
        // This produces the exact right answer, but it is way too slow!
        // The GvL algorithm is order mn.  This is order mn^2.
        Matrix<T> minv(rowsize(),colsize());
        if (isSingular()) {
            return normInf / TMV_Epsilon<RT>();
        } else {
            makeInverse(minv);
            return normInf * minv.normInf();
        }
    }

    template <class M>
    int QRD<M>::colsize() const
    { return pimpl->istrans ? pimpl->QRx.rowsize() : pimpl->QRx.colsize(); }

    template <class M>
    int QRD<M>::rowsize() const
    { return pimpl->istrans ? pimpl->QRx.colsize() : pimpl->QRx.rowsize(); }


    template <class M, class M2>
    static bool CheckDecomp(
        const QRD<M>& qrd, const BaseMatrix_Calc<M2>& m, std::ostream* fout=0)
    {
        typedef typename M2::real_type RT;
        //bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
        bool printmat = fout;
        if (printmat) {
            *fout << "QR:\n";
            if (qrd.isTrans()) *fout << m.transpose() << std::endl;
            else *fout << m << std::endl;
            *fout << "Q = "<<qrd.getQ()<<std::endl;
            *fout << "R = "<<qrd.getR()<<std::endl;
        }
        typename M::copy_type qr = qrd.getQ()*qrd.getR();
        if (printmat) {
            *fout << "QR = "<<qr<<std::endl;
        }
        RT nm = qrd.isTrans() ? Norm(qr-m.transpose()) : Norm(qr-m);
        nm /= Norm(qrd.getQ())*Norm(qrd.getR());
        if (fout) {
            *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<" <? ";
            *fout << RT(m.colsize())<<"*"<<TMV_Epsilon<RT>();
            *fout << " = "<<RT(m.colsize())*TMV_Epsilon<RT>()<<std::endl;
        }
        return nm < RT(m.colsize())*TMV_Epsilon<RT>();
    }


} // namespace mv


#endif
