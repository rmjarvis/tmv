

//---------------------------------------------------------------------------
//
// This file contains the code for doing division using 
// QRP Decomposition.
//
// In a QR Decomposition, if R is singular, then back-substitution will 
// fail in solving R x = Qt b.  However, there is a trick that still 
// gives us a valid least squares solution to A x = b.  
// This is to use column pivoting to obtain a decomposition:
// A = Q [ R11 R12 ] P
//       [  0   0  ]
// where R11 is upper triangular. 
//
// With this decomposition, 
// Q R P x = b
// R P x = Qt b
// Let z = P x (ie. we will solve for z first, then x = z / P = P * z)
// and c = Qt b
// Then R z = c
// Since R is singular, there is no exact solution, but for the least-squares
// problem, we just want to minimize |R z - c|.
// If z = [ z1 ] and c = [ c1 ] (with the same dimensional split as R11 and R12)
//        [ z2 ]         [ c2 ]
// then R z - c = [ R11 z1 + R12 z2 - c1 ]
//                [        -c2           ]
// The minimum norm will occur for any z2 if z1 = R11^-1 (c1 - R12 z2)
// So a solution can be found with z2 = 0.
// This solution then has minimum |A x - b| (ie. the least squares
// solution).  
//
// For dividing from the other side with m > n, we can always find a 
// possible solution (assuming R is non-singular):
//
// xt Q R P = bt 
// xt Q R = bt Pt
// xt Q = bt Pt R^-1  (This is done by front-substitution)
// xt = bt Pt R^-1 Qt
//
// This solution does satisfy the equation xt A = bt, but
// it will not be the solution with minimum |x|.  
// Use SVD for the min |x| solution.
//
// If R is singular, we of course have a problem with the front-substitution
// step.  This time, the QRP decomposition does not work quite as well.
// Take the front-substitution step (the Q and P steps don't affect 
// this minimization), and write R as above:
//
// zt R = ct
// [ z1t z2t ] [ R11  R12 ] = [ c1t c2t ]
//             [  0    0  ]
// [ (z1t R11) (z1t R12) ] = [ c1t c2t]
// |zt R - ct|^2 = |z1t R11 - c1t|^2 + |z1t R12 - c2t|^2
// We can set z2t = 0, since it is arbitrary.
// But it is not so clear what the correct z1t is in this case.
// We take z1t = c1t R11^-1, as an approximate solution, but point out 
// that this is not really correct. You should use SVD instead for 
// the correct least squares solution in this case.
//
// You can control whether QRP does a strict reordering so that the 
// diagonal elements of R are in decreasing order of absolute value 
// (which I call Strict QRP), or whether they are just reordered well
// enough to put the correct zeros at the end (which I call Loose QRP) using
// the function tmv::UseStrictQRP();
// The default value is false, which is faster, but possibly less accurate 
// for some matrices.
// To turn the strict algorithm back off, use tmv::UseStrictQRP(false);
//


#ifndef TMV_QRPD_H
#define TMV_QRPD_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_BaseVector.h"
#include "TMV_Divider.h"
#include "TMV_PackedQ.h"
#include "TMV_Array.h"
#include "TMV_Permutation.h"

namespace tmv {

    // In TMV_QRPDecompose.h
    template <class M, class V>
    static inline void QR_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta);
    template <class M1, class M2>
    static inline void QR_Decompose(
        BaseMatrix_Rec_Mutable<M1>& Q, BaseMatrix_Tri_Mutable<M2>& R);
    template <class M>
    static inline void QR_Decompose(BaseMatrix_Rec_Mutable<M>& m);

    // In TMV_QRInverse.h
    template <class M1, class V1, class M2>
    static inline void QR_Inverse(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& minv);
    template <class M1, class V1, class M2>
    static inline void QR_InverseATA(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& ata);

    // In TMV_QRDiv.h
    template <class M1, class V1, class M2, class M3>
    static inline void QR_Solve(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <class M1, class V1, class V2, class V3>
    static inline void QR_Solve(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <class M1, class V1, class M2, class M3>
    static inline void QR_SolveTranspose(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <class M1, class V1, class V2, class V3>
    static inline void QR_SolveTranspose(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <class M1, class V1, class M2>
    static inline void QR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class V1, class V2>
    static inline void QR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2);
    template <class M1, class V1, class M2>
    static inline void QR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class V1, class V2>
    static inline void QR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2);


    // 
    // Store the global variable for StrictQRP as a singleton to make
    // is thread safe and able to be done inline.
    //

    class QRP_StrictSingleton 
    {
    public:
        // Technically, I think this isn't thread safe, but I'd be 
        // pretty shocked if people were having multiple threads
        // call this funtion at the same time.
        static TMV_INLINE bool& inst() { 
            static bool strict;
            return strict;
        }
    private:
        QRP_StrictSingleton();
    };

    static TMV_INLINE void UseStrictQRP(bool newstrict=true)
    { QRP_StrictSingleton::inst() = newstrict; }



    // The point of the Impl class here is to implement the transfer of 
    // ownership copy semantics.
    // It also differentiates between small and non-small implementations.
    template <bool small, class M>
    struct QRPD_Impl;

    template <class M>
    class QRPD 
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
                M::_colsize != TMV_UNKNOWN && M::_rowsize != TMV_UNKNOWN
                && M::_colsize <= 32 && M::_rowsize <= 32 ) };

        typedef typename QRPD_Impl<small,M>::qrx_type qrx_type;
        typedef typename QRPD_Impl<small,M>::beta_type beta_type;

        typedef typename qrx_type::const_view_type getqr_type;
        typedef PackedQ<qrx_type,beta_type> getq_type;
        typedef typename qrx_type::const_uppertri_type getr_type;
        typedef const beta_type& getbeta_type;
        typedef const Permutation& getp_type;

        //
        // Constructors
        //

        // Sets up the internal storage and does the decomposition.
        template <class M2>
        QRPD(const BaseMatrix<M2>& A, bool _inplace=false);

        // The copy constructor has transfer of ownership semantics.
        // This way an QRPD object can be returned by value, and the 
        // copy is cheap.  I don't think there is any reason to use
        // a more sophisticated technique like shared_ptr or something
        // similar, since there are no non-const methods.
        QRPD(const QRPD<M>& rhs);

        // Clean up the internal storage
        ~QRPD();


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
        getp_type getP() const;

        bool preferInPlace() const { return false; }

    private :

        // mutable so the normal copy constructor with the argument
        // const QRPD<M>& can release the memory.
        mutable std::auto_ptr<QRPD_Impl<small,M> > pimpl;

        size_t colsize() const;
        size_t rowsize() const;

        // op= not allowed.
        QRPD<M>& operator=(const QRPD<M>&);
    };


    template <class T>
    class InstQRPD :
        public QRPD<Matrix<T,ColMajor> >,
        public Divider<T>
    {
    public :
        typedef QRPD<Matrix<T,ColMajor> > base;
        typedef typename base::RT RT;
        typedef typename base::CT CT;
        typedef typename base::FT FT;
        typedef typename base::ZFT ZFT;

        // Sets up the internal storage and does the decomposition.
        template <int C>
        InstQRPD(const ConstMatrixView<T,C>& A, bool _inplace=false);
        InstQRPD(const InstQRPD<T>& rhs);
        ~InstQRPD();

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
        InstQRPD<T>& operator=(const InstQRPD<T>&);
    };

    template <bool isvalid, bool istrans>
    struct QRPHelper;

    template <>
    struct QRPHelper<true,false>
    {
        template <class M1, class V1, class M2, class M3>
        static TMV_INLINE void solve(
            const M1& QRx, const V1& beta, const Permutation& P, int N1,
            const M2& m2, M3& m3)
        { QR_Solve(QRx,beta,&P,N1,m2,m3); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void solveInPlace(
            const M1& QRx, const V1& beta, const Permutation& P, int N1, M2& m2)
        { QR_SolveInPlace(QRx,beta,&P,N1,m2); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverse(
            const M1& QRx, const V1& beta, const Permutation& P, int N1, M2& m2)
        { QR_Inverse(QRx,beta,&P,N1,m2); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverseATA(
            const M1& QRx, const V1& beta, const Permutation& P, int N1, M2& m2)
        { QR_InverseATA(QRx,beta,&P,N1,m2); }
    };
    template <>
    struct QRPHelper<true,true>
    {
        template <class M1, class V1, class M2, class M3>
        static TMV_INLINE void solve(
            const M1& QRx, const V1& beta, const Permutation& P,
            int N1, const M2& m2, M3& m3)
        { QR_SolveTranspose(QRx,beta,&P,N1,m2,m3); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void solveInPlace(
            const M1& QRx, const V1& beta, const Permutation& P, int N1, M2& m2)
        { QR_SolveTransposeInPlace(QRx,beta,&P,N1,m2); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverse(
            const M1& QRx, const V1& beta, const Permutation& P, int N1, M2& m2)
        { 
            typename M2::transpose_type m2t = m2.transpose();
            QR_Inverse(QRx,beta,&P,N1,m2t);
        }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverseATA(
            const M1& QRx, const V1& beta, const Permutation& P, int N1, M2& m2)
        { QR_InverseATA(QRx,beta,&P,N1,m2); }
    };
    template <bool istrans>
    struct QRPHelper<false,istrans>
    {
        template <class M1, class V1, class M2, class M3>
        static TMV_INLINE void solve(
            const M1& , const V1& , const Permutation& P, int N1,
            const M2& , M3& ) 
        { TMVAssert(false && "Calling invalid QRPHelper::solve\n"); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void solveInPlace(
            const M1& , const V1& , const Permutation& P, int N1, M2& )
        { TMVAssert(false && "Calling invalid QRPHelper::solveInPlace\n"); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverse(
            const M1& , const V1& , const Permutation& P, int N1, M2& )
        { TMVAssert(false && "Calling invalid QRPHelper::makeInverse\n"); }
        template <class M1, class V1, class M2>
        static TMV_INLINE void makeInverseATA(
            const M1& , const V1& , const Permutation& P, int N1, M2& ) 
        { TMVAssert(false && "Calling invalid QRPHelper::makeInverseATA\n"); }
    };

    template <class M>
    struct QRPD_Impl<true,M> 
    // small = true, so cs,rs both known 
    {
        typedef typename M::real_type RT;
        typedef typename M::value_type T;
        enum { istrans = int(M::_colsize) < int(M::_rowsize) };
        enum { cs = istrans ? int(M::_rowsize) : int(M::_colsize) };
        enum { rs = istrans ? int(M::_colsize) : int(M::_rowsize) };
        typedef typename MCopyHelper<T,Rec,M::_colsize,M::_rowsize,istrans>::type Mc;
        typedef typename TypeSelect< istrans ,
                typename Mc::transpose_type ,
                typename Mc::view_type >::type qrx_type;
        typedef typename VCopyHelper<RT,rs>::type beta_type;

        template <class M2>
        QRPD_Impl(const BaseMatrix<M2>& A, bool ) : 
            QRx(Maybe<istrans>::transposeview(SmallQRx) ), P(rs), N1(rs)
        {
            TMVStaticAssert(M::_colsize != TMV_UNKNOWN);
            TMVStaticAssert(M::_rowsize != TMV_UNKNOWN);
            TMVAssert(A.colsize() == istrans ? int(rs) : int(cs));
            TMVAssert(A.rowsize() == istrans ? int(cs) : int(rs));
            //std::cout<<"QRPD_Impl small\n";
            //std::cout<<"istrans = "<<istrans<<std::endl;
            //std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            //std::cout<<"A = "<<A<<std::endl;
            //std::cout<<"SmallQRx = "<<SmallQRx<<std::endl;
            //std::cout<<"QRx = "<<QRx<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            A.newAssignTo(SmallQRx);
            //std::cout<<"SmallQRx => "<<SmallQRx<<std::endl;
            //std::cout<<"QRx => "<<QRx<<std::endl;
            QRP_Decompose(QRx,beta,P,QRP_StrictSingleton::inst());
            //std::cout<<"SmallQRx => "<<SmallQRx<<std::endl;
            //std::cout<<"QRx => "<<QRx<<std::endl;
            //std::cout<<"beta => "<<beta<<std::endl;
            while (N1 > 0 && QRx.cref(N1-1,N1-1) == T(0)) --N1;
        }
        template <class M2, class M3>
        void solve(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            QRPHelper<isvalid,istrans>::solve(QRx,beta,P,N1,m2,m3);
        }
        template <class M2, class M3>
        void solveTranspose(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            QRPHelper<isvalid,!istrans>::solve(QRx,beta,P,N1,m2,m3);
        }
        template <class M2>
        void solveInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRPHelper<isvalid,istrans>::solveInPlace(QRx,beta,P,N1,m2);
        }
        template <class M2>
        void solveTransposeInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRPHelper<isvalid,!istrans>::solveInPlace(QRx,beta,P,N1,m2);
        }
        template <class M2>
        void makeInverse(M2& minv)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRPHelper<isvalid,istrans>::makeInverse(QRx,beta,P,N1,minv);
        }
        template <class M2>
        void makeInverseATA(M2& ata)
        { 
            const bool isvalid = M::isreal || M2::iscomplex;
            QRPHelper<isvalid,istrans>::makeInverseATA(QRx,beta,P,N1,ata);
        }

        Mc SmallQRx;
        qrx_type QRx;
        beta_type beta;
        Permutation P;
        int N1;
    };
    
    template <class M>
    struct QRPD_Impl<false,M>
    {
        typedef typename M::value_type T;
        typedef typename M::real_type RT;
        enum { cs1 = M::_colsize };
        enum { rs1 = M::_rowsize };
        enum { knownsizes = cs1 != TMV_UNKNOWN && rs1 != TMV_UNKNOWN };
        enum { istrans1 = knownsizes && cs1 < int(rs1) };
        enum { cs = IntTraits2<cs1,rs1>::max };
        enum { rs = IntTraits2<cs1,rs1>::min };
        typedef typename MViewHelper<T,Rec,cs,rs,1,TMV_UNKNOWN>::type qrx_type;
        typedef Vector<RT> beta_type;

        template <class M2>
        QRPD_Impl(const BaseMatrix_Rec<M2>& A, bool _inplace) :
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
                  int(istrans ? A.rowsize() : A.colsize()) ) // stepj
            ),
            beta(istrans ? A.colsize() : A.rowsize()),
            P(beta.size()), N1(beta.size())
            {
                //std::cout<<"QRD_Impl A = "<<TMV_Text(A)<<std::endl;
                if (!inplace) {
                    if (istrans) {
                        typename qrx_type::transpose_type QRxt = QRx.transpose();
                        Maybe<!knownsizes||istrans1>::newAssignTo(A,QRxt);
                    } else {
                        Maybe<!knownsizes||!istrans1>::newAssignTo(A,QRx);
                    }
                } else {
                    Maybe<M2::_conj>::conjself(QRx);
                }
                QRP_Decompose(QRx,beta,P,QRP_StrictSingleton::inst());
                //std::cout<<"After QRP_Decompose"<<std::endl;
                //std::cout<<"N1 = "<<N1<<std::endl;
                while (N1 > 0 && QRx.cref(N1-1,N1-1) == T(0)) {
                    //std::cout<<"N1 = "<<N1<<std::endl;
                    //std::cout<<"QRx(N1-1,N1-1) = "<<QRx.cref(N1-1,N1-1)<<std::endl;
                    --N1;
                }
                //std::cout<<"After N1 loop"<<std::endl;
            }

        // If A is not a BaseMatrix_Rec, can't do it in place.
        template <class M2>
        QRPD_Impl(const BaseMatrix<M2>& A, bool _inplace) :
            istrans(A.colsize() < A.rowsize()), inplace(false),
            Aptr( A.rowsize()*A.colsize() ),
            QRx(
                Aptr.get() , // ptr
                istrans ? A.rowsize() : A.colsize() ,  // colsize
                istrans ? A.colsize() : A.rowsize() ,  // rowsize
                1 , // stepi
                int(istrans ? A.rowsize() : A.colsize()) // stepj
            ),
            beta(istrans ? A.colsize() : A.rowsize()),
            P(beta.size()), N1(beta.size())
            {
                //std::cout<<"QRD_Impl non-Rec A = "<<TMV_Text(A)<<std::endl;
                if (istrans) {
                    typename qrx_type::transpose_type QRxt = QRx.transpose();
                    Maybe<!knownsizes||istrans1>::newAssignTo(A,QRxt);
                } else {
                    Maybe<!knownsizes||!istrans1>::newAssignTo(A,QRx);
                }
                QRP_Decompose(QRx,beta,P,QRP_StrictSingleton::inst());
                //std::cout<<"After QRP_Decompose"<<std::endl;
                //std::cout<<"N1 = "<<N1<<std::endl;
                while (N1 > 0 && QRx.cref(N1-1,N1-1) == T(0)) {
                    //std::cout<<"N1 = "<<N1<<std::endl;
                    //std::cout<<"QRx(N1-1,N1-1) = "<<QRx.cref(N1-1,N1-1)<<std::endl;
                    --N1;
                }
                //std::cout<<"After N1 loop"<<std::endl;
            }

        template <class M2, class M3>
        void solve(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRPHelper<isvalid1,true>::solve(QRx,beta,P,N1,m2,m3);
            else
                QRPHelper<isvalid2,false>::solve(QRx,beta,P,N1,m2,m3);
        }
        template <class M2, class M3>
        void solveTranspose(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRPHelper<isvalid1,false>::solve(QRx,beta,P,N1,m2,m3);
            else
                QRPHelper<isvalid2,true>::solve(QRx,beta,P,N1,m2,m3);
        }
        template <class M2>
        void solveInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRPHelper<isvalid1,true>::solveInPlace(QRx,beta,P,N1,m2);
            else
                QRPHelper<isvalid2,false>::solveInPlace(QRx,beta,P,N1,m2);
        }
        template <class M2>
        void solveTransposeInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRPHelper<isvalid1,false>::solveInPlace(QRx,beta,P,N1,m2);
            else
                QRPHelper<isvalid2,true>::solveInPlace(QRx,beta,P,N1,m2);
        }
        template <class M2>
        void makeInverse(M2& minv)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                QRPHelper<isvalid1,true>::makeInverse(QRx,beta,P,N1,minv);
            else
                QRPHelper<isvalid2,false>::makeInverse(QRx,beta,P,N1,minv);
        }
        template <class M2>
        void makeInverseATA(M2& ata)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            // For this one, it doesn't actually matter if istrans = true
            // So just use trans = false arbitrarily.
            QRPHelper<isvalid,false>::makeInverseATA(QRx,beta,P,N1,ata);
        }

        const bool istrans;
        const bool inplace;
        AlignedArray<typename M::value_type> Aptr;
        qrx_type QRx;
        beta_type beta;
        Permutation P;
        int N1;
    };

    template <class M> template <class M2>
    QRPD<M>::QRPD(const BaseMatrix<M2>& A, bool inplace) :
        pimpl(new QRPD_Impl<small,M>(A.mat(),inplace)) {}

    template <class M>
    QRPD<M>::QRPD(const QRPD<M>& rhs) : pimpl(rhs.pimpl.release()) {}

    template <class M>
    QRPD<M>::~QRPD() {}

    template <class M> template <class M2, class M3>
    void QRPD<M>::doSolve(
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
    void QRPD<M>::doSolve(
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<V3::_size,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(v3.size() == rowsize());
        pimpl->solve(v2.vec(),v3.vec());
    }

    template <class M> template <class M2>
    void QRPD<M>::doSolveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == colsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveInPlace(m2.mat());
    }

    template <class M> template <class V2>
    void QRPD<M>::doSolveInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveInPlace(v2.vec());
    }

    template <class M> template <class M2, class M3>
    void QRPD<M>::doSolveTranspose(
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
    void QRPD<M>::doSolveTranspose(
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_rowsize>::same));
        TMVStaticAssert((Sizes<V3::_size,M::_colsize>::same));
        TMVAssert(v2.size() == rowsize());
        TMVAssert(v3.size() == colsize());
        pimpl->solveTranspose(v2.vec(),v3.vec());
    }

    template <class M> template <class M2>
    void QRPD<M>::doSolveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == rowsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveTransposeInPlace(m2.mat());
    }

    template <class M> template <class V2>
    void QRPD<M>::doSolveTransposeInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveTransposeInPlace(v2.vec());
    }

    template <class M>
    typename M::value_type QRPD<M>::det() const
    { 
        return getR().det() * 
            typename M::real_type(CalculateDetQ(pimpl->beta)*getP().det());
    }

    template <class M>
    typename M::float_type QRPD<M>::logDet(typename M::zfloat_type* sign) const
    {
        typename M::float_type ret = getR().logDet(sign);
        if (sign) *sign *= 
            typename M::float_type(CalculateDetQ(pimpl->beta)) * getP().det();
        return ret;
    }

    template <class M>
    bool QRPD<M>::isSingular() const 
    { return getR().isSingular(); }

    template <class M> template <class M2>
    void QRPD<M>::doMakeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const
    {
        TMVStaticAssert((Sizes<M::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M::_rowsize,M2::_colsize>::same));
        TMVAssert(minv.colsize() == rowsize());
        TMVAssert(minv.rowsize() == colsize());
        pimpl->makeInverse(minv.mat());
    }

    template <class M> template <class M2>
    void QRPD<M>::doMakeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M2::_rowsize>::same));
        TMVAssert(ata.rowsize() == rowsize());
        TMVAssert(ata.colsize() == rowsize());
        pimpl->makeInverseATA(ata.mat());
    }

    template <class M>
    bool QRPD<M>::isTrans() const 
    { return pimpl->istrans; }

    template <class M>
    typename QRPD<M>::getq_type QRPD<M>::getQ() const 
    { return typename QRPD<M>::getq_type(pimpl->QRx,pimpl->beta); }

    template <class M>
    typename QRPD<M>::getp_type QRPD<M>::getP() const 
    { return pimpl->P; }

    template <class M>
    typename QRPD<M>::getr_type QRPD<M>::getR() const 
    { return pimpl->QRx.upperTri(); }

    template <class M>
    typename QRPD<M>::getqr_type QRPD<M>::getQR() const 
    { return pimpl->QRx; }

    template <class M>
    typename QRPD<M>::getbeta_type QRPD<M>::getBeta() const 
    { return pimpl->beta; }

    template <class M>
    typename M::real_type QRPD<M>::condition(RT normInf) const 
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
    size_t QRPD<M>::colsize() const
    { return pimpl->istrans ? pimpl->QRx.rowsize() : pimpl->QRx.colsize(); }

    template <class M>
    size_t QRPD<M>::rowsize() const
    { return pimpl->istrans ? pimpl->QRx.colsize() : pimpl->QRx.rowsize(); }


    template <class M, class M2>
    static bool CheckDecomp(
        const QRPD<M>& qrpd, const BaseMatrix_Calc<M2>& m, std::ostream* fout=0)
    {
        typedef typename M2::real_type RT;
        //bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
        bool printmat = fout;
        if (printmat) {
            *fout << "QRP:\n";
            if (qrpd.isTrans()) *fout << m.transpose() << std::endl;
            else *fout << m << std::endl;
            *fout << "Q = "<<qrpd.getQ()<<std::endl;
            *fout << "R = "<<qrpd.getR()<<std::endl;
            *fout << "P = "<<qrpd.getP()<<std::endl;
            *fout << "  or by interchanges: ";
            for(int i=0;i<int(qrpd.getP().size());i++)
                *fout<<(qrpd.getP().getValues())[i]<<" ";
            *fout<<std::endl;
        }
        typename M::copy_type qrp = qrpd.getQ()*qrpd.getR()*qrpd.getP();
        if (printmat) {
            *fout << "QRP = "<<qrp<<std::endl;
        }
        RT nm = qrpd.isTrans() ? Norm(qrp-m.transpose()) : Norm(qrp-m);
        nm /= Norm(qrpd.getQ())*Norm(qrpd.getR());
        if (fout) {
            *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<" <? ";
            *fout << RT(m.colsize())<<"*"<<TMV_Epsilon<RT>();
            *fout << " = "<<RT(m.colsize())*TMV_Epsilon<RT>()<<std::endl;
        }
        return nm < RT(m.colsize())*TMV_Epsilon<RT>();
    }


} // namespace mv


#endif
