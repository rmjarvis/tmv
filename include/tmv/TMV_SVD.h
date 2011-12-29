


#ifndef TMV_SVD_H
#define TMV_SVD_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Diag.h"
#include "TMV_BaseVector.h"
#include "TMV_Divider.h"
#include "TMV_Array.h"

namespace tmv {

    // In TMV_SVDecompose.h
    template <class Mu, class Ms, class Mv>
    inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U, BaseMatrix_Diag_Mutable<Ms>& S,
        BaseMatrix_Rec_Mutable<Mv>& V, 
        typename Mu::zfloat_type& signdet, typename Mu::float_type& logdet,
        bool StoreU);

    // In TMV_SVInverse.h
    template <class M1u, class M1s, class M1v, class M2>
    inline void SV_Inverse(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, 
        BaseMatrix_Rec_Mutable<M2>& minv);
    template <class M1u, class M1s, class M1v, class M2>
    inline void SV_InverseATA(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, BaseMatrix_Rec_Mutable<M2>& ata);

    // In TMV_SVDiv.h
    template <class M1u, class M1s, class M1v, class M2, class M3>
    inline void SV_Solve(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <class M1u, class M1s, class M1v, class V2, class V3>
    inline void SV_Solve(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3);




    // The point of the Impl class here is to implement the transfer of 
    // ownership copy semantics.
    // It also differentiates between small and non-small implementations.
    template <bool small, class M>
    struct SVD_Impl;

    template <class M>
    class SVD 
    {
    public :

        typedef typename M::value_type T;
        typedef typename M::real_type RT;
        typedef typename M::complex_type CT;
        typedef typename M::float_type FT;
        typedef typename M::zfloat_type ZFT;

        // This next bit finds the storage type to use for the u matrix
        // regardless of what kind of matrix M is.  e.g. this should
        // work even if M is a TriMatrix or a BandMatrix, etc.
        enum { cs = IntTraits2<M::_colsize,M::_rowsize>::max };
        enum { rs = IntTraits2<M::_colsize,M::_rowsize>::min };

        enum { small = (
                M::_colsize != TMV_UNKNOWN && M::_rowsize != TMV_UNKNOWN
                && M::_colsize <= 32 && M::_rowsize <= 32 ) };

        typedef typename SVD_Impl<small,M>::getu_type getu_type;
        typedef typename SVD_Impl<small,M>::gets_type gets_type;
        typedef typename SVD_Impl<small,M>::getv_type getv_type;

        //
        // Constructors
        //

        // Sets up the internal storage and does the decomposition.
        template <class M2>
        SVD(const BaseMatrix<M2>& A, bool _inplace=false);

        // The copy constructor has transfer of ownership semantics.
        // This way an SVD object can be returned by value, and the 
        // copy is cheap.  I don't think there is any reason to use
        // a more sophisticated technique like shared_ptr or something
        // similar, since there are no non-const methods.
        SVD(const SVD<M>& rhs);

        // Clean up the internal storage
        ~SVD();


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
        // Condition (kappa_2) and Norm2
        //

        RT norm2() const;
        RT condition() const;

        //
        // Access Decomposition
        //

        getu_type getU() const;
        gets_type getS() const;
        getv_type getV() const;

        bool preferInPlace() const { return false; }

        // 
        // Determine which (if any) S values to zero out
        //

        void thresh(RT toler, std::ostream* debugout=0) const;
        void top(int neigen, std::ostream* debugout=0) const;
        int getKMax() const;

    private :

        // mutable so the normal copy constructor with the argument
        // const SVD<M>& can release the memory.
        mutable std::auto_ptr<SVD_Impl<small,M> > pimpl;

        int colsize() const;
        int rowsize() const;

        // op= not allowed.
        SVD<M>& operator=(const SVD<M>&);
    };


    template <class T>
    class InstSVD :
        public SVD<Matrix<T,ColMajor> >,
        public Divider<T>
    {
    public :
        typedef SVD<Matrix<T,ColMajor> > base;
        typedef typename base::RT RT;
        typedef typename base::CT CT;
        typedef typename base::FT FT;
        typedef typename base::ZFT ZFT;

        // Sets up the internal storage and does the decomposition.
        template <int C>
        InstSVD(const ConstMatrixView<T,C>& A, bool _inplace=false);
        InstSVD(const InstSVD<T>& rhs);
        ~InstSVD();

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
        
        RT condition() const;
        inline RT condition(RT ) const
        { return condition(); }
        bool preferInPlace() const;

    private :
        // op= not allowed.
        InstSVD<T>& operator=(const InstSVD<T>&);
    };

    template <bool isvalid, bool istrans>
    struct SVHelper;

    template <>
    struct SVHelper<true,false>
    {
        template <class M1u, class M1s, class M1v, class M2, class M3>
        static TMV_INLINE void solve(
            const M1u& U, const M1s& S, const M1v& V, int kmax,
            const M2& m2, M3& m3)
        {
            SV_Solve(
                U.colRange(0,kmax),S.subDiagMatrix(0,kmax),
                V.rowRange(0,kmax),m2,m3); 
        }
        template <class M1u, class M1s, class M1v, class M2>
        static TMV_INLINE void solveInPlace(
            const M1u& U, const M1s& S, const M1v& V, int kmax, M2& m2)
        { 
            SV_Solve(
                U.colRange(0,kmax),S.subDiagMatrix(0,kmax),
                V.rowRange(0,kmax),m2.copy(),m2); 
        }
        template <class M1u, class M1s, class M1v, class M2>
        static TMV_INLINE void makeInverse(
            const M1u& U, const M1s& S, const M1v& V, int kmax, M2& minv)
        { 
            SV_Inverse(
                U.colRange(0,kmax),S.subDiagMatrix(0,kmax),
                V.rowRange(0,kmax),minv); 
        }
        template <class M1u, class M1s, class M1v, class M2>
        static TMV_INLINE void makeInverseATA(
            const M1u& U, const M1s& S, const M1v& V, int kmax, M2& ata)
        { 
            SV_InverseATA(
                U.colRange(0,kmax),S.subDiagMatrix(0,kmax),
                V.rowRange(0,kmax),ata); 
        }
    };
    template <>
    struct SVHelper<true,true>
    {
        template <class M1u, class M1s, class M1v, class M2, class M3>
        static TMV_INLINE void solve(
            const M1u& U, const M1s& S, const M1v& V, int kmax,
            const M2& m2, M3& m3)
        {
            SV_Solve(
                V.rowRange(0,kmax).transpose(),S.subDiagMatrix(0,kmax),
                U.colRange(0,kmax).transpose(),m2,m3); 
        }
        template <class M1u, class M1s, class M1v, class M2>
        static TMV_INLINE void solveInPlace(
            const M1u& U, const M1s& S, const M1v& V, int kmax, M2& m2)
        {
            SV_Solve(
                V.rowRange(0,kmax).transpose(),S.subDiagMatrix(0,kmax),
                U.colRange(0,kmax).transpose(),m2.copy(),m2); 
        }
        template <class M1u, class M1s, class M1v, class M2>
        static TMV_INLINE void makeInverse(
            const M1u& U, const M1s& S, const M1v& V, int kmax, M2& minv)
        { 
            SV_Inverse(
                V.rowRange(0,kmax).transpose(),S.subDiagMatrix(0,kmax),
                U.colRange(0,kmax).transpose(),minv); 
        }
        template <class M1u, class M1s, class M1v, class M2>
        static TMV_INLINE void makeInverseATA(
            const M1u& U, const M1s& S, const M1v& V, int kmax, M2& ata)
        { 
            SV_InverseATA(
                U.colRange(0,kmax),S.subDiagMatrix(0,kmax),V.rowRange(0,kmax),
                ata); 
        }
    };
    template <bool istrans>
    struct SVHelper<false,istrans>
    {
        template <class M1u, class M1s, class M1v, class M2, class M3>
        static TMV_INLINE void solve(
            const M1u& , const M1s& , const M1v& , int , const M2& , M3& )
        { TMVAssert(false && "Calling invalid SVHelper::solve\n"); }
        template <class M1u, class M1s, class M1v, class M2>
        static TMV_INLINE void solveInPlace(
            const M1u& , const M1s& , const M1v& , int , M2& )
        { TMVAssert(false && "Calling invalid SVHelper::solveInPlace\n"); }
        template <class M1u, class M1s, class M1v, class M2>
        static TMV_INLINE void makeInverse(
            const M1u& , const M1s& , const M1v& , int , M2& )
        { TMVAssert(false && "Calling invalid SVHelper::makeInverse\n"); }
        template <class M1u, class M1s, class M1v, class M2>
        static TMV_INLINE void makeInverseATA(
            const M1u& , const M1s& , const M1v& , int , M2& )
        { TMVAssert(false && "Calling invalid SVHelper::makeInverseATA\n"); }
    };

    template <class M>
    struct SVD_Impl<true,M> 
    // small = true, so cs,rs both known 
    {
        typedef typename M::real_type RT;
        typedef typename M::value_type T;
        typedef typename M::float_type FT;
        typedef typename M::zfloat_type ZT;
        enum { cs1 = M::_colsize };
        enum { rs1 = M::_rowsize };
        enum { istrans = int(cs1) < int(rs1) };
        enum { cs = IntTraits2<cs1,rs1>::max };
        enum { rs = IntTraits2<cs1,rs1>::min };
        enum { A = (istrans ? RowMajor : ColMajor) | NoAlias };
        typedef typename MCopyHelper<T,Rec,cs1,rs1,A>::type Mu;
        typedef typename TypeSelect< istrans ,
                typename Mu::transpose_type ,
                typename Mu::view_type >::type ux_type;
        typedef typename MCopyHelper<T,Diag,rs,rs>::type sx_type;
        typedef typename MCopyHelper<T,Rec,rs,rs,RowMajor|NoAlias>::type vx_type;

        enum { csu = istrans ? int(rs) : int(cs) };
        enum { rsu = rs };
        enum { csv = rs };
        enum { rsv = istrans ? int(cs) : int(rs) };

        typedef typename MViewHelper<T,Rec,csu,rsu,1,csu>::type getu_type;
        typedef typename sx_type::const_view_type gets_type;
        typedef typename MViewHelper<T,Rec,csv,rsv,rsv,1>::type getv_type;

        template <class M2>
        SVD_Impl(const BaseMatrix<M2>& A, bool ) : 
            Ux(Maybe<istrans>::transposeview(SmallUx) ), kmax(0),
            signdet(1), logdet(0)
        {
            TMVStaticAssert(M::_colsize != TMV_UNKNOWN);
            TMVStaticAssert(M::_rowsize != TMV_UNKNOWN);
            TMVAssert(A.colsize() == istrans ? int(rs) : int(cs));
            TMVAssert(A.rowsize() == istrans ? int(cs) : int(rs));
            SmallUx = A;
            SV_Decompose(Ux,Sx,Vx,signdet,logdet,true);
        }
        template <class M2, class M3>
        void solve(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            SVHelper<isvalid,istrans>::solve(Ux,Sx,Vx,kmax,m2,m3);
        }
        template <class M2, class M3>
        void solveTranspose(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            SVHelper<isvalid,!istrans>::solve(Ux,Sx,Vx,kmax,m2,m3);
        }
        template <class M2>
        void solveInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            SVHelper<isvalid,istrans>::solveInPlace(Ux,Sx,Vx,kmax,m2);
        }
        template <class M2>
        void solveTransposeInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            SVHelper<isvalid,!istrans>::solveInPlace(Ux,Sx,Vx,kmax,m2);
        }
        template <class M2>
        void makeInverse(M2& minv)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            SVHelper<isvalid,istrans>::makeInverse(Ux,Sx,Vx,kmax,minv);
        }
        template <class M2>
        void makeInverseATA(M2& ata)
        { 
            const bool isvalid = M::isreal || M2::iscomplex;
            SVHelper<isvalid,istrans>::makeInverseATA(Ux,Sx,Vx,kmax,ata);
        }

        Mu SmallUx;
        ux_type Ux;
        sx_type Sx;
        vx_type Vx;
        ZT signdet;
        FT logdet;
        int kmax;
    };
    
    template <class M>
    struct SVD_Impl<false,M>
    {
        typedef typename M::value_type T;
        typedef typename M::real_type RT;
        typedef typename M::float_type FT;
        typedef typename M::zfloat_type ZT;
        enum { cs1 = M::_colsize };
        enum { rs1 = M::_rowsize };
        enum { knownsizes = cs1 != TMV_UNKNOWN && rs1 != TMV_UNKNOWN };
        enum { istrans1 = knownsizes && cs1 < int(rs1) };
        enum { cs = IntTraits2<cs1,rs1>::max };
        enum { rs = IntTraits2<cs1,rs1>::min };
        typedef typename MViewHelper<T,Rec,cs,rs,1,TMV_UNKNOWN,ColMajor|NoAlias>::type ux_type;
        typedef DiagMatrix<RT> sx_type;
        typedef Matrix<T,RowMajor|NoDivider|NoAlias> vx_type;

        enum { csu = (
                knownsizes ? ( istrans1 ? int(rs) : int(cs) ) :
                int(TMV_UNKNOWN) ) };
        enum { rsu = knownsizes ? int(rs) : int(TMV_UNKNOWN) };
        enum { csv = knownsizes ? int(rs) : int(TMV_UNKNOWN) };
        enum { rsv = (
                knownsizes ? ( istrans1 ? int(cs) : int(rs) ) :
                int(TMV_UNKNOWN) ) };

        typedef typename MViewHelper<T,Rec,csu,rsu,1,csu>::type getu_type;
        typedef typename sx_type::const_view_type gets_type;
        typedef typename MViewHelper<T,Rec,csv,rsv,rsv,1>::type getv_type;

        template <class M2>
        SVD_Impl(const BaseMatrix_Rec<M2>& A, bool _inplace) :
            // if A is short, need to transpose
            istrans(knownsizes ? istrans1 : A.colsize() < A.rowsize()),
            // inplace only if it works with a ColMajor Ux object
            inplace( _inplace && 
                     ((A.iscm() && !istrans) || (A.isrm() && istrans)) ),
            // Aptr is the pointer to new storage if any
            Aptr( inplace ? 0 : A.rowsize()*A.colsize() ),
            // Ux views this memory as the U matrix
            Ux(
                inplace ? A.nonConst().ptr() : Aptr.get() , // ptr
                istrans ? A.rowsize() : A.colsize() ,  // colsize
                istrans ? A.colsize() : A.rowsize() ,  // rowsize
                inplace ? (istrans ? A.stepi() : A.stepj()) : 1 , // stepi
                ( inplace ? (istrans ? A.stepj() : A.stepi()) : 
                  (istrans ? A.rowsize() : A.colsize()) ) // stepj
            ),
            Sx(Ux.rowsize()), Vx(Ux.rowsize(),Ux.rowsize()), kmax(0),
            signdet(1), logdet(0)
            {
                if (!inplace) {
                    if (istrans) {
                        typename ux_type::transpose_type Uxt = Ux.transpose();
                        Maybe<!knownsizes||istrans1>::assignTo(A,Uxt);
                    } else {
                        Maybe<!knownsizes||!istrans1>::assignTo(A,Ux);
                    }
                } else {
                    Maybe<M2::_conj>::conjself(Ux);
                }
                SV_Decompose(Ux,Sx,Vx,signdet,logdet,true);
            }

        // If A is not a BaseMatrix_Rec, can't do it in place.
        template <class M2>
        SVD_Impl(const BaseMatrix<M2>& A, bool _inplace) :
            istrans(A.colsize() < A.rowsize()), inplace(false),
            Aptr( A.rowsize()*A.colsize() ),
            Ux(
                Aptr.get() , // ptr
                istrans ? A.rowsize() : A.colsize() ,  // colsize
                istrans ? A.colsize() : A.rowsize() ,  // rowsize
                1 , // stepi
                (istrans ? A.rowsize() : A.colsize()) // stepj
            ),
            Sx(Ux.rowsize()), Vx(Ux.rowsize(),Ux.rowsize()), kmax(0), 
            signdet(1), logdet(0)
            {
                if (istrans) {
                    typename ux_type::transpose_type Uxt = Ux.transpose();
                    Maybe<!knownsizes||istrans1>::assignTo(A,Uxt);
                } else {
                    Maybe<!knownsizes||!istrans1>::assignTo(A,Ux);
                }
                SV_Decompose(Ux,Sx,Vx,signdet,logdet,true);
            }

        template <class M2, class M3>
        void solve(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                SVHelper<isvalid1,true>::solve(Ux,Sx,Vx,kmax,m2,m3);
            else
                SVHelper<isvalid2,false>::solve(Ux,Sx,Vx,kmax,m2,m3);
        }
        template <class M2, class M3>
        void solveTranspose(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                SVHelper<isvalid1,false>::solve(Ux,Sx,Vx,kmax,m2,m3);
            else
                SVHelper<isvalid2,true>::solve(Ux,Sx,Vx,kmax,m2,m3);
        }
        template <class M2>
        void solveInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                SVHelper<isvalid1,true>::solveInPlace(Ux,Sx,Vx,kmax,m2);
            else
                SVHelper<isvalid2,false>::solveInPlace(Ux,Sx,Vx,kmax,m2);
        }
        template <class M2>
        void solveTransposeInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                SVHelper<isvalid1,false>::solveInPlace(Ux,Sx,Vx,kmax,m2);
            else
                SVHelper<isvalid2,true>::solveInPlace(Ux,Sx,Vx,kmax,m2);
        }
        template <class M2>
        void makeInverse(M2& minv)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            const bool isvalid1 = isvalid && (!knownsizes || istrans1);
            const bool isvalid2 = isvalid && (!knownsizes || !istrans1);
            if (istrans)
                SVHelper<isvalid1,true>::makeInverse(Ux,Sx,Vx,kmax,minv);
            else
                SVHelper<isvalid2,false>::makeInverse(Ux,Sx,Vx,kmax,minv);
        }
        template <class M2>
        void makeInverseATA(M2& ata)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            // For this one, it doesn't actually matter if istrans = true
            // So just use trans = false arbitrarily.
            SVHelper<isvalid,false>::makeInverseATA(Ux,Sx,Vx,kmax,ata);
        }

        const bool istrans;
        const bool inplace;
        AlignedArray<typename M::value_type> Aptr;
        ux_type Ux;
        sx_type Sx;
        vx_type Vx;
        int kmax;
        ZT signdet;
        FT logdet;
    };

    template <class M> template <class M2>
    SVD<M>::SVD(const BaseMatrix<M2>& A, bool inplace) :
        pimpl(new SVD_Impl<small,M>(A.mat(),inplace))
    { thresh(TMV_Epsilon<T>()); }

    template <class M>
    SVD<M>::SVD(const SVD<M>& rhs) : pimpl(rhs.pimpl.release()) {}

    template <class M>
    SVD<M>::~SVD() {}

    template <class M> template <class M2, class M3>
    void SVD<M>::doSolve(
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
    void SVD<M>::doSolve(
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<V3::_size,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(v3.size() == rowsize());
        pimpl->solve(v2.vec(),v3.vec());
    }

    template <class M> template <class M2>
    void SVD<M>::doSolveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == colsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveInPlace(m2.mat());
    }

    template <class M> template <class V2>
    void SVD<M>::doSolveInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveInPlace(v2.vec());
    }

    template <class M> template <class M2, class M3>
    void SVD<M>::doSolveTranspose(
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
    void SVD<M>::doSolveTranspose(
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_rowsize>::same));
        TMVStaticAssert((Sizes<V3::_size,M::_colsize>::same));
        TMVAssert(v2.size() == rowsize());
        TMVAssert(v3.size() == colsize());
        pimpl->solveTranspose(v2.vec(),v3.vec());
    }

    template <class M> template <class M2>
    void SVD<M>::doSolveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == rowsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveTransposeInPlace(m2.mat());
    }

    template <class M> template <class V2>
    void SVD<M>::doSolveTransposeInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(colsize() == rowsize());
        pimpl->solveTransposeInPlace(v2.vec());
    }

    template <class M>
    typename M::value_type SVD<M>::det() const
    { 
        typedef typename M::value_type T;
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);
    }

    template <class M>
    typename M::float_type SVD<M>::logDet(typename M::zfloat_type* sign) const
    {
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;
    }                  

    template <class M>
    bool SVD<M>::isSingular() const 
    { return getS().isSingular(); }

    template <class M> template <class M2>
    void SVD<M>::doMakeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const
    {
        TMVStaticAssert((Sizes<M::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M::_rowsize,M2::_colsize>::same));
        TMVAssert(minv.colsize() == rowsize());
        TMVAssert(minv.rowsize() == colsize());
        pimpl->makeInverse(minv.mat());
    }

    template <class M> template <class M2>
    void SVD<M>::doMakeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M2::_rowsize>::same));
        TMVAssert(ata.rowsize() == rowsize());
        TMVAssert(ata.colsize() == rowsize());
        pimpl->makeInverseATA(ata.mat());
    }

    template <class M>
    typename SVD<M>::getu_type SVD<M>::getU() const 
    { return pimpl->istrans ? pimpl->Vx.transpose() : pimpl->Ux.view(); }

    template <class M>
    typename SVD<M>::gets_type SVD<M>::getS() const 
    { return pimpl->Sx.view(); }

    template <class M>
    typename SVD<M>::getv_type SVD<M>::getV() const 
    { return pimpl->istrans ? pimpl->Ux.transpose() : pimpl->Vx.view(); }

    template <class M>
    typename M::real_type SVD<M>::norm2() const 
    { return pimpl->Sx.cref(0); }

    template <class M>
    typename M::real_type SVD<M>::condition() const 
    { return pimpl->Sx.cref(0) / pimpl->Sx.cref(pimpl->Sx.size()-1); }

    template <class M>
    int SVD<M>::colsize() const
    { return pimpl->istrans ? pimpl->Ux.rowsize() : pimpl->Ux.colsize(); }

    template <class M>
    int SVD<M>::rowsize() const
    { return pimpl->istrans ? pimpl->Ux.colsize() : pimpl->Ux.rowsize(); }

    template <class M>
    void SVD<M>::thresh(RT toler, std::ostream* debugout) const
    {
        if (pimpl->Sx.size() > 0) {
            TMVAssert(toler < RT(1));
            RT thresh = pimpl->Sx(0)*toler;
            for(pimpl->kmax=pimpl->Sx.size();
                pimpl->kmax>0 && pimpl->Sx(pimpl->kmax-1)<=thresh;
                --pimpl->kmax);
            if(debugout) {
                (*debugout)<<"S.diag() = "<<pimpl->Sx.diag()<<std::endl;
                (*debugout)<<"Smax = "<<pimpl->Sx(0)<<
                    ", thresh = "<<thresh<<std::endl;
                (*debugout)<<"kmax = "<<pimpl->kmax;
                (*debugout)<<" (S.size = "<<pimpl->Sx.size()<<")"<<std::endl;
            }
        }
    }

    template <class M>
    void SVD<M>::top(int neigen, std::ostream* debugout) const
    {
        TMVAssert(neigen > 0);
        if (neigen < pimpl->Sx.size()) pimpl->kmax = neigen;
        else pimpl->kmax = pimpl->Sx.size();
        if(debugout) {
            (*debugout)<<"S.diag() = "<<pimpl->Sx.diag()<<std::endl;
            (*debugout)<<"kmax = "<<pimpl->kmax;
            (*debugout)<<" (S.size = "<<pimpl->Sx.size()<<")"<<std::endl;
        }
    }

    template <class M>
    int SVD<M>::getKMax() const
    { return pimpl->kmax; }

    template <class M, class M2>
    static bool CheckDecomp(
        const SVD<M>& svd, const BaseMatrix_Calc<M2>& m, std::ostream* fout=0)
    {
        typedef typename M2::real_type RT;
        //bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
        bool printmat = fout;
        if (printmat) {
            *fout << "SV:\n";
            *fout << m << std::endl;
            *fout << "U = "<<svd.getU()<<std::endl;
            *fout << "S = "<<svd.getS()<<std::endl;
            *fout << "V = "<<svd.getV()<<std::endl;
        }
        typename M::copy_type usv = svd.getU()*svd.getS()*svd.getV();
        if (printmat) {
            *fout << "USV = "<<usv<<std::endl;
        }
        RT nm = Norm(usv-m);
        nm /= Norm(svd.getU())*Norm(svd.getS())*Norm(svd.getV());
        if (fout) {
            *fout << "Norm(M-USV)/Norm(USV) = "<<nm<<" <? ";
            *fout << RT(m.colsize())<<"*"<<TMV_Epsilon<RT>();
            *fout << " = "<<RT(m.colsize())*TMV_Epsilon<RT>()<<std::endl;
        }
        return nm < RT(m.colsize())*TMV_Epsilon<RT>();
    }


} // namespace mv


#endif
