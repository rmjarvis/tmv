

//---------------------------------------------------------------------------
//
// This file contains the driver code for doing division using 
// LU Decomposition for band matrices.
//
// The basics of LU decomposition for band matrices are the same as
// for regular matrices.  However there are a few wrinkles about doing 
// it efficiently.
//
// The main difference as far as this file is concerned is that
// L can be stored in a lower band matrix with m.nlo() subdiagonals.
// However, U needs m.nlo() + m.nhi() superdiagonals for its storage.
// So the decomposition cannot be done in place.


#ifndef TMV_BandLU_H
#define TMV_BandLU_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_Permutation.h"
#include "TMV_Divider.h"

#ifdef PRINTALGO_BandLU
#include <iostream>
#endif

namespace tmv {

    // In TMV_BandLUDecompose.h
    template <class M>
    inline void BandLU_Decompose(
        BaseMatrix_Band_Mutable<M>& m, Permutation& P);

    // In TMV_BandLUInverse.h
    template <class M1>
    inline void BandLU_Inverse(
        BaseMatrix_Band_Mutable<M1>& m1, const Permutation& P);
    template <class M1, class M2>
    inline void BandLU_InverseATA(
        const BaseMatrix_Band<M1>& m1, const Permutation& P,
        const bool trans, BaseMatrix_Band_Mutable<M2>& m2);

    // In TMV_BandLUDiv.h
    template <class M1, class M2>
    inline void BandLU_SolveInPlace(
        const BaseMatrix_Band<M1>& m1, const Permutation& P,
        BaseMatrix_Band_Mutable<M2>& m2);
    template <class M1, class V2>
    inline void BandLU_SolveInPlace(
        const BaseMatrix_Band<M1>& m1, const Permutation& P,
        BaseVector_Mutable<V2>& v2);
    template <class M1, class M2>
    inline void BandLU_SolveTransposeInPlace(
        const BaseMatrix_Band<M1>& m1, const Permutation& P,
        BaseMatrix_Band_Mutable<M2>& m2);
    template <class M1, class V2>
    inline void BandLU_SolveTransposeInPlace(
        const BaseMatrix_Band<M1>& m1, const Permutation& P,
        BaseVector_Mutable<V2>& v2);


    // The point of the Impl class here is to implement the transfer of
    // ownership copy semantics.
    // It also differentiates between small and non-small implementations.
    template <bool small, class M>
    struct BandLUD_Impl;

    template <class M>
    class BandLUD 
    {
    public :

        typedef typename M::value_type T;
        typedef typename M::real_type RT;
        typedef typename M::complex_type CT;
        typedef typename M::float_type FT;
        typedef typename M::zfloat_type ZFT;

        // This next bit finds the storage type to use for the lu matrix
        // regardless of what kind of matrix M is.  e.g. this should
        // work even if M is a TriMatrix or a BandMatrix, etc.
        enum { cs = M::_colsize };
        enum { rs = M::_rowsize };

        enum { small = (
                M::_colsize != Unknown && M::_rowsize != Unknown &&
                M::_nlo != Unknown && M::_nhi != Unknown &&
                M::_colsize <= 8 && M::_rowsize <= 8 ) };

        typedef typename BandLUD_Impl<small,M>::lux_type lux_type;
        typedef const lux_type& getlu_type;
        typedef LowerTriMatrix<T,UnitDiag> getl_type;
        typedef typename lux_type::const_upperband_type getu_type;
        typedef Permutation getp_type;

        //
        // Constructors
        //

        // Sets up the internal storage and does the decomposition.
        template <class M2>
        BandLUD(const BaseMatrix<M2>& A, bool _inplace=false);

        // The copy constructor has transfer of ownership semantics.
        // This way an BandLUD object can be returned by value, and the 
        // copy is cheap.  I don't think there is any reason to use
        // a more sophisticated technique like shared_ptr or something
        // similar, since there are no non-const methods.
        BandLUD(const BandLUD<M>& rhs);

        // Clean up the internal storage
        ~BandLUD();


        //
        // Perform the division in place
        //
        template <class M2>
        void solveInPlace(BaseMatrix_Mutable<M2>& m2) const;
        template <class V2>
        void solveInPlace(BaseVector_Mutable<V2>& v2) const;

        template <class M2>
        void solveTransposeInPlace(BaseMatrix_Mutable<M2>& m2) const;
        template <class V2>
        void solveTransposeInPlace(BaseVector_Mutable<V2>& v2) const;

        //
        // Next the not-in-place division 
        // For BandLU, we just copy m1->m3 and do the in-place version.
        // 
        
        template <class M1, class M2>
        void solve(
            const BaseMatrix<M1>& m1, BaseMatrix_Mutable<M2>& m2) const
        { solveInPlace(m2=m1); }
        template <class V1, class V2>
        void solve(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const
        { solveInPlace(v2=v1); }

        template <class M1, class M2>
        void solveTranspose(
            const BaseMatrix<M1>& m1, BaseMatrix_Mutable<M2>& m2) const
        { solveTransposeInPlace(m2=m1); }
        template <class V1, class V2>
        void solveTranspose(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const
        { solveTransposeInPlace(v2=v1); }


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
        void makeInverse(BaseMatrix_Mutable<M2>& minv) const;


        //
        // InverseATA
        //
        
        template <class M2>
        void makeInverseATA(BaseMatrix_Mutable<M2>& ata) const;


        // 
        // Condition (kappa_inf)
        //

        RT condition(RT normInf) const;


        //
        // Access Decomposition
        //

        bool isTrans() const;
        getl_type getL() const;
        getu_type getU() const;
        getlu_type getLU() const;
        getp_type getP() const;

        bool preferInPlace() const { return true; }

    private :

        // mutable so the normal copy constructor with the argument
        // const BandLUD<M>& can relase the memory.
        mutable std::auto_ptr<BandLUD_Impl<small,M> > pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

        // op= not allowed.
        BandLUD<M>& operator=(const BandLUD<M>&);
    };
    
    template <class T>
    class InstBandLUD :
        public BandLUD<BandMatrix<T,ColMajor> >,
        public Divider<T>
    {
    public :
        typedef BandLUD<BandMatrix<T,ColMajor> > base;
        typedef typename base::RT RT;
        typedef typename base::CT CT;
        typedef typename base::FT FT;
        typedef typename base::ZFT ZFT;

        // Sets up the internal storage and does the decomposition.
        template <int C>
        InstBandLUD(const ConstBandMatrixView<T,C>& A, bool _inplace=false);
        InstBandLUD(const InstBandLUD<T>& rhs);
        ~InstBandLUD();

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
        InstBandLUD<T>& operator=(const InstBandLUD<T>&);
    };
    
    // Now the instantiation and definition.
    template <bool isvalid, bool istrans>
    struct BandLUHelper;

    template <>
    struct BandLUHelper<true,false>
    {
        template <class M1, class M2, class M3>
        static inline void solve(
            const M1& LUx, const Permutation& P, const M2& m2, M3& m3)
        { BandLU_SolveInPlace(LUx,P,m3=m2); }
        template <class M1, class M2>
        static inline void solveInPlace(
            const M1& LUx, const Permutation& P, M2& m2)
        { BandLU_SolveInPlace(LUx,P,m2); }
        template <class M1, class M2>
        static inline void makeInverse(const M1& LUx, const Permutation& P, M2& m2)
        {
            // This one might be same storage if LUx was done in place.
            // e.g. { A.divideInPlace(); A = A.inverse(); }
            // So go ahead and check. (i.e. Don't use noAlias().)
            // TODO: Add this behavior to test suite.
            m2 = LUx;
            BandLU_Inverse(m2,P);
        }
    };
    template <>
    struct BandLUHelper<true,true>
    {
        template <class M1, class M2, class M3>
        static inline void solve(
            const M1& LUx, const Permutation& P, const M2& m2, M3& m3)
        { BandLU_SolveTransposeInPlace(LUx,P,m3=m2); }
        template <class M1, class M2>
        static inline void solveInPlace(
            const M1& LUx, const Permutation& P, M2& m2)
        { BandLU_SolveTransposeInPlace(LUx,P,m2); }
        template <class M1, class M2>
        static inline void makeInverse(
            const M1& LUx, const Permutation& P, M2& m2)
        {
            typename M2::transpose_type m2t = m2.transpose();
            Copy(LUx,m2t);
            BandLU_Inverse(m2t,P);
        }
    };
    template <bool istrans>
    struct BandLUHelper<false,istrans>
    {
        template <class M1, class M2, class M3>
        static inline void solve(
            const M1& , const Permutation& , const M2&, M3& ) {}
        template <class M1, class M2>
        static inline void solveInPlace(
            const M1& , const Permutation& , M2& ) {}
        template <class M1, class M2>
        static inline void makeInverse(
            const M1& , const Permutation& , M2& ) {}
    };

    template <class M>
    struct BandLUD_Impl<true,M>
    {
        enum { size = M::_colsize };
        enum { lo1 = M::_nlo };
        enum { hi1 = M::_nhi };
        enum { istridiag = lo1 == 1 && hi1 == 1 };
        enum { istrans = (ptrdiff_t(hi1) < ptrdiff_t(lo1) ||
                          (ptrdiff_t(hi1) == ptrdiff_t(lo1) && M::_rowmajor)) };
        enum { lo = (istrans ? (lo1+hi1) : lo1) };
        enum { hi = (istrans ? hi1 : (lo1+hi1)) };
        typedef typename M::value_type T;
        enum { A = NoAlias | (istridiag ? DiagMajor : istrans ? RowMajor : ColMajor) };
        typedef typename BCopyHelper<T,Band,size,size,lo,hi,A>::type Mc;
        typedef typename TypeSelect< istrans ,
                typename Mc::transpose_type ,
                typename Mc::view_type>::type lux_type;

        template <class M2>
        BandLUD_Impl(const BaseMatrix<M2>& A, bool ) : 
            LUx( Maybe<istrans>::transposeview(SmallLUx) ),
            P(size)
        {
            TMVStaticAssert(M::_colsize != Unknown);
            TMVStaticAssert(M::_rowsize != Unknown);
            TMVStaticAssert(M::_nlo != Unknown);
            TMVStaticAssert(M::_nhi != Unknown);
            TMVStaticAssert(M::_colsize == M::_rowsize);
            TMVStaticAssert(M::_colsize == M2::_colsize);
            TMVStaticAssert(M::_rowsize == M2::_rowsize);
            TMVStaticAssert(lux_type::_colmajor);
            SmallLUx.setZero();
            BandMatrixViewOf(SmallLUx,lo,hi) = A;
            BandLU_Decompose(LUx,P);
        }
        template <class M2, class M3>
        void solve(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            BandLUHelper<isvalid,istrans>::solve(LUx,P,m2,m3); 
        }
        template <class M2, class M3>
        void solveTranspose(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            BandLUHelper<isvalid,!istrans>::solve(LUx,P,m2,m3); 
        }
        template <class M2>
        void solveInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            BandLUHelper<isvalid,istrans>::solveInPlace(LUx,P,m2); 
        }
        template <class M2>
        void solveTransposeInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            BandLUHelper<isvalid,!istrans>::solveInPlace(LUx,P,m2); 
        }
        template <class M2>
        void makeInverse(M2& minv)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            BandLUHelper<isvalid,istrans>::makeInverse(LUx,P,minv);
        }
        template <class M2>
        void makeInverseATA(M2& ata)
        { BandLU_InverseATA(LUx,P,istrans,ata); }

        Mc SmallLUx;
        lux_type LUx;
        Permutation P;
    };
    
    template <class M>
    struct BandLUD_Impl<false,M>
    {
        enum { size = M::_colsize };
        enum { lo1 = M::_nlo };
        enum { hi1 = M::_nhi };
        enum { knownsizes = lo1 != Unknown && hi1 != Unknown };
        enum { istridiag1 = lo1 == 1 && hi1 == 1 };
        enum { istrans1 = knownsizes && (
                ptrdiff_t(hi1) < ptrdiff_t(lo1) || 
                (ptrdiff_t(hi1) == ptrdiff_t(lo1) && M::_rowmajor)) };
        enum { lo = istrans1 ? (lo1+hi1) : lo1 };
        enum { hi = istrans1 ? hi1 : (lo1+hi1) };
        typedef typename M::value_type T;
        enum { A = NoAlias | (istridiag1 ? DiagMajor : istrans1 ? RowMajor : ColMajor) };
        typedef typename BCopyHelper<T,Band,size,size,lo,hi,A>::type Mc;
        typedef typename TypeSelect< istrans1 ,
                typename Mc::transpose_type::noalias_type ,
                typename Mc::noalias_type>::type lux_type;

#define NEWLO TMV_MIN(A.nlo(),A.nhi())
#define NEWHI TMV_MIN(A.nlo()+A.nhi(),A.colsize()-1)
#define APTR1 (inplace ? 0 : \
               BandStorageLength(ColMajor,A.colsize(),A.colsize(),NEWLO,NEWHI))
#define TRID (A.nlo() == 1 && A.nhi() == 1)
#define APTR (inplace ? A.nonConst().ptr() : Aptr1.get())

#define LUX (istrans1 ? \
             (inplace ? \
              BandMatrixView<T>(A.nonConst().ptr(),A.colsize(),A.colsize(),\
                                A.nhi(),NEWHI,A.stepj(),A.stepi(),A.diagstep(),\
                                A.ct() \
                                TMV_FIRSTLAST1(A.nonConst()._first,\
                                               A.nonConst()._last) ) : \
              BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nhi(), \
                               NEWHI, TRID ? DiagMajor : ColMajor)) : \
             (inplace ? \
              BandMatrixView<T>(A.nonConst().ptr(),A.colsize(),\
                                A.colsize(),A.nlo(),NEWHI,\
                                A.stepi(),A.stepj(),A.diagstep(),\
                                A.ct() \
                                TMV_FIRSTLAST1(A.nonConst()._first,\
                                               A.nonConst()._last) ) : \
              BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nlo(), \
                               NEWHI, TRID ? DiagMajor : ColMajor)))

        template <class M2>
        BandLUD_Impl(const BaseMatrix_Band<M2>& A, bool _inplace) :
            istrans(A.nlo() < A.nhi() || (A.nlo() == A.nhi() && A.isrm())),
            istridiag(A.nlo() == 1 && A.nhi() == 1),
            inplace(NEWLO == 0 ||
                    (_inplace && (
                        (A.isrm() && istrans && A.stepi() >= A.nlo()+2*A.nhi()) || 
                        (A.iscm() && !istrans && A.stepj() >= 2*A.nlo()+A.nhi()) || 
                        (A.isdm() && istridiag)
                    ))
            ),
            // Aptr1 is the pointer to new storage if any
            Aptr1(APTR1),
            // Aptr is the pointer to the (0,0) element (which isn't necessarily Aptr1)
            Aptr(APTR),
            // LUx views this memory as the LU matrix
            LUx(
                inplace ? A.nonConst().ptr() : (
                    Aptr + (istridiag ? A.colsize() : 0) ), // ptr
                // A is square, so no need to check istrans for the right
                // values of rowsize,colsize here.
                A.colsize() ,  // colsize
                A.colsize() ,  // rowsize
                istrans ? A.nhi() : A.nlo(),  // nlo
                TMV_MIN(A.nlo()+A.nhi(),A.colsize()-1) , // nhi
                istridiag ? 1-ptrdiff_t(A.colsize()) : 1 ,  // stepi
                istridiag ? ptrdiff_t(A.colsize()) : 1 ,  // stepj
                // Here we do need to check istrans for the right step.
                ( inplace ? (istrans ? A.stepi() : A.stepj()) :
                  A.rowsize() ) // stepj
            ),
            // allocate memory for the permutation
            P(A.rowsize())
            {
                TMVStaticAssert(lux_type::_colmajor);
                TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
                TMVStaticAssert((Sizes<M2::_colsize,M2::_rowsize>::same));
                TMVAssert(A.colsize() == A.rowsize());
                if (!inplace) {
                    if (istrans) LUx = A.transpose();
                    else LUx = A;
                } else {
                    Maybe<M2::_conj>::conjself(LUx);
                }
                BandLU_Decompose(LUx,P);
#ifdef PRINTALGO_BandLU
                std::cout<<"BandLUD_Impl (not small) constructor\n";
                std::cout<<"this = "<<this<<std::endl;
                std::cout<<"istrans = "<<istrans<<std::endl;
                std::cout<<"inplace = "<<inplace<<std::endl;
                std::cout<<"Aptr = "<<Aptr<<"  "<<(inplace ? 0 : A.rowsize()*A.rowsize())<<std::endl;
                std::cout<<"A = "<<A<<std::endl;
                std::cout<<"LUx = "<<LUx<<std::endl;
                std::cout<<"P = "<<P<<std::endl;
#endif
            }

        template <class M2, class M3>
        void solve(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            if (istrans)
                BandLUHelper<isvalid,true>::solve(LUx,P,m2,m3); 
            else
                BandLUHelper<isvalid,false>::solve(LUx,P,m2,m3); 
        }
        template <class M2, class M3>
        void solveTranspose(const M2& m2, M3& m3)
        {
            const bool isvalid = (M::isreal && M2::isreal) || M3::iscomplex;
            if (istrans)
                BandLUHelper<isvalid,false>::solve(LUx,P,m2,m3); 
            else
                BandLUHelper<isvalid,true>::solve(LUx,P,m2,m3); 
        }
        template <class M2>
        void solveInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            if (istrans)
                BandLUHelper<isvalid,true>::solveInPlace(LUx,P,m2); 
            else
                BandLUHelper<isvalid,false>::solveInPlace(LUx,P,m2); 
        }
        template <class M2>
        void solveTransposeInPlace(M2& m2)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            if (istrans)
                BandLUHelper<isvalid,false>::solveInPlace(LUx,P,m2); 
            else
                BandLUHelper<isvalid,true>::solveInPlace(LUx,P,m2); 
        }
        template <class M2>
        void makeInverse(M2& minv)
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            if (istrans) 
                BandLUHelper<isvalid,true>::makeInverse(LUx,P,minv);
            else
                BandLUHelper<isvalid,false>::makeInverse(LUx,P,minv);
        }
        template <class M2>
        void makeInverseATA(M2& ata)
        { BandLU_InverseATA(LUx,P,istrans,ata); }


        const bool istrans,istridiag,inplace;
        AlignedArray<typename M::value_type> Aptr1;
        typename M::value_type* Aptr;
        lux_type LUx;
        Permutation P;
    };

    template <class M> template <class M2>
    BandLUD<M>::BandLUD(const BaseMatrix<M2>& A, bool inplace) :
        pimpl(new BandLUD_Impl<small,M>(A.mat(),inplace)) 
    {
#ifdef PRINTALGO_BandLU
        std::cout<<"Done BandLUD constructor for "<<this<<std::endl;
        std::cout<<"pimpl = "<<pimpl.get()<<std::endl;
        std::cout<<"A = "<<TMV_Text(A)<<"  "<<A.mat().cptr()<<std::endl;
        std::cout<<"inplace = "<<inplace<<std::endl;
        std::cout<<"small = "<<small<<std::endl;
        if (pimpl.get())
            std::cout<<"LUx = "<<TMV_Text(pimpl->LUx)<<"  "<<pimpl->LUx.cptr()<<std::endl;
#endif
    }

    template <class M>
    BandLUD<M>::BandLUD(const BandLUD<M>& rhs) : pimpl(rhs.pimpl.release()) 
    {
#ifdef PRINTALGO_BandLU
        std::cout<<"BandLUD copy constructor for "<<this<<std::endl;
        std::cout<<"pimpl = "<<pimpl.get()<<std::endl;
        if (pimpl.get())
            std::cout<<"LUx = "<<TMV_Text(pimpl->LUx)<<"  "<<pimpl->LUx.cptr()<<std::endl;
#endif
    }


    template <class M>
    BandLUD<M>::~BandLUD()
    {
#ifdef PRINTALGO_BandLU
        std::cout<<"BandLUD destructor for "<<this<<std::endl;
        std::cout<<"pimpl = "<<pimpl.get()<<std::endl;
        if (pimpl.get())
            std::cout<<"LUx = "<<TMV_Text(pimpl->LUx)<<"  "<<pimpl->LUx.cptr()<<std::endl;
#endif
    }

    template <class M> template <class M2>
    void BandLUD<M>::solveInPlace(BaseMatrix_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVAssert(m2.colsize() == colsize());
        pimpl->solveInPlace(m2.mat());
    }

    template <class M> template <class M2>
    void BandLUD<M>::solveTransposeInPlace(BaseMatrix_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == rowsize());
        pimpl->solveTransposeInPlace(m2.mat());
    }

    template <class M> template <class V2>
    void BandLUD<M>::solveInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVAssert(v2.size() == colsize());
        pimpl->solveInPlace(v2.vec());
    }

    template <class M> template <class V2>
    void BandLUD<M>::solveTransposeInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_rowsize>::same));
        TMVAssert(v2.size() == rowsize());
        pimpl->solveTransposeInPlace(v2.vec());
    }

    template <class M>
    typename M::value_type BandLUD<M>::det() const
    { return typename M::real_type(getP().det()) * getU().det(); }

    template <class M>
    typename M::float_type BandLUD<M>::logDet(typename M::zfloat_type* sign) const
    {
        typename M::float_type ret = getU().logDet(sign);
        if (sign) *sign *= typename M::float_type(getP().det());
        return ret;
    }                  

    template <class M>
    bool BandLUD<M>::isSingular() const 
    { return getU().isSingular(); }

    template <class M> template <class M2>
    void BandLUD<M>::makeInverse(BaseMatrix_Mutable<M2>& minv) const
    {
        TMVStaticAssert((Sizes<M::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M::_rowsize,M2::_colsize>::same));
        TMVAssert(minv.colsize() == rowsize());
        TMVAssert(minv.rowsize() == colsize());
        pimpl->makeInverse(minv.mat());
    }

    template <class M> template <class M2>
    void BandLUD<M>::makeInverseATA(BaseMatrix_Mutable<M2>& ata) const
    {
        TMVAssert(ata.rowsize() == rowsize());
        TMVAssert(ata.colsize() == rowsize());
        pimpl->makeInverseATA(ata.mat());
    }

    template <class M>
    bool BandLUD<M>::isTrans() const 
    { return pimpl->istrans; }

    template <class M>
    typename BandLUD<M>::getl_type BandLUD<M>::getL() const 
    {

        return BandLUD<M>::getl_type(pimpl->LUx.unitLowerTri()); 
    }

    template <class M>
    typename BandLUD<M>::getu_type BandLUD<M>::getU() const 
    { return pimpl->LUx.upperBand(); }

    template <class M>
    typename BandLUD<M>::getlu_type BandLUD<M>::getLU() const 
    {
#ifdef PRINTALGO_BandLU
        std::cout<<"BandLU<M>::getLU\n";
        std::cout<<"this = "<<this<<std::endl;
        std::cout<<"pimpl = "<<pimpl.get()<<std::endl;
        std::cout<<"pimpl->LUx.cptr() = "<<pimpl->LUx.cptr()<<std::endl;
        std::cout<<"LUx = "<<pimpl->LUx<<std::endl;
#endif
        return pimpl->LUx; 
    }

    template <class M>
    typename BandLUD<M>::getp_type BandLUD<M>::getP() const 
    { return pimpl->P; }

    template <class M>
    typename M::real_type BandLUD<M>::condition(RT normInf) const 
    {
        // FIXME: This is a placeholder until I write the real function.
        // Make sure to do this before releasing the code!
        // See page 129 of Golub and van Loan.
        //
        // This produces the exact right answer, but it is way too slow!
        // The GvL algorithm is order n^2.  This is order n^3.
        if (isSingular()) {
            return normInf / TMV_Epsilon<RT>();
        } else {
            Matrix<T> minv(rowsize(),colsize());
            makeInverse(minv);
            return normInf * minv.normInf();
        }
    }

    template <class M>
    ptrdiff_t BandLUD<M>::colsize() const
    { return pimpl->LUx.colsize(); }

    template <class M>
    ptrdiff_t BandLUD<M>::rowsize() const
    { return pimpl->LUx.rowsize(); }


    template <class M, class M2> 
    static bool CheckDecomp(
        const BandLUD<M>& lud, const BaseMatrix_Calc<M2>& m,
        std::ostream* fout=0) 
    {
        typedef typename M2::real_type RT;
        //bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
        bool printmat = fout;
        if (printmat) {
            *fout << "LU:\n";
            if (lud.isTrans()) *fout << m.transpose() << std::endl;
            else *fout << m << std::endl;
            *fout << "L = "<<lud.getL()<<std::endl;
            *fout << "U = "<<lud.getU()<<std::endl;
        }
        typename M::copy_type lu = lud.getL()*lud.getU();
        if (printmat) {
            *fout << "LU = "<<lu<<std::endl;
        }
        //lu = getP() * lu;
        lud.getP().applyOnLeft(lu);
        if (printmat) {
            *fout << "PLU = "<<lu<<std::endl;
        }
        RT nm = lud.isTrans() ? Norm(lu-m.transpose()) : Norm(lu-m);
        nm /= Norm(lud.getL())*Norm(lud.getU());
        if (fout) {
            *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<" <? ";
            *fout << RT(m.colsize())<<"*"<<TMV_Epsilon<RT>();
            *fout << " = "<<RT(m.colsize())*TMV_Epsilon<RT>()<<std::endl;
        }
        return nm < RT(m.colsize())*TMV_Epsilon<RT>();
    }

} // namespace tmv

#endif
