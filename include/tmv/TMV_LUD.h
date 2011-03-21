///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------
//
// This file contains the driver code for doing division using 
// LU Decomposition.
//
// The basic idea of an LU decomposition is that any 
// square matrix A can be decomposed into a lower triangular
// matrix times an upper triangular matrix.
//
// For stability reasons, we actually decompose a permutation
// of A instead, so:
//
// A = P L U
//
// Only one of L or U needs a non-unit diagonal, so we choose L to 
// have unit diagonal, and U to have the non-unit diagonal.
//
// This means that we can store L and U both in a square matrix
// the same size as A, with L being the elements below the diagonal
// and U being the elements above and including the diagonal.
//
// The determinant of A can be calculated easily from the LU
// decomposition:
//
// det(A) = det(P) * det(L) * det(U)
// det(A) = (+-1) * 1 * det(U)
// As we calculate the decomposition, we keep track of whether
// det(P) is +-1 
// The determinant of U is just the product of the diagonal elements.
// So the determinant of A is just det(P) times the diagonal elements 
// of U.
//


#ifndef TMV_LU_H
#define TMV_LU_H

#include "TMV_Divider.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Matrix.h"
#include "TMV_Permutation.h"
#include "TMV_CopyM.h"

#include "TMV_LUDiv.h"
#include "TMV_LUDecompose.h"
#include "TMV_LUInverse.h"


namespace tmv {

    // The point of the Impl class here is to implement the transfer of
    // ownership copy semantics.
    template <bool small, class M>
    struct LUD_Impl;

    template <class M>
    class LUD 
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
                M::_colsize != UNKNOWN && M::_rowsize != UNKNOWN &&
                M::_colsize <= 8 && M::_rowsize <= 8 ) };

        typedef typename LUD_Impl<small,M>::lux_type lux_type;
        typedef typename lux_type::const_view_type getlu_type;
        typedef typename lux_type::const_unit_lowertri_type getl_type;
        typedef typename lux_type::const_uppertri_type getu_type;
        typedef Permutation getp_type;

        //
        // Constructors
        //

        // Sets up the internal storage and does the decomposition.
        template <class M2>
        LUD(const BaseMatrix<M2>& A, bool _inplace=false);

        // The copy constructor has transfer of ownership semantics.
        // This way an LUD object can be returned by value, and the 
        // copy is cheap.  I don't think there is any reason to use
        // a more sophisticated technique like shared_ptr or something
        // similar, since there are no non-const methods.
        LUD(const LUD<M>& rhs);

        // Clean up the internal storage
        ~LUD();


        //
        // Perform the division in place
        //
        template <class M2>
        void solveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class V2>
        void solveInPlace(BaseVector_Mutable<V2>& v2) const;

        template <class M2>
        void solveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class V2>
        void solveTransposeInPlace(BaseVector_Mutable<V2>& v2) const;

        //
        // Next the not-in-place division 
        // For LU, we just copy m1->m3 and do the in-place version.
        // 
        
        template <class M1, class M2>
        void solve(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const
        { solveInPlace(m2=m1); }
        template <class V1, class V2>
        void solve(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const
        { solveInPlace(v2=v1); }

        template <class M1, class M2>
        void solveTranspose(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const
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
        void makeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const;


        //
        // InverseATA
        //
        
        template <class M2>
        void makeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const;


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
        // const LUD<M>& can relase the memory.
        mutable std::auto_ptr<LUD_Impl<small,M> > pimpl;

        size_t colsize() const;
        size_t rowsize() const;

        // op= not allowed.
        LUD<M>& operator=(const LUD<M>&);
    };
    
    template <class T>
    class InstLUD :
        public LUD<Matrix<T,ColMajor> >,
        public Divider<T>
    {
    public :
        typedef LUD<Matrix<T,ColMajor> > base;
        typedef typename base::RT RT;
        typedef typename base::CT CT;
        typedef typename base::FT FT;
        typedef typename base::ZFT ZFT;

        // Sets up the internal storage and does the decomposition.
        template <class M2>
        InstLUD(const BaseMatrix<M2>& A, bool _inplace=false);
        InstLUD(const InstLUD<T>& rhs);
        ~InstLUD();

        // These are the virtual functions from the Divider base class.
        void doSolveInPlace(MatrixView<RT> m2) const;
        void doSolveInPlace(MatrixView<CT> m2) const;
        void doSolveInPlace(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const;
        void doSolveInPlace(VectorView<RT> v2) const;
        void doSolveInPlace(VectorView<CT> v2) const;
        void doSolveInPlace(VectorView<CT,UNKNOWN,true> v2) const;

        void doSolveTransposeInPlace(MatrixView<RT> m2) const;
        void doSolveTransposeInPlace(MatrixView<CT> m2) const;
        void doSolveTransposeInPlace(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const;
        void doSolveTransposeInPlace(VectorView<RT> v2) const;
        void doSolveTransposeInPlace(VectorView<CT> v2) const;
        void doSolveTransposeInPlace(
            VectorView<CT,UNKNOWN,true> v2) const;

        void doSolve(
            const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const;
        void doSolve(
            const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const;
        void doSolve(
            const ConstMatrixView<RT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const;
        void doSolve(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const;
        void doSolve(
            const ConstMatrixView<CT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const;
        void doSolve(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
            MatrixView<CT> m2) const;
        void doSolve(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const;
        void doSolve(
            const ConstVectorView<RT>& v1, VectorView<RT> v2) const;
        void doSolve(
            const ConstVectorView<RT>& v1, VectorView<CT> v2) const;
        void doSolve(
            const ConstVectorView<RT>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const;
        void doSolve(
            const ConstVectorView<CT>& v1, VectorView<CT> v2) const;
        void doSolve(
            const ConstVectorView<CT>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const;
        void doSolve(
            const ConstVectorView<CT,UNKNOWN,true>& v1, 
            VectorView<CT> v2) const;
        void doSolve(
            const ConstVectorView<CT,UNKNOWN,true>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const;

        void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<RT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<CT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
            MatrixView<CT> m2) const;
        void doSolveTranspose(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const;
        void doSolveTranspose(
            const ConstVectorView<RT>& v1, VectorView<RT> v2) const;
        void doSolveTranspose(
            const ConstVectorView<RT>& v1, VectorView<CT> v2) const;
        void doSolveTranspose(
            const ConstVectorView<RT>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const;
        void doSolveTranspose(
            const ConstVectorView<CT>& v1, VectorView<CT> v2) const;
        void doSolveTranspose(
            const ConstVectorView<CT>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const;
        void doSolveTranspose(
            const ConstVectorView<CT,UNKNOWN,true>& v1, 
            VectorView<CT> v2) const;
        void doSolveTranspose(
            const ConstVectorView<CT,UNKNOWN,true>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const;

        T det() const;
        FT logDet(ZFT* sign) const;
        bool isSingular() const;

        void doMakeInverse(MatrixView<RT> minv) const;
        void doMakeInverse(MatrixView<CT> minv) const;
        void doMakeInverse(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> minv) const;

        void doMakeInverseATA(MatrixView<RT> ata) const;
        void doMakeInverseATA(MatrixView<CT> ata) const;
        void doMakeInverseATA(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> ata) const;
        
        RT condition(RT normInf) const;
        bool preferInPlace() const;

    private :
        // op= not allowed.
        Divider<T>& operator=(const Divider<T>&);
    };
    
    // Now the instantiation and definition.
    template <bool isvalid>
    struct LUHelper;

    template <>
    struct LUHelper<true>
    {
        template <class M1, class M2>
        static void makeInverse(
            const M1& LUx, const Permutation& P, bool istrans, M2& m2)
        {
            // This one might be same storage if LUx was done in place.
            // e.g. { A.divideInPlace(); A = A.inverse(); }
            // So go ahead and check. (i.e. Don't use NoAliasCopy.)
            if (istrans) {
                typename M2::transpose_type m2t = m2.transpose();
                Copy(LUx,m2t);
                LU_Inverse(m2t,P);
            } else {
                Copy(LUx,m2);
                LU_Inverse(m2,P);
            }
        }
        template <class M1, class M2>
        static void solveInPlace(
            const M1& LUx, const Permutation& P, bool istrans, M2& m2)
        {
            if (istrans) LU_SolveTransposeInPlace(LUx,P,m2);
            else LU_SolveInPlace(LUx,P,m2); 
        }
    };
    template <>
    struct LUHelper<false>
    {
        template <class M1, class M2>
        static void makeInverse(const M1& , const Permutation& , bool, M2& ) {}
        template <class M1, class M2>
        static void solveInPlace(const M1& , const Permutation& , bool, M2& ) {}
    };

    template <class M>
    struct LUD_Impl<true,M>
    {
        enum { istrans = M::_rowmajor };
        enum { size = M::_colsize };
        typedef typename TypeSelect< istrans ,
                typename M::transpose_type ,
                typename M::view_type>::type lux_type;

        template <class M2>
        LUD_Impl(const BaseMatrix<M2>& A, bool ) : 
            LUx( Maybe<istrans>::transposeview(SmallLUx) ),
            P(size)
        {
            TMVStaticAssert(M::_colsize != UNKNOWN);
            TMVStaticAssert(M::_rowsize != UNKNOWN);
            TMVStaticAssert(M::_colsize == int(M::_rowsize));
            TMVStaticAssert(M::_colsize == int(M2::_colsize));
            TMVStaticAssert(M::_rowsize == int(M2::_rowsize));
            TMVStaticAssert(lux_type::_colmajor);
            A.newAssignTo(SmallLUx);
        }

        M SmallLUx;
        lux_type LUx;
        Permutation P;
    };
    
    template <class M>
    struct LUD_Impl<false,M>
    {
        enum { istrans1 = M::_rowmajor };
        typedef typename TypeSelect< istrans1 ,
                typename M::transpose_type ,
                typename M::view_type>::type lux_type;

        template <class M2>
        LUD_Impl(const BaseMatrix_Rec<M2>& A, bool _inplace) :
            // if A is rm, copy it to the transpose of LUx (which is cm)
            istrans(A.isrm()),
            // inplace only if matrix is rowmajor or colmajor
            inplace((A.iscm() || A.isrm()) && _inplace),
            // Aptr is the pointer to new storage if any
            Aptr( inplace ? 0 : A.rowsize()*A.rowsize() ),
            // LUx views this memory as the LU matrix
            LUx(
                inplace ? A.nonConst().ptr() : Aptr.get() , // ptr
                // A is square, so no need to check istrans for the right
                // values of rowsize,colsize here.
                A.rowsize() ,  // colsize
                A.rowsize() ,  // rowsize
                1 ,  // stepi
                // Here we do need to check istrans for the right step.
                ( inplace ? (istrans ? A.stepi() : A.stepj()) :
                  int(A.rowsize()) ) // stepj
            ),
            // allocate memory for the permutation
            P(A.rowsize())
            {
                TMVStaticAssert(lux_type::_colmajor);
                TMVAssert(A.colsize() == A.rowsize());
                if (!inplace) {
                    if (istrans) NoAliasCopy(A.transpose(),LUx); 
                    else NoAliasCopy(A,LUx); 
                } else {
                    Maybe<M2::_conj>::conjself(LUx);
                }
            }

        // If A is not a BaseMatrix_Rec, can't do it in place.
        template <class M2>
        LUD_Impl(const BaseMatrix<M2>& A, bool ) :
            istrans(false), inplace(false),
            Aptr( A.rowsize()*A.rowsize() ),
            LUx(Aptr.get(),A.rowsize(), A.rowsize(),1,A.rowsize()),
            P(A.rowsize())
        {
            TMVStaticAssert(lux_type::_colmajor);
            TMVAssert(A.colsize() == A.rowsize());
            A.newAssignTo(LUx);
        }

        const bool istrans,inplace;
        AlignedArray<typename M::value_type> Aptr;
        lux_type LUx;
        Permutation P;
    };

    template <class M> template <class M2>
    LUD<M>::LUD(const BaseMatrix<M2>& A, bool inplace) :
        pimpl(new LUD_Impl<small,M>(A,inplace)) 
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(A.mat().isSquare());
        LU_Decompose(pimpl->LUx,pimpl->P);
    }

    template <class M>
    LUD<M>::LUD(const LUD<M>& rhs) : pimpl(rhs.pimpl.release()) {}

    template <class M>
    LUD<M>::~LUD() {}

    template <class M> template <class M2>
    void LUD<M>::solveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVAssert(m2.colsize() == colsize());
        const bool isvalid = M::isreal || M2::iscomplex;
        LUHelper<isvalid>::solveInPlace(pimpl->LUx,pimpl->P,pimpl->istrans,m2);
    }

    template <class M> template <class M2>
    void LUD<M>::solveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == rowsize());
        const bool isvalid = M::isreal || M2::iscomplex;
        LUHelper<isvalid>::solveInPlace(pimpl->LUx,pimpl->P,!pimpl->istrans,m2);
    }

    template <class M> template <class V2>
    void LUD<M>::solveInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVAssert(v2.size() == colsize());
        const bool isvalid = M::isreal || V2::iscomplex;
        LUHelper<isvalid>::solveInPlace(pimpl->LUx,pimpl->P,pimpl->istrans,v2);
    }

    template <class M> template <class V2>
    void LUD<M>::solveTransposeInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_rowsize>::same));
        TMVAssert(v2.size() == rowsize());
        const bool isvalid = M::isreal || V2::iscomplex;
        LUHelper<isvalid>::solveInPlace(pimpl->LUx,pimpl->P,!pimpl->istrans,v2);
    }

    template <class M>
    typename M::value_type LUD<M>::det() const
    { return typename M::real_type(getP().det()) * getU().det(); }

    template <class M>
    typename M::float_type LUD<M>::logDet(typename M::zfloat_type* sign) const
    {
        typename M::float_type ret = getU().logDet(sign);
        if (sign) *sign *= getP().det();
        return ret;
    }                  

    template <class M>
    bool LUD<M>::isSingular() const 
    { return getU().isSingular(); }

    template <class M> template <class M2>
    void LUD<M>::makeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const
    {
        TMVAssert(minv.colsize() == rowsize());
        TMVAssert(minv.rowsize() == colsize());
        const bool isvalid = M::isreal || M2::iscomplex;
        LUHelper<isvalid>::makeInverse(pimpl->LUx,pimpl->P,pimpl->istrans,minv);
    }

    template <class M> template <class M2>
    void LUD<M>::makeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const
    {
        TMVAssert(ata.rowsize() == rowsize());
        TMVAssert(ata.colsize() == rowsize());
        LU_InverseATA(pimpl->LUx,pimpl->P,pimpl->istrans,ata);
    }

    template <class M>
    bool LUD<M>::isTrans() const 
    { return pimpl->istrans; }

    template <class M>
    typename LUD<M>::getl_type LUD<M>::getL() const 
    { return pimpl->LUx.unitLowerTri(); }

    template <class M>
    typename LUD<M>::getu_type LUD<M>::getU() const 
    { return pimpl->LUx.upperTri(); }

    template <class M>
    typename LUD<M>::getlu_type LUD<M>::getLU() const 
    { return pimpl->LUx; }

    template <class M>
    typename LUD<M>::getp_type LUD<M>::getP() const 
    { return pimpl->P; }

    template <class M>
    typename M::real_type LUD<M>::condition(RT normInf) const 
    {
        // FIXME: This is a placeholder until I write the real function.
        // Make sure to do this before releasing the code!
        // See page 129 of Golub and van Loan.
        //
        // This produces the exact right answer, but it is way too slow!
        // The GvL algorithm is order n^2.  This is order n^3.
        Matrix<T> minv(rowsize(),colsize());
        makeInverse(minv);
        return normInf * minv.normInf();
    }

    template <class M>
    size_t LUD<M>::colsize() const
    { return pimpl->LUx.colsize(); }

    template <class M>
    size_t LUD<M>::rowsize() const
    { return pimpl->LUx.rowsize(); }

    template <class M, class M2> 
    static bool CheckDecomp(
        const LUD<M>& lud, const BaseMatrix_Calc<M2>& m, std::ostream* fout=0) 
    {
        typedef typename M2::real_type RT;
        bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
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
        RT kappa = lud.condition(m.normInf());
        if (fout) {
            *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<" <? ";
            *fout << kappa<<"*"<<RT(m.colsize())<<"*"<<TMV_Epsilon<RT>();
            *fout << " = "<<kappa*RT(m.colsize())*TMV_Epsilon<RT>()<<std::endl;
        }
        return nm < kappa*RT(m.colsize())*TMV_Epsilon<RT>();
    }

} // namespace tmv

#endif
