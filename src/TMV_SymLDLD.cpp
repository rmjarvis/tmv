///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "tmv/TMV_SymLDLD.h"
#include "TMV_SymLDLDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "TMV_SymSquare.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_PermutationArith.h"
#include <ostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T>
    struct SymLDLDiv<T>::SymLDLDiv_Impl
    {
    public :
        SymLDLDiv_Impl(const GenSymMatrix<T>& m, bool inplace);

        const bool inplace;
        AlignedArray<T> Aptr1;
        SymMatrixView<T> LLx;
        Vector<T> xD;
        Permutation P;
        mutable RT logdet;
        mutable T signdet;

        const GenSymMatrix<T>& A0;
    };

#define APTR1 (inplace ? 0 : (A.size()*A.size()))
#define LLX \
    (inplace ? \
     (A.uplo()==Upper ? \
      (A.isherm() ? A.nonConst().adjoint() : A.nonConst().transpose()) : \
      A.nonConst()) : \
     SymMatrixView<T>(Aptr1.get(),A.size(),1,A.size(),A.isherm()?Herm:Sym,Lower,NonConj\
                      TMV_FIRSTLAST1(Aptr1.get(),Aptr1.get()+A.size()*A.size())))

    template <class T>
    SymLDLDiv<T>::SymLDLDiv_Impl::SymLDLDiv_Impl(
        const GenSymMatrix<T>& A, bool _inplace
    ) :
        inplace(_inplace && (A.isrm() || A.iscm())), 
        Aptr1(APTR1), LLx(LLX), xD(A.size()-1),
        P(A.colsize()), logdet(0), signdet(1), A0(A) {}

#undef APTR1
#undef LLX

    template <class T>
    SymLDLDiv<T>::SymLDLDiv(const GenSymMatrix<T>& A, bool inplace) :
        pimpl(new SymLDLDiv_Impl(A,inplace))
    {
        if (inplace) TMVAssert(A == pimpl->LLx); 
        else pimpl->LLx = A;

        LDL_Decompose(pimpl->LLx,pimpl->xD.view(),pimpl->P,
                      pimpl->logdet,pimpl->signdet);
    }

    template <class T>
    SymLDLDiv<T>::~SymLDLDiv() {}

    template <class T> template <class T1> 
    void SymLDLDiv<T>::doLDivEq(MatrixView<T1> m) const
    { LDL_LDivEq(pimpl->LLx,pimpl->xD,pimpl->P.getValues(),m); }

    template <class T> template <class T1> 
    void SymLDLDiv<T>::doRDivEq(MatrixView<T1> m) const
    { LDL_RDivEq(pimpl->LLx,pimpl->xD,pimpl->P.getValues(),m); }

    template <class T> template <class T1, class T2> 
    void SymLDLDiv<T>::doLDiv(
        const GenMatrix<T1>& m1, MatrixView<T2> m0) const
    { LDL_LDivEq(pimpl->LLx,pimpl->xD,pimpl->P.getValues(),m0=m1); }

    template <class T> template <class T1, class T2> 
    void SymLDLDiv<T>::doRDiv(
        const GenMatrix<T1>& m1, MatrixView<T2> m0) const
    { LDL_RDivEq(pimpl->LLx,pimpl->xD,pimpl->P.getValues(),m0=m1); }

    template <class T>
    T SymLDLDiv<T>::det() const 
    {
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }

    template <class T>
    RT SymLDLDiv<T>::logDet(T* sign) const 
    { 
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;
    }

    template <class T> template <class T1> 
    void SymLDLDiv<T>::doMakeInverse(SymMatrixView<T1> sinv) const
    { 
        TMVAssert(isReal(T()) || issym() == sinv.issym());
        TMVAssert(isReal(T()) || isherm() == sinv.isherm());
        LDL_Inverse(pimpl->LLx,pimpl->xD,pimpl->P.getValues(),sinv); 
    }

    template <class T> template <class T1> 
    void SymLDLDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    {
        if (pimpl->LLx.isherm()) {
            if (isComplex(T1())) minv.diag().imagPart().setZero();
            doMakeInverse(HermMatrixViewOf(minv,Lower));
            if (minv.colsize() > 1)
                minv.upperTri().offDiag() =
                    minv.lowerTri().offDiag().adjoint();
        } else {
            doMakeInverse(SymMatrixViewOf(minv,Lower));
            if (minv.colsize() > 1)
                minv.upperTri().offDiag() =
                    minv.lowerTri().offDiag().transpose();
        }
    }

    template <class T>
    void SymLDLDiv<T>::doMakeInverseATA(MatrixView<T> ata) const
    {
        TMVAssert(ata.colsize() == pimpl->LLx.size());
        TMVAssert(ata.rowsize() == pimpl->LLx.size());

        if (pimpl->LLx.isherm()) {
            if (isComplex(T())) ata.diag().imagPart().setZero();
            doMakeInverse(HermMatrixViewOf(ata,Lower));
            SymSquare<true>(ata);
        } else {
            doMakeInverse(SymMatrixViewOf(ata,Lower));
            SymSquare<false>(ata);
        }
    }

    template <class T>
    bool SymLDLDiv<T>::isSingular() const 
    {
        return pimpl->signdet == T(0) ||
            ( pimpl->LLx.diag().minAbs2Element() <=
              TMV_Epsilon<T>() * pimpl->LLx.diag().maxAbs2Element() ); 
    }

    template <class T>
    const ConstLowerTriMatrixView<T> SymLDLDiv<T>::getL() const
    { return pimpl->LLx.lowerTri().viewAsUnitDiag(); }

    template <class T>
    const BandMatrix<T> SymLDLDiv<T>::getD() const 
    { 
        BandMatrix<T> temp(pimpl->LLx.size(),pimpl->LLx.size(),1,1);
        temp.diag() = pimpl->LLx.diag();
        temp.diag(-1) = pimpl->xD;
        temp.diag(1) = (pimpl->LLx.isherm() ? pimpl->xD.conjugate() : 
                        pimpl->xD.view());
        return temp;
    }

    template <class T>
    const Permutation& SymLDLDiv<T>::getP() const 
    { return pimpl->P; }

    template <class T>
    const GenSymMatrix<T>& SymLDLDiv<T>::getLL() const 
    { return pimpl->LLx; }

    template <class T>
    const GenVector<T>& SymLDLDiv<T>::getxD() const 
    { return pimpl->xD; }

    template <class T>
    bool SymLDLDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "SymLDLDiv:\n";
            *fout << "M = "<<mm<<std::endl;
            *fout << "L = "<<getL()<<std::endl;
            *fout << "D = "<<getD()<<std::endl;
            *fout << "P = "<<getP()<<std::endl;
            *fout << "  or by interchanges: ";
            for(ptrdiff_t i=0;i<getP().size();i++)
                *fout<<(getP().getValues())[i]<<" ";
            *fout <<std::endl;
        }
        Matrix<T> lu = getP()*getL()*getD()*
            (pimpl->LLx.isherm()?getL().adjoint():getL().transpose())*
            getP().transpose();
        RT nm = Norm(lu-mm);
        nm /= TMV_SQR(Norm(getL()))*Norm(getD());
        if (fout) {
            *fout << "LDLt = "<<lu<<std::endl;
            *fout << "Norm(M-LDLt)/Norm(LDLt) = "<<nm<<std::endl;
        }
        return nm < mm.doCondition()*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T>
    ptrdiff_t SymLDLDiv<T>::colsize() const
    { return pimpl->LLx.size(); }

    template <class T>
    ptrdiff_t SymLDLDiv<T>::rowsize() const
    { return pimpl->LLx.size(); }

    template <class T>
    bool SymLDLDiv<T>::isherm() const
    { return pimpl->LLx.isherm(); }

    template <class T>
    bool SymLDLDiv<T>::issym() const
    { return pimpl->LLx.issym(); }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymLDLD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


