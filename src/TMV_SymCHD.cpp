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



#include "tmv/TMV_SymCHD.h"
#include "TMV_SymCHDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "TMV_SymSquare.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include <ostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T>
    struct HermCHDiv<T>::HermCHDiv_Impl
    {
    public :
        HermCHDiv_Impl(const GenSymMatrix<T>& m, bool inplace);

        const bool inplace;
        AlignedArray<T> Aptr1;
        SymMatrixView<T> LLx;
        mutable bool zerodet;
        mutable RT logdet;
        mutable bool donedet;
    };

#define APTR1 (inplace ? 0 : (A.size()*A.size()))
#define LLX \
    (inplace ? (A.uplo()==Upper ? A.nonConst().adjoint() : A.nonConst()) : \
     SymMatrixView<T>(Aptr1.get(),A.size(),1,A.size(),Herm,Lower,NonConj \
                      TMV_FIRSTLAST1(Aptr1.get(),Aptr1.get()+A.size()*A.size())))

    template <class T>
    HermCHDiv<T>::HermCHDiv_Impl::HermCHDiv_Impl(
        const GenSymMatrix<T>& A, bool _inplace) :
        inplace(_inplace && (A.iscm() || A.isrm())), 
        Aptr1(APTR1), LLx(LLX), 
        zerodet(false), logdet(0), donedet(false) {}

#undef APTR1
#undef LLX

    template <class T>
    HermCHDiv<T>::HermCHDiv(const GenSymMatrix<T>& A, bool inplace) :
        pimpl(new HermCHDiv_Impl(A,inplace))
    {
        TMVAssert(isReal(T()) || A.isherm());
        if (inplace) TMVAssert(A==pimpl->LLx); 
        else pimpl->LLx = A;
#ifndef NOTHROW
        try {
#endif
            CH_Decompose(pimpl->LLx);
#ifndef NOTHROW
        } catch (NonPosDef) {
            if (inplace) throw NonPosDefHermMatrix<T>(pimpl->LLx);
            else throw NonPosDefHermMatrix2<T>(pimpl->LLx,A);
        }
#endif
    }

    template <class T>
    HermCHDiv<T>::~HermCHDiv() {}

    template <class T> template <class T1> 
    void HermCHDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(pimpl->LLx.size() == m.colsize());
        CH_LDivEq(pimpl->LLx,m);
    }

    template <class T> template <class T1> 
    void HermCHDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(pimpl->LLx.size() == m.rowsize());
        CH_RDivEq(pimpl->LLx,m);
    }

    template <class T> template <class T1, class T2> 
    void HermCHDiv<T>::doLDiv(
        const GenMatrix<T1>& m1, MatrixView<T2> m0) const
    {
        TMVAssert(m1.colsize() == m0.colsize());
        TMVAssert(m1.rowsize() == m0.rowsize());
        TMVAssert(pimpl->LLx.size() == m1.colsize());
        CH_LDivEq(pimpl->LLx,m0=m1);
    }

    template <class T> template <class T1, class T2> 
    void HermCHDiv<T>::doRDiv(
        const GenMatrix<T1>& m1, MatrixView<T2> m0) const
    {
        TMVAssert(m1.colsize() == m0.colsize());
        TMVAssert(m1.rowsize() == m0.rowsize());
        TMVAssert(pimpl->LLx.size() == m1.rowsize());
        CH_RDivEq(pimpl->LLx,m0=m1);
    }

    template <class T>
    T HermCHDiv<T>::det() const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->LLx.diag()).logDet(&s);
            pimpl->logdet *= RT(2);
            pimpl->zerodet = s == T(0);
            pimpl->donedet = true;
        }         
        if (pimpl->zerodet) return T(0);
        else return TMV_EXP(pimpl->logdet);  
    }                  

    template <class T>
    RT HermCHDiv<T>::logDet(T* sign) const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->LLx.diag()).logDet(&s);
            pimpl->logdet *= RT(2);
            pimpl->zerodet = s == T(0);
            pimpl->donedet = true;
        }
        if (sign) {
            if (pimpl->zerodet) *sign = T(0);
            else *sign = T(1);
        }
        return pimpl->logdet;  
    }                  

    template <class T> template <class T1> 
    void HermCHDiv<T>::doMakeInverse(SymMatrixView<T1> sinv) const
    {
        TMVAssert(sinv.size() == pimpl->LLx.size());
        TMVAssert(sinv.isherm());
        CH_Inverse(pimpl->LLx,sinv);
    }

    template <class T> template <class T1> 
    void HermCHDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    {
        TMVAssert(minv.colsize() == pimpl->LLx.size());
        TMVAssert(minv.rowsize() == pimpl->LLx.size());

        if (isComplex(T1())) minv.diag().imagPart().setZero();
        doMakeInverse(HermMatrixViewOf(minv,Lower));
        if (minv.colsize() > 1)
            minv.upperTri().offDiag() = minv.lowerTri().offDiag().adjoint();
    }

    template <class T>
    void HermCHDiv<T>::doMakeInverseATA(MatrixView<T> ata) const
    {
        // ata = (At A)^-1 = A^-1 (A^-1)t
        //     = A^-1 A^-1
        if (isComplex(T())) ata.diag().imagPart().setZero();
        SymMatrixView<T> hermata = HermMatrixViewOf(ata,Lower);
        doMakeInverse(hermata);
        SymSquare<true>(ata);
    }

    template <class T>
    bool HermCHDiv<T>::isSingular() const 
    {
        return pimpl->LLx.diag().minAbs2Element() <=
            TMV_Epsilon<T>() * pimpl->LLx.diag().maxAbs2Element(); 
    }

    template <class T> 
    const ConstLowerTriMatrixView<T> HermCHDiv<T>::getL() const 
    { return pimpl->LLx.lowerTri(); }

    template <class T>
    const GenSymMatrix<T>& HermCHDiv<T>::getLL() const 
    { return pimpl->LLx; }

    template <class T>
    bool HermCHDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "HermCHDiv:\n";
            *fout << "M = "<<mm<<std::endl;
            *fout << "L = "<<getL()<<std::endl;
        }
        Matrix<T> lu = getL()*getL().adjoint();
        RT nm = Norm(lu-mm);
        nm /= TMV_SQR(Norm(getL()));
        if (fout) {
            *fout << "LLt = "<<lu<<std::endl;
            *fout << "Norm(M-LLt)/Norm(LLt) = "<<nm<<std::endl;
        }
        return nm < mm.doCondition()*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T>
    ptrdiff_t HermCHDiv<T>::colsize() const
    { return pimpl->LLx.size(); }

    template <class T>
    ptrdiff_t HermCHDiv<T>::rowsize() const
    { return pimpl->LLx.size(); }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymCHD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv
