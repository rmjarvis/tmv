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



#include "tmv/TMV_SymBandCHD.h"
#include "TMV_SymBandCHDiv.h"
#include "tmv/TMV_SymBandMatrix.h"
#include "TMV_SymSquare.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include <iostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    struct HermBandCHDiv<T>::HermBandCHDiv_Impl
    {
    public :
        HermBandCHDiv_Impl(const GenSymBandMatrix<T>& m, bool inplace);

        const bool inplace;
        AlignedArray<T> Aptr1;
        T* Aptr;
        SymBandMatrixView<T> LLx;
        mutable bool zerodet;
        mutable RT logdet;
        mutable bool donedet;
    };

#define APTR1 (inplace ? 0 : \
               BandStorageLength(ColMajor,A.size(),A.size(),A.nlo(),0))
#define TRID (A.nlo() == 1)
#define APTR (inplace ? A.nonConst().ptr() : Aptr1.get())
#define LLX \
    (inplace ? (A.uplo()==Upper ? A.nonConst().adjoint() : A.nonConst()) : \
     HermBandMatrixViewOf(Aptr,A.size(),A.nlo(),Lower, \
                          (TRID ? DiagMajor : ColMajor)))

    template <class T> 
    HermBandCHDiv<T>::HermBandCHDiv_Impl::HermBandCHDiv_Impl(
        const GenSymBandMatrix<T>& A, bool _inplace) :
        inplace( ( _inplace && 
                   ( ((A.iscm() || A.isrm()) && A.nlo()>1) || 
                     (A.isdm() && TRID)
                   ) ) || A.nlo()==0 ),
        Aptr1(APTR1), Aptr(APTR), LLx(LLX),
        zerodet(false), logdet(1), donedet(false) 
    {}

#undef APTR
#undef APTR1
#undef LLX

    template <class T> 
    HermBandCHDiv<T>::HermBandCHDiv(
        const GenSymBandMatrix<T>& A, bool inplace) :
        pimpl(new HermBandCHDiv_Impl(A,inplace))
    {
        TMVAssert(isReal(T()) || A.isherm());
        if (inplace) TMVAssert(A==pimpl->LLx); 
        else pimpl->LLx = A;
#ifndef NOTHROW
        try {
#endif
            if (A.nlo() > 1)
                CH_Decompose(pimpl->LLx);
            else if (A.nlo() == 1)
                LDL_Decompose(pimpl->LLx);
            else {
                if (A.realPart().diag().minElement() <= RT(0)) {
#ifdef NOTHROW
                    std::cerr<<"Non Posdef HermBandMatrix found\n";
                    exit(1); 
#else
                    throw NonPosDef();
#endif
                }
            }
#ifndef NOTHROW
        } catch (NonPosDef) {
            if (inplace) throw NonPosDefHermBandMatrix<T>(pimpl->LLx);
            else throw NonPosDefHermBandMatrix2<T>(pimpl->LLx,A);
        }
#endif
        TMVAssert(pimpl->LLx.isHermOK());
    }

    template <class T> 
    HermBandCHDiv<T>::~HermBandCHDiv() {}

    template <class T> template <class T1> 
    void HermBandCHDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(pimpl->LLx.size() == m.colsize());
        if (pimpl->LLx.nlo() > 1)
            CH_LDivEq(pimpl->LLx,m);
        else if (pimpl->LLx.nlo() == 1)
            LDL_LDivEq(pimpl->LLx,m);
        else
            DiagMatrixViewOf(pimpl->LLx.diag()).LDivEq(m);
    }

    template <class T> template <class T1> 
    void HermBandCHDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(pimpl->LLx.size() == m.rowsize());
        if (pimpl->LLx.nlo() > 1)
            CH_RDivEq(pimpl->LLx,m);
        else if (pimpl->LLx.nlo() == 1)
            LDL_RDivEq(pimpl->LLx,m);
        else
            DiagMatrixViewOf(pimpl->LLx.diag()).RDivEq(m);
    }

    template <class T> template <class T1, class T2> 
    void HermBandCHDiv<T>::doLDiv(
        const GenMatrix<T1>& m1, MatrixView<T2> m0) const
    {
        TMVAssert(m1.colsize() == m0.colsize());
        TMVAssert(m1.rowsize() == m0.rowsize());
        TMVAssert(pimpl->LLx.size() == m1.colsize());
        if (pimpl->LLx.nlo() > 1)
            CH_LDivEq(pimpl->LLx,m0=m1);
        else if (pimpl->LLx.nlo() == 1)
            LDL_LDivEq(pimpl->LLx,m0=m1);
        else
            DiagMatrixViewOf(pimpl->LLx.diag()).LDiv(m1,m0);
    }

    template <class T> template <class T1, class T2> 
    void HermBandCHDiv<T>::doRDiv(
        const GenMatrix<T1>& m1, MatrixView<T2> m0) const
    {
        TMVAssert(m1.colsize() == m0.colsize());
        TMVAssert(m1.rowsize() == m0.rowsize());
        TMVAssert(pimpl->LLx.size() == m1.rowsize());
        if (pimpl->LLx.nlo() > 1)
            CH_RDivEq(pimpl->LLx,m0=m1);
        else if (pimpl->LLx.nlo() == 1)
            LDL_RDivEq(pimpl->LLx,m0=m1);
        else
            DiagMatrixViewOf(pimpl->LLx.diag()).RDiv(m1,m0);
    }

    template <class T> 
    T HermBandCHDiv<T>::det() const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->LLx.diag()).logDet(&s);
            if (pimpl->LLx.nlo() > 1) { pimpl->logdet *= RT(2); }
            pimpl->zerodet = s == T(0);
            pimpl->donedet = true;
        }         
        if (pimpl->zerodet) return T(0);
        else return TMV_EXP(pimpl->logdet);  
    }                  

    template <class T> 
    RT HermBandCHDiv<T>::logDet(T* sign) const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->LLx.diag()).logDet(&s);
            if (pimpl->LLx.nlo() > 1) { pimpl->logdet *= RT(2); }
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
    void HermBandCHDiv<T>::doMakeInverse(SymMatrixView<T1> sinv) const
    {
        TMVAssert(sinv.size() == pimpl->LLx.size());
        TMVAssert(sinv.isherm());
        if (pimpl->LLx.nlo() > 1)
            CH_Inverse(pimpl->LLx,sinv);
        else if (pimpl->LLx.nlo() == 1)
            LDL_Inverse(pimpl->LLx,sinv);
        else
            sinv = DiagMatrixViewOf(pimpl->LLx.diag().realPart()).inverse();
    }

    template <class T> template <class T1> 
    void HermBandCHDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    {
        TMVAssert(minv.colsize() == pimpl->LLx.size());
        TMVAssert(minv.rowsize() == pimpl->LLx.size());

        if (isComplex(T1())) minv.diag().imagPart().setZero();
        doMakeInverse(HermMatrixViewOf(minv,Lower));
        if (minv.colsize() > 1)
            minv.upperTri().offDiag() =
                minv.lowerTri().offDiag().adjoint();
    }

    template <class T> 
    void HermBandCHDiv<T>::doMakeInverseATA(MatrixView<T> ata) const
    {
        // ata = (At A)^-1 = A^-1 (A^-1)t
        //     = A^-1 A^-1
        if (isComplex(T())) ata.diag().imagPart().setZero();
        SymMatrixView<T> hermata = HermMatrixViewOf(ata,Lower);
        doMakeInverse(hermata);
        SymSquare<true>(ata);
    }

    template <class T> 
    bool HermBandCHDiv<T>::isSingular() const 
    {
        return pimpl->LLx.diag().minAbs2Element() <=
            TMV_Epsilon<T>() * pimpl->LLx.diag().maxAbs2Element(); 
    }

    template <class T> 
    const BandMatrix<T> HermBandCHDiv<T>::getL() const 
    { 
        BandMatrix<T> L = pimpl->LLx.lowerBand();
        if (pimpl->LLx.nlo() <= 1) L.diag().setAllTo(T(1));
        return L;
    }

    template <class T> 
    const DiagMatrix<T> HermBandCHDiv<T>::getD() const
    { 
        if (pimpl->LLx.nlo() <= 1) return DiagMatrix<T>(pimpl->LLx.diag());
        else return DiagMatrix<T>(pimpl->LLx.size(),T(1));
    }

    template <class T> 
    const GenSymBandMatrix<T>& HermBandCHDiv<T>::getLL() const 
    { return pimpl->LLx; }

    template <class T> 
    bool HermBandCHDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "HermBandCHDiv:\n";
            *fout << "M = "<<mm<<std::endl;
            *fout << "L = "<<getL()<<std::endl;
            *fout << "D = "<<getD()<<std::endl;
        }
        BandMatrix<T> lu = getL()*getD()*getL().adjoint();
        RT nm = Norm(lu-mm);
        nm /= TMV_SQR(Norm(getL())) * Norm(getD());
        if (fout) {
            *fout << "LDLt = "<<lu<<std::endl;
            *fout << "M-LDLt = "<<(mm-lu)<<std::endl;
            *fout << "Norm(M-LDLt)/Norm(LDLt) = "<<nm<<std::endl;
        }
        return nm < mm.doCondition()*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T> 
    ptrdiff_t HermBandCHDiv<T>::colsize() const
    { return pimpl->LLx.size(); }

    template <class T> 
    ptrdiff_t HermBandCHDiv<T>::rowsize() const
    { return pimpl->LLx.size(); }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymBandCHD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


