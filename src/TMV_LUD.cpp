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


#include "tmv/TMV_LUD.h"
#include "TMV_LUDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_PermutationArith.h"
#include <ostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    struct LUDiv<T>::LUDiv_Impl
    {
    public :
        LUDiv_Impl(const GenMatrix<T>& m, bool inplace);

        const bool istrans;
        const bool inplace;
        AlignedArray<T> Aptr1;
        T* Aptr;
        MatrixView<T> LUx;
        Permutation P;
        mutable RT logdet;
        mutable T signdet;
        mutable bool donedet;

        const GenMatrix<T>& A0;
    };

#define APTR1 (inplace ? 0 : (A.colsize()*A.rowsize()))
#define APTR (inplace ? A.nonConst().ptr() : Aptr1.get())
#define LUX (istrans ? \
             (inplace ? A.nonConst().transpose() : \
              MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) ): \
             (inplace ? A.nonConst().view() : \
              MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)))

    template <class T> 
    LUDiv<T>::LUDiv_Impl::LUDiv_Impl(
        const GenMatrix<T>& A, bool _inplace) :
        istrans(A.isrm()), inplace(_inplace && (A.iscm() || A.isrm())),
        Aptr1(APTR1), Aptr(APTR), LUx(LUX), P(A.colsize()),
        logdet(0), signdet(1), donedet(false), A0(A) {}

#undef LUX
#undef APTR

    template <class T> 
    LUDiv<T>::LUDiv(const GenMatrix<T>& A, bool inplace) :
        pimpl(new LUDiv_Impl(A,inplace)) 
    {
        TMVAssert(A.isSquare());
        if (pimpl->istrans) {
            if (pimpl->inplace) TMVAssert(A.transpose() == pimpl->LUx);
            else pimpl->LUx = A.transpose();
        } else {
            if (inplace) TMVAssert(A == pimpl->LUx);
            else pimpl->LUx = A;
        }
        LU_Decompose(pimpl->LUx,pimpl->P);
    }

    template <class T> 
    LUDiv<T>::~LUDiv() {}

    template <class T> template <class T1> 
    void LUDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(pimpl->LUx.colsize() == m.colsize());
        if (pimpl->istrans) 
            LU_RDivEq(pimpl->LUx,pimpl->P.getValues(),m.transpose());
        else 
            LU_LDivEq(pimpl->LUx,pimpl->P.getValues(),m);
    }

    template <class T> template <class T1> 
    void LUDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(pimpl->LUx.colsize() == m.rowsize());
        if (pimpl->istrans) 
            LU_LDivEq(pimpl->LUx,pimpl->P.getValues(),m.transpose());
        else 
            LU_RDivEq(pimpl->LUx,pimpl->P.getValues(),m);
    }

    template <class T> template <class T1, class T2> 
    void LUDiv<T>::doLDiv(
        const GenMatrix<T1>& m1, MatrixView<T2> m0) const
    {
        TMVAssert(m1.colsize() == m0.colsize());
        TMVAssert(m1.rowsize() == m0.rowsize());
        TMVAssert(pimpl->LUx.colsize() == m1.colsize());
        doLDivEq(m0=m1);
    }

    template <class T> template <class T1, class T2> 
    void LUDiv<T>::doRDiv(
        const GenMatrix<T1>& m1, MatrixView<T2> m0) const
    {
        TMVAssert(m1.colsize() == m0.colsize());
        TMVAssert(m1.rowsize() == m0.rowsize());
        TMVAssert(pimpl->LUx.colsize() == m1.rowsize());
        doRDivEq(m0=m1);
    }

    template <class T> 
    T LUDiv<T>::det() const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->LUx.diag()).logDet(&s);
            pimpl->signdet = RT(pimpl->P.det()) * s;
            pimpl->donedet = true;
        }         
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }                  

    template <class T> 
    RT LUDiv<T>::logDet(T* sign) const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->LUx.diag()).logDet(&s);
            pimpl->signdet = RT(pimpl->P.det()) * s;
            pimpl->donedet = true;
        }
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;  
    }                  

    template <class T> template <class T1> 
    void LUDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    {
        TMVAssert(minv.colsize() == pimpl->LUx.colsize());
        TMVAssert(minv.rowsize() == pimpl->LUx.colsize());
        // m = P L U
        // m^-1 = U^-1 L^-1 Pt
        if (pimpl->istrans) 
            LU_Inverse(pimpl->LUx,pimpl->P.getValues(),minv.transpose());
        else 
            LU_Inverse(pimpl->LUx,pimpl->P.getValues(),minv);
    }

    template <class T> 
    void LUDiv<T>::doMakeInverseATA(MatrixView<T> ata) const
    {
        TMVAssert(ata.colsize() == pimpl->LUx.colsize());
        TMVAssert(ata.rowsize() == pimpl->LUx.colsize());
        // (At A)^-1 = A^-1 (A^-1)t
        // = (U^-1 L^-1 Pt) (P L^-1t U^-1t)
        // = U^-1 L^-1 L^-1t U^-1t
        //
        // if PLU is really AT, then
        // A^-1 = P L^-1T U^-1T
        // (At A)^-1 = P L^-1T U^-1T U^-1* L^-1* Pt

        LowerTriMatrixView<T> L = pimpl->LUx.lowerTri(UnitDiag);
        UpperTriMatrixView<T> U = pimpl->LUx.upperTri();

        if (pimpl->istrans) {
            UpperTriMatrixView<T> uinv = ata.upperTri();
            uinv = U.inverse();
            ata = uinv.transpose() * uinv.conjugate();
            ata /= L.transpose();
            ata %= L.conjugate();
            ata.reversePermuteCols(pimpl->P.getValues());
            ata.reversePermuteRows(pimpl->P.getValues());
        } else {
            LowerTriMatrixView<T> linv = ata.lowerTri(UnitDiag);
            linv = L.inverse();
            ata = linv * linv.adjoint();
            ata /= U;
            ata %= U.adjoint();
        }
    }

    template <class T> 
    bool LUDiv<T>::isSingular() const 
    {
        return pimpl->LUx.diag().minAbs2Element() <=
            TMV_Epsilon<T>() * pimpl->LUx.diag().maxAbs2Element(); 
    }

    template <class T> 
    bool LUDiv<T>::isTrans() const 
    { return pimpl->istrans; }

    template <class T> 
    ConstLowerTriMatrixView<T> LUDiv<T>::getL() const 
    { return pimpl->LUx.lowerTri(UnitDiag); }

    template <class T> 
    ConstUpperTriMatrixView<T> LUDiv<T>::getU() const 
    { return pimpl->LUx.upperTri(NonUnitDiag); }

    template <class T> 
    const GenMatrix<T>& LUDiv<T>::getLU() const 
    { return pimpl->LUx; }

    template <class T> 
    const Permutation& LUDiv<T>::getP() const 
    { return pimpl->P; }

    template <class T> 
    bool LUDiv<T>::checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "LUDiv:\n";
            *fout << "M = "<<
                (pimpl->istrans?mm.transpose():mm.view())<<std::endl;
            *fout << "L = "<<getL()<<std::endl;
            *fout << "U = "<<getU()<<std::endl;
            *fout << "P = "<<getP()<<std::endl;
            *fout << "  or by interchanges: ";
            for(ptrdiff_t i=0;i<getP().size();i++)
                *fout<<(getP().getValues())[i]<<" ";
        }
        Matrix<T> lu = getP()*getL()*getU();
        RT nm = Norm(lu-(pimpl->istrans ? mm.transpose() : mm.view()));
        nm /= Norm(getL())*Norm(getU());
        if (fout) {
            *fout << "PLU = "<<lu<<std::endl;
            *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<std::endl;
        }
        return nm < mm.doCondition()*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T> 
    ptrdiff_t LUDiv<T>::colsize() const
    { return pimpl->LUx.colsize(); }

    template <class T> 
    ptrdiff_t LUDiv<T>::rowsize() const
    { return pimpl->LUx.rowsize(); }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_LUD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


