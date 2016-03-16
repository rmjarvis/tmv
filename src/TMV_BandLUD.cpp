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


//#include <iostream>

#include "tmv/TMV_BandLUD.h"
#include "TMV_BandLUDiv.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Permutation.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_PermutationArith.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
#include <ostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    struct BandLUDiv<T>::BandLUDiv_Impl
    {
    public :
        BandLUDiv_Impl(const GenBandMatrix<T>& A, bool _inplace);
        BandLUDiv_Impl(const AssignableToBandMatrix<T>& A);

        const bool istrans;
        const bool inplace;
        AlignedArray<T> Aptr1;
        T* Aptr;
        BandMatrixView<T> LUx;
        Permutation P;
        mutable RT logdet;
        mutable T signdet;
        mutable bool donedet;
    };

#define NEWLO TMV_MIN(A.nlo(),A.nhi())
#define NEWHI TMV_MIN(A.nlo()+A.nhi(),A.colsize()-1)
#define APTR1 (inplace ? 0 : \
               BandStorageLength(ColMajor,A.colsize(),A.colsize(),NEWLO,NEWHI))
#define TRID (A.nlo() == 1 && A.nhi() == 1)
#define APTR (inplace ? A.nonConst().ptr() : Aptr1.get())

#define LUX (istrans ? \
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

    template <class T> 
    BandLUDiv<T>::BandLUDiv_Impl::BandLUDiv_Impl(
        const GenBandMatrix<T>& A, bool _inplace) :
        istrans(A.nhi()<A.nlo() || (A.nhi()==A.nlo() && A.isrm())),
        inplace(NEWLO == 0 || 
                (_inplace && 
                 ((A.isrm() && istrans) || (A.iscm() && !istrans) || 
                  (A.isdm() && TRID)))),
        Aptr1(APTR1), Aptr(APTR), LUx(LUX),
        P(A.colsize()), logdet(0), signdet(1), donedet(false) 
    {
        //std::cout<<"BandLUDiv_Impl constructor\n";
        //std::cout<<"A = "<<TMV_Text(A)<<" = "<<A<<std::endl;
        //LUx.setZero();
        //std::cout<<"LUx = "<<TMV_Text(LUx)<<" = "<<LUx<<std::endl;
    }

#undef LUX
#undef APTR
#undef APTR1
#undef NEWLO
#undef NEWHI

    template <class T> 
    BandLUDiv<T>::BandLUDiv(const GenBandMatrix<T>& A, bool inplace) :
        pimpl(new BandLUDiv_Impl(A,inplace))
    {
        TMVAssert(A.isSquare());
        if (inplace) {
            // For inplace decomposition, make sure the original band matrix
            // has room for the extra upper diagonals...
            // if iscm stepj >= (2*A.nlo()+A.nhi())
            // if isdm extra diags appear at end, so can't really check
            // if isrm stepi >= (2*A.nhi()+A.nlo())
            if (A.iscm()) {
                TMVAssert(!pimpl->istrans);
                TMVAssert(A.stepj() >= TMV_MIN(A.colsize(),2*A.nlo()+A.nhi()));
                TMVAssert(pimpl->LUx.diagRange(-A.nlo(),A.nhi()+1) == A);
            } else if (A.isrm()) {
                TMVAssert(pimpl->istrans);
                TMVAssert(A.stepi() >= TMV_MIN(A.colsize(),2*A.nhi()+A.nlo()));
                TMVAssert(
                    pimpl->LUx.diagRange(-A.nhi(),A.nlo()+1).transpose() == A);
            } else {
                TMVAssert(A.isdm());
                if (pimpl->istrans)
                    TMVAssert(
                        pimpl->LUx.diagRange(-A.nhi(),A.nlo()+1).transpose() == A);
                else
                    TMVAssert(pimpl->LUx.diagRange(-A.nlo(),A.nhi()+1) == A);
            }
        } else {
            if (pimpl->istrans) 
                BandMatrixViewOf(pimpl->LUx,A.nhi(),A.nlo()) =
                    A.transpose();
            else BandMatrixViewOf(pimpl->LUx,A.nlo(),A.nhi()) = A;
        }
        //std::cout<<"LUx => "<<pimpl->LUx<<std::endl;

        if (pimpl->LUx.nlo() > 0) {
            ptrdiff_t Anhi = pimpl->istrans ? A.nlo() : A.nhi();
            if (Anhi < pimpl->LUx.nhi())
                pimpl->LUx.diagRange(Anhi+1,pimpl->LUx.nhi()+1).setZero();
            //std::cout<<"LUx => "<<pimpl->LUx<<std::endl;
            LU_Decompose(pimpl->LUx,pimpl->P,Anhi);
            //std::cout<<"LUx => "<<pimpl->LUx<<std::endl;
        }
    }

#define NEWLO TMV_MIN(A.nlo(),A.nhi())
#define NEWHI TMV_MIN(A.nlo()+A.nhi(),A.colsize()-1)
#define APTR1 BandStorageLength(ColMajor,A.colsize(),A.colsize(),NEWLO,NEWHI)
#define APTR Aptr1.get()

#define LUX \
    BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),NEWLO,NEWHI, \
                     TRID ? DiagMajor : ColMajor)

    template <class T> 
    BandLUDiv<T>::BandLUDiv_Impl::BandLUDiv_Impl(
        const AssignableToBandMatrix<T>& A) :
        istrans(A.nhi()<A.nlo()), inplace(false),
        Aptr1(APTR1), Aptr(APTR), LUx(LUX),
        P(A.colsize()), logdet(0), signdet(1), donedet(false) 
    {
        //std::cout<<"BandLUDivImpl Assignable constructor\n";
        //std::cout<<"A = "<<TMV_Text(A)<<" = "<<BandMatrix<T>(A)<<std::endl;
        //std::cout<<"LUx = "<<TMV_Text(LUx)<<" = "<<LUx<<std::endl;
        //std::cout<<"Aptr1 = "<<Aptr1.get()<<std::endl;
        //std::cout<<"Aptr = "<<Aptr<<std::endl;
        //std::cout<<"len = "<<APTR1<<std::endl;
        //std::cout<<"LUx.stepi = "<<LUx.stepi()<<std::endl;
        //std::cout<<"LUx.stepj = "<<LUx.stepj()<<std::endl;
        //std::cout<<"LUx.first = "<<LUx._first<<std::endl;
        //std::cout<<"LUx.last = "<<LUx._last<<std::endl;
        //std::cout<<"&LUx(0,0) = "<<LUx.cptr()<<std::endl;
        //std::cout<<"&LUx("<<NEWLO<<",0) = "<<LUx.cptr()+NEWLO*LUx.stepi()<<std::endl;
        //std::cout<<"&LUx(0,"<<NEWHI<<") = "<<LUx.cptr()+NEWHI*LUx.stepj()<<std::endl;
    }

#undef LUX
#undef TRID
#undef APTR
#undef APTR1
#undef NEWLO
#undef NEWHI

    template <class T> 
    BandLUDiv<T>::BandLUDiv(
        const AssignableToBandMatrix<T>& A) : pimpl(new BandLUDiv_Impl(A))
    {
        TMVAssert(A.isSquare());
        if (pimpl->istrans) 
            BandMatrixViewOf(pimpl->LUx,A.nhi(),A.nlo()).transpose() = A;
        else 
            BandMatrixViewOf(pimpl->LUx,A.nlo(),A.nhi()) = A;
        //std::cout<<"LUx => "<<pimpl->LUx<<std::endl;

        if (pimpl->LUx.nlo() > 0) {
            ptrdiff_t Anhi = pimpl->istrans ? A.nlo() : A.nhi();
            if (Anhi < pimpl->LUx.nhi())
                pimpl->LUx.diagRange(Anhi+1,pimpl->LUx.nhi()+1).setZero();
            //std::cout<<"LUx => "<<pimpl->LUx<<std::endl;
            LU_Decompose(pimpl->LUx,pimpl->P,Anhi);
            //std::cout<<"LUx => "<<pimpl->LUx<<std::endl;
        }
    }

    template <class T> 
    BandLUDiv<T>::~BandLUDiv() {}

    template <class T> template <class T1> 
    void BandLUDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        if (pimpl->istrans) 
            LU_RDivEq(pimpl->LUx,pimpl->P.getValues(),m.transpose());
        else LU_LDivEq(pimpl->LUx,pimpl->P.getValues(),m);
    }

    template <class T> template <class T1> 
    void BandLUDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        if (pimpl->istrans) 
            LU_LDivEq(pimpl->LUx,pimpl->P.getValues(),m.transpose());
        else LU_RDivEq(pimpl->LUx,pimpl->P.getValues(),m);
    }

    template <class T> template <class T1, class T2> 
    void BandLUDiv<T>::doLDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        if (pimpl->istrans) 
            LU_RDivEq(pimpl->LUx,pimpl->P.getValues(),(x=m).transpose());
        else LU_LDivEq(pimpl->LUx,pimpl->P.getValues(),x=m);
    }

    template <class T> template <class T1, class T2> 
    void BandLUDiv<T>::doRDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        if (pimpl->istrans) 
            LU_LDivEq(pimpl->LUx,pimpl->P.getValues(),(x=m).transpose());
        else LU_RDivEq(pimpl->LUx,pimpl->P.getValues(),x=m);
    }

    template <class T> T BandLUDiv<T>::det() const
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
    RT BandLUDiv<T>::logDet(T* sign) const
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
    void BandLUDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    {
        if (pimpl->istrans)
            LU_Inverse(pimpl->LUx,pimpl->P.getValues(),minv.transpose());
        else
            LU_Inverse(pimpl->LUx,pimpl->P.getValues(),minv);
    }

    template <class T> 
    void BandLUDiv<T>::doMakeInverseATA(MatrixView<T> ata) const
    {
        // See corresponding routine in TMV_LUD.cpp
        if (pimpl->istrans) {
            UpperTriMatrixView<T> uinv = ata.upperTri();
            uinv = getU();
            TriInverse(uinv,pimpl->LUx.nhi());
            ata = uinv.transpose() * uinv.conjugate();
            // ata /= L.transpose()
            // ata = LT^-1 ata
            // ataT = ataT L^-1
            LU_PackedPL_RDivEq(pimpl->LUx,pimpl->P.getValues(),ata.transpose());
            // ata %= L.conjugate()
            // ata* = ata* L^-1
            LU_PackedPL_RDivEq(pimpl->LUx,pimpl->P.getValues(),ata.conjugate());
        } else {
            LowerTriMatrixView<T> linv = ata.lowerTri(UnitDiag);
            // linv = L.inverse()
            ata.setToIdentity();
            LU_PackedPL_LDivEq(pimpl->LUx,pimpl->P.getValues(),ata);
            ata = linv * linv.adjoint();
            ConstBandMatrixView<T> U = getU();
            // ata /= U;
            TriLDivEq(U,ata,NonUnitDiag);
            // ata %= U.adjoint();
            // ata = ata Ut^-1
            // atat = U^-1 atat
            TriLDivEq(U,ata.adjoint(),NonUnitDiag);
        }
    }

    template <class T> 
    bool BandLUDiv<T>::isSingular() const 
    { 
        return pimpl->LUx.diag().minAbs2Element() <=
            TMV_Epsilon<T>() * pimpl->LUx.diag().maxAbs2Element(); 
    }

    template <class T> 
    bool BandLUDiv<T>::isTrans() const 
    { return pimpl->istrans; }

    template <class T> 
    ConstBandMatrixView<T> BandLUDiv<T>::getU() const
    { return BandMatrixViewOf(pimpl->LUx,0,pimpl->LUx.nhi()); }

    template <class T> 
    LowerTriMatrix<T,UnitDiag> BandLUDiv<T>::getL() const
    {
        LowerTriMatrix<T,UnitDiag> L(pimpl->LUx.colsize());
        LU_PackedPL_Unpack(pimpl->LUx,pimpl->P.getValues(),L.view());
        return L;
    }

    template <class T> 
    const GenBandMatrix<T>& BandLUDiv<T>::getLU() const 
    { return pimpl->LUx; }

    template <class T> 
    const Permutation& BandLUDiv<T>::getP() const
    { return pimpl->P; }

    template <class T> 
    bool BandLUDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "BandLUDiv:\n";
            *fout << "M = "<<
                (pimpl->istrans?mm.transpose():mm.view())<<std::endl;
            *fout << "L = "<<getL()<<std::endl;
            *fout << "U = "<<getU()<<std::endl;
        }
        Matrix<T> lu = getP() * getL() * getU();
        RT nm = Norm(lu-(pimpl->istrans ? mm.transpose() : mm.view()));
        nm /= Norm(getL())*Norm(getU());
        if (fout) {
            *fout << "PLU = "<<lu<<std::endl;
            *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<std::endl;
        }
        return nm < mm.doCondition()*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T> 
    ptrdiff_t BandLUDiv<T>::colsize() const
    { return pimpl->LUx.colsize(); }

    template <class T> 
    ptrdiff_t BandLUDiv<T>::rowsize() const
    { return pimpl->LUx.rowsize(); }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandLUD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


