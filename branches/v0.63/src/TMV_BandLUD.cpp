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



#include "tmv/TMV_BandLUD.h"
#include "TMV_BandLUDiv.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_TriMatrixArith.h"
#include <ostream>

namespace tmv {

    template <class T> 
    struct BandLUDiv<T>::BandLUDiv_Impl
    {
    public :
        BandLUDiv_Impl(const GenBandMatrix<T>& A, bool _inplace);
        BandLUDiv_Impl(const AssignableToBandMatrix<T>& A);

        const bool istrans;
        const bool inplace;
        auto_array<T> Aptr1;
        T* Aptr;
        BandMatrixView<T> LUx;
        auto_array<int> P;
        mutable TMV_RealType(T) logdet;
        mutable T signdet;
        mutable bool donedet;
    };

#define NEWLO TMV_MIN(A.nlo(),A.nhi())
#define NEWHI TMV_MIN(A.nlo()+A.nhi(),int(A.colsize())-1)
#define APTR1 (inplace ? 0 : \
               new T[BandStorageLength(ColMajor,A.colsize(),A.colsize(), \
                                       NEWLO,NEWHI)])
#define APTR (inplace ? A.nonConst().ptr() : Aptr1.get())

#define LUX (istrans ? \
             (inplace ? \
              BandMatrixView<T>(A.nonConst().ptr(),A.colsize(),A.colsize(),\
                                A.nhi(),NEWHI,A.stepj(),A.stepi(),A.diagstep(),\
                                TMV_TransOf(A.stor()),A.ct() \
                                TMV_FIRSTLAST1(A.nonConst().first,\
                                               A.nonConst().last) ) : \
              BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nhi(), \
                               NEWHI, (A.nlo() == 1 && A.nhi() == 1) ? \
                               DiagMajor : ColMajor)) : \
             (inplace ? \
              BandMatrixView<T>(A.nonConst().ptr(),A.colsize(),\
                                A.colsize(),A.nlo(),NEWHI,\
                                A.stepi(),A.stepj(),A.diagstep(),\
                                A.stor(),A.ct() \
                                TMV_FIRSTLAST1(A.nonConst().first,\
                                               A.nonConst().last) ) : \
              BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nlo(), \
                               NEWHI, (A.nlo() == 1 && A.nhi() == 1) ?\
                               DiagMajor : ColMajor)))

    template <class T> 
    BandLUDiv<T>::BandLUDiv_Impl::BandLUDiv_Impl(
        const GenBandMatrix<T>& A, bool _inplace
    ) :
        istrans(A.nhi()<A.nlo() || (A.nhi()==A.nlo() && A.isrm())),
        inplace(NEWLO == 0 || 
                (_inplace && 
                 ((A.isrm() && istrans) || (A.iscm() && !istrans) || 
                  (A.isdm() && A.nlo()==1 && A.nhi()==1)))),
        Aptr1(APTR1), Aptr(APTR), LUx(LUX),
        P(new int[A.colsize()]), logdet(0), signdet(1), donedet(false) {}

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
                TMVAssert(A.stepj() >= TMV_MIN(int(A.colsize()),2*A.nlo()+A.nhi()));
                TMVAssert(pimpl->LUx.diagRange(-A.nlo(),A.nhi()+1) == A);
            } else if (A.isrm()) {
                TMVAssert(pimpl->istrans);
                TMVAssert(A.stepi() >= TMV_MIN(int(A.colsize()),2*A.nhi()+A.nlo()));
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

        if (pimpl->LUx.nlo() > 0) {
            int Anhi = pimpl->istrans ? A.nlo() : A.nhi();
            if (Anhi < pimpl->LUx.nhi())
                pimpl->LUx.diagRange(Anhi+1,pimpl->LUx.nhi()+1).zero();
            LU_Decompose(pimpl->LUx,pimpl->P.get(),pimpl->signdet,Anhi);
        } else {
            pimpl->P.get()[0] = 1; // A tag to indicate that P is not set yet.
        }
    }

#define NEWLO TMV_MIN(A.nlo(),A.nhi())
#define NEWHI TMV_MIN(A.nlo()+A.nhi(),int(A.colsize())-1)
#define APTR1 \
    new T[BandStorageLength(ColMajor,A.colsize(),A.colsize(),NEWLO,NEWHI)]
#define APTR Aptr1.get()

#define LUX \
    BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),NEWLO,NEWHI, \
                     (A.nlo()==1 && A.nhi()==1) ? DiagMajor : ColMajor)

    template <class T> 
    BandLUDiv<T>::BandLUDiv_Impl::BandLUDiv_Impl(
        const AssignableToBandMatrix<T>& A) :
        istrans(A.nhi()<A.nlo()), inplace(false),
        Aptr1(APTR1), Aptr(APTR), LUx(LUX),
        P(new int[A.colsize()]), logdet(0), signdet(1), donedet(false) {}

#undef LUX
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

        if (pimpl->LUx.nlo() > 0) {
            int Anhi = pimpl->istrans ? A.nlo() : A.nhi();
            if (Anhi < pimpl->LUx.nhi())
                pimpl->LUx.diagRange(Anhi+1,pimpl->LUx.nhi()+1).zero();
            LU_Decompose(pimpl->LUx,pimpl->P.get(),pimpl->signdet,Anhi);
        } else {
            pimpl->P.get()[0] = 1; // A tag to indicate that P is not set yet.
        }
    }

    template <class T> 
    BandLUDiv<T>::~BandLUDiv() {}

    template <class T> template <class T1> 
    void BandLUDiv<T>::doLDivEq(const MatrixView<T1>& m) const
    {
        if (pimpl->istrans) LU_RDivEq(pimpl->LUx,pimpl->P.get(),m.transpose());
        else LU_LDivEq(pimpl->LUx,pimpl->P.get(),m);
    }

    template <class T> template <class T1> 
    void BandLUDiv<T>::doRDivEq(const MatrixView<T1>& m) const
    {
        if (pimpl->istrans) LU_LDivEq(pimpl->LUx,pimpl->P.get(),m.transpose());
        else LU_RDivEq(pimpl->LUx,pimpl->P.get(),m);
    }

    template <class T> template <class T1, class T2> 
    void BandLUDiv<T>::doLDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        if (pimpl->istrans) 
            LU_RDivEq(pimpl->LUx,pimpl->P.get(),(x=m).transpose());
        else LU_LDivEq(pimpl->LUx,pimpl->P.get(),x=m);
    }

    template <class T> template <class T1, class T2> 
    void BandLUDiv<T>::doRDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        if (pimpl->istrans) 
            LU_LDivEq(pimpl->LUx,pimpl->P.get(),(x=m).transpose());
        else LU_RDivEq(pimpl->LUx,pimpl->P.get(),x=m);
    }

    template <class T> T BandLUDiv<T>::det() const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->LUx.diag()).logDet(&s);
            pimpl->signdet *= s;
            pimpl->donedet = true;
        }         
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }                  

    template <class T> 
    TMV_RealType(T) BandLUDiv<T>::logDet(T* sign) const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->LUx.diag()).logDet(&s);
            pimpl->signdet *= s;
            pimpl->donedet = true;
        }
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;  
    }                  

    template <class T> template <class T1> 
    void BandLUDiv<T>::doMakeInverse(const MatrixView<T1>& minv) const
    {
        if (pimpl->istrans)
            LU_Inverse(pimpl->LUx,pimpl->P.get(),minv.transpose());
        else
            LU_Inverse(pimpl->LUx,pimpl->P.get(),minv);
    }

    template <class T> 
    void BandLUDiv<T>::doMakeInverseATA(const MatrixView<T>& ata) const
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
            LU_PackedPL_RDivEq(pimpl->LUx,pimpl->P.get(),ata.transpose());
            // ata %= L.conjugate()
            // ata* = ata* L^-1
            LU_PackedPL_RDivEq(pimpl->LUx,pimpl->P.get(),ata.conjugate());
        } else {
            LowerTriMatrixView<T> linv = ata.lowerTri(UnitDiag);
            // linv = L.inverse()
            ata.setToIdentity();
            LU_PackedPL_LDivEq(pimpl->LUx,pimpl->P.get(),ata);
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
    { return det() == T(0); }

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
        LU_PackedPL_Unpack(pimpl->LUx,pimpl->P.get(),L.view());
        return L;
    }

    template <class T> 
    const GenBandMatrix<T>& BandLUDiv<T>::getLU() const 
    { return pimpl->LUx; }

    template <class T> 
    const int* BandLUDiv<T>::getP() const
    {
        if (pimpl->LUx.nlo() == 0 && pimpl->LUx.colsize() > 0 && 
            pimpl->P.get()[0] == 1) {
            // Then P hasn't been set up yet.
            const int N = pimpl->LUx.colsize();
            for(int i=0;i<N;++i) pimpl->P.get()[i] = i;
        }
        return pimpl->P.get();
    }

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
        Matrix<T> lu = getL() * Matrix<T>(getU());
        lu.reversePermuteRows(getP());
        TMV_RealType(T) nm = Norm(
            lu-(pimpl->istrans ? mm.transpose() : mm.view()));
        nm /= Norm(getL())*Norm(getU());
        if (fout) {
            *fout << "PLU = "<<lu<<std::endl;
            *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<std::endl;
        }
        return nm < mm.doCondition()*mm.colsize()*TMV_Epsilon<T>();
    }

    template <class T> 
    size_t BandLUDiv<T>::colsize() const
    { return pimpl->LUx.colsize(); }

    template <class T> 
    size_t BandLUDiv<T>::rowsize() const
    { return pimpl->LUx.rowsize(); }

#define InstFile "TMV_BandLUD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


