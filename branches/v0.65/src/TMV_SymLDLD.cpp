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
        T* Aptr;
        SymMatrixView<T> LLx;
        Vector<T> xD;
        Permutation P;
        mutable RT logdet;
        mutable T signdet;
    };

#define APTR1 (inplace ? 0 : (A.size()*A.size()))
#define APTR (inplace ? A.nonConst().ptr() : Aptr1.get())
#define LLX \
    (inplace ? \
     (A.uplo()==Upper ? \
      (A.isherm() ? A.nonConst().adjoint() : A.nonConst().transpose()) : \
      A.nonConst()) : \
     (A.isherm() ? \
      HermMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.lowerTri())) : \
      SymMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.lowerTri()))))

    template <class T>
    SymLDLDiv<T>::SymLDLDiv_Impl::SymLDLDiv_Impl(
        const GenSymMatrix<T>& A, bool _inplace
    ) :
        inplace(_inplace && (A.isrm() || A.iscm())), 
        Aptr1(APTR1), Aptr(APTR), LLx(LLX), xD(A.size()-1),
        P(A.colsize()), logdet(0), signdet(1) {}

#undef APTR1
#undef APTR
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
    void SymLDLDiv<T>::doLDivEq(const MatrixView<T1>& m) const
    { LDL_LDivEq(pimpl->LLx,pimpl->xD,pimpl->P.getValues(),m); }

    template <class T> template <class T1> 
    void SymLDLDiv<T>::doRDivEq(const MatrixView<T1>& m) const
    { LDL_RDivEq(pimpl->LLx,pimpl->xD,pimpl->P.getValues(),m); }

    template <class T> template <class T1, class T2> 
    void SymLDLDiv<T>::doLDiv(
        const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
    { LDL_LDivEq(pimpl->LLx,pimpl->xD,pimpl->P.getValues(),m0=m1); }

    template <class T> template <class T1, class T2> 
    void SymLDLDiv<T>::doRDiv(
        const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
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
    void SymLDLDiv<T>::doMakeInverse(const SymMatrixView<T1>& sinv) const
    { 
        TMVAssert(isReal(T()) || issym() == sinv.issym());
        TMVAssert(isReal(T()) || isherm() == sinv.isherm());
        LDL_Inverse(pimpl->LLx,pimpl->xD,pimpl->P.getValues(),sinv); 
    }

    template <class T> template <class T1> 
    void SymLDLDiv<T>::doMakeInverse(const MatrixView<T1>& minv) const
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
    void SymLDLDiv<T>::doMakeInverseATA(const MatrixView<T>& ata) const
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
            for(int i=0;i<int(getP().size());i++)
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
    size_t SymLDLDiv<T>::colsize() const
    { return pimpl->LLx.size(); }

    template <class T>
    size_t SymLDLDiv<T>::rowsize() const
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


