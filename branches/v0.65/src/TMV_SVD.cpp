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



#include "tmv/TMV_SVD.h"
#include "TMV_SVDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include <ostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    struct SVDiv<T>::SVDiv_Impl
    {
    public :
        SVDiv_Impl(const GenMatrix<T>& m, bool inplace); 

        const bool istrans;
        const bool inplace;
        AlignedArray<T> Aptr1;
        T* Aptr;
        MatrixView<T> U;
        DiagMatrix<RT> S;
        Matrix<T,ColMajor> V;
        RT logdet;
        T signdet;
        mutable int kmax;
    };

#define APTR1 (inplace ? 0 : (A.colsize()*A.rowsize()))
#define APTR (inplace ? A.nonConst().ptr() : Aptr1.get())
#define UX (istrans ? \
            (inplace ? A.nonConst().transpose() : \
             MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor)) : \
            (inplace ? A.nonConst().view() : \
             MatrixViewOf(Aptr,A.colsize(),A.rowsize(), \
                          A.isSquare() ? BaseStorOf(A) : ColMajor)))

    template <class T> 
    SVDiv<T>::SVDiv_Impl::SVDiv_Impl(const GenMatrix<T>& A, bool _inplace) :
        istrans(A.colsize() < A.rowsize()), 
        inplace(_inplace && (A.isrm() || A.iscm())), Aptr1(APTR1), Aptr(APTR),
        U(UX), S(U.rowsize()), V(U.rowsize(),U.rowsize()), 
        logdet(0), signdet(1), kmax(0) {}

#undef UX
#undef APTR

    template <class T> 
    SVDiv<T>::SVDiv(const GenMatrix<T>& A, bool inplace) :
        pimpl(new SVDiv_Impl(A,inplace))
    {
        if (pimpl->istrans) {
            if (inplace) TMVAssert(A.transpose() == pimpl->U); 
            else pimpl->U = A.transpose();
        } else {
            if (inplace) TMVAssert(A == pimpl->U); 
            else pimpl->U = A;
        }

        SV_Decompose<T>(pimpl->U,pimpl->S.view(),pimpl->V.view(),
                        pimpl->logdet, pimpl->signdet,true);

        // Set kmax for actual 0 elements (to within machine precision).
        // Any further cut in the number of singular values to use
        // should be done by the user.
        thresh(TMV_Epsilon<T>());
    }

    template <class T> 
    SVDiv<T>::~SVDiv() {}

    template <class T> template <class T1> 
    void SVDiv<T>::doLDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(m.colsize() == colsize());
        TMVAssert(m.colsize() == rowsize());
        if (pimpl->istrans) 
            SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
                    m.transpose(),m.transpose());
        else 
            SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
    }

    template <class T> template <class T1> 
    void SVDiv<T>::doRDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(m.rowsize() == colsize());
        TMVAssert(m.rowsize() == rowsize());

        if (pimpl->istrans) 
            SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
                    m.transpose(),m.transpose());
        else 
            SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
    }

    template <class T> template <class T1, class T2> 
    void SVDiv<T>::doLDiv(const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        if (pimpl->istrans) 
            SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
                   m.transpose(),x.transpose());
        else 
            SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
    }

    template <class T> template <class T1, class T2> 
    void SVDiv<T>::doRDiv(const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(x.rowsize() == colsize());

        if (pimpl->istrans) 
            SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
                   m.transpose(),x.transpose());
        else 
            SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
    }

    template <class T> 
    T SVDiv<T>::det() const 
    { 
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }

    template <class T> 
    RT SVDiv<T>::logDet(T* sign) const 
    { 
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;
    }

    template <class T> template <class T1> 
    void SVDiv<T>::doMakeInverse(const MatrixView<T1>& minv) const
    { 
        if (pimpl->istrans) {
            // A^-1 = (Vt S^-1 Ut)T = U* S^-1 V*
            Matrix<T,ColMajor> SinvV =
                pimpl->V.conjugate().rowRange(0,pimpl->kmax) /
                pimpl->S.subDiagMatrix(0,pimpl->kmax);
            minv = pimpl->U.conjugate().colRange(0,pimpl->kmax) * SinvV;
        } else {
            // A^-1 = Vt S^-1 Ut
            Matrix<T,ColMajor> SinvUt =
                pimpl->U.adjoint().rowRange(0,pimpl->kmax) /
                pimpl->S.subDiagMatrix(0,pimpl->kmax);
            minv = pimpl->V.adjoint().colRange(0,pimpl->kmax) * SinvUt;
        }
    }

    template <class T> 
    void SVDiv<T>::doMakeInverseATA(const MatrixView<T>& minv) const
    {
        // A = U S V
        // At = Vt S Ut
        // AtA = Vt S^2 V
        // (AtA)^-1 = Vt S^-2 V
        //
        // if istrans:
        // AT = U S V
        // At = U* S V*
        // AAt = VT S^2 V*
        // (AAt)^-1 = VT S^-2 V*
        //
        Matrix<T,ColMajor> SinvV = pimpl->V.rowRange(0,pimpl->kmax) /
            pimpl->S.subDiagMatrix(0,pimpl->kmax);
        if (pimpl->istrans)
            minv = SinvV.transpose() * SinvV.conjugate();
        else
            minv = SinvV.adjoint() * SinvV;
    }

    template <class T> 
    bool SVDiv<T>::isSingular() const
    { return pimpl->kmax < int(pimpl->S.size()); }

    template <class T> 
    RT SVDiv<T>::norm2() const 
    { return pimpl->S.size() > 0 ? pimpl->S(0) : RT(0) ; }

    template <class T> 
    RT SVDiv<T>::condition() const 
    { 
        return pimpl->S.size() > 0 ?
            pimpl->S(0)/pimpl->S(pimpl->S.size()-1) : RT(1); 
    }

    template <class T> 
    void SVDiv<T>::thresh(RT toler, std::ostream* debugout) const
    {
        if (pimpl->S.size() == 0) pimpl->kmax = 0;
        else {
            TMVAssert(toler < RT(1));
            RT thresh = pimpl->S(0)*toler;
            for(pimpl->kmax=pimpl->S.size(); 
                pimpl->kmax>0 && pimpl->S(pimpl->kmax-1)<=thresh;
                --pimpl->kmax);
            if(debugout) {
                (*debugout)<<"S = "<<pimpl->S<<std::endl;
                (*debugout)<<"Smax = "<<pimpl->S(0)<<
                    ", thresh = "<<thresh<<std::endl;
                (*debugout)<<"kmax = "<<pimpl->kmax;
                (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
            }
        }
    }

    template <class T> 
    void SVDiv<T>::top(int neigen, std::ostream* debugout) const
    {
        TMVAssert(neigen > 0);
        if (neigen < int(pimpl->S.size())) pimpl->kmax = neigen;
        else pimpl->kmax = pimpl->S.size();
        if(debugout) {
            (*debugout)<<"S = "<<pimpl->S<<std::endl;
            (*debugout)<<"kmax = "<<pimpl->kmax;
            (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
        }
    }

    template <class T> 
    int SVDiv<T>::getKMax() const
    { return pimpl->kmax; }

    template <class T> 
    ConstMatrixView<T> SVDiv<T>::getU() const
    {
        if (pimpl->istrans) return pimpl->V.transpose(); 
        else return pimpl->U.view(); 
    }
    template <class T> 
    ConstDiagMatrixView<RT> SVDiv<T>::getS() const
    { return pimpl->S.view(); }

    template <class T> 
    ConstMatrixView<T> SVDiv<T>::getV() const
    {
        if (pimpl->istrans) return pimpl->U.transpose(); 
        else return pimpl->V.view(); 
    }

    template <class T> 
    bool SVDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "SVDiv:\n";
            *fout << "M = "<<mm<<std::endl;
            *fout << "U = "<<getU()<<std::endl;
            *fout << "S = "<<getS()<<std::endl;
            *fout << "V = "<<getV()<<std::endl;
        }
        Matrix<T> usv = getU()*getS()*getV();
        RT nm = Norm(usv-mm);
        nm /= Norm(getU())*Norm(getS())*Norm(getV());
        RT cond = getS()(0) / getS()(getKMax()-1);
        if (fout) {
            *fout << "USV = "<<usv<<std::endl;
            *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
            *fout <<"  "<<cond<<" * "<<TMV_Epsilon<T>()<<std::endl;
        }
        return nm < cond*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T> 
    size_t SVDiv<T>::colsize() const
    { return pimpl->istrans ? pimpl->U.rowsize() : pimpl->U.colsize(); }

    template <class T> 
    size_t SVDiv<T>::rowsize() const
    { return pimpl->istrans ? pimpl->U.colsize() : pimpl->U.rowsize(); }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


