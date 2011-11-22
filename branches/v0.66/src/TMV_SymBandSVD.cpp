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



#include "tmv/TMV_SymBandSVD.h"
#include "TMV_SymBandSVDiv.h"
#include "tmv/TMV_SymBandMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include "TMV_SymSVDiv.h"
#include <ostream>
#include <iostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T>
    struct HermBandSVDiv<T>::HermBandSVDiv_Impl
    {
    public :
        HermBandSVDiv_Impl(const GenSymBandMatrix<T>& m) :
            U(m.size(),m.size()), S(m.size()), logdet(0), signdet(1), kmax(0) {}

        Matrix<T,ColMajor> U;
        DiagMatrix<RT> S;
        RT logdet;
        T signdet;
        mutable int kmax;
    }; // HermBandSVDiv

    template <class T>
    HermBandSVDiv<T>::HermBandSVDiv(const GenSymBandMatrix<T>& A) :
        pimpl(new HermBandSVDiv_Impl(A))
    {
        TMVAssert(A.isherm());

        SV_Decompose<T>(A, pimpl->U.view(), pimpl->S.view(), 0,
                       pimpl->logdet, pimpl->signdet);
        thresh(TMV_Epsilon<T>());
    }

    template <class T>
    HermBandSVDiv<T>::~HermBandSVDiv() {}

    template <class T> template <class T1>
    void HermBandSVDiv<T>::doLDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(m.colsize() == colsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,m);
    }

    template <class T> template <class T1>
    void HermBandSVDiv<T>::doRDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(m.rowsize() == rowsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,m);
    }

    template <class T> template <class T1, class T2> 
    void HermBandSVDiv<T>::doLDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,x);
    }

    template <class T> template <class T1, class T2> 
    void HermBandSVDiv<T>::doRDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(x.rowsize() == colsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,x);
    }

    template <class T>
    T HermBandSVDiv<T>::det() const 
    { 
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }

    template <class T>
    RT HermBandSVDiv<T>::logDet(T* sign) const 
    { 
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;
    }

    template <class T> template <class T1>
    void HermBandSVDiv<T>::doMakeInverse(const SymMatrixView<T1>& sinv) const
    {
        TMVAssert(sinv.size() == pimpl->S.size());
        TMVAssert(sinv.isherm());
        CallHermSV_Inverse(T(),pimpl->U,pimpl->S,pimpl->kmax,sinv);
    }

    template <class T> template <class T1>
    void HermBandSVDiv<T>::doMakeInverse(const MatrixView<T1>& minv) const
    { 
        // A^-1 = U S^-1 Ut
        CallHermSV_Inverse(
            T(),pimpl->U,pimpl->S,pimpl->kmax,HermMatrixViewOf(minv,Upper));
        if (pimpl->S.size() > 1)
            minv.lowerTri().offDiag() = minv.upperTri().offDiag().adjoint();
    }

    template <class T> 
    void HermBandSVDiv<T>::doMakeInverseATA(const MatrixView<T>& minv) const
    {
        // A = U S Ut
        // At = U S Ut
        // AtA = U S^2 Ut
        // (AtA)^-1 = U S^-2 Ut
        //
        Matrix<T,RowMajor> SinvUt = 
            pimpl->U.adjoint().rowRange(0,pimpl->kmax) /
            pimpl->S.subDiagMatrix(0,pimpl->kmax);
        minv = SinvUt.adjoint() * SinvUt;
    }

    template <class T> 
    bool HermBandSVDiv<T>::isSingular() const 
    { return pimpl->kmax < int(pimpl->S.size()); }

    template <class T>
    RT HermBandSVDiv<T>::norm2() const 
    { return pimpl->S.size() > 0 ? TMV_ABS(pimpl->S(0)) : RT(0); }

    template <class T>
    RT HermBandSVDiv<T>::condition() const 
    { 
        return pimpl->S.size() > 0 ? 
            TMV_ABS(pimpl->S(0)/pimpl->S(pimpl->S.size()-1)) : RT(1);
    }


    template <class T>
    void HermBandSVDiv<T>::thresh(RT toler, std::ostream* debugout) const
    {
        if (pimpl->S.size() == 0) pimpl->kmax = 0;
        else {
            TMVAssert(toler < RT(1));
            RT thresh = TMV_ABS(pimpl->S(0))*toler;
            for(pimpl->kmax=pimpl->S.size(); 
                pimpl->kmax>0 && TMV_ABS(pimpl->S(pimpl->kmax-1))<=thresh;
                --pimpl->kmax);
            if(debugout) {
                (*debugout)<<"S = "<<pimpl->S<<std::endl;
                (*debugout)<<"Smax = "<<TMV_ABS(pimpl->S(0));
                (*debugout)<<", thresh = "<<thresh<<std::endl;
                (*debugout)<<"kmax = "<<pimpl->kmax;
                (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
            }
        }
    }

    template <class T>
    void HermBandSVDiv<T>::top(int neigen, std::ostream* debugout) const
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
    int HermBandSVDiv<T>::getKMax() const 
    { return pimpl->kmax; }

    template <class T>
    ConstMatrixView<T> HermBandSVDiv<T>::getU() const
    { return pimpl->U.view(); }

    template <class T>
    DiagMatrix<RT> HermBandSVDiv<T>::getS() const
    {
        DiagMatrix<RT> temp = pimpl->S;
        const size_t N = pimpl->S.size();
        for(size_t i=0;i<N;i++) 
            if (temp(i) < RT(0)) temp(i) = -temp(i);
        return temp;
    }

    template <class T>
    Matrix<T> HermBandSVDiv<T>::getV() const
    {
        Matrix<T> temp = pimpl->U.adjoint();
        const size_t N = pimpl->S.size();
        for(size_t i=0;i<N;i++) 
            if (pimpl->S(i) < 0) temp.row(i) *= RT(-1);
        return temp;
    }

    template <class T> 
    bool HermBandSVDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "HermBandSVDiv:\n";
            *fout << "M = "<<mm<<std::endl;
            *fout << "U = "<<getU()<<std::endl;
            *fout << "S = "<<getS()<<std::endl;
            *fout << "V = "<<getV()<<std::endl;
        }
        Matrix<T> usv = getU()*getS()*getV();
        RT nm = Norm(usv-mm);
        nm /= Norm(getU())*Norm(getS())*Norm(getV());
        RT cond = condition();
        if (fout) {
            *fout << "USV = "<<usv<<std::endl;
            *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
            *fout <<"  "<<cond<<" * "<<TMV_Epsilon<T>()<<std::endl;
        }
        return nm < cond*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T>
    size_t HermBandSVDiv<T>::colsize() const
    { return pimpl->S.size(); }

    template <class T>
    size_t HermBandSVDiv<T>::rowsize() const
    { return pimpl->S.size(); }

    template <class T>
    struct SymBandSVDiv<T>::SymBandSVDiv_Impl
    {
    public :
        SymBandSVDiv_Impl(const GenSymBandMatrix<T>& m) :
            U(m.size(),m.size()), S(m.size()), V(m.size(),m.size()), 
            logdet(0), signdet(1), kmax(0) {}

        Matrix<T,ColMajor> U;
        DiagMatrix<RT> S;
        Matrix<T,ColMajor> V;
        RT logdet;
        T signdet;
        mutable int kmax;
    }; // SymBandSVDiv

    template <class T>
    SymBandSVDiv<T>::SymBandSVDiv(const GenSymBandMatrix<T>& A) :
        pimpl(new SymBandSVDiv_Impl(A))
    {
        TMVAssert(isComplex(T()));

        SV_Decompose<T>(A,pimpl->U.view(),pimpl->S.view(),
                       pimpl->V.view(),pimpl->logdet,pimpl->signdet);
        thresh(TMV_Epsilon<T>());
    }

    template <class T>
    SymBandSVDiv<T>::~SymBandSVDiv() {}

    template <class T> template <class T1>
    void SymBandSVDiv<T>::doLDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(m.colsize() == pimpl->U.colsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
    }

    template <class T> template <class T1>
    void SymBandSVDiv<T>::doRDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(m.rowsize() == pimpl->U.rowsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
    }

    template <class T> template <class T1, class T2> 
    void SymBandSVDiv<T>::doLDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
    }

    template <class T> template <class T1, class T2> 
    void SymBandSVDiv<T>::doRDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(x.rowsize() == colsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
    }

    template <class T>
    T SymBandSVDiv<T>::det() const 
    { 
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }

    template <class T>
    RT SymBandSVDiv<T>::logDet(T* sign) const 
    { 
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;
    }

    template <class T> template <class T1>
    void SymBandSVDiv<T>::doMakeInverse(const SymMatrixView<T1>& sinv) const
    {
        TMVAssert(sinv.size() == pimpl->S.size());
        TMVAssert(sinv.issym());
        CallSymSV_Inverse(T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,sinv);
    }

    template <class T> template <class T1>
    void SymBandSVDiv<T>::doMakeInverse(const MatrixView<T1>& minv) const
    { 
        CallSymSV_Inverse(
            T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
            SymMatrixViewOf(minv,Upper));
        if (pimpl->S.size() > 1)
            minv.lowerTri().offDiag() = minv.upperTri().offDiag().transpose();
    }

    template <class T>
    void SymBandSVDiv<T>::doMakeInverseATA(const MatrixView<T>& minv) const
    {
        // A = U S V
        // At = Vt S Ut
        // AtA = Vt S^2 V
        // (AtA)^-1 = Vt S^-2 V
        //
        Matrix<T,RowMajor> SinvV =
            pimpl->V.rowRange(0,pimpl->kmax) /
            pimpl->S.subDiagMatrix(0,pimpl->kmax);
        minv = SinvV.adjoint() * SinvV;
    }

    template <class T> 
    bool SymBandSVDiv<T>::isSingular() const 
    { return pimpl->kmax < int(pimpl->S.size()); }

    template <class T>
    RT SymBandSVDiv<T>::norm2() const 
    { return pimpl->S.size() > 0 ? pimpl->S(0) : RT(0); }

    template <class T>
    RT SymBandSVDiv<T>::condition() const 
    { 
        return pimpl->S.size() > 0 ? 
            pimpl->S(0)/pimpl->S(pimpl->S.size()-1) : RT(1);
    }

    template <class T>
    void SymBandSVDiv<T>::thresh(RT toler, std::ostream* debugout) const
    {
        TMVAssert(toler < RT(1));
        RT thresh = pimpl->S(0)*toler;
        for(pimpl->kmax=pimpl->S.size(); 
            pimpl->kmax>0 && TMV_ABS(pimpl->S(pimpl->kmax-1))<=thresh; 
            --pimpl->kmax);
        if(debugout) {
            (*debugout)<<"S = "<<pimpl->S<<std::endl;
            (*debugout)<<"Smax = "<<pimpl->S(0)<<
                ", thresh = "<<thresh<<std::endl;
            (*debugout)<<"kmax = "<<pimpl->kmax;
            (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
        }
    }

    template <class T>
    void SymBandSVDiv<T>::top(int neigen,std::ostream* debugout) const
    {
        TMVAssert(neigen > 0);
        TMVAssert(neigen <= int(pimpl->S.size()));
        pimpl->kmax = neigen;
        if(debugout) {
            (*debugout)<<"S = "<<pimpl->S<<std::endl;
            (*debugout)<<"kmax = "<<pimpl->kmax;
            (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
        }
    }

    template <class T>
    int SymBandSVDiv<T>::getKMax() const 
    { return pimpl->kmax; }

    template <class T>
    ConstMatrixView<T> SymBandSVDiv<T>::getU() const
    { return pimpl->U.view(); }

    template <class T>
    ConstDiagMatrixView<RT> SymBandSVDiv<T>::getS() const
    { return pimpl->S.view(); }

    template <class T>
    ConstMatrixView<T> SymBandSVDiv<T>::getV() const
    { return pimpl->V.view(); }

    template <class T> 
    bool SymBandSVDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "SymBandSVDiv:\n";
            *fout << "M = "<<mm<<std::endl;
            *fout << "U = "<<getU()<<std::endl;
            *fout << "S = "<<getS()<<std::endl;
            *fout << "V = "<<getV()<<std::endl;
        }
        Matrix<T> usv = getU()*getS()*getV();
        RT nm = Norm(usv-mm);
        nm /= Norm(getU())*Norm(getS())*Norm(getV());
        RT cond = condition();
        if (fout) {
            *fout << "USV = "<<usv<<std::endl;
            *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
            *fout <<"  "<<cond<<" * "<<TMV_Epsilon<T>()<<std::endl;
        }
        return nm < cond*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T>
    size_t SymBandSVDiv<T>::colsize() const
    { return pimpl->S.size(); }

    template <class T>
    size_t SymBandSVDiv<T>::rowsize() const
    { return pimpl->S.size(); }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymBandSVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

