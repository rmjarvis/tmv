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
        mutable ptrdiff_t kmax;
    }; // HermBandSVDiv

    template <class T>
    HermBandSVDiv<T>::HermBandSVDiv(const GenSymBandMatrix<T>& A) :
        pimpl(new HermBandSVDiv_Impl(A))
    {
        TMVAssert(A.isherm());

        MatrixView<T> Vt(0,0,0,1,1,NonConj);
        SV_Decompose<T>(A, pimpl->U.view(), pimpl->S.view(), Vt,
                       pimpl->logdet, pimpl->signdet);
        thresh(TMV_Epsilon<T>());
    }

    template <class T>
    HermBandSVDiv<T>::~HermBandSVDiv() {}

    template <class T> template <class T1>
    void HermBandSVDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.colsize() == colsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,m);
    }

    template <class T> template <class T1>
    void HermBandSVDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.rowsize() == rowsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,m);
    }

    template <class T> template <class T1, class T2> 
    void HermBandSVDiv<T>::doLDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,x);
    }

    template <class T> template <class T1, class T2> 
    void HermBandSVDiv<T>::doRDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
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
    void HermBandSVDiv<T>::doMakeInverse(SymMatrixView<T1> sinv) const
    {
        TMVAssert(sinv.size() == pimpl->S.size());
        TMVAssert(sinv.isherm());
        CallHermSV_Inverse(T(),pimpl->U,pimpl->S,pimpl->kmax,sinv);
    }

    template <class T> template <class T1>
    void HermBandSVDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    { 
        // A^-1 = U S^-1 Ut
        if (isComplex(T1())) minv.diag().imagPart().setZero();
        CallHermSV_Inverse(
            T(),pimpl->U,pimpl->S,pimpl->kmax,HermMatrixViewOf(minv,Upper));
        if (pimpl->S.size() > 1)
            minv.lowerTri().offDiag() = minv.upperTri().offDiag().adjoint();
    }

    template <class T> 
    void HermBandSVDiv<T>::doMakeInverseATA(MatrixView<T> minv) const
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
    { return pimpl->kmax < pimpl->S.size(); }

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
    void HermBandSVDiv<T>::top(ptrdiff_t neigen, std::ostream* debugout) const
    {
        TMVAssert(neigen > 0);
        if (neigen < pimpl->S.size()) pimpl->kmax = neigen;
        else pimpl->kmax = pimpl->S.size();
        if(debugout) {
            (*debugout)<<"S = "<<pimpl->S<<std::endl;
            (*debugout)<<"kmax = "<<pimpl->kmax;
            (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
        }
    }

    template <class T>
    ptrdiff_t HermBandSVDiv<T>::getKMax() const 
    { return pimpl->kmax; }

    template <class T>
    ConstMatrixView<T> HermBandSVDiv<T>::getU() const
    { return pimpl->U.view(); }

    template <class T>
    DiagMatrix<RT> HermBandSVDiv<T>::getS() const
    {
        DiagMatrix<RT> temp = pimpl->S;
        const ptrdiff_t N = pimpl->S.size();
        for(ptrdiff_t i=0;i<N;i++) 
            if (temp(i) < RT(0)) temp(i) = -temp(i);
        return temp;
    }

    template <class T>
    Matrix<T> HermBandSVDiv<T>::getVt() const
    {
        Matrix<T> temp = pimpl->U.adjoint();
        const ptrdiff_t N = pimpl->S.size();
        for(ptrdiff_t i=0;i<N;i++) 
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
            *fout << "Vt = "<<getVt()<<std::endl;
        }
        Matrix<T> usv = getU()*getS()*getVt();
        RT nm = Norm(usv-mm);
        nm /= Norm(getU())*Norm(getS())*Norm(getVt());
        RT cond = condition();
        if (fout) {
            *fout << "USVt = "<<usv<<std::endl;
            *fout << "Norm(M-USVt)/Norm(USVt) = "<<nm;
            *fout <<"  "<<cond<<" * "<<TMV_Epsilon<T>()<<std::endl;
        }
        return nm < cond*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T>
    ptrdiff_t HermBandSVDiv<T>::colsize() const
    { return pimpl->S.size(); }

    template <class T>
    ptrdiff_t HermBandSVDiv<T>::rowsize() const
    { return pimpl->S.size(); }

    template <class T>
    struct SymBandSVDiv<T>::SymBandSVDiv_Impl
    {
    public :
        SymBandSVDiv_Impl(const GenSymBandMatrix<T>& m) :
            U(m.size(),m.size()), S(m.size()), Vt(m.size(),m.size()), 
            logdet(0), signdet(1), kmax(0) {}

        Matrix<T,ColMajor> U;
        DiagMatrix<RT> S;
        Matrix<T,ColMajor> Vt;
        RT logdet;
        T signdet;
        mutable ptrdiff_t kmax;
    }; // SymBandSVDiv

    template <class T>
    SymBandSVDiv<T>::SymBandSVDiv(const GenSymBandMatrix<T>& A) :
        pimpl(new SymBandSVDiv_Impl(A))
    {
        TMVAssert(isComplex(T()));

        SV_Decompose<T>(A,pimpl->U.view(),pimpl->S.view(),
                       pimpl->Vt.view(),pimpl->logdet,pimpl->signdet);
        thresh(TMV_Epsilon<T>());
    }

    template <class T>
    SymBandSVDiv<T>::~SymBandSVDiv() {}

    template <class T> template <class T1>
    void SymBandSVDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.colsize() == pimpl->U.colsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,m);
    }

    template <class T> template <class T1>
    void SymBandSVDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.rowsize() == pimpl->U.rowsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,m);
    }

    template <class T> template <class T1, class T2> 
    void SymBandSVDiv<T>::doLDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,x);
    }

    template <class T> template <class T1, class T2> 
    void SymBandSVDiv<T>::doRDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(x.rowsize() == colsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,x);
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
    void SymBandSVDiv<T>::doMakeInverse(SymMatrixView<T1> sinv) const
    {
        TMVAssert(sinv.size() == pimpl->S.size());
        TMVAssert(sinv.issym());
        CallSymSV_Inverse(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,sinv);
    }

    template <class T> template <class T1>
    void SymBandSVDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    { 
        CallSymSV_Inverse(
            T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,
            SymMatrixViewOf(minv,Upper));
        if (pimpl->S.size() > 1)
            minv.lowerTri().offDiag() = minv.upperTri().offDiag().transpose();
    }

    template <class T>
    void SymBandSVDiv<T>::doMakeInverseATA(MatrixView<T> minv) const
    {
        // A = U S Vt
        // At = V S Ut
        // AtA = V S^2 Vt
        // (AtA)^-1 = V S^-2 Vt
        //
        Matrix<T,RowMajor> SinvVt =
            pimpl->Vt.rowRange(0,pimpl->kmax) /
            pimpl->S.subDiagMatrix(0,pimpl->kmax);
        minv = SinvVt.adjoint() * SinvVt;
    }

    template <class T> 
    bool SymBandSVDiv<T>::isSingular() const 
    { return pimpl->kmax < pimpl->S.size(); }

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
    void SymBandSVDiv<T>::top(ptrdiff_t neigen,std::ostream* debugout) const
    {
        TMVAssert(neigen > 0);
        TMVAssert(neigen <= pimpl->S.size());
        pimpl->kmax = neigen;
        if(debugout) {
            (*debugout)<<"S = "<<pimpl->S<<std::endl;
            (*debugout)<<"kmax = "<<pimpl->kmax;
            (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
        }
    }

    template <class T>
    ptrdiff_t SymBandSVDiv<T>::getKMax() const 
    { return pimpl->kmax; }

    template <class T>
    ConstMatrixView<T> SymBandSVDiv<T>::getU() const
    { return pimpl->U.view(); }

    template <class T>
    ConstDiagMatrixView<RT> SymBandSVDiv<T>::getS() const
    { return pimpl->S.view(); }

    template <class T>
    ConstMatrixView<T> SymBandSVDiv<T>::getVt() const
    { return pimpl->Vt.view(); }

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
            *fout << "Vt = "<<getVt()<<std::endl;
        }
        Matrix<T> usv = getU()*getS()*getVt();
        RT nm = Norm(usv-mm);
        nm /= Norm(getU())*Norm(getS())*Norm(getVt());
        RT cond = condition();
        if (fout) {
            *fout << "USVt = "<<usv<<std::endl;
            *fout << "Norm(M-USVt)/Norm(USVt) = "<<nm;
            *fout <<"  "<<cond<<" * "<<TMV_Epsilon<T>()<<std::endl;
        }
        return nm < cond*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T>
    ptrdiff_t SymBandSVDiv<T>::colsize() const
    { return pimpl->S.size(); }

    template <class T>
    ptrdiff_t SymBandSVDiv<T>::rowsize() const
    { return pimpl->S.size(); }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymBandSVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


