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



#include "tmv/TMV_SymSVD.h"
#include "TMV_SymSVDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include <ostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T>
    struct HermSVDiv<T>::HermSVDiv_Impl
    {
    public :
        HermSVDiv_Impl(const GenSymMatrix<T>& A, bool inplace);

        const bool inplace;
        AlignedArray<T> Uptr1;
        T* Uptr;
        MatrixView<T> U;
        DiagMatrix<RT> S;
        mutable RT logdet;
        mutable RT signdet;
        mutable bool donedet;
        mutable ptrdiff_t kmax;
    }; // HermSVDiv

#define UPTR1 (inplace ? 0 : (A.size()*A.size()))
#define UX (inplace ? \
            MatrixView<T>(A.nonConst().ptr(),A.size(),A.size(),\
                          A.stepi(),A.stepj(),NonConj) : \
            MatrixView<T>(Uptr1.get(),A.size(),A.size(),1,A.size(),NonConj))

    template <class T>
    HermSVDiv<T>::HermSVDiv_Impl::HermSVDiv_Impl(
        const GenSymMatrix<T>& A, bool _inplace) :
        inplace(_inplace && (A.iscm() || A.isrm())),
        Uptr1(UPTR1), U(UX), S(A.size()), 
        logdet(0), signdet(1), donedet(false), kmax(0) {}

    template <class T>
    HermSVDiv<T>::HermSVDiv(const GenSymMatrix<T>& A, bool inplace) :
        pimpl(new HermSVDiv_Impl(A,inplace))
    {
        TMVAssert(A.isherm());
        pimpl->U.lowerTri() = A.lowerTri();
        HermSV_Decompose<T>(pimpl->U.view(),pimpl->S.view());

        thresh(TMV_Epsilon<T>());
    }

    template <class T>
    HermSVDiv<T>::~HermSVDiv() {}

    template <class T> template <class T1> 
    void HermSVDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.colsize() == colsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,m);
    }

    template <class T> template <class T1> 
    void HermSVDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.rowsize() == rowsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,m);
    }

    template <class T> template <class T1, class T2> 
    void HermSVDiv<T>::doLDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,x);
    }

    template <class T> template <class T1, class T2> 
    void HermSVDiv<T>::doRDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(x.rowsize() == colsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->U.adjoint(),pimpl->kmax,m,x);
    }

    template <class T>
    T HermSVDiv<T>::det() const
    {
        if (!pimpl->donedet) {
            pimpl->logdet = pimpl->S.logDet(&pimpl->signdet);
            pimpl->donedet = true;
        }         
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }                  

    template <class T>
    RT HermSVDiv<T>::logDet(T* sign) const
    {
        if (!pimpl->donedet) {
            pimpl->logdet = pimpl->S.logDet(&pimpl->signdet);
            pimpl->donedet = true;
        }
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;  
    }                  
    template <class T> template <class T1> 
    void HermSVDiv<T>::doMakeInverse(SymMatrixView<T1> sinv) const
    {
        TMVAssert(sinv.size() == pimpl->S.size());
        TMVAssert(sinv.isherm());
        CallHermSV_Inverse(T(),pimpl->U,pimpl->S,pimpl->kmax,sinv);
    }

    template <class T> template <class T1>  
    void HermSVDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    { 
        // A^-1 = U S^-1 Ut
        if (isComplex(T1())) minv.diag().imagPart().setZero();
        CallHermSV_Inverse(
            T(),pimpl->U,pimpl->S,pimpl->kmax,HermMatrixViewOf(minv,Upper));
        if (pimpl->S.size() > 1)
            minv.lowerTri().offDiag() = minv.upperTri().offDiag().adjoint();
    }

    template <class T>
    void HermSVDiv<T>::doMakeInverseATA(MatrixView<T> minv) const
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
    bool HermSVDiv<T>::isSingular() const 
    { return pimpl->kmax < pimpl->S.size(); }

    template <class T>
    RT HermSVDiv<T>::norm2() const 
    { return pimpl->S.size() > 0 ? TMV_ABS(pimpl->S(0)) : RT(0); }

    template <class T>
    RT HermSVDiv<T>::condition() const 
    { 
        return pimpl->S.size() > 0 ? 
            TMV_ABS(pimpl->S(0)/pimpl->S(pimpl->S.size()-1)) : RT(1);
    }

    template <class T>
    void HermSVDiv<T>::thresh(RT toler, std::ostream* debugout) const
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
    void HermSVDiv<T>::top(ptrdiff_t neigen, std::ostream* debugout) const
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
    ptrdiff_t HermSVDiv<T>::getKMax() const 
    { return pimpl->kmax; }

    template <class T>
    ConstMatrixView<T> HermSVDiv<T>::getU() const
    { return pimpl->U.view(); }

    template <class T>
    DiagMatrix<RT> HermSVDiv<T>::getS() const
    { 
        DiagMatrix<RT> temp = pimpl->S;
        const ptrdiff_t N = pimpl->S.size();
        for(ptrdiff_t i=0;i<N;i++) 
            if (temp(i) < 0) temp(i) = -temp(i);
        return temp;
    }

    template <class T>
    Matrix<T> HermSVDiv<T>::getVt() const
    {
        Matrix<T> temp = pimpl->U.adjoint();
        const ptrdiff_t N = pimpl->S.size();
        for(ptrdiff_t i=0;i<N;i++) 
            if (pimpl->S(i) < 0) temp.row(i) *= RT(-1);
        return temp;
    }

    template <class T>
    bool HermSVDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "HermSVDiv:\n";
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
    ptrdiff_t HermSVDiv<T>::colsize() const
    { return pimpl->S.size(); }

    template <class T>
    ptrdiff_t HermSVDiv<T>::rowsize() const
    { return pimpl->S.size(); }

    template <class T>
    struct SymSVDiv<T>::SymSVDiv_Impl
    {
    public :
        SymSVDiv_Impl(const GenSymMatrix<T>& A, bool _inplace);

        const bool inplace;
        AlignedArray<T> Uptr1;
        MatrixView<T> U;
        DiagMatrix<RT> S;
        Matrix<T,ColMajor> Vt;
        RT logdet;
        T signdet;
        mutable ptrdiff_t kmax;
    }; // SymSVDiv

    template <class T>
    SymSVDiv<T>::SymSVDiv_Impl::SymSVDiv_Impl(
        const GenSymMatrix<T>& A, bool _inplace) :
        inplace(_inplace && (A.iscm() || A.isrm())),
        Uptr1(UPTR1), U(UX), S(A.size()), Vt(A.size(),A.size()),
        logdet(0), signdet(1), kmax(0) {}

#undef UPTR1
#undef UX


    template <class T>
    SymSVDiv<T>::SymSVDiv(const GenSymMatrix<T>& A, bool inplace) :
        pimpl(new SymSVDiv_Impl(A,inplace))
    {
        TMVAssert(isComplex(T()));
        pimpl->U.lowerTri() = A.lowerTri();
        SymSV_Decompose<T>(
            pimpl->U.view(),pimpl->S.view(),
            pimpl->Vt.view(),pimpl->logdet,pimpl->signdet);

        thresh(TMV_Epsilon<T>());
    }

    template <class T>
    SymSVDiv<T>::~SymSVDiv() {}

    template <class T> template <class T1> 
    void SymSVDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.colsize() == pimpl->U.colsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,m);
    }

    template <class T> template <class T1> 
    void SymSVDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.rowsize() == pimpl->U.rowsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,m);
    }

    template <class T> template <class T1, class T2> 
    void SymSVDiv<T>::doLDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,x);
    }

    template <class T> template <class T1, class T2> 
    void SymSVDiv<T>::doRDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(x.rowsize() == colsize());
        CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,x);
    }

    template <class T>
    T SymSVDiv<T>::det() const 
    { 
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }

    template <class T>
    RT SymSVDiv<T>::logDet(T* sign) const 
    { 
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;
    }

    template <class T> template <class T1> 
    void SymSVDiv<T>::doMakeInverse(SymMatrixView<T1> sinv) const
    {
        TMVAssert(sinv.size() == pimpl->S.size());
        TMVAssert(sinv.issym());
        CallSymSV_Inverse(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,sinv);
    }

    template <class T> template <class T1>  
    void SymSVDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    { 
        CallSymSV_Inverse(
            T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,
            SymMatrixViewOf(minv,Upper));
        if (pimpl->S.size() > 1)
            minv.lowerTri().offDiag() = minv.upperTri().offDiag().transpose();
    }

    template <class T>
    void SymSVDiv<T>::doMakeInverseATA(MatrixView<T> minv) const
    {
        // A = U S Vt
        // At = V S Ut
        // AtA = V S^2 Vt
        // (AtA)^-1 = V S^-2 Vt
        //
        Matrix<T,RowMajor> SinvVt = pimpl->Vt.rowRange(0,pimpl->kmax) /
            pimpl->S.subDiagMatrix(0,pimpl->kmax);
        minv = SinvVt.adjoint() * SinvVt;
    }

    template <class T>
    bool SymSVDiv<T>::isSingular() const 
    { return pimpl->kmax < pimpl->S.size(); }

    template <class T>
    RT SymSVDiv<T>::norm2() const 
    { return pimpl->S.size() > 0 ? pimpl->S(0) : RT(0); }

    template <class T>
    RT SymSVDiv<T>::condition() const 
    { 
        return pimpl->S.size() > 0 ? 
            pimpl->S(0)/pimpl->S(pimpl->S.size()-1) : RT(1);
    }

    template <class T>
    void SymSVDiv<T>::thresh(RT toler, std::ostream* debugout) const
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
    void SymSVDiv<T>::top(ptrdiff_t neigen, std::ostream* debugout) const
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
    ptrdiff_t SymSVDiv<T>::getKMax() const 
    { return pimpl->kmax; }

    template <class T>
    ConstMatrixView<T> SymSVDiv<T>::getU() const
    { return pimpl->U.view(); }

    template <class T>
    ConstDiagMatrixView<RT> SymSVDiv<T>::getS() const
    { return pimpl->S.view(); }

    template <class T>
    ConstMatrixView<T> SymSVDiv<T>::getVt() const
    { return pimpl->Vt.view(); }

    template <class T>
    bool SymSVDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "SymSVDiv:\n";
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
    ptrdiff_t SymSVDiv<T>::colsize() const
    { return pimpl->S.size(); }

    template <class T>
    ptrdiff_t SymSVDiv<T>::rowsize() const
    { return pimpl->S.size(); }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymSVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


