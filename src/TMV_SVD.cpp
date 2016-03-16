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
        Matrix<T,ColMajor> Vt;
        RT logdet;
        T signdet;
        mutable ptrdiff_t kmax;
    };

#define APTR1 (inplace ? 0 : (A.colsize()*A.rowsize()))
#define APTR (inplace ? A.nonConst().ptr() : Aptr1.get())
#define UX (istrans ? \
            (inplace ? A.nonConst().transpose() : \
             MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor)) : \
            (inplace ? A.nonConst().view() : \
             MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)))

    template <class T> 
    SVDiv<T>::SVDiv_Impl::SVDiv_Impl(const GenMatrix<T>& A, bool _inplace) :
        istrans(A.colsize() < A.rowsize()), 
        inplace(_inplace && (A.isrm() || A.iscm())), Aptr1(APTR1), Aptr(APTR),
        U(UX), S(U.rowsize()), Vt(U.rowsize(),U.rowsize()), 
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

        SV_Decompose<T>(pimpl->U,pimpl->S.view(),pimpl->Vt.view(),
                        pimpl->logdet, pimpl->signdet,true);

        // Set kmax for actual 0 elements (to within machine precision).
        // Any further cut in the number of singular values to use
        // should be done by the user.
        thresh(TMV_Epsilon<T>());
    }

    template <class T> 
    SVDiv<T>::~SVDiv() {}

    template <class T> template <class T1> 
    void SVDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.colsize() == colsize());
        TMVAssert(m.colsize() == rowsize());
        if (pimpl->istrans) 
            CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,
                    m.transpose(),m.transpose());
        else 
            CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,m);
    }

    template <class T> template <class T1> 
    void SVDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(m.rowsize() == colsize());
        TMVAssert(m.rowsize() == rowsize());

        if (pimpl->istrans) 
            CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,
                    m.transpose(),m.transpose());
        else 
            CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,m);
    }

    template <class T> template <class T1, class T2> 
    void SVDiv<T>::doLDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        if (pimpl->istrans) 
            CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,
                   m.transpose(),x.transpose());
        else 
            CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,x);
    }

    template <class T> template <class T1, class T2> 
    void SVDiv<T>::doRDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(x.rowsize() == colsize());

        if (pimpl->istrans) 
            CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,
                   m.transpose(),x.transpose());
        else 
            CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->Vt,pimpl->kmax,m,x);
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
    void SVDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    { 
        if (pimpl->istrans) {
            // A^-1 = (V S^-1 Ut)T = U* S^-1 Vt*
            Matrix<T,ColMajor> SinvVt =
                pimpl->Vt.conjugate().rowRange(0,pimpl->kmax) /
                pimpl->S.subDiagMatrix(0,pimpl->kmax);
            minv = pimpl->U.conjugate().colRange(0,pimpl->kmax) * SinvVt;
        } else {
            // A^-1 = V S^-1 Ut
            Matrix<T,ColMajor> SinvUt =
                pimpl->U.adjoint().rowRange(0,pimpl->kmax) /
                pimpl->S.subDiagMatrix(0,pimpl->kmax);
            minv = pimpl->Vt.adjoint().colRange(0,pimpl->kmax) * SinvUt;
        }
    }

    template <class T> 
    void SVDiv<T>::doMakeInverseATA(MatrixView<T> minv) const
    {
        // A = U S Vt
        // At = V S Ut
        // AtA = V S^2 Vt
        // (AtA)^-1 = V S^-2 Vt
        //
        // if istrans:
        // AT = U S Vt
        // At = U* S VT
        // AAt = V* S^2 VT
        // (AAt)^-1 = V* S^-2 VT
        //
        Matrix<T,ColMajor> SinvVt = pimpl->Vt.rowRange(0,pimpl->kmax) /
            pimpl->S.subDiagMatrix(0,pimpl->kmax);
        if (pimpl->istrans)
            minv = SinvVt.transpose() * SinvVt.conjugate();
        else
            minv = SinvVt.adjoint() * SinvVt;
    }

    template <class T> 
    bool SVDiv<T>::isSingular() const
    { return pimpl->kmax < pimpl->S.size(); }

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
    void SVDiv<T>::top(ptrdiff_t neigen, std::ostream* debugout) const
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
    ptrdiff_t SVDiv<T>::getKMax() const
    { return pimpl->kmax; }

    template <class T> 
    ConstMatrixView<T> SVDiv<T>::getU() const
    {
        if (pimpl->istrans) return pimpl->Vt.transpose(); 
        else return pimpl->U.view(); 
    }
    template <class T> 
    ConstDiagMatrixView<RT> SVDiv<T>::getS() const
    { return pimpl->S.view(); }

    template <class T> 
    ConstMatrixView<T> SVDiv<T>::getVt() const
    {
        if (pimpl->istrans) return pimpl->U.transpose(); 
        else return pimpl->Vt.view(); 
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
            *fout << "Vt = "<<getVt()<<std::endl;
        }
        Matrix<T> usv = getU()*getS()*getVt();
        RT nm = Norm(usv-mm);
        nm /= Norm(getU())*Norm(getS())*Norm(getVt());
        RT cond = getS()(0) / getS()(getKMax()-1);
        if (fout) {
            *fout << "USVt = "<<usv<<std::endl;
            *fout << "Norm(M-USVt)/Norm(USVt) = "<<nm;
            *fout <<"  "<<cond<<" * "<<TMV_Epsilon<T>()<<std::endl;
        }
        return nm < cond*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T> 
    ptrdiff_t SVDiv<T>::colsize() const
    { return pimpl->istrans ? pimpl->U.rowsize() : pimpl->U.colsize(); }

    template <class T> 
    ptrdiff_t SVDiv<T>::rowsize() const
    { return pimpl->istrans ? pimpl->U.colsize() : pimpl->U.rowsize(); }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


