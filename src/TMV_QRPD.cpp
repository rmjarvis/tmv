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



#include "tmv/TMV_QRPD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "TMV_QRDiv.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_PackedQ.h"
#include "tmv/TMV_PermutationArith.h"
#include <ostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    struct QRPDiv<T>::QRPDiv_Impl
    {
    public :
        QRPDiv_Impl(const GenMatrix<T>& m, bool inplace);

        const bool istrans;
        const bool inplace;
        AlignedArray<T> Aptr1;
        T* Aptr;
        MatrixView<T> QRx;
        Vector<T> beta;
        Permutation P;
        mutable RT logdet;
        mutable T signdet;
        mutable bool donedet;
        ptrdiff_t N1;
    };

#define APTR1 (inplace ? 0 : (A.colsize()*A.rowsize()))
#define APTR (inplace ? A.nonConst().ptr() : Aptr1.get())
#define QRX (istrans ? \
             (inplace ? A.nonConst().transpose() : \
              MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor)) : \
             (inplace ? A.nonConst().view() : \
              MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)))

    template <class T> 
    QRPDiv<T>::QRPDiv_Impl::QRPDiv_Impl(
        const GenMatrix<T>& A, bool _inplace) :
        istrans(A.colsize()<A.rowsize()),
        inplace(_inplace && (A.iscm() || A.isrm())), 
        Aptr1(APTR1), Aptr(APTR), QRx(QRX), 
        beta(QRx.rowsize()), P(beta.size()),
        logdet(0), signdet(1), donedet(false), N1(beta.size()) {}

#undef QRX
#undef APTR

    template <class T> 
    QRPDiv<T>::QRPDiv(const GenMatrix<T>& A, bool inplace) :
        pimpl(new QRPDiv_Impl(A,inplace)) 
    {
        if (pimpl->istrans) {
            if (inplace) TMVAssert(A.transpose() == pimpl->QRx); 
            else pimpl->QRx = A.transpose();
        } else {
            if (inplace) TMVAssert(A == pimpl->QRx); 
            else pimpl->QRx = A;
        }
        QRP_Decompose(
            pimpl->QRx,pimpl->beta.view(),pimpl->P,pimpl->signdet,
            QRP_IsStrict());
        ptrdiff_t sd = pimpl->QRx.diag().step();
        TMVAssert(pimpl->N1>0);
        const T* endofdiag = pimpl->QRx.diag().cptr() + (pimpl->N1-1)*sd;
        while(pimpl->N1>0 && *endofdiag==T(0)) { 
            --pimpl->N1;
            endofdiag -= sd; 
        }
    }

    template <class T> 
    QRPDiv<T>::~QRPDiv() {}

    template <class T> template <class T1> 
    void QRPDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(pimpl->QRx.isSquare());
        TMVAssert(m.colsize() == pimpl->QRx.colsize());
        if (pimpl->istrans) 
            QR_LDivEq(
                pimpl->QRx,pimpl->beta,pimpl->P.getValues(),m.transpose(),
                pimpl->N1);
        else 
            QR_LDivEq(pimpl->QRx,pimpl->beta,pimpl->P.getValues(),m,pimpl->N1);
    }

    template <class T> template <class T1> 
    void QRPDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(pimpl->QRx.isSquare());
        TMVAssert(m.rowsize() == pimpl->QRx.rowsize());
        if (pimpl->istrans) 
            QR_RDivEq(pimpl->QRx,pimpl->beta,pimpl->P.getValues(),m.transpose(),
                      pimpl->N1);
        else 
            QR_RDivEq(pimpl->QRx,pimpl->beta,pimpl->P.getValues(),m,pimpl->N1);
    }

    template <class T> template <class T1, class T2> 
    void QRPDiv<T>::doLDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == (pimpl->istrans ? pimpl->QRx.rowsize() :
                                  pimpl->QRx.colsize()));
        TMVAssert(x.colsize() == (pimpl->istrans ? pimpl->QRx.colsize() : 
                                  pimpl->QRx.rowsize()));
        if (pimpl->istrans) 
            QR_RDiv(pimpl->QRx,pimpl->beta,pimpl->P.getValues(),
                    m.transpose(),x.transpose(),pimpl->N1);
        else 
            QR_LDiv(pimpl->QRx,pimpl->beta,pimpl->P.getValues(),m,x,pimpl->N1);
    }

    template <class T> template <class T1, class T2> 
    void QRPDiv<T>::doRDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == (pimpl->istrans ? pimpl->QRx.colsize() :
                                  pimpl->QRx.rowsize()));
        TMVAssert(x.rowsize() == (pimpl->istrans ? pimpl->QRx.rowsize() :
                                  pimpl->QRx.colsize()));
        if (pimpl->istrans) QR_LDiv(pimpl->QRx,pimpl->beta,pimpl->P.getValues(),
                                    m.transpose(),x.transpose(),pimpl->N1);
        else QR_RDiv(pimpl->QRx,pimpl->beta,pimpl->P.getValues(),m,x,pimpl->N1);
    }


    template <class T> 
    T QRPDiv<T>::det() const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->QRx.diag()).logDet(&s);
            pimpl->signdet *= s;
            pimpl->donedet = true;
        }         
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }                  

    template <class T> 
    RT QRPDiv<T>::logDet(T* sign) const
    {
        if (!pimpl->donedet) {
            T s;
            pimpl->logdet = DiagMatrixViewOf(pimpl->QRx.diag()).logDet(&s);
            pimpl->signdet *= s;
            pimpl->donedet = true;
        }
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;  
    }                  

    template <class T> template <class T1> 
    void QRPDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    {
        if (pimpl->istrans)
            QR_Inverse(pimpl->QRx,pimpl->beta,pimpl->P.getValues(),
                       minv.transpose(),pimpl->N1);
        else
            QR_Inverse(pimpl->QRx,pimpl->beta,pimpl->P.getValues(),
                       minv,pimpl->N1);
    }

    template <class T> 
    void QRPDiv<T>::doMakeInverseATA(MatrixView<T> ata) const
    {
        // At A = Pt Rt R P
        // (At A)^-1 = (Pt Rt R P)^-1 = Pt R^-1 R^-1t P
        UpperTriMatrix<T> rinv = pimpl->QRx.upperTri().inverse();
        ata = rinv * rinv.adjoint();
        ata.reversePermuteRows(pimpl->P.getValues());
        ata.reversePermuteCols(pimpl->P.getValues());
    }

    template <class T> 
    bool QRPDiv<T>::isSingular() const
    { 
        return pimpl->QRx.diag().minAbs2Element() <=
            TMV_Epsilon<T>() * pimpl->QRx.diag().maxAbs2Element(); 
    }

    template <class T> 
    bool QRPDiv<T>::isTrans() const
    { return pimpl->istrans; }

    template <class T> 
    PackedQ<T> QRPDiv<T>::getQ() const
    { return PackedQ<T>(pimpl->QRx,pimpl->beta); }

    template <class T> 
    ConstUpperTriMatrixView<T> QRPDiv<T>::getR() const
    { return pimpl->QRx.upperTri(); }

    template <class T> 
    const GenMatrix<T>& QRPDiv<T>::getQRx() const 
    { return pimpl->QRx; }

    template <class T> 
    const GenVector<T>& QRPDiv<T>::getBeta() const 
    { return pimpl->beta; }

    template <class T> 
    const Permutation& QRPDiv<T>::getP() const 
    { return pimpl->P; }

    template <class T> 
    bool QRPDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
        if (printmat) {
            *fout << "QRPDiv:\n";
            *fout << "M = "<<
                (pimpl->istrans?mm.transpose():mm.view())<<std::endl;
            *fout << "Q = "<<getQ()<<std::endl;
            *fout << "R = "<<getR()<<std::endl;
            *fout << "P = "<<getP()<<std::endl;
            *fout << "  or by interchanges: ";
            for(ptrdiff_t i=0;i<getP().size();i++)
                *fout<<(getP().getValues())[i]<<" ";
            *fout<<std::endl;
        }
        Matrix<T> qr = getQ()*getR()*getP();
        RT nm = Norm(qr-(pimpl->istrans ? mm.transpose() : mm.view()));
        nm /= Norm(getQ())*Norm(getR());
        if (printmat) {
            *fout << "QRP = "<<qr<<std::endl;
        }
        if (fout) {
            *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<"  ";
            *fout <<pimpl->QRx.rowsize()*TMV_Epsilon<T>()<<std::endl;
        }
        return nm < mm.doCondition()*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T> 
    ptrdiff_t QRPDiv<T>::colsize() const
    { return pimpl->istrans ? pimpl->QRx.rowsize() : pimpl->QRx.colsize(); }

    template <class T> 
    ptrdiff_t QRPDiv<T>::rowsize() const
    { return pimpl->istrans ? pimpl->QRx.colsize() : pimpl->QRx.rowsize(); }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_QRPD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


