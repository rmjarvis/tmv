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



#include "tmv/TMV_BandQRD.h"
#include "TMV_BandQRDiv.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Householder.h"
#include "TMV_BandLUDiv.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include <ostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    struct BandQRDiv<T>::BandQRDiv_Impl
    {
    public :
        BandQRDiv_Impl(const GenBandMatrix<T>& A, bool _inplace);

        const bool istrans;
        const bool inplace;
        AlignedArray<T> Aptr1;
        T* Aptr;
        BandMatrixView<T> QRx;
        Vector<T> Qbeta;
        mutable RT logdet;
        mutable T signdet;
        mutable bool donedet;
    };

#define NEWLO (istrans ? A.nhi() : A.nlo())
#define NEWHI TMV_MIN(A.nlo()+A.nhi(),(istrans?A.colsize():A.rowsize())-1)
#define APTR1 inplace ? 0 : \
    BandStorageLength(ColMajor, istrans ? A.rowsize() : A.colsize(), \
                      istrans ? A.colsize() : A.rowsize(), NEWLO, NEWHI)
#define APTR inplace ? A.nonConst().ptr() : Aptr1.get()
#define QRX (istrans ? \
             (inplace ? A.nonConst().transpose() : \
              BandMatrixViewOf(Aptr,A.rowsize(),A.colsize(),A.nhi(), \
                               NEWHI, ColMajor)) : \
             (inplace ? A.nonConst().view() : \
              BandMatrixViewOf(Aptr,A.colsize(),A.rowsize(),A.nlo(),  \
                               NEWHI, ColMajor)))

    template <class T> 
    BandQRDiv<T>::BandQRDiv_Impl::BandQRDiv_Impl(
        const GenBandMatrix<T>& A, bool _inplace) :
        istrans(A.colsize()<A.rowsize() || (A.isSquare() && A.nhi()<A.nlo())
                || (A.isSquare() && A.nhi()==A.nlo() && A.isrm())),
        inplace(_inplace || NEWLO == 0),
        Aptr1(APTR1), Aptr(APTR), QRx(QRX), Qbeta(QRx.rowsize()),
        logdet(0), signdet(1), donedet(false) {}

#undef NEWLO
#undef NEWHI
#undef APTR
#undef QRX

    template <class T> 
    BandQRDiv<T>::BandQRDiv(
        const GenBandMatrix<T>& A, bool inplace) :
        pimpl(new BandQRDiv_Impl(A,inplace))
    {
        if (inplace) {
            // For inplace decomposition, make sure the original band matrix
            // has room for the extra upper diagonals...
            // if isrm stepi >= (2*A.nlo()+A.nhi())
            // if iscm stepj >= (2*A.nlo()+A.nhi())
            // if isdm extra diags appear at end, so can't really check
            TMVAssert(!pimpl->QRx.isrm() || 
                      pimpl->QRx.stepi()>=pimpl->QRx.nlo()+pimpl->QRx.nhi());
            TMVAssert(!pimpl->QRx.iscm() || 
                      pimpl->QRx.stepj()>=pimpl->QRx.nlo()+pimpl->QRx.nhi());
            TMVAssert(pimpl->QRx == (
                    pimpl->istrans ? A.transpose() : A.view()));
            if (pimpl->QRx.nlo() > 0)
                pimpl->QRx.diagRange(pimpl->QRx.nhi()-pimpl->QRx.nlo()+1,
                                     pimpl->QRx.nhi()+1).setZero();
        } else {
            if (pimpl->istrans) pimpl->QRx = A.transpose();
            else pimpl->QRx = A;
        }
        TMVAssert(pimpl->QRx.colsize() == 
                  (pimpl->istrans?A.rowsize():A.colsize()));
        TMVAssert(pimpl->QRx.rowsize() ==
                  (pimpl->istrans?A.colsize():A.rowsize()));
        if (pimpl->QRx.nlo() > 0) 
            QR_Decompose(pimpl->QRx,pimpl->Qbeta.view(),pimpl->signdet);
    }

    template <class T> 
    BandQRDiv<T>::~BandQRDiv() {}

    template <class T> template <class T1> 
    void BandQRDiv<T>::doLDivEq(MatrixView<T1> m) const
    {
        if (pimpl->istrans)
            QR_RDivEq(pimpl->QRx,pimpl->Qbeta,m.transpose());
        else 
            QR_LDivEq(pimpl->QRx,pimpl->Qbeta,m);
    }

    template <class T> template <class T1> 
    void BandQRDiv<T>::doRDivEq(MatrixView<T1> m) const
    {
        if (pimpl->istrans) QR_LDivEq(pimpl->QRx,pimpl->Qbeta,m.transpose());
        else QR_RDivEq(pimpl->QRx,pimpl->Qbeta,m);
    }

    template <class T> template <class T1, class T2> 
    void BandQRDiv<T>::doLDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        if (pimpl->istrans) 
            QR_RDiv(pimpl->QRx,pimpl->Qbeta,m.transpose(),x.transpose());
        else 
            QR_LDiv(pimpl->QRx,pimpl->Qbeta,m,x);
    }

    template <class T> template <class T1, class T2> 
    void BandQRDiv<T>::doRDiv(
        const GenMatrix<T1>& m, MatrixView<T2> x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(x.rowsize() == colsize());
        if (pimpl->istrans) 
            QR_LDiv(pimpl->QRx,pimpl->Qbeta,m.transpose(),x.transpose());
        else 
            QR_RDiv(pimpl->QRx,pimpl->Qbeta,m,x);
    }

    template <class T> 
    T BandQRDiv<T>::det() const
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
    RT BandQRDiv<T>::logDet(T* sign) const
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
    void BandQRDiv<T>::doMakeInverse(MatrixView<T1> minv) const
    {
        MatrixView<T1> minv2 = pimpl->istrans ? minv.transpose() : minv;
        TMVAssert(pimpl->QRx.colsize() >= pimpl->QRx.rowsize());
        TMVAssert(minv2.colsize() == pimpl->QRx.rowsize());
        TMVAssert(minv2.rowsize() == pimpl->QRx.colsize());
        QR_Inverse(pimpl->QRx,pimpl->Qbeta,minv2);
    }

    template <class T> 
    void BandQRDiv<T>::doMakeInverseATA(MatrixView<T> minv) const
    {
        TMVAssert(minv.colsize() == pimpl->QRx.rowsize());
        TMVAssert(minv.rowsize() == pimpl->QRx.rowsize());
        TMVAssert(pimpl->QRx.colsize() >= pimpl->QRx.rowsize());
        const ptrdiff_t N = pimpl->QRx.rowsize();
        // At A = Rt R
        // (At A)^-1 = (Rt R)^-1 = R^-1 Rt^-1
        UpperTriMatrixView<T> Rinv = minv.upperTri();
        Rinv = pimpl->QRx.subBandMatrix(0,N,0,N,0,pimpl->QRx.nhi());
        TriInverse(Rinv,pimpl->QRx.nhi());
        minv = Rinv * Rinv.transpose();
    }

    template <class T> 
    bool BandQRDiv<T>::isSingular() const
    {
        return pimpl->QRx.diag().minAbs2Element() <=
            TMV_Epsilon<T>() * pimpl->QRx.diag().maxAbs2Element(); 
    }

    template <class T> 
    bool BandQRDiv<T>::isTrans() const
    { return pimpl->istrans; }

    template <class T> 
    void GetQFromBandQR(
        MatrixView<T> Q, const GenVector<T>& Qbeta, const ptrdiff_t nlo) 
    {
        // Extract the Q matrix from a combined QRx matrix
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Qbeta.size() == Q.rowsize());
        const ptrdiff_t M = Q.colsize();
        const ptrdiff_t N = Q.rowsize();
        Q.upperTri().setZero();
        for(ptrdiff_t j=N-1;j>=0;j--) {
            if (j+nlo+1 > M)
                HouseholderUnpack(Q.subMatrix(j,M,j,N),Qbeta(j));
            else
                HouseholderUnpack(Q.subMatrix(j,j+nlo+1,j,N),Qbeta(j));
        }
    }

    template <class T> 
    Matrix<T> BandQRDiv<T>::getQ() const 
    {
        // Extract the Q matrix from a combined QRx matrix
        TMVAssert(pimpl->QRx.colsize() >= pimpl->QRx.rowsize());
        Matrix<T> temp(pimpl->QRx.colsize(),pimpl->QRx.rowsize());
        temp = pimpl->QRx;
        if (pimpl->QRx.nlo() == 0) temp.setToIdentity();
        else GetQFromBandQR(temp.view(),pimpl->Qbeta.view(),pimpl->QRx.nlo());
        return temp;
    }

    template <class T> 
    ConstBandMatrixView<T> BandQRDiv<T>::getR() const
    { 
        return pimpl->QRx.subBandMatrix(
            0,pimpl->QRx.rowsize(),0,pimpl->QRx.rowsize(),0,pimpl->QRx.nhi()); 
    }
    template <class T> 
    const GenBandMatrix<T>& BandQRDiv<T>::getQR() const 
    { return pimpl->QRx; }
    template <class T> 
    const GenVector<T>& BandQRDiv<T>::getQBeta() const 
    { return pimpl->Qbeta; }

    template <class T> 
    bool BandQRDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        Matrix<T> Q = getQ();
        if (fout) {
            *fout <<"BandQRDiv:\n";
            *fout << "M = "<<
                (pimpl->istrans?mm.transpose():mm.view())<<std::endl;
            *fout << "Q = "<<Q<<std::endl;
            *fout << "R = "<<getR()<<std::endl;
        }
        Matrix<T> qr = Q*getR();
        RT nm = Norm(qr-(pimpl->istrans ? mm.transpose() : mm.view()));
        nm /= Norm(Q)*Norm(getR());
        if (fout) {
            *fout << "QR = "<<qr<<std::endl;
            *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<std::endl;
        }
        return nm < mm.doCondition()*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T> 
    ptrdiff_t BandQRDiv<T>::colsize() const
    { return pimpl->istrans ? pimpl->QRx.rowsize() : pimpl->QRx.colsize(); }

    template <class T> 
    ptrdiff_t BandQRDiv<T>::rowsize() const
    { return pimpl->istrans ? pimpl->QRx.colsize() : pimpl->QRx.rowsize(); }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandQRD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv
