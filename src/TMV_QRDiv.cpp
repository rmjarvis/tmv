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

#include "TMV_QRDiv.h"
#include "tmv/TMV_QRD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"

namespace tmv {

    //
    // LDiv
    //

    template <class T, class T1, class T2> 
    void QR_LDiv(
        const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const ptrdiff_t* P,
        const GenMatrix<T2>& m, MatrixView<T> x, ptrdiff_t N1)
    {
        TMVAssert(QRx.colsize() >= QRx.rowsize());
        TMVAssert(beta.size() == QRx.rowsize());
        TMVAssert(m.colsize() == QRx.colsize());
        TMVAssert(x.colsize() == QRx.rowsize());
        TMVAssert(x.rowsize() == m.rowsize());
        TMVAssert(QRx.isrm() || QRx.iscm());
        TMVAssert(QRx.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        // Solve Q R P x = m
        // where Q and R are stored in QRx, and beta are the beta

        //std::cout<<"Start QR_LDiv"<<std::endl;
        //std::cout<<"QRx = "<<QRx<<std::endl;
        //std::cout<<"beta = "<<beta<<std::endl;
        //std::cout<<"m = "<<m<<std::endl;
        //std::cout<<"Norm(QRx) = "<<Norm(QRx)<<std::endl;
        //std::cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
        //std::cout<<"Norm(m) = "<<Norm(m)<<std::endl;

        // First Solve Q y = m
        if (QRx.isSquare()) {
            x = m;
            Q_LDivEq(QRx,beta,x);
            //std::cout<<"x = m/Q = "<<x<<std::endl;
        } else {
            if (m.isrm()) {
                Matrix<T,RowMajor> m1 = m;
                // Q is Q1 [ I ]
                //         [ 0 ]
                // where Q1 is the part of Q that is stored in QRx and beta
                // m1 = Q^-1 m1
                Q_LDivEq(QRx,beta,m1.view());
                // y = [ I 0 ] m1
                x = m1.rowRange(0,x.colsize()); // x = y here
            } else {
                Matrix<T,ColMajor> m1 = m;
                Q_LDivEq(QRx,beta,m1.view());
                x = m1.rowRange(0,x.colsize()); // x = y here
            }
            //std::cout<<"x = m/Q (non-sq) = "<<x<<std::endl;
        }

        // Now solve R z = y
        x.rowRange(N1,x.colsize()).setZero();
        //std::cout<<"x(N1:N).zero = "<<x<<std::endl;
        //std::cout<<"1 Norm(x) = "<<Norm(x)<<std::endl;

        //x.rowRange(0,N1) /= QRx.upperTri().subTriMatrix(0,N1);
        QRx.upperTri().subTriMatrix(0,N1).LDivEq(x.rowRange(0,N1));
        //std::cout<<"x /= R = "<<x<<std::endl;
        //std::cout<<"2 Norm(x) = "<<Norm(x)<<std::endl;

        // Finally P x = z
        if (P) x.reversePermuteRows(P);
        //std::cout<<"x /= P = "<<x<<std::endl;
        //std::cout<<"3 Norm(x) = "<<Norm(x)<<std::endl;
    }

    //
    // LDivEq
    //

    template <class T, class T1> 
    void QR_LDivEq(
        const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const ptrdiff_t* P,
        MatrixView<T> m, ptrdiff_t N1)
    {
        TMVAssert(QRx.colsize() == QRx.rowsize());
        TMVAssert(beta.size() == QRx.rowsize());
        TMVAssert(m.colsize() == QRx.colsize());
        TMVAssert(QRx.isrm() || QRx.iscm());
        TMVAssert(QRx.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        // Solves Q R P x = m in place (m <- x)
        Q_LDivEq(QRx,beta,m);
        m.rowRange(N1,m.colsize()).setZero();
        //m.rowRange(0,N1) /= QRx.upperTri().subTriMatrix(0,N1);
        QRx.upperTri().subTriMatrix(0,N1).LDivEq(m.rowRange(0,N1));
        if (P) m.reversePermuteRows(P);
    }


    //
    // RDiv
    //

    template <class T, class T1, class T2> 
    void QR_RDiv(
        const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const ptrdiff_t* P,
        const GenMatrix<T2>& m, MatrixView<T> x, ptrdiff_t N1)
    {
        TMVAssert(QRx.colsize() >= QRx.rowsize());
        TMVAssert(beta.size() == QRx.rowsize());
        TMVAssert(x.rowsize() == QRx.colsize());
        TMVAssert(m.rowsize() == QRx.rowsize());
        TMVAssert(x.colsize() == m.colsize());
        TMVAssert(QRx.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        // Solve x Q R P = m
        // where Q and R are stored in QRx, and beta are the beta
        //std::cout<<"Start QR_RDiv"<<std::endl;
        //std::cout<<"Norm(QRx) = "<<Norm(QRx)<<std::endl;
        //std::cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
        //std::cout<<"Norm(m) = "<<Norm(m)<<std::endl;

        // First solve y P = m
        x.colRange(0,m.rowsize()) = m;
        if (P) x.colRange(0,m.rowsize()).permuteCols(P);

        // Next solve z R = y by forward substitution
        x.colRange(N1,x.rowsize()).setZero();
        //std::cout<<"1 Norm(x) = "<<Norm(x)<<std::endl;

        //x.colRange(0,N1) %= QRx.upperTri().subTriMatrix(0,N1);
        QRx.upperTri().subTriMatrix(0,N1).RDivEq(x.colRange(0,N1));
        //std::cout<<"2 Norm(x) = "<<Norm(x)<<std::endl;

        // Finally solve x Q = z
        // Q = Q1 [ I ]
        //        [ 0 ]
        // where Q1 is the part of Q that is stored in QRx and beta
        // We've already dealt with the first part by zeroing out the 
        // right columns of x.
        Q_RDivEq(QRx,beta,x);
        //std::cout<<"3 Norm(x) = "<<Norm(x)<<std::endl;
    }

    //
    // RDivEq
    //

    template <class T, class T1> 
    void QR_RDivEq(
        const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const ptrdiff_t* P,
        MatrixView<T> m, ptrdiff_t N1)
    {
        TMVAssert(QRx.colsize() == QRx.rowsize());
        TMVAssert(beta.size() == QRx.rowsize());
        TMVAssert(m.rowsize() == QRx.colsize());
        TMVAssert(QRx.isrm() || QRx.iscm());
        TMVAssert(QRx.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        // Solve x Q R P = m in place (m <- x)

        if (P) m.permuteCols(P);
        m.colRange(N1,m.rowsize()).setZero();
        //m.colRange(0,N1) %= QRx.upperTri().subTriMatrix(0,N1);
        QRx.upperTri().subTriMatrix(0,N1).RDivEq(m.colRange(0,N1));
        Q_RDivEq(QRx,beta,m);
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_QRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


