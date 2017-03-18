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


//#define XDEBUG


#include "TMV_BandQRDiv.h"
#include "tmv/TMV_BandMatrix.h"
#include "TMV_Householder.h"
#include "TMV_BandLUDiv.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // Packed BandQ - LDivEq/RDivEq
    //

    template <class T, class T1> 
    void Q_LDivEq(
        const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == Qbeta.size());
        TMVAssert(Q.colsize() == m.colsize());

        if (Q.nlo() > 0) {
            ptrdiff_t i2 = Q.nlo()+1;
            const ptrdiff_t M = Q.colsize();
            const ptrdiff_t N = Q.rowsize();
            for(ptrdiff_t j=0,i1=1;j<N;++j,++i1) {
                if (Qbeta(j) != T1(0)) 
                    HouseholderLMult(Q.col(j,i1,i2),Qbeta(j),m.rowRange(j,i2));
                if (i2<M) ++i2;
            }
        }
    }

    template <class T, class T1> 
    void Q_RDivEq(
        const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Qbeta.size() == Q.rowsize());
        TMVAssert(m.rowsize() == Q.colsize());

        if (Q.nlo() > 0) {
            const ptrdiff_t M = Q.colsize();
            const ptrdiff_t N = Q.rowsize();
            ptrdiff_t i1 = N;
            ptrdiff_t i2 = Q.isSquare() ? N : TMV_MIN(N+Q.nlo(),M);
            ptrdiff_t k=Q.isSquare() ? Q.nlo() : N+Q.nlo()-i2;
            for(ptrdiff_t j=N-1;i1>0;--j,--i1) {
                if (Qbeta(j) != T1(0)) 
                    HouseholderLMult(Q.col(j,i1,i2).conjugate(),Qbeta(j),
                                      m.colRange(j,i2).transpose());
                if (k>0) --k; else --i2;
            }
        }
    }

    //
    // LDiv
    //

    template <class T, class T1, class T2> 
    void QR_LDiv(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta,
        const GenMatrix<T2>& m, const MatrixView<T>& x)
    {
        TMVAssert(QRx.colsize() >= QRx.rowsize());
        TMVAssert(Qbeta.size() == QRx.rowsize());
        TMVAssert(m.colsize() == QRx.colsize());
        TMVAssert(x.colsize() == QRx.rowsize());
        TMVAssert(x.rowsize() == m.rowsize());

#ifdef XDEBUG
        cout<<"Start QR_LDiv\n";
        cout<<"QRx = "<<QRx<<endl;
        cout<<"Qbeta = "<<Qbeta<<endl;
        cout<<"m = "<<m<<endl;
        //cout<<"x = "<<x<<endl;
        Matrix<T1> QR2 = QRx;
        GetQFromBandQR(QR2.view(),Qbeta,QRx.nlo());
        QR2 *= QRx.diagRange(0,QRx.nhi()+1);
        QR2.divideUsing(QR);
        Matrix<T> x2 = m / QR2;
#endif
        const ptrdiff_t N = QRx.rowsize();

        if (QRx.isSquare()) {
            x = m;
            Q_LDivEq(QRx,Qbeta,x);
        } else if (QRx.nlo() > 0) {
            if (m.isrm()) {
                Matrix<T,RowMajor> m1 = m;
                Q_LDivEq(QRx,Qbeta,m1.view());
                x = m1.rowRange(0,N);
            } else {
                Matrix<T,ColMajor> m1 = m;
                Q_LDivEq(QRx,Qbeta,m1.view());
                x = m1.rowRange(0,N);
            }
        } else {
            x = m.rowRange(0,N);
        }

        TriLDivEq(QRx.subBandMatrix(0,N,0,N,0,QRx.nhi()),x,NonUnitDiag);

#ifdef XDEBUG
        cout<<"x => "<<x<<endl;
        cout<<"Norm(x2-x) = "<<Norm(x2-x)<<endl;
        if (!(Norm(x2-x) <= 0.001*Norm(x))) {
            cerr<<"QR_LDiv: \n";
            cerr<<"m = "<<TMV_Text(m)<<"  "<<m<<endl;
            cerr<<"x = "<<TMV_Text(x)<<endl;
            cerr<<"-> x = "<<x<<endl;
            cerr<<"x2 = "<<x2<<endl;
            cerr<<"QR = "<<QR2<<endl;
            cerr<<"QR x = "<<QR2*x<<endl;
            cerr<<"QR x2 = "<<QR2*x2<<endl;
            abort();
        }
#endif
    }

    //
    // LDivEq
    //

    template <class T, class T1> 
    void QR_LDivEq(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
        const MatrixView<T>& m)
    {
        TMVAssert(QRx.colsize() == QRx.rowsize());
        TMVAssert(Qbeta.size() == QRx.rowsize());
        TMVAssert(m.colsize() == QRx.colsize());

#ifdef XDEBUG
        cout<<"Start QR_LDivEq\n";
        cout<<"QRx = "<<QRx<<endl;
        cout<<"Qbeta = "<<Qbeta<<endl;
        cout<<"m = "<<m<<endl;
        Matrix<T> m0(m);
        Matrix<T1> QR2 = QRx;
        GetQFromBandQR(QR2.view(),Qbeta,QRx.nlo());
        QR2 *= QRx.diagRange(0,QRx.nhi()+1);
        QR2.divideUsing(QR);
        Matrix<T> m2 = m / QR2;
#endif
        const ptrdiff_t N = QRx.rowsize();

        Q_LDivEq(QRx,Qbeta,m);
        TriLDivEq(QRx.subBandMatrix(0,N,0,N,0,QRx.nhi()),m,NonUnitDiag);

#ifdef XDEBUG
        cout<<"m => "<<m<<endl;
        cout<<"Norm(m2-m) = "<<Norm(m2-m)<<endl;
        if (!(Norm(m2-m) <= 0.001*Norm(m))) {
            cerr<<"QR_LDivEq: \n";
            cerr<<"m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"-> m = "<<m<<endl;
            cerr<<"m2 = "<<m2<<endl;
            cerr<<"QR = "<<QR2<<endl;
            cerr<<"QR m = "<<QR2*m<<endl;
            cerr<<"QR m2 = "<<QR2*m2<<endl;
            abort();
        }
#endif
    }

    //
    // RDiv
    //

    template <class T, class T1, class T2> 
    void QR_RDiv(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
        const GenMatrix<T2>& m, const MatrixView<T>& x)
    {
        TMVAssert(QRx.colsize() >= QRx.rowsize());
        TMVAssert(Qbeta.size() == QRx.rowsize());
        TMVAssert(x.rowsize() == QRx.colsize());
        TMVAssert(m.rowsize() == QRx.rowsize());
        TMVAssert(x.colsize() == m.colsize());

#ifdef XDEBUG
        cout<<"Start QR_RDiv\n";
        cout<<"QRx = "<<QRx<<endl;
        cout<<"Qbeta = "<<Qbeta<<endl;
        cout<<"m = "<<m<<endl;
        //cout<<"x = "<<x<<endl;
        Matrix<T1> QR2 = QRx;
        GetQFromBandQR(QR2.view(),Qbeta,QRx.nlo());
        QR2 *= QRx.diagRange(0,QRx.nhi()+1);
        QR2.divideUsing(QR);
        Matrix<T> x2 = m % QR2;
#endif
        const ptrdiff_t M = QRx.colsize();
        const ptrdiff_t N = QRx.rowsize();

        x.colRange(N,M).setZero();
        x.colRange(0,N) = m;
        TriLDivEq(QRx.subBandMatrix(0,N,0,N,0,QRx.nhi()).transpose(),
                  x.colRange(0,N).transpose(),NonUnitDiag);

        Q_RDivEq(QRx,Qbeta,x);

#ifdef XDEBUG
        cout<<"x => "<<x<<endl;
        cout<<"Norm(x2-x) = "<<Norm(x2-x)<<endl;
        if (!(Norm(x2-x) <= 0.001*Norm(x))) {
            cerr<<"QR_RDiv: \n";
            cerr<<"m = "<<TMV_Text(m)<<"  "<<m<<endl;
            cerr<<"x = "<<TMV_Text(x)<<endl;
            cerr<<"-> x = "<<x<<endl;
            cerr<<"x2 = "<<x2<<endl;
            cerr<<"QR = "<<QR2<<endl;
            cerr<<"x QR = "<<x*QR2<<endl;
            cerr<<"x2 QR = "<<x2*QR2<<endl;
            abort();
        }
#endif
    }

    //
    // RDivEq
    //

    template <class T, class T1> 
    void QR_RDivEq(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
        const MatrixView<T>& m)
    {
        TMVAssert(QRx.colsize() == QRx.rowsize());
        TMVAssert(Qbeta.size() == QRx.rowsize());
        TMVAssert(m.rowsize() == QRx.colsize());

#ifdef XDEBUG
        cout<<"Start QR_RDivEq\n";
        cout<<"QRx = "<<QRx<<endl;
        cout<<"Qbeta = "<<Qbeta<<endl;
        cout<<"m = "<<m<<endl;
        Matrix<T> m0(m);
        Matrix<T1> QR2 = QRx;
        GetQFromBandQR(QR2.view(),Qbeta,QRx.nlo());
        QR2 *= QRx.diagRange(0,QRx.nhi()+1);
        QR2.divideUsing(QR);
        Matrix<T> m2 = m % QR2;
#endif

        const ptrdiff_t N = QRx.rowsize();

        // Solve x Q R = m in place (m <- x)
        TriLDivEq(QRx.subBandMatrix(0,N,0,N,0,QRx.nhi()).transpose(),
                  m.transpose(),NonUnitDiag);
        Q_RDivEq(QRx,Qbeta,m);

#ifdef XDEBUG
        cout<<"m => "<<m<<endl;
        cout<<"Norm(m2-m) = "<<Norm(m2-m)<<endl;
        if (!(Norm(m2-m) <= 0.001*Norm(m))) {
            cerr<<"QR_RDivEq: \n";
            cerr<<"m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"-> m = "<<m<<endl;
            cerr<<"m2 = "<<m2<<endl;
            cerr<<"QR = "<<QR2<<endl;
            cerr<<"m QR = "<<m*QR2<<endl;
            cerr<<"m2 QR = "<<m2*QR2<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandQRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


