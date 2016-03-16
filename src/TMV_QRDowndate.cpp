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


//#define XDEBUG


#include "tmv/TMV_QRD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_Householder.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef NOTHROW
#include <iostream>
#endif
#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define QR_BLOCKSIZE TMV_BLOCKSIZE
#else
#define QR_BLOCKSIZE 64
#endif

    //
    // QR Downdate
    //

#ifndef NOTHROW
    template <class T> 
    class BadQRDowndate : public NonPosDef
    {
    public:
        UpperTriMatrix<T> R;
        Matrix<T> A;

        BadQRDowndate(
            const GenUpperTriMatrix<T>& _R, const GenMatrix<T>& _A) :
            NonPosDef("QR Downdate."), R(_R), A(_A) {}
        BadQRDowndate(const BadQRDowndate<T>& rhs) :
            R(rhs.R), A(rhs.A) {}
        ~BadQRDowndate() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV NonPosDef: QR Downdate found that the resulting "<<std::endl;
            os<<"down-dated RtR is not positive definite. "<<std::endl;
            os<<"(and hence the down date is impossible)"<<std::endl;
            os<<"The partially downdated matrix is \n"<<R<<std::endl;
            os<<"The matrix attempting to be down-dated was \n"<<A<<std::endl;
        }
    };
#endif

    template <class T> 
    static void NonBlockQRDowndate(UpperTriMatrixView<T> R, MatrixView<T> A)
    {
        TMVAssert(A.rowsize() == R.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(R.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(R.dt() == NonUnitDiag);
        // Given that A0 = Q0 R0
        // Given that [ A0 ] = Q1 R1
        //            [ A  ] 
        // Find R0 so that A0 = Q0 R0
        // Input R is R1, output is R0

        const ptrdiff_t N = A.rowsize();

        T* Rdiag = R.ptr();
        const ptrdiff_t ds = R.stepi()+R.stepj();

        for(ptrdiff_t j=0;j<N;++j,Rdiag+=ds) {
            // Apply the Householder Reflection for this column
            VectorView<T> v = A.col(j);
            T beta;
            if (!(HouseholderUnReflect(*Rdiag,v,beta))) {
#ifdef NOTHROW
                std::cerr<<"Bad QR Downdate\n";
                exit(1); 
#else
                throw BadQRDowndate<T>(R,A);
#endif
            }
            TMVAssert(beta != T(1));
            VectorView<T> m0 = R.row(j,j+1,N);
            MatrixView<T> mx = A.colRange(j+1,N);

            // m0' = m0 - beta m0 - beta vtmx
            // m0 = (m0' + beta btmx)/(1-beta)
            Vector<T> bvtmx = beta*v.conjugate()*mx;
            m0 += bvtmx;
            m0 /= T(1)-beta;

            // mx' = mx - beta v (m0 + vtmx)
            bvtmx += beta*m0;
            mx -= v ^ bvtmx;
        }
    }

    template <class T> 
    static void RecursiveQRDowndate(
        UpperTriMatrixView<T> R, MatrixView<T> A,
        UpperTriMatrixView<T> Z, bool makeZ)
    {
        TMVAssert(A.rowsize() == R.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(R.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(R.dt() == NonUnitDiag);

        const ptrdiff_t N = A.rowsize();

        if (N==1) {
            T b;
            if (!(HouseholderUnReflect(*R.ptr(),A.col(0),b))) {
#ifdef NOTHROW
                std::cerr<<"Bad QR Downdate\n";
                exit(1); 
#else
                throw BadQRDowndate<T>(R,A);

#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(Z.ptr() >= Z._first);
            TMVAssert(Z.ptr() < Z._last);
#endif
            *Z.ptr() = TMV_CONJ(b);
        } else if (N==2) {
            T* R00 = R.ptr();
            T* R01 = R00+R.stepj();
            T* R11 = R01+R.stepi();
            T* Z00 = Z.ptr();
            T* Z01 = Z00+Z.stepj();
            T* Z11 = Z01+Z.stepi();
            T b0;
            if (!(HouseholderUnReflect(*R00,A.col(0),b0))) {
#ifdef NOTHROW
                std::cerr<<"Bad QR Downdate\n"; 
                exit(1); 
#else
                throw BadQRDowndate<T>(R,A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(Z00 >= Z._first);
            TMVAssert(Z00 < Z._last);
#endif
            *Z00 = TMV_CONJ(b0);
            if (b0 != T(0)) {
                TMVAssert(b0 != T(1));
                T vtmx = A.col(0).conjugate() * A.col(1);
#ifdef TMVFLDEBUG
                TMVAssert(R01 >= R._first);
                TMVAssert(R01 < R._last);
#endif
                *R01 = (*R01 + b0*vtmx)/(T(1)-b0);
                A.col(1) -= b0*(vtmx + *R01) * A.col(0);
            }

            T b1;
            if (!(HouseholderUnReflect(*R11,A.col(1),b1))) {
#ifdef NOTHROW
                std::cerr<<"Bad QR Downdate\n";
                exit(1); 
#else
                throw BadQRDowndate<T>(R,A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(Z11 >= Z._first);
            TMVAssert(Z11 < Z._last);
#endif
            *Z11 = TMV_CONJ(b1);

            if (makeZ) {
                T vtmx = A.col(0).conjugate() * A.col(1);
#ifdef TMVFLDEBUG
                TMVAssert(Z01 >= Z._first);
                TMVAssert(Z01 < Z._last);
#endif
                *Z01 = -TMV_CONJ(b0*b1)*vtmx;
            }
        } else {
            ptrdiff_t j1 = N/2;

            UpperTriMatrixView<T> R1 = R.subTriMatrix(0,j1);
            MatrixView<T> Rx = R.subMatrix(0,j1,j1,N);
            UpperTriMatrixView<T> R2 = R.subTriMatrix(j1,N);

            MatrixView<T> A1 = A.colRange(0,j1);
            MatrixView<T> A2 = A.colRange(j1,N);

            UpperTriMatrixView<T> Z1 = Z.subTriMatrix(0,j1);
            MatrixView<T> Zx = Z.subMatrix(0,j1,j1,N);
            UpperTriMatrixView<T> Z2 = Z.subTriMatrix(j1,N);

#ifndef NOTHROW
            try {
#endif
                RecursiveQRDowndate(R1,A1,Z1,true);
#ifndef NOTHROW
            } catch (BadQRDowndate<T>) {
                throw BadQRDowndate<T>(R,A);
            }
#endif

            Zx = A1.adjoint() * A2;
            Zx = Z1.adjoint() * Zx;

            Rx += Zx;
            LowerTriMatrix<T> ImZt = T(1)-Z1.adjoint();
            Rx /= ImZt;

            Zx += Z1.adjoint() * Rx;
            A2 -= A1 * Zx;

#ifndef NOTHROW
            try {
#endif
                RecursiveQRDowndate(R2,A2,Z2,makeZ);
#ifndef NOTHROW
            } catch (BadQRDowndate<T>) {
                throw BadQRDowndate<T>(R,A);
            }
#endif

            if (makeZ) {
                Zx = A1.adjoint() * A2; 
                Zx = -Z1*Zx;
                Zx *= Z2;
            }
        }
    }

    template <class T> 
    static void BlockQRDowndate(UpperTriMatrixView<T> R, MatrixView<T> A)
    {
        TMVAssert(A.rowsize() == R.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(R.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(R.dt() == NonUnitDiag);

        const ptrdiff_t N = A.rowsize();

        UpperTriMatrix<T,NonUnitDiag|ColMajor> BaseZ(
            TMV_MIN(QR_BLOCKSIZE,int(N)));
        for(ptrdiff_t j1=0;j1<N;) {
            ptrdiff_t j2 = TMV_MIN(N,j1+QR_BLOCKSIZE);
            MatrixView<T> A1 = A.colRange(j1,j2);
            UpperTriMatrixView<T> R1 = R.subTriMatrix(j1,j2);
            UpperTriMatrixView<T> Z = BaseZ.subTriMatrix(0,j2-j1);

#ifndef NOTHROW
            try {
#endif
                RecursiveQRDowndate(R1,A1,Z,j2<N);
#ifndef NOTHROW
            } catch (BadQRDowndate<T>) {
                throw BadQRDowndate<T>(R,A);
            }
#endif

            if (j2 < N) {
                // m0' = m0 - Zt(Ytmx+m0)
                // m0' + ZtYtmx = (I-Zt) m0;
                MatrixView<T> m0 = R.subMatrix(j1,j2,j2,N);
                MatrixView<T> Y = A.colRange(j1,j2);
                MatrixView<T> mx = A.colRange(j2,N);

                Matrix<T,ColMajor> ZtYtm = Y.adjoint() * mx;
                ZtYtm = Z.adjoint() * ZtYtm;

                m0 += ZtYtm;
                LowerTriMatrix<T> ImZt = T(1)-Z.adjoint();
                m0 /= ImZt;

                ZtYtm += Z.adjoint() * m0;
                mx -= Y * ZtYtm;
            }
            j1 = j2;
        }
    }

    template <class T> 
    void QR_Downdate(UpperTriMatrixView<T> R, MatrixView<T> A)
    {
        TMVAssert(A.rowsize() == R.size());
        TMVAssert(R.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(R.dt() == NonUnitDiag);

#ifdef XDEBUG
        Matrix<T> R0(R);
        Matrix<T> A0(A);
        UpperTriMatrix<T> R2(R);
        Matrix<T> A2(A);
        NonBlockQRDowndate(R2.view(),A2.view());
#endif

        if (A.rowsize() > 0) {
            if (A.rowsize() > QR_BLOCKSIZE)
                BlockQRDowndate(R,A);
            else {
                UpperTriMatrix<T,NonUnitDiag|ColMajor> Z(A.rowsize());
                RecursiveQRDowndate(R,A,Z.view(),false);
            }
        }

#ifdef XDEBUG
        if (Norm(R2-R) > 1.e-5*Norm(A0)*Norm(R0)) {
            cerr<<"QR_Downdate\n";
            cerr<<"R0 = "<<TMV_Text(R)<<"  "<<R0<<endl;
            cerr<<"A0 = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"R -> "<<R<<endl;
            cerr<<"NonBlock R -> "<<R2<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_QRDowndate.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


