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
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Householder.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

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
    // QR Update
    //

    template <class T> 
    static void NonBlockQRUpdate(UpperTriMatrixView<T> R, MatrixView<T> A)
    {
        TMVAssert(A.rowsize() == R.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(R.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(R.dt() == NonUnitDiag);
        // Given that A0 = Q0 R0
        // Find R1, so that [ A0 ] = Q1 R1
        //                  [ A  ] 
        // Input R is R0, output is R1

        const ptrdiff_t N = A.rowsize();

        T* Rdiag = R.ptr();
        const ptrdiff_t ds = R.stepi()+R.stepj();
        T det(0);

        for(ptrdiff_t j=0;j<N;++j,Rdiag+=ds) {
            // Apply the Householder Reflection for this column
            VectorView<T> v = A.col(j);
            T beta = HouseholderReflect(*Rdiag,v,det);
            if (beta != T(0))
                HouseholderLMult(v,beta,R.row(j,j+1,N),A.colRange(j+1,N));
        }
    }

    template <class T> 
    static void RecursiveQRUpdate(
        UpperTriMatrixView<T> R, MatrixView<T> A,
        UpperTriMatrixView<T> Z, bool makeZ)
    {
        TMVAssert(A.rowsize() == R.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(R.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(R.dt() == NonUnitDiag);

        const ptrdiff_t N = A.rowsize();
        T det(0);

        TMVAssert(!R.isconj());
        TMVAssert(!Z.isconj());

        if (N==1) {
            T b = HouseholderReflect(*R.ptr(),A.col(0),det);
#ifdef TMVFLDEBUG
            TMVAssert(Z.ptr() >= Z._first);
            TMVAssert(Z.ptr() < Z._last);
#endif
            *Z.ptr() = TMV_CONJ(b);
        } else if (N==2) {
            T* R00 = R.ptr(); // = R(0,0)
            T* R01 = R00 + R.stepj(); // = R(0,1)
            T* R11 = R01 + R.stepi(); // = R(1,1)
            T* Z00 = Z.ptr(); // = Z(0,0)
            T* Z01 = Z00 + Z.stepj(); // = Z(0,1)
            T* Z11 = Z01 + Z.stepi(); // = Z(1,1)
            T b0 = HouseholderReflect(*R00,A.col(0),det);
            if (b0 != T(0)) {
                T temp = b0*(A.col(0).conjugate()*A.col(1) + *R01);
#ifdef TMVFLDEBUG
                TMVAssert(R01 >= R._first);
                TMVAssert(R01 < R._last);
#endif
                *R01 -= temp;
                A.col(1) -= temp * A.col(0);
            }
#ifdef TMVFLDEBUG
            TMVAssert(Z00 >= Z._first);
            TMVAssert(Z00 < Z._last);
            TMVAssert(Z11 >= Z._first);
            TMVAssert(Z11 < Z._last);
#endif
            *Z00 = TMV_CONJ(b0);
            T b1 = HouseholderReflect(*R11,A.col(1),det);
            *Z11 = TMV_CONJ(b1);

            if (makeZ) {
                T temp = A.col(0).conjugate()*A.col(1);
#ifdef TMVFLDEBUG
                TMVAssert(Z01 >= Z._first);
                TMVAssert(Z01 < Z._last);
#endif
                *Z01 = -(*Z00 * *Z11)*temp;
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

            RecursiveQRUpdate(R1,A1,Z1,true);

            // Zx is a temporary here - it happens to be the right shape.
            Zx = A1.adjoint() * A2; 
            Zx += Rx;
            Zx = Z1.adjoint()*Zx;
            Rx -= Zx;
            A2 -= A1 * Zx;

            RecursiveQRUpdate(R2,A2,Z2,makeZ);

            if (makeZ) {
                Zx = A1.adjoint() * A2; 
                Zx = -Z1*Zx;
                Zx *= Z2;
            }
        }
    }

    template <class T> 
    static void BlockQRUpdate(UpperTriMatrixView<T> R, MatrixView<T> A)
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

            RecursiveQRUpdate(R1,A1,Z,j2<N);

            if (j2 < N) {
                Matrix<T,ColMajor> ZtYtm =
                    A.colRange(j1,j2).adjoint() * A.colRange(j2,N);
                ZtYtm += R.subMatrix(j1,j2,j2,N);
                ZtYtm = Z.adjoint() * ZtYtm;
                R.subMatrix(j1,j2,j2,N) -= ZtYtm;
                A.colRange(j2,N) -= A.colRange(j1,j2) * ZtYtm;
            }
            j1 = j2;
        }
    }

    template <class T> 
    void QR_Update(UpperTriMatrixView<T> R, MatrixView<T> A)
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
        NonBlockQRUpdate(R2.view(),A2.view());
#endif
        if (A.rowsize() > 0) {
            if (A.rowsize() > QR_BLOCKSIZE)
                BlockQRUpdate(R,A);
            else {
                UpperTriMatrix<T,NonUnitDiag|ColMajor> Z(A.rowsize());
                RecursiveQRUpdate(R,A,Z.view(),false);
            }
        }
#ifdef XDEBUG
        if (Norm(R2-R) > 1.e-5*Norm(R0)*Norm(A0)) {
            cerr<<"QR_Update\n";
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

#define InstFile "TMV_QRUpdate.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


