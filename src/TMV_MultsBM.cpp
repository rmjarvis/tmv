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


#include "tmv/TMV_SymBandMatrixArithFunc.h"
#include "tmv/TMV_SymBandMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_BandMatrixArith.h"

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE TMV_BLOCKSIZE
#else
#define SYM_MM_BLOCKSIZE 64
#endif

    //
    // MultMM
    //

    template <bool add, class T, class Ta, class Tb> 
    static void DoMultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));

        if (add) C += alpha * A.lowerBand() * B;
        else C = alpha * A.lowerBand() * B;

        const ptrdiff_t N = A.size();
        if (N > 1 && A.nlo() > 0) {
            C.rowRange(0,N-1) += A.upperBandOff() * B.rowRange(1,N);
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void FullTempMultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        if (C.isrm()) {
            Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
            DoMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        } else {
            Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
            DoMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void BlockTempMultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        const ptrdiff_t N = C.rowsize();
        for(ptrdiff_t j=0;j<N;) {
            ptrdiff_t j2 = TMV_MIN(N,j+SYM_MM_BLOCKSIZE);
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                if (C.isrm()) {
                    Matrix<Tb,RowMajor> B2 = TMV_REAL(alpha) * B.colRange(j,j2);
                    DoMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                } else {
                    Matrix<Tb,ColMajor> B2 = TMV_REAL(alpha) * B.colRange(j,j2);
                    DoMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                }
            } else {
                if (C.isrm()) {
                    Matrix<T,RowMajor> B2 = alpha * B.colRange(j,j2);
                    DoMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                } else {
                    Matrix<T,ColMajor> B2 = alpha * B.colRange(j,j2);
                    DoMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                }
            }
            j = j2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
        // C (+)= alpha * A * B
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
        //cout<<"Start MultMM: alpha = "<<alpha<<endl;
        //cout<<"A = "<<A.cptr()<<"  "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"B = "<<B.cptr()<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<"C = "<<C.cptr()<<"  "<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C0 = C;
        Matrix<T> C2 = alpha*A0*B0;
        if (add) C2 += C0;
#endif

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (alpha == T(0)) {
                if (!add) C.setZero();
            } else if (C.isconj()) {
                MultMM<add>(
                    TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
            } else if (SameStorage(A,C))  {
                FullTempMultMM<add>(alpha,A,B,C);
            } else if (SameStorage(B,C))  {
                if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
                    BlockTempMultMM<add>(alpha,A,B,C);
                else
                    FullTempMultMM<add>(alpha,A,B,C);
            } else {
                BlockTempMultMM<add>(alpha, A, B, C);
            }
        }

#ifdef XDEBUG
        //cout<<"Done: C = "<<C<<endl;
        if (!(Norm(C-C2) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     (add?Norm(C0):TMV_RealType(T)(0))))) {
            cerr<<"MultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultsBM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


