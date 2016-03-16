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
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_SymBandMatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_MM_BLOCKSIZE 64
#define SYM_MM_BLOCKSIZE2 32
#endif

    //
    // MultMM
    //

    template <bool add, class T, class Ta, class Tb> 
    static void DoMultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenBandMatrix<Tb>& B, BandMatrixView<T> C)
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));

        const ptrdiff_t N = A.size();

        if (add) C += alpha * A.lowerBand() * B;
        else C = alpha * A.lowerBand() * B;

        if (N > 1 && A.nlo() > 0) {
            const ptrdiff_t M = C.rowsize();
            if (B.nlo() > 0) {
                C.subBandMatrix(
                    0,N-1,0,M,
                    (C.nlo()==C.colsize()-1)?(C.nlo()-1):C.nlo(),
                    C.nhi()
                ) += 
                    alpha * A.upperBandOff() * 
                    B.subBandMatrix(
                        1,N,0,M,B.nlo()-1,
                        (B.nhi()==B.rowsize()-1)?B.nhi():(B.nhi()+1));
            }
            else  {
                C.subBandMatrix(
                    0,N-1,1,M,(C.nlo()>=C.colsize()-2)?
                    (C.colsize()-2):(C.nlo()+1),C.nhi()-1
                ) += 
                    alpha * A.upperBandOff() *
                    B.subBandMatrix(1,N,1,M,B.nlo(),B.nhi());
            }
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void TempMultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenBandMatrix<Tb>& B, BandMatrixView<T> C)
    {
        if (C.isrm()) {
            BandMatrix<T,RowMajor> C2(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
            DoMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        } else if (C.iscm()) {
            BandMatrix<T,ColMajor> C2(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
            DoMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        } else {
            BandMatrix<T,DiagMajor> C2(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
            DoMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenBandMatrix<Tb>& B, BandMatrixView<T> C)
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
            } else if (SameStorage(A,C) || SameStorage(B,C)) {
                TempMultMM<add>(alpha,A,B,C);
            } else {
                DoMultMM<add>(alpha, A, B, C);
            }
        }

#ifdef XDEBUG
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
            cerr<<"Norm(diff) = "<<Norm(C-C2)<<std::endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A, 
        const GenSymBandMatrix<Tb>& B, BandMatrixView<T> C)
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
            } else {
                if (A.nlo() > B.nlo()) {
                    if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                        BandMatrix<Tb> B2 = TMV_REAL(alpha)*B;
                        MultMM<add>(T(1),A,B2,C);
                    } else {
                        BandMatrix<T> B2 = alpha*B;
                        MultMM<add>(T(1),A,B2,C);
                    }
                } else {
                    if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                        BandMatrix<Ta> A2 = TMV_REAL(alpha)*A;
                        MultMM<add>(
                            T(1),B.transpose(),A2.transpose(),C.transpose());
                    } else {
                        BandMatrix<T> A2 = alpha*A;
                        MultMM<add>(
                            T(1),B.transpose(),A2.transpose(),C.transpose());
                    }
                }
            }
        }

#ifdef XDEBUG
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

#define InstFile "TMV_MultsBB.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


