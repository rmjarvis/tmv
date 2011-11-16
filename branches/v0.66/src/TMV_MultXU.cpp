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


#include "tmv/TMV_TriMatrixArithFunc.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // MultXM
    //

    template <class T, class T1> 
    static void RowMajorMultXM(const T1 alpha, const UpperTriMatrixView<T>& A)
    {
        TMVAssert(A.isrm());
        TMVAssert(!A.isunit());
        TMVAssert(alpha != T1(1));
        TMVAssert(A.size() > 0);
        TMVAssert(A.ct() == NonConj);

        T* Aii = A.ptr();
        const int ds = A.stepi()+1;
        const int N = A.size();

        for(int len=N;len>0;--len,Aii+=ds) {
            // A.row(i,i,N) *= alpha;
            T* Aij = Aii;
            for(int j=len;j>0;--j,++Aij) {
#ifdef TMVFLDEBUG
                TMVAssert(Aij >= A.first);
                TMVAssert(Aij < A.last);
#endif
                *Aij *= alpha;
            }
        }
    }

    template <class T, class T1> 
    static void ColMajorMultXM(const T1 alpha, const UpperTriMatrixView<T>& A)
    {
        TMVAssert(A.iscm());
        TMVAssert(!A.isunit());
        TMVAssert(alpha != T1(1));
        TMVAssert(A.size() > 0);
        TMVAssert(A.ct() == NonConj);

        T* A0j = A.ptr();
        const int Astepj = A.stepj();
        const int N = A.size();

        for(int j=N,len=1;j>0;--j,++len,A0j+=Astepj) {
            // A.col(j,0,j+1) *= alpha;
            T* Aij = A0j;
            for(int i=len;i>0;--i,++Aij) {
#ifdef TMVFLDEBUG
                TMVAssert(Aij >= A.first);
                TMVAssert(Aij < A.last);
#endif
                *Aij *= alpha;
            }
        }
    }

    template <class T> 
    void MultXM(const T alpha, const UpperTriMatrixView<T>& A)
    // A = alpha * A
    {
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T> A2 = alpha * A0;
        //cout<<"MultXM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<" "<<A<<endl;
#endif

        if (A.size() > 0 && alpha != T(1)) {
            TMVAssert(!A.isunit());
            if (A.isconj()) {
                MultXM(TMV_CONJ(alpha),A.conjugate());
            } else if (alpha == T(0)) {
                A.setZero();
            } else if (A.isrm()) {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0))
                    RowMajorMultXM(TMV_REAL(alpha),A);
                else
                    RowMajorMultXM(alpha,A);
            } else if (A.iscm()) {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0))
                    ColMajorMultXM(TMV_REAL(alpha),A);
                else
                    ColMajorMultXM(alpha,A);
            } else {
                const int M = A.colsize();
                const int N = A.rowsize();
                for(int i=0;i<M;++i) 
                    A.row(i,i,N) *= alpha;
            }
        }
#ifdef XDEBUG
        //cout<<"Done MultXM: A = "<<A<<endl;
        if (!(Norm(Matrix<T>(A)-A2) <= 0.001*TMV_ABS(alpha)*Norm(A))) {
            cerr<<"MultXM: alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> A = "<<A<<endl;
            cerr<<"A2 = "<<A2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta> 
    void ElementProd(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const UpperTriMatrixView<T>& B)
    {
        TMVAssert(A.size() == B.size());
        const int N = B.size();
        if (B.isrm()) {
            for(int i=0;i<N;i++)
                ElementProd(alpha,A.row(i,i,N),B.row(i,i,N));
        } else {
            for(int j=0;j<N;j++)
                ElementProd(alpha,A.col(j,0,j+1),B.col(j,0,j+1));
        }
    }

    template <class T, class Ta, class Tb> 
    void AddElementProd(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<T>& C)
    {
        TMVAssert(A.size() == C.size());
        TMVAssert(B.size() == C.size());
        const int N = C.size();
        if (C.isrm()) {
            for(int i=0;i<N;i++)
                AddElementProd(alpha,A.row(i,i,N),B.row(i,i,N),C.row(i,i,N));
        } else {
            for(int j=0;j<N;j++)
                AddElementProd(alpha,A.col(j,0,j+1),B.col(j,0,j+1),
                               C.col(j,0,j+1));
        }
    }

#define InstFile "TMV_MultXU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

