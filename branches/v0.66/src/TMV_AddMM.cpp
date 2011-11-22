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


#include "tmv/TMV_MatrixArithFunc.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_VIt.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // AddMM
    //

    template <bool rm, bool a1, bool ca, class T, class T1, class T2> 
    static void DoRowAddMM(const T1 alpha, const GenMatrix<T2>& A, 
                           const MatrixView<T>& B)
    {
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(alpha != T1(0));
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.ct() == NonConj);
        TMVAssert(!SameStorage(A,B));
        TMVAssert(rm == (A.isrm() && B.isrm()));
        TMVAssert(a1 == (alpha == T1(1)));
        TMVAssert(ca == A.isconj());

        const T2* Arowi = A.cptr();
        T* Browi = B.ptr();
        const int M = A.colsize();
        const int N = A.rowsize();
        const int Asi = A.stepi();
        const int Asj = (rm ? 1 : A.stepj());
        const int Bsi = B.stepi();
        const int Bsj = (rm ? 1 : B.stepj());

        for(int i=M;i>0;--i,Arowi+=Asi,Browi+=Bsi) {
            const T2* Aij = Arowi;
            T* Bij = Browi;
            for(int j=N;j>0;--j,(rm?++Aij:Aij+=Asj),(rm?++Bij:Bij+=Bsj)) {
#ifdef TMVFLDEBUG
                TMVAssert(Bij >= B._first);
                TMVAssert(Bij < B._last);
#endif
                if (a1) *Bij += (ca ? TMV_CONJ(*Aij) : *Aij);
                else *Bij += alpha * (ca ? TMV_CONJ(*Aij) : *Aij);
            }
        }
    }

    template <bool rm, class T, class Ta> static void rowAddMM(
        const T alpha, const GenMatrix<Ta>& A, const MatrixView<T>& B)
    { 
        if (TMV_IMAG(alpha) == TMV_RealType(T)(0))
            if (TMV_REAL(alpha) == TMV_RealType(T)(1))
                if (A.isconj()) DoRowAddMM<rm,true,true>(TMV_REAL(alpha),A,B);
                else DoRowAddMM<rm,true,false>(TMV_REAL(alpha),A,B);
            else
                if (A.isconj()) DoRowAddMM<rm,false,true>(TMV_REAL(alpha),A,B);
                else DoRowAddMM<rm,false,false>(TMV_REAL(alpha),A,B);
        else
            if (A.isconj()) DoRowAddMM<rm,false,true>(alpha,A,B);
            else DoRowAddMM<rm,false,false>(alpha,A,B);
    }

    template <class T, class Ta> static void DoAddMM(
        const T alpha, const GenMatrix<Ta>& A, const MatrixView<T>& B)
    { 
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.ct() == NonConj);
        TMVAssert(!SameStorage(A,B));

        if (A.stor() == B.stor() && A.canLinearize() && B.canLinearize()) {
            TMVAssert(A.stepi() == B.stepi() && A.stepj() == B.stepj());
            B.linearView() += alpha * A.constLinearView();
        } else {
            if (A.isrm() && B.isrm())
                rowAddMM<true>(alpha,A,B); 
            else if (A.iscm() && B.iscm())
                rowAddMM<true>(alpha,A.transpose(),B.transpose()); 
            else if (A.rowsize() > A.colsize())
                rowAddMM<false>(alpha,A,B); 
            else
                rowAddMM<false>(alpha,A.transpose(),B.transpose()); 
        }
    }

    template <class T, class Ta> void AddMM(const T alpha,
              const GenMatrix<Ta>& A, const MatrixView<T>& B)
        // B += alpha * A 
    {
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> B2 = B;
        for(int i=0;i<int(A.colsize());i++)
            for(int j=0;j<int(A.rowsize());j++)
                B2(i,j) += alpha*A(i,j);
        //cout<<"AddMM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<", B = "<<TMV_Text(B)<<"  "<<B<<endl;
#endif

        if (alpha != T(0) && B.colsize() > 0 && B.rowsize() > 0) {
            if (B.isconj()) 
                AddMM(TMV_CONJ(alpha),A.conjugate(),B.conjugate());
            else {
                if (SameStorage(A,B)) {
                    if (B.isrm()) {
                        Matrix<T,RowMajor> A2 = A;
                        DoAddMM(alpha,A2,B);
                    } else {
                        Matrix<T,ColMajor> A2 = A;
                        DoAddMM(alpha,A2,B);
                    }
                } 
                else DoAddMM(alpha,A,B);
            }
        }
#ifdef XDEBUG
        //cout<<"done: B = "<<B<<endl;
        Matrix<T> diff(B.colsize(),B.rowsize());
        for(int i=0;i<int(B.colsize());i++)
            for(int j=0;j<int(B.rowsize());j++)
                diff(i,j) = B(i,j) - B2(i,j);
        if (Norm(diff) > 0.001*(TMV_ABS(alpha)*Norm(A0)+Norm(B0))) {
            cerr<<"AddMM: alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"->B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta, class Tb> void AddMM(const T alpha,
              const GenMatrix<Ta>& A, const T beta, const GenMatrix<Tb>& B,
              const MatrixView<T>& C)
        // C = alpha * A + beta * B
    {
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == C.rowsize());
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> C2(A.colsize(),A.rowsize());
        for(int i=0;i<int(A.colsize());i++)
            for(int j=0;j<int(A.rowsize());j++)
                C2(i,j) = alpha*A(i,j) + beta*B(i,j);
        //cout<<"AddMM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"beta = "<<beta<<", B = "<<TMV_Text(B)<<"  "<<B;
        //cout<<", C = "<<TMV_Text(C)<<"  "<<C<<endl;
#endif

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (SameStorage(A,C)) {
                if (SameStorage(B,C)) {
                    if (A.isrm()) {
                        Matrix<Ta,RowMajor> tempA = A;
                        C = B;
                        C *= beta;
                        AddMM(alpha,tempA,C);
                    } else {
                        Matrix<Ta,ColMajor> tempA = A;
                        C = B;
                        C *= beta;
                        AddMM(alpha,tempA,C);
                    }
                } else {
                    C = A;
                    C *= alpha;
                    AddMM(beta,B,C);
                }
            } else {
                C = B;
                C *= beta;
                AddMM(alpha,A,C);
            }
        }

#ifdef XDEBUG
        //cout<<"Done: C = "<<C<<endl;
        Matrix<T> diff(C.colsize(),C.rowsize());
        for(int i=0;i<int(C.colsize());i++)
            for(int j=0;j<int(C.rowsize());j++)
                diff(i,j) = C(i,j) - C2(i,j);
        if (Norm(diff) >
            0.001*(1.+TMV_ABS(alpha)*Norm(A0)+TMV_ABS(beta)*Norm(B0))) {
            cerr<<"AddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"->C = "<<TMV_Text(C)<<"  "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_AddMM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


