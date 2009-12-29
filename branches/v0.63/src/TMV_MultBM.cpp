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


#include "TMV_Blas.h"
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define MM_BLOCKSIZE TMV_BLOCKSIZE
#else
#define MM_BLOCKSIZE 64
#endif

    //
    // MultMM
    //

    template <bool add, class T, class Ta, class Tb> 
    static void RowMultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha!= T(0));
        TMVAssert(A.rowsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.ct()==NonConj);

        int j1=0;
        int k=A.nlo();
        int j2=A.nhi()+1;
        const int M = A.colsize();
        const int N = A.rowsize();

        for(int i=0;i<M; ++i) {
            // C.row(i) (+)= alpha * A.row(i,j1,j2) * B.rowRange(j1,j2);
            MultMV<add>(alpha,B.rowRange(j1,j2).transpose(),A.row(i,j1,j2),C.row(i));
            if (k>0) --k; else ++j1;
            if (j2<N) ++j2;
            else if (j1==N) {
                if (!add) C.rowRange(i+1,M).zero();
                break;
            }
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void OPMultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha!= T(0));
        TMVAssert(A.rowsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.ct()==NonConj);

        int i1=0;
        int k=A.nhi();
        int i2=A.nlo()+1;
        const int M = A.colsize();
        const int N = A.rowsize();

        if (!add) C.zero();
        for(int j=0;j<N;++j) {
            C.rowRange(i1,i2) += alpha * A.col(j,i1,i2) ^ B.row(j);
            if (k>0) --k; else ++i1;
            if (i2<M) ++i2;
            else if (i1==M) break;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void ColMultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha!= T(0));
        TMVAssert(A.rowsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.ct()==NonConj);

        const int N = B.rowsize();
        for(int j=0;j<N;j++)
            // C.col(j) (+)= alpha * A * B.col(j);
            MultMV<add>(alpha,A,B.col(j),C.col(j));
    }

    template <bool add, class T, class Ta, class Tb> 
    static void NonLapTriDiagMultMM(
        const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.nlo() == 1);
        TMVAssert(A.nhi() == 1);
        TMVAssert(A.isdm());

        const int N = A.diag().size();
        const int Nu = A.rowsize()>A.colsize() ? N : N-1;
        const int Nl = A.rowsize()<A.colsize() ? N : N-1;

        const Ta* di = A.cptr();
        const Ta* dui = A.diag(1).cptr();
        const Ta* dli = A.diag(-1).cptr()-1;

        for(int i=0;i<N;++i,++di,++dui,++dli) {
            if (add) C.row(i) += *di*B.row(i);
            else C.row(i) = *di*B.row(i);
            if (i>0) C.row(i) += *dli*B.row(i-1);
            if (i<Nu) C.row(i) += *dui*B.row(i+1);
        }
        if (Nl == N) {
            if (add) C.row(N) += *dli*B.row(N-1);
            else C.row(N) = *dli*B.row(N-1);
        }
    }

#ifdef ELAP
    template <class T, class Ta, class Tb> 
    static inline void LapTriDiagMultMM(
        const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B, 
        const int beta, const MatrixView<T>& C)
    {
        if (beta == 1) NonLapTriDiagMultMM<true>(A,B,C); 
        else NonLapTriDiagMultMM<false>(A,B,C); 
    }
#ifdef INST_DOUBLE
    template <> 
    void LapTriDiagMultMM(
        const GenBandMatrix<double>& A, const GenMatrix<double>& B,
        const int beta, const MatrixView<double>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(B.ct()==NonConj);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.nlo() == 1);
        TMVAssert(A.nhi() == 1);
        TMVAssert(A.isdm());
        TMVAssert(A.isSquare());
        TMVAssert(B.iscm());
        TMVAssert(C.iscm());

        int n = A.colsize();
        int nrhs = B.rowsize();
        double a(1);
        int ldB = B.stepj();
        double xbeta(beta);
        int ldC = C.stepj();
        LAPNAME(dlagtm) (
            LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
                         LAPV(a),LAPP(A.cptr()+A.stepj()),
                         LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
                         LAPP(B.cptr()),LAPV(ldB),LAPV(xbeta),LAPP(C.ptr()),LAPV(ldC) LAP1);
    }
    template <> 
    void LapTriDiagMultMM(
        const GenBandMatrix<std::complex<double> >& A, 
        const GenMatrix<std::complex<double> >& B,
        const int beta, const MatrixView<std::complex<double> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.nlo() == 1);
        TMVAssert(A.nhi() == 1);
        TMVAssert(A.isdm());
        TMVAssert(A.isSquare());
        TMVAssert(B.iscm());
        TMVAssert(C.iscm());
        TMVAssert(B.isconj() == C.isconj());

        int n = A.colsize();
        int nrhs = B.rowsize();
        double a(1);
        int ldB = B.stepj();
        double xbeta(beta);
        int ldC = C.stepj();
        LAPNAME(zlagtm) (
            LAPCM B.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
                         LAPV(a),LAPP(A.cptr()+A.stepj()),
                         LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
                         LAPP(B.cptr()),LAPV(ldB),LAPV(xbeta),LAPP(C.ptr()),LAPV(ldC) LAP1);
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapTriDiagMultMM(
        const GenBandMatrix<float>& A, const GenMatrix<float>& B,
        const int beta, const MatrixView<float>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(B.ct()==NonConj);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.nlo() == 1);
        TMVAssert(A.nhi() == 1);
        TMVAssert(A.isdm());
        TMVAssert(A.isSquare());
        TMVAssert(B.iscm());
        TMVAssert(C.iscm());

        int n = A.colsize();
        int nrhs = B.rowsize();
        float a(1);
        int ldB = B.stepj();
        float xbeta(beta);
        int ldC = C.stepj();
        LAPNAME(slagtm) (
            LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
                         LAPV(a),LAPP(A.cptr()+A.stepj()),
                         LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
                         LAPP(B.cptr()),LAPV(ldB),LAPV(xbeta),LAPP(C.ptr()),LAPV(ldC) LAP1);
    }
    template <> 
    void LapTriDiagMultMM(
        const GenBandMatrix<std::complex<float> >& A, 
        const GenMatrix<std::complex<float> >& B,
        const int beta, const MatrixView<std::complex<float> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.nlo() == 1);
        TMVAssert(A.nhi() == 1);
        TMVAssert(A.isdm());
        TMVAssert(A.isSquare());
        TMVAssert(B.iscm());
        TMVAssert(C.iscm());
        TMVAssert(B.isconj() == C.isconj());

        int n = A.colsize();
        int nrhs = B.rowsize();
        float a(1);
        int ldB = B.stepj();
        float xbeta(beta);
        int ldC = C.stepj();
        LAPNAME(clagtm) (
            LAPCM B.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
                         LAPV(a),LAPP(A.cptr()+A.stepj()),
                         LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
                         LAPP(B.cptr()),LAPV(ldB),LAPV(xbeta),LAPP(C.ptr()),LAPV(ldC) LAP1);
    }
#endif // FLOAT
#endif // ELAP

    template <bool add, class T, class Ta, class Tb> 
    static inline void doTriDiagMultMM(
        const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
#ifdef ELAP
        if (A.isSquare() && B.iscm() && C.iscm() && B.isconj() == C.isconj())
            LapTriDiagMultMM(A,B,add?1:0,C);
        else
#endif
            NonLapTriDiagMultMM<add>(A,B,C);
    }

    template <bool add, class T, class Ta, class Tb> 
    static void TriDiagMultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        if (alpha == T(1) && A.isdm()) {
            if (A.isconj()) 
                doTriDiagMultMM<add>(A.conjugate(),B.conjugate(),C.conjugate());
            else
                doTriDiagMultMM<add>(A,B,C);
        } else if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
            BandMatrix<Ta,DiagMajor> A1 = TMV_REAL(alpha)*A;
            doTriDiagMultMM<add>(A1,B,C);
        } else {
            BandMatrix<T,DiagMajor> A1 = alpha*A;
            doTriDiagMultMM<add>(A1,B,C);
        }
    }

    // MJ: Put in a recursive block calculation here.  (Also in Band*Band)
    template <bool add, class T, class Ta, class Tb> 
    static void doMultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
#ifdef XDEBUG
        //cout<<"Start doMultMM:\n";
        //cout<<"alpha = "<<alpha<<endl;
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A<<endl;
        //cout<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"  "<<B<<endl;
        //cout<<"add = "<<add<<endl;
        //cout<<"C = "<<TMV_Text(C)<<"  "<<C.cptr();
        //if (add) cout<<"  "<<C;
        //cout<<endl;
        Matrix<Tb> B0 = B;
        Matrix<Ta> A0 = A;
        Matrix<T> C0 = C;
        Matrix<T> C2 = C;
        MultMM<add>(alpha,A0,B0,C2.view());
        //cout<<"C2 = "<<C2<<endl;
#endif
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha!= T(0));
        TMVAssert(A.rowsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.ct()==NonConj);

        if (A.isrm() && C.isrm()) RowMultMM<add>(alpha,A,B,C);
        else if (A.iscm() && B.isrm()) OPMultMM<add>(alpha,A,B,C);
        else if (B.iscm() && C.iscm()) ColMultMM<add>(alpha,A,B,C);
        else if (A.nlo() == 1 && A.nhi() == 1) {
            if (A.isconj())
                TriDiagMultMM<add>(TMV_CONJ(alpha),A.conjugate(),B.conjugate(),
                                   C.conjugate());
            else
                TriDiagMultMM<add>(alpha,A,B,C);
        }
        else if (C.colsize() < C.rowsize()) RowMultMM<add>(alpha,A,B,C);
        else ColMultMM<add>(alpha,A,B,C);
#ifdef XDEBUG
        //cout<<"C -> "<<C<<endl;
        if (Norm(C2-C) > 0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                                (add?Norm(C0):TMV_RealType(T)(0)))) {
            cerr<<"doMultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C.cptr()<<"  "<<C0<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    static void FullTempMultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        if (C.isrm()) {
            Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
            doMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        } else {
            Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
            doMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void BlockTempMultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        const int N = C.rowsize();
        for (int j=0;j<N;) {
            int j2 = TMV_MIN(N,j+MM_BLOCKSIZE);
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                if (C.isrm()) {
                    Matrix<Tb,RowMajor> B2 = TMV_REAL(alpha) * B.colRange(j,j2);
                    doMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                } else  {
                    Matrix<Tb,ColMajor> B2 = TMV_REAL(alpha) * B.colRange(j,j2);
                    doMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                }
            } else {
                if (C.isrm()) {
                    Matrix<T,RowMajor> B2 = alpha * B.colRange(j,j2);
                    doMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                } else  {
                    Matrix<T,ColMajor> B2 = alpha * B.colRange(j,j2);
                    doMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                }
            }
            j=j2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(const T alpha,
              const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
              const MatrixView<T>& C)
        // C (+)= alpha * A * B
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());

#ifdef XDEBUG
        //cout<<"Start MultMM:\n";
        //cout<<"alpha = "<<alpha<<endl;
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A<<endl;
        //cout<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"  "<<B<<endl;
        //cout<<"add = "<<add<<endl;
        //cout<<"C = "<<TMV_Text(C)<<"  "<<C.cptr();
        //if (add) cout<<"  "<<C;
        //cout<<endl;
        Matrix<Tb> B0 = B;
        Matrix<Ta> A0 = A;
        Matrix<T> C0 = C;
        Matrix<T> C2 = C;
        MultMM<add>(alpha,A0,B0,C2.view());
        //cout<<"C2 = "<<C2<<endl;
#endif

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (A.rowsize() == 0 || alpha == T(0)) {
                if (!add) C.zero();
            } else if (A.rowsize() > A.colsize()+A.nhi()) {
                MultMM<add>(alpha,A.colRange(0,A.colsize()+A.nhi()),
                            B.rowRange(0,A.colsize()+A.nhi()),C);
            } else if (A.colsize() > A.rowsize()+A.nlo()) {
                MultMM<add>(alpha,A.rowRange(0,A.rowsize()+A.nlo()),
                            B,C.rowRange(0,A.rowsize()+A.nlo()));
                if (!add) C.rowRange(A.rowsize()+A.nlo(),A.colsize()).zero();
            } else if (C.isconj()) {
                MultMM<add>(TMV_CONJ(alpha),A.conjugate(),B.conjugate(),
                            C.conjugate());
            } else if (SameStorage(A,C)) {
                FullTempMultMM<add>(alpha,A,B,C);
            } else if (SameStorage(B,C)) {
                if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
                    BlockTempMultMM<add>(alpha,A,B,C);
                else 
                    FullTempMultMM<add>(alpha,A,B,C);
            } else {
                doMultMM<add>(alpha, A, B, C);
            }
        }
#ifdef XDEBUG
        //cout<<"C -> "<<C<<endl;
        if (Norm(C2-C) > 0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                                (add?Norm(C0):TMV_RealType(T)(0)))) {
            cerr<<"MultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C.cptr()<<"  "<<C0<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultBM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv
