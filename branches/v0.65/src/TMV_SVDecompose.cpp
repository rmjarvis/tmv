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
#include "TMV_SVDiv.h"
#include "tmv/TMV_SVD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "TMV_QRDiv.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include <iostream>

#ifdef XDEBUG
#define THRESH 1.e-11
#include "tmv/TMV_DiagMatrixArith.h"
#define dbgcout std::cout 
//#define dbgcout if(false) std::cout
using std::cerr;
using std::endl;
#else
#define dbgcout if(false) std::cout
#endif

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    static inline void DoSVDecomposeFromBidiagonal_NZ(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V,
        bool UisI, bool VisI)
    {
        SV_DecomposeFromBidiagonal_DC<T>(U,D,E,V,UisI,VisI); 
        //SV_DecomposeFromBidiagonal_QR<T>(U,D,E,V); 
        //UisI = VisI; if (UisI != VisI) abort();
    }

    template <class T> 
    void DoSVDecomposeFromBidiagonal(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V,
        bool UisI, bool VisI)
    {
        const int N = D.size();

        // D and E may have zeros on entry.  
        // This routine does the trivial deflations to get to subproblems
        // in which D and E are fully non-zero

        // First chop any small elements in D,E
        BidiagonalChopSmallElements(D,E);
        dbgcout<<"After Chop: D = "<<D<<std::endl;
        dbgcout<<"After Chop: E = "<<E<<std::endl;

        // Find sub-problems to solve:
        for(int q = N-1; q>0; ) {
            dbgcout<<"Looking for sub-problem:\n";
            dbgcout<<"q = "<<q<<std::endl;
            if (E(q-1) == T(0)) --q;
            else if (D(q) == T(0)) {
                dbgcout<<"D(q) == 0, so do ZeroLastCol\n";
                // We have the end looking like:
                //   ? ?
                //     ? x
                //       0
                // So we need to find a p where all E(i) with p<=i<q are 
                // non-zero.
                int p = q-1;
                while (p>0 && E(p-1) != T(0)) --p;
                // Now Zero out the last column:
                if (V) BidiagonalZeroLastCol<T>(
                    D.subVector(p,q),E.subVector(p,q),V->rowRange(p,q+1));
                else BidiagonalZeroLastCol<T>(
                    D.subVector(p,q),E.subVector(p,q),0);
                VisI = false;
                --q;
            } else {
                // Find first p before q with either E(p) = 0 or D(p) = 0
                int p=q-1;
                while (p>0 && E(p-1)!=T(0) && D(p)!=T(0)) --p; 
                dbgcout<<"p = "<<p<<std::endl;
                if (D(p) == T(0)) {
                    dbgcout<<"D(p) == 0, so do ZeroFirstRow\n";
                    // We have a block looking like:
                    //   0 x
                    //     x x 
                    //       x x
                    //         x
                    if (U) 
                        BidiagonalZeroFirstRow<T>(
                            U->colRange(p,q+1),
                            D.subVector(p+1,q+1),E.subVector(p,q));
                    else 
                        BidiagonalZeroFirstRow<T>(
                            0, D.subVector(p+1,q+1),E.subVector(p,q));
                    UisI = false;
                    ++p;
                }
                if (q > p) {
                    dbgcout<<"No zeros in D,E:\n";
                    dbgcout<<"D = "<<D.subVector(p,q+1)<<std::endl;
                    dbgcout<<"E = "<<E.subVector(p,q)<<std::endl;
                    if (U)
                        if (V) 
                            DoSVDecomposeFromBidiagonal_NZ<T>(
                                U->colRange(p,q+1),D.subVector(p,q+1),
                                E.subVector(p,q),V->rowRange(p,q+1),
                                UisI && p==0 && q+1==N, VisI && p==0 && q+1==N);
                        else 
                            DoSVDecomposeFromBidiagonal_NZ<T>(
                                U->colRange(p,q+1),D.subVector(p,q+1),
                                E.subVector(p,q),0, 
                                UisI && p==0 && q+1==N, false);
                    else
                        if (V) 
                            DoSVDecomposeFromBidiagonal_NZ<T>(
                                0,D.subVector(p,q+1),
                                E.subVector(p,q),V->rowRange(p,q+1),
                                false, VisI && p==0 && q+1==N);
                        else 
                            DoSVDecomposeFromBidiagonal_NZ<T>(
                                0,D.subVector(p,q+1),E.subVector(p,q),0,
                                false, false);
                }
                q = p;
            }
        }
    }

    template <class T> 
    static void NonLapSVDecomposeFromBidiagonal(
        MVP<T> U, const VectorView<RT>& D, 
        const VectorView<RT>& E, MVP<T> V, bool setUV)
    {
#ifdef XDEBUG
        dbgcout<<"Start Decompose from Bidiag (NonLap):\n";
        if (U) dbgcout<<"U = "<<TMV_Text(*U)<<endl;
        if (V) dbgcout<<"V = "<<TMV_Text(*V)<<endl;
        dbgcout<<"D = "<<TMV_Text(D)<<"  step "<<D.step()<<"  "<<D<<endl;
        dbgcout<<"E = "<<TMV_Text(E)<<"  step "<<E.step()<<"  "<<E<<endl;
        //if (U) dbgcout<<"U = "<<*U<<endl;
        //if (V) dbgcout<<"V = "<<*V<<endl;

        dbgcout<<"setUV = "<<setUV<<endl;
        Matrix<RT> B(D.size(),D.size(),RT(0));
        B.diag() = D;
        B.diag(1) = E;
        Matrix<T> A0(U&&V ? U->colsize() : D.size(),D.size());
        if (U && V && !setUV) A0 = (*U) * B * (*V);
        else A0 = B;
        //dbgcout<<"A0 = "<<A0<<endl;
#endif

        const int N = D.size();

        if (setUV) {
            TMVAssert(U && V);
            U->setToIdentity();
            V->setToIdentity();
        }

        // Before running the normal algorithms, rescale D,E by the maximum
        // value to help avoid overflow and underflow.
        RT scale = TMV_MAX(D.maxAbs2Element(),E.maxAbs2Element());
        dbgcout<<"scale = "<<scale<<std::endl;
        dbgcout<<"1/scale = "<<RT(1)/scale<<std::endl;
        if (TMV_Underflow(scale)) {
            // Hopeless case.  Just zero out D,E and call it done.
            D.setZero();
            E.setZero();
            return;
        }
        D /= scale;
        E /= scale;
        dbgcout<<"After scale: \n";
        dbgcout<<"D = "<<D<<std::endl;
        dbgcout<<"E = "<<E<<std::endl;

        dbgcout<<"Before Call DoSVDecompose: D = "<<D<<std::endl;
        DoSVDecomposeFromBidiagonal<T>(U,D,E,V,setUV,setUV);
        dbgcout<<"After Call DoSVDecompose: D = "<<D<<std::endl;

        // Make all of the singular values positive
        RT* Di = D.ptr();
        for(int i=0;i<N;++i,++Di) if (*Di < 0) {
#ifdef TMVFLDEBUG
            TMVAssert(Di >= D.first);
            TMVAssert(Di < D.last);
#endif
            *Di = -(*Di);
            if (V) V->row(i) = -V->row(i);
        }
        dbgcout<<"After make positive: \n";
        dbgcout<<"D = "<<D<<std::endl;

        // Now A = U * S * V
        // Sort output singular values 
        AlignedArray<int> sortp(N);
        D.sort(sortp.get(),Descend);
        if (U) U->permuteCols(sortp.get());
        if (V) V->permuteRows(sortp.get());

        // Undo the scaling
        D *= scale;
        dbgcout<<"After undo scale: \n";
        dbgcout<<"D = "<<D<<std::endl;

#ifdef XDEBUG
        if (U && V) {
            Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
            if (Norm(A0-AA) > THRESH*Norm(A0)) {
                cerr<<"SV_DecomposeFromBidiagonal: \n";
                cerr<<"input B = "<<B<<endl;
                cerr<<"U => "<<*U<<endl;
                cerr<<"S => "<<D<<endl;
                cerr<<"V => "<<*V<<endl;
                cerr<<"UBV = "<<A0<<endl;
                cerr<<"USV = "<<AA<<endl;
                cerr<<"diff = ";
                (A0-AA).write(cerr,(A0-AA).maxAbsElement()*1.e-3);
                cerr<<endl;
                cerr<<"Norm(diff) = "<<Norm(A0-AA)<<std::endl;
                cerr<<"THRESH*Norm(A0) = "<<THRESH<<"*"<<Norm(A0)<<" = "<<THRESH*Norm(A0)<<std::endl;
                abort();
            }
        }
#endif
    }

#ifdef LAP 
    template <class T> 
    static inline void LapSVDecomposeFromBidiagonal(
        MVP<T> U, const VectorView<RT>& D, 
        const VectorView<RT>& E, MVP<T> V, bool setUV)
    { NonLapSVDecomposeFromBidiagonal<T>(U,D,E,V,setUV); }
#ifdef INST_DOUBLE
    template <> 
    void LapSVDecomposeFromBidiagonal(
        MVP<double> U, const VectorView<double>& D, 
        const VectorView<double>& E, MVP<double> V, bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
            TMVAssert(U->ct()==NonConj);
        }
        if (V) { 
            TMVAssert(V->rowsize() == V->colsize()); 
            TMVAssert(V->rowsize() == D.size()); 
            TMVAssert(V->ct()==NonConj);
            if (U && setUV) TMVAssert(U->stor() == V->stor());
        }
        char u = 'U';
        int n = D.size();
        Vector<double> E1(n);
        E1.subVector(0,n-1) = E;
        if (setUV) {
            char c = 'I';
            TMVAssert(U && V);
            if (U->iscm()) {
                TMVAssert(V->iscm());
                int ldu = U->stepj();
                int ldv = V->stepj();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<double> work(lwork);
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(dbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(U->ptr()),LAPV(ldu),LAPP(V->ptr()),LAPV(ldv),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            } else {
                u = 'L';
                TMVAssert(U->isrm());
                TMVAssert(V->isrm());
                int ldu = U->stepi();
                int ldv = V->stepi();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<double> work(lwork);
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(dbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(V->ptr()),LAPV(ldv),LAPP(U->ptr()),LAPV(ldu),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            }
        } else if (U || V) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n);
            Matrix<double,ColMajor> V1(n,n);
            int ldu = U1.stepj();
            int ldv = V1.stepj();
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<double> work(lwork);
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (U) *U = *U*U1;
            if (V) *V = V1*(*V);
        } else {
            int ldu = n;
            int ldv = n;
            char c = 'N';
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<double> work(lwork);
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results("dbdsdc");
        E = E1.subVector(0,n-1);
    }
    template <> 
    void LapSVDecomposeFromBidiagonal(
        MVP<std::complex<double> > U, const VectorView<double>& D, 
        const VectorView<double>& E, MVP<std::complex<double> > V, 
        bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
            TMVAssert(U->ct()==NonConj);
        }
        if (V) { 
            TMVAssert(V->rowsize() == V->colsize()); 
            TMVAssert(V->rowsize() == D.size()); 
            TMVAssert(V->ct()==NonConj);
        }

        char u = 'U';
        int n = D.size();
        Vector<double> E1(n);
        E1.subVector(0,n-1) = E;
        if (U || V) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n);
            Matrix<double,ColMajor> V1(n,n);
            int ldu = U1.stepj();
            int ldv = V1.stepj();
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<double> work(lwork);
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (setUV) {
                if (U) *U = U1;
                if (V) *V = V1;
            } else {
                if (U) *U = *U*U1;
                if (V) *V = V1*(*V);
            }
        } else {
            TMVAssert(!setUV);
            int ldu = n;
            int ldv = n;
            char c = 'N';
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<double> work(lwork);
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results("dbdsdc");
        E = E1.subVector(0,n-1);
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapSVDecomposeFromBidiagonal(
        MVP<float> U, const VectorView<float>& D, 
        const VectorView<float>& E, MVP<float> V, bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
            TMVAssert(U->ct()==NonConj);
        }
        if (V) { 
            TMVAssert(V->rowsize() == V->colsize()); 
            TMVAssert(V->rowsize() == D.size()); 
            TMVAssert(V->ct()==NonConj);
            if (U && setUV) TMVAssert(U->stor() == V->stor());
        }
        char u = 'U';
        int n = D.size();
        Vector<float> E1(n);
        E1.subVector(0,n-1) = E;
        if (setUV) {
            char c = 'I';
            TMVAssert(U && V);
            if (U->iscm()) {
                TMVAssert(V->iscm());
                int ldu = U->stepj();
                int ldv = V->stepj();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<float> work(lwork);
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(sbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(U->ptr()),LAPV(ldu),LAPP(V->ptr()),LAPV(ldv),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            } else {
                u = 'L';
                TMVAssert(U->isrm());
                TMVAssert(V->isrm());
                int ldu = U->stepi();
                int ldv = V->stepi();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<float> work(lwork);
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(sbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(V->ptr()),LAPV(ldv),LAPP(U->ptr()),LAPV(ldu),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            }
        } else if (U || V) {
            char c = 'I';
            Matrix<float,ColMajor> U1(n,n);
            Matrix<float,ColMajor> V1(n,n);
            int ldu = U1.stepj();
            int ldv = V1.stepj();
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<float> work(lwork);
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(sbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (U) *U = *U*U1;
            if (V) *V = V1*(*V);
        } else {
            int ldu = n;
            int ldv = n;
            char c = 'N';
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<float> work(lwork);
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(sbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results("sbdsdc");
        E = E1.subVector(0,n-1);
    }
    template <> 
    void LapSVDecomposeFromBidiagonal(
        MVP<std::complex<float> > U, const VectorView<float>& D, 
        const VectorView<float>& E, MVP<std::complex<float> > V, 
        bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
            TMVAssert(U->ct()==NonConj);
        }
        if (V) { 
            TMVAssert(V->rowsize() == V->colsize()); 
            TMVAssert(V->rowsize() == D.size()); 
            TMVAssert(V->ct()==NonConj);
        }

        char u = 'U';
        int n = D.size();
        Vector<float> E1(n);
        E1.subVector(0,n-1) = E;
        if (U || V) {
            char c = 'I';
            Matrix<float,ColMajor> U1(n,n);
            Matrix<float,ColMajor> V1(n,n);
            int ldu = U1.stepj();
            int ldv = V1.stepj();
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<float> work(lwork);
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(sbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (setUV) {
                TMVAssert(U && V);
                *U = U1;
                *V = V1;
            } else {
                if (U) *U = *U*U1;
                if (V) *V = V1*(*V);
            }
        } else {
            TMVAssert(!setUV);
            int ldu = n;
            int ldv = n;
            char c = 'N';
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<float> work(lwork);
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(sbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results("sbdsdc");
        E = E1.subVector(0,n-1);
    }
#endif // FLOAT
#endif // LAP

    template <class T> 
    void SV_DecomposeFromBidiagonal(
        MVP<T> U, const VectorView<RT>& D, 
        const VectorView<RT>& E, MVP<T> V, bool setUV)
    {
#ifdef XDEBUG
        Matrix<RT> B(D.size(),D.size(),RT(0));
        B.diag() = D;
        B.diag(1) = E;
        Matrix<T> A0(U&&V ? U->colsize() : D.size(),D.size());
        if (U && V && !setUV) A0 = (*U) * B * (*V);
        else A0 = B;
        Matrix<T> U0(U ? U->colsize() : D.size(),D.size());
        if (U) U0 = (*U);
        Matrix<T> V0(D.size(),D.size());
        if (V) V0 = (*V);
        Vector<RT> D0 = D;
        Vector<RT> E0 = E;
#endif
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
            TMVAssert(U->ct()==NonConj);
        }
        if (V) { 
            TMVAssert(V->rowsize() == V->colsize()); 
            TMVAssert(V->rowsize() == D.size()); 
            TMVAssert(V->ct()==NonConj);
        }
        TMVAssert((!U || U->iscm() || U->isrm()));
        TMVAssert((!V || V->iscm() || V->isrm()));
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(!U || !V || !setUV || U->stor()==V->stor());

        if (D.size() > 0) {
#ifdef LAP
            LapSVDecomposeFromBidiagonal(U,D,E,V,setUV);
            //RT Dmax = MaxAbsElement(D);
            //D.clip(TMV_Epsilon<T>()*Dmax);
#else 
            NonLapSVDecomposeFromBidiagonal(U,D,E,V,setUV);
#endif
        }
#ifdef XDEBUG
        if (U && V) {
            Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
            dbgcout<<"U = "<<*U<<std::endl;
            dbgcout<<"D = "<<D<<std::endl;
            dbgcout<<"V = "<<*V<<std::endl;
            dbgcout<<"AA = "<<AA<<std::endl;
            dbgcout<<"A0 = "<<A0<<std::endl;
            dbgcout<<"A0-AA = "<<Matrix<T>(A0-AA).clip(1.e-3)<<std::endl;
            dbgcout<<"SVDecomposeFromBidiag: Norm(A0-AA) = "<<Norm(A0-AA)<<std::endl;
            dbgcout<<"cf "<<THRESH*Norm(A0)<<std::endl;
            if (Norm(A0-AA) > THRESH*Norm(A0)) {
                cerr<<"SV_DecomposeFromBidiagonal: \n";
                cerr<<"input B = "<<B<<endl;
                cerr<<"UBV = "<<A0<<endl;
                cerr<<"USV = "<<AA<<endl;
                cerr<<"diff = ";
                (A0-AA).write(cerr,(A0-AA).maxAbsElement()*1.e-3);
                cerr<<endl;
                cerr<<"U = "<<*U<<endl;
                cerr<<"S = "<<D<<endl;
                cerr<<"V = "<<*V<<endl;
#ifdef LAP
                Matrix<T,ColMajor> U2 = U0;
                Matrix<T,ColMajor> V2 = V0;
                Vector<RT> D2 = D0;
                Vector<RT> E2 = E0;
                NonLapSVDecomposeFromBidiagonal<T>(
                    U2.view(),D2.view(),E2.view(),V2.view(),setUV);
                cerr<<"U = "<<*U<<endl;
                cerr<<"NonLap U: "<<U2<<endl;
                cerr<<"Diff U = "<<Matrix<T>(*U-U2).clip(1.e-3);
                cerr<<"V = "<<*V<<endl;
                cerr<<"NonLap V: "<<V2<<endl;
                cerr<<"Diff V = "<<Matrix<T>(*V-V2).clip(1.e-3);
                cerr<<"D = "<<D<<endl;
                cerr<<"NonLap D: = "<<D2<<endl;
                cerr<<"Diff D = "<<Vector<T>(D-D2).clip(1.e-3);
                cerr<<"E = "<<E<<endl;
                cerr<<"NonLap: = "<<E2<<endl;
                cerr<<"Norm(U-U2) = "<<Norm(*U-U2)<<std::endl;
                cerr<<"Norm(V-V2) = "<<Norm(*V-V2)<<std::endl;
                cerr<<"Norm(D-D2) = "<<Norm(D-D2)<<std::endl;
#endif
                abort();
            }
        }
#endif
    }

    //
    // Main SVD Drivers
    //

    template <class T> 
    void SV_Decompose(
        const MatrixView<T>& U, const DiagMatrixView<RT>& S, 
        MVP<T> V, RT& logdet, T& signdet, bool StoreU)
    {
#ifdef XDEBUG
        Matrix<T> A0(U);
        dbgcout<<"SVDecompose:\n";
        //dbgcout<<"A0 = "<<A0<<std::endl;
        dbgcout<<"StoreU = "<<StoreU<<std::endl;
#endif
        // Decompose A (input as U) into U S V
        // where S is a diagonal real matrix, and U,V are unitary matrices.
        // A,U are M x N (M >= N)
        // S,V are N x N
        // The determinant is returned in logdet, signdet.

        TMVAssert(U.rowsize() <= U.colsize());
        TMVAssert(U.ct() == NonConj);
        if (V) {
            TMVAssert(V->ct() == NonConj);
            TMVAssert(V->colsize() == U.rowsize());
            TMVAssert(V->rowsize() == U.rowsize());
        }
        TMVAssert(S.size() == U.rowsize());
        TMVAssert(U.iscm() || U.isrm());

        const int M = U.colsize();
        const int N = U.rowsize();
        if (N == 0) return;

        // If M is much larger than N (technically M > 5/3 N),
        // then it is quicker to start by doing a QR decomposition 
        // and then do SVD on the square R matrix.  
        // Thus, the final U of the SVD is Q (from the QR decomp)
        // times U from R's SVD.
        if (M > 5*N/3) {
            if (StoreU) {
                Matrix<T,ColMajor> R(N,N);
                R.lowerTri().offDiag().setZero();
                QR_Decompose(U,R.upperTri(),signdet);
                SV_Decompose(R.view(),S,V,logdet,signdet,StoreU);
                // Now R is a Unitary Matrix U'.  Need to multiply U by U'
                U = U*R;
            } else {
                Vector<T> Qbeta(N);
                QR_Decompose(U,Qbeta.view(),signdet);
                if (N > 1) U.rowRange(0,N).lowerTri().offDiag().setZero();
                SV_Decompose(U.rowRange(0,N),S,V,logdet,signdet,StoreU);
            }
        } else {
            // First we reduce A to bidiagonal form: A = U * B * V
            // using a series of Householder transformations.
            // The diagonal of the Bidiagonal Matrix B is stored in D.
            // The superdiagonal is stored in E.
            Vector<RT> E(N-1);
            Vector<T> Ubeta(N);
            Vector<T> Vbeta(N-1);
            dbgcout<<"Before Bidiagonalize:\n";
            dbgcout<<"U.maxAbs = "<<U.maxAbsElement()<<std::endl;
            Bidiagonalize(
                U,Ubeta.view(),Vbeta.view(),S.diag(),E.view(),signdet);
            dbgcout<<"After Bidiagonalize:\n";
            dbgcout<<"U.maxAbs = "<<U.maxAbsElement()<<std::endl;
            dbgcout<<"D.maxAbs = "<<S.maxAbsElement()<<std::endl;
            dbgcout<<"E.maxAbs = "<<E.maxAbsElement()<<std::endl;
            dbgcout<<"D = "<<S.diag()<<std::endl;
            dbgcout<<"E = "<<E<<std::endl;
            // The determinant of B is just the product of the diagonal elements:
            if (signdet != T(0)) {
                RT s;
                logdet += S.logDet(&s);
                signdet *= s;
            }

            // Now UV stores Householder vectors for U in lower diagonal columns 
            // (HLi) and Householder vectors for V in upper diagonal rows (HRi)
            // The Householder matrices for U are actually the adjoints of the 
            // matrices that bidiagonalize A, and for V are the transposes:
            // U = HLn-1t ... HL1t HL0t A HR0T HR1T ... HRn-2T
            // Using the fact that H Ht = I, we get A = U B V with:
            // U = HL0 ... HLn-1 
            if (V) {
                V->row(0).makeBasis(0);
                V->rowRange(1,N) = U.rowRange(0,N-1);
                V->col(0,1,N).setZero();
                GetQFromQR(V->subMatrix(1,N,1,N).transpose(),Vbeta);
            }
            if (StoreU) {
                GetQFromQR(U,Ubeta);
            }
            dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;

            if (StoreU) SV_DecomposeFromBidiagonal<T>(U,S.diag(),E.view(),V);
            else SV_DecomposeFromBidiagonal<T>(0,S.diag(),E.view(),V);
            dbgcout<<"After DecomposeFromBidiag: Norm(UtU-1) = "<<
                Norm(U.adjoint()*U-T(1))<<std::endl;

        }
#ifdef XDEBUG
        if (StoreU && V && S.size()>0) {
            Matrix<T> A2 = U * S * (*V);
            dbgcout<<"SVDecompose: Norm(A0-A2) = "<<Norm(A0-A2)<<std::endl;
            dbgcout<<"cf "<<THRESH*Norm(U)*Norm(S)*Norm(*V)<<std::endl;
            if (Norm(A0-A2) > THRESH * Norm(U) * Norm(S) * Norm(*V)) {
                cerr<<"SV_Decompose:\n";
                cerr<<"A = "<<A0<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"S = "<<S<<endl;
                cerr<<"V = "<<*V<<endl;
                cerr<<"USV = "<<A2<<endl;
                abort();
            }
        }
#endif
    }

    //
    // Driver routines:
    //
    template <class T> 
    void SV_Decompose(
        const MatrixView<T>& U, const DiagMatrixView<RT>& SS, 
        const MatrixView<T>& V, bool StoreU)
    {
        TMVAssert(U.colsize() >= U.rowsize());
        TMVAssert(SS.size() == U.rowsize());
        TMVAssert(V.colsize() == U.rowsize());
        TMVAssert(V.rowsize() == U.rowsize());
        if (U.isconj()) {
            if (V.isconj()) {
                SV_Decompose(U.conjugate(),SS,V.conjugate(),StoreU);
            } else {
                SV_Decompose(U.conjugate(),SS,V,StoreU);
                V.conjugateSelf();
            }
        } else {
            if (V.isconj()) {
                SV_Decompose(U,SS,V.conjugate(),StoreU);
                V.conjugateSelf();
            } else {
                RT ld(0);
                T d(0);
                SV_Decompose<T>(U,SS,V,ld,d,StoreU); 
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const MatrixView<T>& U, const DiagMatrixView<RT>& SS, bool StoreU)
    {
        TMVAssert(U.colsize() >= U.rowsize());
        TMVAssert(SS.size() == U.rowsize());
        if (U.isconj()) {
            SV_Decompose(U.conjugate(),SS,StoreU);
        } else {
            RT ld(0);
            T d(0);
            SV_Decompose<T>(U,SS,0,ld,d,StoreU); 
        }
    }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


