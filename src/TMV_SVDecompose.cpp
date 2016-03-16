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
#define THRESH 1.e-15
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
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E, 
        MatrixView<T> Vt, bool UisI, bool VisI)
    {
        SV_DecomposeFromBidiagonal_DC<T>(U,D,E,Vt,UisI,VisI); 
        //SV_DecomposeFromBidiagonal_QR<T>(U,D,E,Vt); 
        //UisI = VisI; if (UisI != VisI) abort();
    }

    template <class T> 
    void DoSVDecomposeFromBidiagonal(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<T> Vt, bool UisI, bool VisI)
    {
        const ptrdiff_t N = D.size();

        // D and E may have zeros on entry.  
        // This routine does the trivial deflations to get to subproblems
        // in which D and E are fully non-zero

        // First chop any small elements in D,E
        BidiagonalChopSmallElements(D,E);
        dbgcout<<"After Chop: D = "<<D<<std::endl;
        dbgcout<<"After Chop: E = "<<E<<std::endl;

        // Find sub-problems to solve:
        for(ptrdiff_t q = N-1; q>0; ) {
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
                ptrdiff_t p = q-1;
                while (p>0 && !(E(p-1) == T(0))) --p;
                // Now Zero out the last column:
                if (Vt.cptr()) BidiagonalZeroLastCol<T>(
                    D.subVector(p,q),E.subVector(p,q),Vt.rowRange(p,q+1));
                else BidiagonalZeroLastCol<T>(
                    D.subVector(p,q),E.subVector(p,q),Vt);
                VisI = false;
                --q;
            } else {
                // Find first p before q with either E(p) = 0 or D(p) = 0
                ptrdiff_t p=q-1;
                while (p>0 && !(E(p-1)==T(0)) && !(D(p)==T(0))) --p; 
                dbgcout<<"p = "<<p<<std::endl;
                if (D(p) == T(0)) {
                    dbgcout<<"D(p) == 0, so do ZeroFirstRow\n";
                    // We have a block looking like:
                    //   0 x
                    //     x x 
                    //       x x
                    //         x
                    if (U.cptr()) 
                        BidiagonalZeroFirstRow<T>(
                            U.colRange(p,q+1),
                            D.subVector(p+1,q+1),E.subVector(p,q));
                    else 
                        BidiagonalZeroFirstRow<T>(
                            U, D.subVector(p+1,q+1),E.subVector(p,q));
                    UisI = false;
                    ++p;
                }
                if (q > p) {
                    dbgcout<<"No zeros in D,E:\n";
                    dbgcout<<"D = "<<D.subVector(p,q+1)<<std::endl;
                    dbgcout<<"E = "<<E.subVector(p,q)<<std::endl;
                    if (U.cptr())
                        if (Vt.cptr()) 
                            DoSVDecomposeFromBidiagonal_NZ<T>(
                                U.colRange(p,q+1),D.subVector(p,q+1),
                                E.subVector(p,q),Vt.rowRange(p,q+1),
                                UisI && p==0 && q+1==N, VisI && p==0 && q+1==N);
                        else 
                            DoSVDecomposeFromBidiagonal_NZ<T>(
                                U.colRange(p,q+1),D.subVector(p,q+1),
                                E.subVector(p,q),Vt, 
                                UisI && p==0 && q+1==N, false);
                    else
                        if (Vt.cptr()) 
                            DoSVDecomposeFromBidiagonal_NZ<T>(
                                U,D.subVector(p,q+1),
                                E.subVector(p,q),Vt.rowRange(p,q+1),
                                false, VisI && p==0 && q+1==N);
                        else 
                            DoSVDecomposeFromBidiagonal_NZ<T>(
                                U,D.subVector(p,q+1),E.subVector(p,q),Vt,
                                false, false);
                }
                q = p;
            }
        }
    }

    template <class T> 
    static void NonLapSVDecomposeFromBidiagonal(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<T> Vt, bool setUV)
    {
#ifdef XDEBUG
        dbgcout<<"Start Decompose from Bidiag (NonLap):\n";
        if (U.cptr()) dbgcout<<"U = "<<TMV_Text(U)<<endl;
        if (Vt.cptr()) dbgcout<<"V = "<<TMV_Text(Vt)<<endl;
        dbgcout<<"D = "<<TMV_Text(D)<<"  step "<<D.step()<<"  "<<D<<endl;
        dbgcout<<"E = "<<TMV_Text(E)<<"  step "<<E.step()<<"  "<<E<<endl;
        //if (U.cptr()) dbgcout<<"U = "<<U<<endl;
        //if (Vt.cptr()) dbgcout<<"Vt = "<<Vt<<endl;

        dbgcout<<"setUV = "<<setUV<<endl;
        Matrix<RT> B(D.size(),D.size(),RT(0));
        B.diag() = D;
        B.diag(1) = E;
        Matrix<T> A0(U.cptr()&&Vt.cptr() ? U.colsize() : D.size(),D.size());
        if (U.cptr() && Vt.cptr() && !setUV) A0 = U * B * Vt;
        else A0 = B;
        //dbgcout<<"A0 = "<<A0<<endl;
#endif

        const ptrdiff_t N = D.size();

        if (setUV) {
            TMVAssert(U.cptr() && Vt.cptr());
            U.setToIdentity();
            Vt.setToIdentity();
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
        DoSVDecomposeFromBidiagonal<T>(U,D,E,Vt,setUV,setUV);
        dbgcout<<"After Call DoSVDecompose: D = "<<D<<std::endl;

        // Make all of the singular values positive
        RT* Di = D.ptr();
        for(ptrdiff_t i=0;i<N;++i,++Di) if (*Di < 0) {
#ifdef TMVFLDEBUG
            TMVAssert(Di >= D._first);
            TMVAssert(Di < D._last);
#endif
            *Di = -(*Di);
            if (Vt.cptr()) Vt.row(i) = -Vt.row(i);
        }
        dbgcout<<"After make positive: \n";
        dbgcout<<"D = "<<D<<std::endl;

        // Now A = U * S * Vt
        // Sort output singular values 
        AlignedArray<ptrdiff_t> sortp(N);
        D.sort(sortp.get(),Descend);
        if (U.cptr()) U.permuteCols(sortp.get());
        if (Vt.cptr()) Vt.permuteRows(sortp.get());

        // Undo the scaling
        D *= scale;
        dbgcout<<"After undo scale: \n";
        dbgcout<<"D = "<<D<<std::endl;

#ifdef XDEBUG
        if (U.cptr() && Vt.cptr()) {
            Matrix<T> AA = U * DiagMatrixViewOf(D) * Vt;
            if (!(Norm(A0-AA) < THRESH*Norm(A0))) {
                cerr<<"SV_DecomposeFromBidiagonal: \n";
                cerr<<"input B = "<<B<<endl;
                cerr<<"U => "<<U<<endl;
                cerr<<"S => "<<D<<endl;
                cerr<<"Vt => "<<Vt<<endl;
                cerr<<"UBVt = "<<A0<<endl;
                cerr<<"USVt = "<<AA<<endl;
                cerr<<"diff = "<<ThreshIO((A0-AA).maxAbsElement()*1.e-3)<<
                    (A0-AA)<<endl;
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
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<T> Vt, bool setUV)
    { NonLapSVDecomposeFromBidiagonal<T>(U,D,E,Vt,setUV); }
#ifdef INST_DOUBLE
    template <> 
    void LapSVDecomposeFromBidiagonal(
        MatrixView<double> U, VectorView<double> D, VectorView<double> E,
        MatrixView<double> Vt, bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
            TMVAssert(U.ct()==NonConj);
        }
        if (Vt.cptr()) { 
            TMVAssert(Vt.rowsize() == Vt.colsize()); 
            TMVAssert(Vt.rowsize() == D.size()); 
            TMVAssert(Vt.ct()==NonConj);
        }

        char u = 'U';
        int n = D.size();
        Vector<double> E1(n);
        E1.subVector(0,n-1) = E;
        E1[n-1] = 0.;
        int Lap_info=0;
        if (setUV) {
            char c = 'I';
            TMVAssert(U.cptr() && Vt.cptr());
            //std::cout<<"setUV\n";
            //std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"Vt = "<<Vt<<std::endl;
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
            if (U.iscm()) {
                TMVAssert(Vt.iscm());
                int ldu = U.stepj();
                int ldv = Vt.stepj();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<double> work(lwork);
                VectorViewOf(work.get(),lwork).setZero();
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(dbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(U.ptr()),LAPV(ldu),LAPP(Vt.ptr()),LAPV(ldv),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            } else {
                u = 'L';
                TMVAssert(U.isrm());
                TMVAssert(Vt.isrm());
                int ldu = U.stepi();
                int ldv = Vt.stepi();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<double> work(lwork);
                VectorViewOf(work.get(),lwork).setZero();
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(dbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(Vt.ptr()),LAPV(ldv),LAPP(U.ptr()),LAPV(ldu),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            }
        } else if (U.cptr() || Vt.cptr()) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n,0.);
            Matrix<double,ColMajor> Vt1(n,n,0.);
            int ldu = U1.stepj();
            int ldv = Vt1.stepj();
            //std::cout<<"U || Vt\n";
            //if (U.cptr()) std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"U1 = "<<U1<<std::endl;
            //if (Vt.cptr()) std::cout<<"Vt = "<<Vt<<std::endl;
            //std::cout<<"Vt1 = "<<Vt1<<std::endl;
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(Vt1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (U.cptr()) U = U*U1;
            if (Vt.cptr()) Vt = Vt1*Vt;
        } else {
            //std::cout<<"!(U || Vt)\n";
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
            int ldu = n;
            int ldv = n;
            char c = 'N';
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results(Lap_info,"dbdsdc");
        E = E1.subVector(0,n-1);
        //std::cout<<"Done: D => "<<D<<std::endl;
        //std::cout<<"E1 => "<<E1<<std::endl;
        //std::cout<<"E => "<<E<<std::endl;
    }
    template <> 
    void LapSVDecomposeFromBidiagonal(
        MatrixView<std::complex<double> > U,
        VectorView<double> D, VectorView<double> E,
        MatrixView<std::complex<double> > Vt, bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
            TMVAssert(U.ct()==NonConj);
        }
        if (Vt.cptr()) { 
            TMVAssert(Vt.rowsize() == Vt.colsize()); 
            TMVAssert(Vt.rowsize() == D.size()); 
            TMVAssert(Vt.ct()==NonConj);
        }

        char u = 'U';
        int n = D.size();
        Vector<double> E1(n);
        E1.subVector(0,n-1) = E;
        E1[n-1] = 0.;
        int Lap_info=0;
        if (U.cptr() || Vt.cptr()) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n,0.);
            Matrix<double,ColMajor> Vt1(n,n,0.);
            int ldu = U1.stepj();
            int ldv = Vt1.stepj();
            //std::cout<<"U || Vt\n";
            //if (U.cptr()) std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"U1 = "<<U1<<std::endl;
            //if (Vt.cptr()) std::cout<<"Vt = "<<Vt<<std::endl;
            //std::cout<<"Vt1 = "<<Vt1<<std::endl;
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(Vt1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (setUV) {
                if (U.cptr()) U = U1;
                if (Vt.cptr()) Vt = Vt1;
            } else {
                if (U.cptr()) U = U*U1;
                if (Vt.cptr()) Vt = Vt1*Vt;
            }
        } else {
            TMVAssert(!setUV);
            int ldu = n;
            int ldv = n;
            char c = 'N';
            //std::cout<<"!(U || Vt)\n";
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results(Lap_info,"dbdsdc");
        E = E1.subVector(0,n-1);
        //std::cout<<"Done: D => "<<D<<std::endl;
        //std::cout<<"E1 => "<<E1<<std::endl;
        //std::cout<<"E => "<<E<<std::endl;
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapSVDecomposeFromBidiagonal(
        MatrixView<float> U, VectorView<float> D, VectorView<float> E,
        MatrixView<float> Vt, bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
            TMVAssert(U.ct()==NonConj);
        }
        if (Vt.cptr()) { 
            TMVAssert(Vt.rowsize() == Vt.colsize()); 
            TMVAssert(Vt.rowsize() == D.size()); 
            TMVAssert(Vt.ct()==NonConj);
        }

        char u = 'U';
        int n = D.size();
        Vector<float> E1(n);
        E1.subVector(0,n-1) = E;
        E1[n-1] = 0.;
        int Lap_info=0;
        if (setUV) {
            char c = 'I';
            TMVAssert(U.cptr() && Vt.cptr());
            if (U.iscm()) {
                TMVAssert(Vt.iscm());
                int ldu = U.stepj();
                int ldv = Vt.stepj();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<float> work(lwork);
                VectorViewOf(work.get(),lwork).setZero();
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(sbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(U.ptr()),LAPV(ldu),LAPP(Vt.ptr()),LAPV(ldv),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            } else {
                u = 'L';
                TMVAssert(U.isrm());
                TMVAssert(Vt.isrm());
                int ldu = U.stepi();
                int ldv = Vt.stepi();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<float> work(lwork);
                VectorViewOf(work.get(),lwork).setZero();
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(sbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(Vt.ptr()),LAPV(ldv),LAPP(U.ptr()),LAPV(ldu),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            }
        } else if (U.cptr() || Vt.cptr()) {
            char c = 'I';
            Matrix<float,ColMajor> U1(n,n,0.F);
            Matrix<float,ColMajor> Vt1(n,n,0.F);
            int ldu = U1.stepj();
            int ldv = Vt1.stepj();
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(sbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(Vt1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (U.cptr()) U = U*U1;
            if (Vt.cptr()) Vt = Vt1*Vt;
        } else {
            int ldu = n;
            int ldv = n;
            char c = 'N';
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(sbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results(Lap_info,"sbdsdc");
        E = E1.subVector(0,n-1);
    }
    template <> 
    void LapSVDecomposeFromBidiagonal(
        MatrixView<std::complex<float> > U,
        VectorView<float> D, VectorView<float> E,
        MatrixView<std::complex<float> > Vt, bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
            TMVAssert(U.ct()==NonConj);
        }
        if (Vt.cptr()) { 
            TMVAssert(Vt.rowsize() == Vt.colsize()); 
            TMVAssert(Vt.rowsize() == D.size()); 
            TMVAssert(Vt.ct()==NonConj);
        }

        char u = 'U';
        int n = D.size();
        Vector<float> E1(n);
        E1.subVector(0,n-1) = E;
        E1[n-1] = 0.;
        int Lap_info=0;
        if (U.cptr() || Vt.cptr()) {
            char c = 'I';
            Matrix<float,ColMajor> U1(n,n,0.F);
            Matrix<float,ColMajor> Vt1(n,n,0.F);
            int ldu = U1.stepj();
            int ldv = Vt1.stepj();
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(sbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(Vt1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (setUV) {
                TMVAssert(U.cptr() && Vt.cptr());
                U = U1;
                Vt = Vt1;
            } else {
                if (U.cptr()) U = U*U1;
                if (Vt.cptr()) Vt = Vt1*Vt;
            }
        } else {
            TMVAssert(!setUV);
            int ldu = n;
            int ldv = n;
            char c = 'N';
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(sbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results(Lap_info,"sbdsdc");
        E = E1.subVector(0,n-1);
    }
#endif // FLOAT
#endif // LAP

    template <class T> 
    void SV_DecomposeFromBidiagonal(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<T> Vt, bool setUV)
    {
#ifdef XDEBUG
        Matrix<RT> B(D.size(),D.size(),RT(0));
        B.diag() = D;
        B.diag(1) = E;
        Matrix<T> A0(U.cptr()&&Vt.cptr() ? U.colsize() : D.size(),D.size());
        if (U.cptr() && Vt.cptr() && !setUV) A0 = U * B * Vt;
        else A0 = B;
        Matrix<T> U0(U.cptr() ? U.colsize() : D.size(),D.size());
        if (U.cptr()) U0 = U;
        Matrix<T> Vt0(D.size(),D.size());
        if (Vt.cptr()) Vt0 = Vt;
        Vector<RT> D0 = D;
        Vector<RT> E0 = E;
#endif
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
            TMVAssert(U.ct()==NonConj);
        }
        if (Vt.cptr()) { 
            TMVAssert(Vt.rowsize() == Vt.colsize()); 
            TMVAssert(Vt.rowsize() == D.size()); 
            TMVAssert(Vt.ct()==NonConj);
        }
        TMVAssert((!U.cptr() || U.iscm() || U.isrm()));
        TMVAssert((!Vt.cptr() || Vt.iscm() || Vt.isrm()));
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);

        if (D.size() > 0) {
#ifdef LAP
            LapSVDecomposeFromBidiagonal(U,D,E,Vt,setUV);
            //RT Dmax = MaxAbsElement(D);
            //D.clip(TMV_Epsilon<T>()*Dmax);
#else 
            NonLapSVDecomposeFromBidiagonal(U,D,E,Vt,setUV);
#endif
        }
#ifdef XDEBUG
        if (U.cptr() && Vt.cptr()) {
            Matrix<T> AA = U * DiagMatrixViewOf(D) * Vt;
            dbgcout<<"U = "<<U<<std::endl;
            dbgcout<<"D = "<<D<<std::endl;
            dbgcout<<"Vt = "<<Vt<<std::endl;
            dbgcout<<"AA = "<<AA<<std::endl;
            dbgcout<<"A0 = "<<A0<<std::endl;
            dbgcout<<"A0-AA = "<<Matrix<T>(A0-AA).clip(1.e-3)<<std::endl;
            dbgcout<<"SVDecomposeFromBidiag: Norm(A0-AA) = "<<Norm(A0-AA)<<std::endl;
            dbgcout<<"cf "<<THRESH*Norm(A0)<<std::endl;
            if (!(Norm(A0-AA) < THRESH*Norm(A0))) {
                cerr<<"SV_DecomposeFromBidiagonal: \n";
                cerr<<"input B = "<<B<<endl;
                cerr<<"UBVt = "<<A0<<endl;
                cerr<<"USVt = "<<AA<<endl;
                cerr<<"diff = "<<ThreshIO((A0-AA).maxAbsElement()*1.e-3)<<
                    (A0-AA)<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"S = "<<D<<endl;
                cerr<<"Vt = "<<Vt<<endl;
#ifdef LAP
                Matrix<T,ColMajor> U2 = U0;
                Matrix<T,ColMajor> Vt2 = Vt0;
                Vector<RT> D2 = D0;
                Vector<RT> E2 = E0;
                NonLapSVDecomposeFromBidiagonal<T>(
                    U2.view(),D2.view(),E2.view(),Vt2.view(),setUV);
                cerr<<"U = "<<U<<endl;
                cerr<<"NonLap U: "<<U2<<endl;
                cerr<<"Diff U = "<<ThreshIO(1.e-3)<<Matrix<T>(U-U2);
                cerr<<"Vt = "<<Vt<<endl;
                cerr<<"NonLap Vt: "<<Vt2<<endl;
                cerr<<"Diff Vt = "<<ThreshIO(1.e-3)<<Matrix<T>(Vt-Vt2);
                cerr<<"D = "<<D<<endl;
                cerr<<"NonLap D: = "<<D2<<endl;
                cerr<<"Diff D = "<<Vector<T>(D-D2).clip(1.e-3);
                cerr<<"E = "<<E<<endl;
                cerr<<"NonLap: = "<<E2<<endl;
                cerr<<"Norm(U-U2) = "<<Norm(U-U2)<<std::endl;
                cerr<<"Norm(Vt-Vt2) = "<<Norm(Vt-Vt2)<<std::endl;
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
        MatrixView<T> U, DiagMatrixView<RT> S, 
        MatrixView<T> Vt, RT& logdet, T& signdet, bool StoreU)
    {
#ifdef XDEBUG
        Matrix<T> A0(U);
        dbgcout<<"SVDecompose:\n";
        //dbgcout<<"A0 = "<<A0<<std::endl;
        dbgcout<<"StoreU = "<<StoreU<<std::endl;
#endif
        // Decompose A (input as U) into U S Vt
        // where S is a diagonal real matrix, and U,Vt are unitary matrices.
        // A,U are M x N (M >= N)
        // S,Vt are N x N
        // The determinant is returned in logdet, signdet.

        TMVAssert(U.rowsize() <= U.colsize());
        TMVAssert(U.ct() == NonConj);
        if (Vt.cptr()) {
            TMVAssert(Vt.ct() == NonConj);
            TMVAssert(Vt.colsize() == U.rowsize());
            TMVAssert(Vt.rowsize() == U.rowsize());
        }
        TMVAssert(S.size() == U.rowsize());
        TMVAssert(U.iscm() || U.isrm());

        const ptrdiff_t M = U.colsize();
        const ptrdiff_t N = U.rowsize();
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
                SV_Decompose(R.view(),S,Vt,logdet,signdet,StoreU);
                // Now R is a Unitary Matrix U'.  Need to multiply U by U'
                U = U*R;
            } else {
                Vector<T> Qbeta(N);
                QR_Decompose(U,Qbeta.view(),signdet);
                if (N > 1) U.rowRange(0,N).lowerTri().offDiag().setZero();
                SV_Decompose(U.rowRange(0,N),S,Vt,logdet,signdet,StoreU);
            }
        } else {
            // First we reduce A to bidiagonal form: A = U * B * Vt
            // using a series of Householder transformations.
            // The diagonal of the Bidiagonal Matrix B is stored in D.
            // The superdiagonal is stored in E.
            Vector<RT> E(N-1);
            Vector<T> Ubeta(N);
            Vector<T> Vtbeta(N-1);
            dbgcout<<"Before Bidiagonalize:\n";
            dbgcout<<"U.maxAbs = "<<U.maxAbsElement()<<std::endl;
            Bidiagonalize(
                U,Ubeta.view(),Vtbeta.view(),S.diag(),E.view(),signdet);
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

            // Now U stores Householder vectors for U in lower diagonal columns 
            // (HLi) and Householder vectors for Vt in upper diagonal rows (HRi)
            // The Householder matrices for U are actually the adjoints of the 
            // matrices that bidiagonalize A, and for Vt are the transposes:
            // U = HLn-1t ... HL1t HL0t A HR0T HR1T ... HRn-2T
            // Using the fact that H Ht = I, we get A = U B Vt with:
            // U = HL0 ... HLn-1 
            if (Vt.cptr()) {
                Vt.row(0).makeBasis(0);
                Vt.rowRange(1,N) = U.rowRange(0,N-1);
                Vt.col(0,1,N).setZero();
                GetQFromQR(Vt.subMatrix(1,N,1,N).transpose(),Vtbeta);
                dbgcout<<"Vt => "<<Vt<<std::endl;
                dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<std::endl;
                dbgcout<<"Norm(VVt-1) = "<<Norm(Vt.adjoint()*Vt-T(1))<<std::endl;
            }
            if (StoreU) {
                GetQFromQR(U,Ubeta);
                dbgcout<<"U => "<<U<<std::endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                SV_DecomposeFromBidiagonal<T>(U,S.diag(),E.view(),Vt);
                dbgcout<<"After DecomposeFromBidiag: Norm(UtU-1) = "<<
                    Norm(U.adjoint()*U-T(1))<<std::endl;
            } else {
                MatrixView<T> U2(0,0,0,1,1,NonConj);
                SV_DecomposeFromBidiagonal<T>(U2,S.diag(),E.view(),Vt);
            }

        }
#ifdef XDEBUG
        dbgcout<<"Done SVDecompose\n";
        dbgcout<<"S = "<<S.diag()<<std::endl;
        if (StoreU && Vt.cptr() && S.size()>0) {
            Matrix<T> A2 = U * S * Vt;
            dbgcout<<"SVDecompose: Norm(A0-A2) = "<<Norm(A0-A2)<<std::endl;
            dbgcout<<"cf "<<THRESH*Norm(U)*Norm(S)*Norm(Vt)<<std::endl;
            if (!(Norm(A0-A2) < THRESH * Norm(U) * Norm(S) * Norm(Vt))) {
                cerr<<"SV_Decompose:\n";
                cerr<<"A = "<<A0<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"S = "<<S.diag()<<endl;
                cerr<<"Vt = "<<Vt<<endl;
                cerr<<"USVt = "<<A2<<endl;
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
        MatrixView<T> U, DiagMatrixView<RT> SS, MatrixView<T> Vt, bool StoreU)
    {
        TMVAssert(U.colsize() >= U.rowsize());
        TMVAssert(SS.size() == U.rowsize());
        TMVAssert(Vt.colsize() == U.rowsize());
        TMVAssert(Vt.rowsize() == U.rowsize());
        if (U.isconj()) {
            if (Vt.isconj()) {
                SV_Decompose(U.conjugate(),SS,Vt.conjugate(),StoreU);
            } else {
                SV_Decompose(U.conjugate(),SS,Vt,StoreU);
                Vt.conjugateSelf();
            }
        } else {
            if (Vt.isconj()) {
                SV_Decompose(U,SS,Vt.conjugate(),StoreU);
                Vt.conjugateSelf();
            } else {
                RT ld(0);
                T d(0);
                SV_Decompose<T>(U,SS,Vt,ld,d,StoreU); 
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        MatrixView<T> U, DiagMatrixView<RT> SS, bool StoreU)
    {
        TMVAssert(U.colsize() >= U.rowsize());
        TMVAssert(SS.size() == U.rowsize());
        if (U.isconj()) {
            SV_Decompose(U.conjugate(),SS,StoreU);
        } else {
            RT ld(0);
            T d(0);
            MatrixView<T> Vt(0,0,0,1,1,NonConj);
            SV_Decompose<T>(U,SS,Vt,ld,d,StoreU); 
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


