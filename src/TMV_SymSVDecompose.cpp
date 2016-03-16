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

#include <sstream>
#include "TMV_Blas.h"
#include "TMV_SymSVDiv.h"
#include "tmv/TMV_SymSVD.h"
#include "tmv/TMV_SymCHD.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_BandMatrix.h"
#include "TMV_QRDiv.h"
#include "TMV_BandSVDiv.h"
#include "tmv/TMV_BandSVD.h"
#include "tmv/TMV_SVD.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include "TMV_IsNaN.h"

#ifdef XDEBUG
#define THRESH 1.e-5
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    static void DoEigenFromTridiagonal_NZ(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E)
    { 
        EigenFromTridiagonal_DC(U,D,E,NonConj);
        //EigenFromTridiagonal_QR(U,D,E); 
    }

    template <class T> 
    static void NonLapEigenFromTridiagonal(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E)
    {
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.isSquare());
            TMVAssert(U.rowsize() == D.size());
        }

        const ptrdiff_t N = D.size();

        // E may have zeros on entry.
        // This routine finds subproblems in which E is fully non-zero.

        // First chop small elements of D,E
        HermTridiagonalChopSmallElements(D,E);

        // Find sub-problems to solve:
        for(ptrdiff_t q = N-1; q>0; ) {
            if (E(q-1) == T(0)) --q;
            else {
                ptrdiff_t p=q-1;
                while (p > 0 && (E(p-1) != T(0))) --p; 
                // Set p such that E(p-1) = 0 and 
                // all E(i) with p<=i<q are non-zero.
                if (U.cptr()) {
                    DoEigenFromTridiagonal_NZ<T>(
                        U.colRange(p,q+1),D.subVector(p,q+1),E.subVector(p,q));
                } else {
                    DoEigenFromTridiagonal_NZ<T>(
                        U,D.subVector(p,q+1),E.subVector(p,q));
                }
                q = p;
            }
        }
    }

#ifdef LAP 
    template <class T> 
    static inline void LapEigenFromTridiagonal(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E)
    { NonLapEigenFromTridiagonal(U,D,E); }

#ifdef INST_DOUBLE
    static void AltLapEigenFromTridiagonal(
        MatrixView<double> U, VectorView<double> D, VectorView<double> E)
    {
        //std::cout<<"AltLapEigen double\n";
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        int n = D.size();
        if (U.cptr() && U.iscm() && U.colsize() == U.stepj()) {
            char c = 'V';
            int ldu = U.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
            int lgn = int(ceil(log(double(n))/log(2.)));
#ifdef NOWORKQUERY
            int lwork = 1+3*n+2*n*lgn+4*n*n;
            // The official recommendation only has 3*n*n at the end,
            // but this leads to seg faults.  I think this is a bug
            // in the LAPACK distribution.  And the above value actually
            // makes sense with respect to the zstedc recommendation for
            // workspace which has n*n in complex<double> and 
            // 1+3*n+2*n*lgn+3*n*n in the double workspace.
            // So the above value is the sum of these two as one might
            // expect, so I think that's probably right.
            int liwork = 6*(1+n)+5*n*lgn;
            AlignedArray<double> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            //std::cout<<"After dstedc work query:\n";
            //std::cout<<"lwork = "<<lwork<<", liwork = "<<liwork<<std::endl;
            //std::cout<<"cf: "<<1+3*n+2*n*lgn+4*n*n<<"  "<<6*(1+n)+5*n*lgn<<std::endl;
            // Even worse -- some implementations do not return the correct
            // work size when using a work query.
            // Add n^2 conditionally on the return being this value in 
            // case an implementation has a totally different work size,
            // in which case, we'll assume that they have it right.
            if (lwork == 1+3*n+2*n*lgn+3*n*n) lwork += n*n;
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            //std::cout<<"Before dstedc:\n";
            //std::cout<<"c = "<<c<<", n = "<<n<<std::endl;
            //std::cout<<"lwork = "<<lwork<<", liwork = "<<liwork<<std::endl;
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            //std::cout<<"After dstedc: lap_info = "<<Lap_info<<std::endl;
            //std::cout<<"work[0] = "<<work[0]<<", iwork[0] = "<<iwork[0]<<std::endl;
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"dstedc");
#endif
        } else if (U.cptr()) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n);
            int ldu = U1.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1+(4+n)*n;
            int liwork = 3+5*n;
            AlignedArray<double> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            //std::cout<<"Before dstedc:\n";
            //std::cout<<"c = "<<c<<", n = "<<n<<std::endl;
            //std::cout<<"lwork = "<<lwork<<", liwork = "<<liwork<<std::endl;
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            //std::cout<<"After dstedc: lap_info = "<<Lap_info<<std::endl;
            //std::cout<<"work[0] = "<<work[0]<<", iwork[0] = "<<iwork[0]<<std::endl;
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"dstedc");
#endif
            U *= U1;
        } else {
            char c = 'N';
            int ldu = n;
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1;
            int liwork = 1;
            AlignedArray<double> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            //std::cout<<"Before dstedc:\n";
            //std::cout<<"c = "<<c<<", n = "<<n<<std::endl;
            //std::cout<<"lwork = "<<lwork<<", liwork = "<<liwork<<std::endl;
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"dstedc");
#endif
        }
    }
    static void AltLapEigenFromTridiagonal(
        MatrixView<std::complex<double> > U, VectorView<double> D, 
        VectorView<double> E)
    {
        //std::cout<<"AltLapEigen complex double\n";
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        int n = D.size();
        if (U.cptr() && U.iscm() && U.colsize() == U.stepj()) {
            char c = 'V';
            int ldu = U.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lgn = int(ceil(log(double(n))/log(2.)));
            int lwork = n*n;
            int lrwork = 1+3*n+2*n*lgn+3*n*n;
            int liwork = 6*(1+n)+5*n*lgn;
            AlignedArray<std::complex<double> > work(lwork);
            AlignedArray<double> rwork(lrwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int lrwork = -1;
            int liwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            AlignedArray<double> rwork(1);
            rwork.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(zstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(rwork.get()) LAPVWK(lrwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(real(work[0]));
            lrwork = int(rwork[0]);
            liwork = iwork[0];
            //std::cout<<"After zstedc work query:\n";
            //std::cout<<"lwork = "<<lwork<<", lrwork = "<<lrwork<<", liwork = "<<liwork<<std::endl;
            //int lgn = int(ceil(log(double(n))/log(2.)));
            //std::cout<<"cf: "<<n*n<<"  "<<1+3*n+2*n*lgn+3*n*n<<"  "<<6*(1+n)+5*n*lgn<<std::endl;
            work.resize(lwork);
            rwork.resize(lrwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
            VectorViewOf(rwork.get(),lrwork).setZero();
#endif
            //std::cout<<"Before zstedc:\n";
            //std::cout<<"c = "<<c<<", n = "<<n<<std::endl;
            //std::cout<<"lwork = "<<lwork<<", lrwork = "<<lrwork<<", liwork = "<<liwork<<std::endl;
            LAPNAME(zstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(rwork.get()) LAPVWK(lrwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"zstedc");
#else
            LAP_Results(Lap_info,int(real(work[0])),n,n,lwork,"zstedc");
#endif
        } else if (U.cptr()) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n);
            int ldu = U1.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1+(4+n)*n;
            int liwork = 3+5*n;
            AlignedArray<double> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            //std::cout<<"Before dstedc:\n";
            //std::cout<<"c = "<<c<<", n = "<<n<<std::endl;
            //std::cout<<"lwork = "<<lwork<<", liwork = "<<liwork<<std::endl;
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"dstedc");
#endif
            U *= U1;
        } else {
            char c = 'N';
            int ldu = n;
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1;
            int liwork = 1;
            AlignedArray<double> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            //std::cout<<"Before dstedc:\n";
            //std::cout<<"c = "<<c<<", n = "<<n<<std::endl;
            //std::cout<<"lwork = "<<lwork<<", liwork = "<<liwork<<std::endl;
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"dstedc");
#endif
        }
    }
#endif
#ifdef INST_FLOAT
    static void AltLapEigenFromTridiagonal(
        MatrixView<float> U, VectorView<float> D, VectorView<float> E)
    {
        //std::cout<<"AltLapEigen float\n";
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        int n = D.size();
        if (U.cptr() && U.iscm() && U.colsize() == U.stepj()) {
            char c = 'V';
            int ldu = U.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
            int lgn = int(ceil(log(float(n))/log(2.)));
#ifdef NOWORKQUERY
            int lwork = 1+3*n+2*n*lgn+4*n*n;
            int liwork = 6*(1+n)+5*n*lgn;
            AlignedArray<float> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            if (lwork == 1+3*n+2*n*lgn+3*n*n) lwork += n*n;
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"sstedc");
#endif
        } else if (U.cptr()) {
            char c = 'I';
            Matrix<float,ColMajor> U1(n,n);
            int ldu = U1.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1+(4+n)*n;
            int liwork = 3+5*n;
            AlignedArray<float> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"sstedc");
#endif
            U *= U1;
        } else {
            char c = 'N';
            int ldu = n;
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1;
            int liwork = 1;
            AlignedArray<float> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"sstedc");
#endif
        }
    }
    static void AltLapEigenFromTridiagonal(
        MatrixView<std::complex<float> > U, VectorView<float> D, 
        VectorView<float> E)
    {
        //std::cout<<"AltLapEigen complex float\n";
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        int n = D.size();
        if (U.cptr() && U.iscm() && U.colsize() == U.stepj()) {
            char c = 'V';
            int ldu = U.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lgn = int(ceil(log(float(n))/log(2.)));
            int lwork = n*n;
            int lrwork = 1+3*n+2*n*lgn+3*n*n;
            int liwork = 6*(1+n)+5*n*lgn;
            AlignedArray<std::complex<float> > work(lwork);
            AlignedArray<float> rwork(lrwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int lrwork = -1;
            int liwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.;
            AlignedArray<float> rwork(1);
            rwork.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(cstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(rwork.get()) LAPVWK(lrwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(real(work[0]));
            lrwork = int(rwork[0]);
            liwork = iwork[0];
            work.resize(lwork);
            rwork.resize(lrwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
            VectorViewOf(rwork.get(),lrwork).setZero();
#endif
            LAPNAME(cstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(rwork.get()) LAPVWK(lrwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"cstedc");
#else
            LAP_Results(Lap_info,int(real(work[0])),n,n,lwork,"cstedc");
#endif
        } else if (U.cptr()) {
            char c = 'I';
            Matrix<float,ColMajor> U1(n,n);
            int ldu = U1.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1+(4+n)*n;
            int liwork = 3+5*n;
            AlignedArray<float> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"sstedc");
#endif
            U *= U1;
        } else {
            char c = 'N';
            int ldu = n;
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1;
            int liwork = 1;
            AlignedArray<float> work(lwork);
            AlignedArray<int> iwork(liwork);
#else
            int lwork = -1;
            int liwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.;
            AlignedArray<int> iwork(1);
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.resize(lwork);
            iwork.resize(liwork);
#endif
            VectorViewOf(work.get(),lwork).setZero();
#endif
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sstedc");
#else
            LAP_Results(Lap_info,int(work[0]),n,n,lwork,"sstedc");
#endif
        }
    }
#endif
#ifdef NOSTEGR
#ifdef INST_DOUBLE
    template <> 
    inline void LapEigenFromTridiagonal(
        MatrixView<double> U, VectorView<double> D, VectorView<double> E)
    { AltLapEigenFromTridiagonal(U,D,E); }
    template <> 
    inline void LapEigenFromTridiagonal(
        MatrixView<std::complex<double> > U, VectorView<double> D, 
        VectorView<double> E)
    { AltLapEigenFromTridiagonal(U,D,E); }
#endif
#ifdef INST_FLOAT
    template <> 
    inline void LapEigenFromTridiagonal(
        MatrixView<float> U, VectorView<float> D, VectorView<float> E)
    { AltLapEigenFromTridiagonal(U,D,E); }
    template <> 
    inline void LapEigenFromTridiagonal(
        MatrixView<std::complex<float> > U, VectorView<float> D, 
        VectorView<float> E)
    { AltLapEigenFromTridiagonal(U,D,E); }
#endif
#else // normal stegr implementaion
#ifdef INST_DOUBLE
    template <> 
    void LapEigenFromTridiagonal(
        MatrixView<double> U, VectorView<double> D, VectorView<double> E)
    {
        //std::cout<<"Regular LapEigen double\n";
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        int n = D.size();
        if (U.cptr()) {
            char c1 = 'V';
            char c2 = 'A';
            double junk = 0;
            int ijunk = 0;
            double tol = 0;
            int neigen = 0;
            Vector<double> Din = D;
            Vector<double> E1(n);
            E1.subVector(0,n-1) = E;  E1(n-1) = 0.;
            Vector<double> Dout(n);
            Matrix<double,ColMajor> U1(n,n);
            int ldu = U1.stepj();
            AlignedArray<int> isuppz(2*n);
            int Lap_info=0;
#ifndef LAPNOWORK
            int lwork = 18*n;
            int liwork = 10*n;
            AlignedArray<double> work(lwork);
            AlignedArray<int> iwork(liwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            //std::cout<<"U size,step = "<<U.colsize()<<" "<<U.rowsize()<<" ";
            //std::cout<<U.stepi()<<" "<<U.stepj()<<" ";
            //std::cout<<"U1 size,step = "<<U1.colsize()<<" "<<U1.rowsize()<<" ";
            //std::cout<<U1.stepi()<<" "<<U1.stepj()<<" ";
            //std::cout<<"dstegr:\n";
            //std::cout<<c1<<"  "<<c2<<"  "<<n<<std::endl;
            //std::cout<<Din.ptr()<<"  "<<E1.ptr()<<"  "<<neigen<<std::endl;
            //std::cout<<Dout.ptr()<<"  "<<U1.ptr()<<"  "<<ldu<<std::endl;
            //std::cout<<isuppz.get()<<"  "<<work.get()<<"  "<<lwork<<std::endl;
            //std::cout<<iwork.get()<<"  "<<liwork<<"  "<<Lap_info<<std::endl;
            LAPNAME(dstegr) (
                LAPCM LAPV(c1),LAPV(c2),LAPV(n),
                LAPP(Din.ptr()),LAPP(E1.ptr()),
                LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk),LAPV(tol),
                LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
                LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
            //cout<<"After dstegr\n";
            //cout<<"lapinfo = "<<Lap_info<<endl;
            //cout<<"E(n-1) = "<<E1(n-1)<<endl;
            //cout<<"neigen = "<<neigen<<endl;
            //cout<<"Norm(UtU-1) = "<<Norm(U1.adjoint()*U1-1.)<<endl;
            //cout<<"UtU.diag() = "<<(U1.adjoint()*U1).diag()<<endl;
            
            // The dstegr algorithm has a bug wherein it sometimes fails
            // to finish, in which case the output values are incorrect.
            // It informs that this has happened via the info variable.
            // Specifically, it sets info = 2 in this case.
            // We test for all info > 0 just to be sure.
            // When this happens, we call dstedc instead.
            //
            // Also, sometimes, the U matrix comes back with a column
            // of all zeros.  No error is indicated by Lap_info in 
            // these (rare) cases, so I'm not sure what is going on,
            // but I check for it here too.
            //
            // In addition to the above problem, the dstegr routine 
            // seems to not be very careful about nan issues.  
            // So we also check for nan's in U1, and call dstedc if 
            // there are any.
            double nantest = U1.linearView().sumElements();
            int badcol = -1;
            for(int j=0;j<n;++j) {
                if (U1.col(j).sumElements() < TMV_Epsilon<double>()) badcol = j;
            }
            //cout<<"nantest = "<<nantest<<endl;
            //cout<<"badcol = "<<badcol<<endl;
            if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n ||
                badcol >= 0 || isNaN(nantest)) {
                std::ostringstream ss;
                ss << "Error in LAPACK function dstegr: ";
                if (Lap_info > 0) 
                    ss << "Returned info = "<<Lap_info<<".  ";
                else if (E1(n-1) > 0.F) 
                    ss << "Returned non-zero E(n) = "<<E1(n-1)<<".  ";
                else if (neigen < n) 
                    ss << "Only found "<<neigen<<" / "<<n<<" eigenvectors.  ";
                else if (badcol >= 0) 
                    ss << "U found to be not unitary: U.col("<<badcol<<") = "<<
                        U1.col(badcol)<<".  ";
                else 
                    ss <<  "NaN found in eigenvector matrix.  ";
                ss << "Calling dstedc instead.";
                TMV_Warning(ss.str());
                return AltLapEigenFromTridiagonal(U,D,E);
            }
            LAP_Results(Lap_info,"dstegr");
            D = Dout;
            U *= U1;
        } else {
            return AltLapEigenFromTridiagonal(U,D,E);
        }
    }
    template <> 
    void LapEigenFromTridiagonal(
        MatrixView<std::complex<double> > U, VectorView<double> D, 
        VectorView<double> E)
    {
        std::cout<<"Regular LapEigen complex double\n";
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        int n = D.size();
        if (U.cptr()) {
            char c1 = 'V';
            char c2 = 'A';
            double junk = 0;
            int ijunk = 0;
            double tol = 0; // This is automatically bumped up to the minimum
            int neigen;
            Vector<double> Din = D;
            Vector<double> E1(n);
            E1.subVector(0,n-1) = E;  E1(n-1) = 0.;
            Vector<double> Dout(n);
            Matrix<double,ColMajor> U1(n,n);
            int ldu = U1.stepj();
            AlignedArray<int> isuppz(2*n);
            int Lap_info=0;
#ifndef LAPNOWORK
            int lwork = 18*n;
            int liwork = 10*n;
            AlignedArray<double> work(lwork);
            AlignedArray<int> iwork(liwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            //std::cout<<"U size,step = "<<U.colsize()<<" "<<U.rowsize()<<" ";
            //std::cout<<U.stepi()<<" "<<U.stepj()<<" ";
            //std::cout<<"U1 size,step = "<<U1.colsize()<<" "<<U1.rowsize()<<" ";
            //std::cout<<U1.stepi()<<" "<<U1.stepj()<<" ";
            //std::cout<<"dstegr:\n";
            //std::cout<<c1<<"  "<<c2<<"  "<<n<<std::endl;
            //std::cout<<Din.ptr()<<"  "<<E1.ptr()<<std::endl;
            //std::cout<<Dout.ptr()<<"  "<<U1.ptr()<<"  "<<ldu<<std::endl;
            LAPNAME(dstegr) (
                LAPCM LAPV(c1),LAPV(c2),LAPV(n),
                LAPP(Din.ptr()),LAPP(E1.ptr()),
                LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk),LAPV(tol),
                LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
                LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
            //cout<<"After dstegr\n";
            //cout<<"lapinfo = "<<Lap_info<<endl;
            //cout<<"E(n-1) = "<<E1(n-1)<<endl;
            //cout<<"neigen = "<<neigen<<endl;
            //cout<<"Norm(UtU-1) = "<<Norm(U1.adjoint()*U1-1.)<<endl;
            //cout<<"UtU.diag() = "<<(U1.adjoint()*U1).diag()<<endl;
            double nantest = U1.linearView().sumElements();
            int badcol = -1;
            for(int j=0;j<n;++j) {
                if (U1.col(j).sumElements() < TMV_Epsilon<double>()) badcol = j;
            }
            //cout<<"nantest = "<<nantest<<endl;
            //cout<<"badcol = "<<badcol<<endl;
            if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n ||
                badcol >= 0 || isNaN(nantest)) {
                std::ostringstream ss;
                ss << "Error in LAPACK function dstegr: ";
                if (Lap_info > 0) 
                    ss << "Returned info = "<<Lap_info<<".  ";
                else if (E1(n-1) > 0.F) 
                    ss << "Returned non-zero E(n) = "<<E1(n-1)<<".  ";
                else if (neigen < n) 
                    ss << "Only found "<<neigen<<" / "<<n<<" eigenvectors.  ";
                else if (badcol >= 0) 
                    ss << "U found to be not unitary: U.col("<<badcol<<") = "<<
                        U1.col(badcol)<<".  ";
                else 
                    ss <<  "NaN found in eigenvector matrix.  ";
                ss << "Calling dstedc instead.";
                TMV_Warning(ss.str());
                return AltLapEigenFromTridiagonal(U,D,E);
            }
            LAP_Results(Lap_info,"dstegr");
            D = Dout;
            U *= U1;
        } else {
            return AltLapEigenFromTridiagonal(U,D,E);
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapEigenFromTridiagonal(
        MatrixView<float> U, VectorView<float> D, VectorView<float> E)
    {
        //std::cout<<"Regular LapEigen float\n";
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        int n = D.size();
        if (U.cptr()) {
            char c1 = 'V';
            char c2 = 'A';
            float junk = 0;
            int ijunk = 0;
            float tol = 0; // This is automatically bumped up to the minimum
            int neigen;
            Vector<float> Din = D;
            Vector<float> E1(n);
            E1.subVector(0,n-1) = E;  E1(n-1) = 0.F;
            Vector<float> Dout(n);
            Matrix<float,ColMajor> U1(n,n);
            int ldu = U1.stepj();
            AlignedArray<int> isuppz(2*n);
            int Lap_info=0;
#ifndef LAPNOWORK
            int lwork = 18*n;
            int liwork = 10*n;
            AlignedArray<float> work(lwork);
            AlignedArray<int> iwork(liwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            LAPNAME(sstegr) (
                LAPCM LAPV(c1),LAPV(c2),LAPV(n),
                LAPP(Din.ptr()),LAPP(E1.ptr()),
                LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk), LAPV(tol),
                LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
                LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork)
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
            //cout<<"After sstegr\n";
            //cout<<"lapinfo = "<<Lap_info<<endl;
            //cout<<"E(n-1) = "<<E1(n-1)<<endl;
            //cout<<"neigen = "<<neigen<<endl;
            //cout<<"Norm(UtU-1) = "<<Norm(U1.adjoint()*U1-1.F)<<endl;
            //cout<<"UtU.diag() = "<<(U1.adjoint()*U1).diag()<<endl;
            float nantest = U1.linearView().sumElements();
            int badcol = -1;
            for(int j=0;j<n;++j) {
                if (U1.col(j).sumElements() < TMV_Epsilon<float>()) badcol = j;
            }
            //cout<<"nantest = "<<nantest<<endl;
            //cout<<"badcol = "<<badcol<<endl;
            if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n || 
                badcol >= 0 || isNaN(nantest)) {
                std::ostringstream ss;
                ss << "Error in LAPACK function sstegr: ";
                if (Lap_info > 0) 
                    ss << "Returned info = "<<Lap_info<<".  ";
                else if (E1(n-1) > 0.F) 
                    ss << "Returned non-zero E(n) = "<<E1(n-1)<<".  ";
                else if (neigen < n) 
                    ss << "Only found "<<neigen<<" / "<<n<<" eigenvectors.  ";
                else if (badcol >= 0) 
                    ss << "U found to be not unitary: U.col("<<badcol<<") = "<<
                        U1.col(badcol)<<".  ";
                else 
                    ss <<  "NaN found in eigenvector matrix.  ";
                ss << "Calling sstedc instead.";
                TMV_Warning(ss.str());
                return AltLapEigenFromTridiagonal(U,D,E);
            }
            LAP_Results(Lap_info,"sstegr");
            D = Dout;
            U *= U1;
        } else {
            return AltLapEigenFromTridiagonal(U,D,E);
        }
    }
    template <> 
    void LapEigenFromTridiagonal(
        MatrixView<std::complex<float> > U, VectorView<float> D, 
        VectorView<float> E)
    {
        //std::cout<<"Regular LapEigen complex float\n";
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        int n = D.size();
        if (U.cptr()) {
            char c1 = 'V';
            char c2 = 'A';
            float junk = 0;
            int ijunk = 0;
            float tol = 0; // This is automatically bumped up to the minimum
            int neigen;
            Vector<float> Din = D;
            Vector<float> E1(n);
            E1.subVector(0,n-1) = E;  E1(n-1) = 0.F;
            Vector<float> Dout(n);
            Matrix<float,ColMajor> U1(n,n);
            int ldu = U1.stepj();
            AlignedArray<int> isuppz(2*n);
            int Lap_info=0;
#ifndef LAPNOWORK
            int lwork = 18*n;
            int liwork = 10*n;
            AlignedArray<float> work(lwork);
            AlignedArray<int> iwork(liwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            LAPNAME(sstegr) (
                LAPCM LAPV(c1),LAPV(c2),LAPV(n),
                LAPP(Din.ptr()),LAPP(E1.ptr()),
                LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk), LAPV(tol),
                LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
                LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
            //cout<<"After sstegr\n";
            //cout<<"lapinfo = "<<Lap_info<<endl;
            //cout<<"E(n-1) = "<<E1(n-1)<<endl;
            //cout<<"neigen = "<<neigen<<endl;
            //cout<<"Norm(UtU-1) = "<<Norm(U1.adjoint()*U1-1.F)<<endl;
            //cout<<"UtU.diag() = "<<(U1.adjoint()*U1).diag()<<endl;
            float nantest = U1.linearView().sumElements();
            int badcol = -1;
            for(int j=0;j<n;++j) {
                if (U1.col(j).sumElements() < TMV_Epsilon<float>()) badcol = j;
            }
            //cout<<"nantest = "<<nantest<<endl;
            //cout<<"badcol = "<<badcol<<endl;
            if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n ||
                badcol >= 0 || isNaN(nantest)) {
                std::ostringstream ss;
                ss << "Error in LAPACK function sstegr: ";
                if (Lap_info > 0) 
                    ss << "Returned info = "<<Lap_info<<".  ";
                else if (E1(n-1) > 0.F) 
                    ss << "Returned non-zero E(n) = "<<E1(n-1)<<".  ";
                else if (neigen < n) 
                    ss << "Only found "<<neigen<<" / "<<n<<" eigenvectors.  ";
                else if (badcol >= 0) 
                    ss << "U found to be not unitary: U.col("<<badcol<<") = "<<
                        U1.col(badcol)<<".  ";
                else 
                    ss <<  "NaN found in eigenvector matrix.  ";
                ss << "Calling sstedc instead.";
                TMV_Warning(ss.str());
                return AltLapEigenFromTridiagonal(U,D,E);
            }
            LAP_Results(Lap_info,"sstegr");
            D = Dout;
            U *= U1;
        } else {
            return AltLapEigenFromTridiagonal(U,D,E);
        }
    }
#endif // FLOAT
#endif // NOSTEGR
#endif // LAP

    template <class T> 
    void EigenFromTridiagonal(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E)
    {
        //std::cout<<"Start EigenFromTridiag"<<std::endl;

        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U.cptr()) {
            TMVAssert(U.colsize() == U.rowsize());
            TMVAssert(U.rowsize() == D.size());
            TMVAssert(U.ct()==NonConj);
        }

        if (D.size() == 0) return;

#ifdef XDEBUG
        std::cout<<"Start EigenFromTridiag\n";
        std::cout<<"D = "<<D<<std::endl;
        std::cout<<"E = "<<E<<std::endl;
        if (U.cptr()) { std::cout<<"U = "<<U<<std::endl; }
        Matrix<T> A0(D.size(),D.size());
        Vector<RT> D0(D);
        Vector<RT> E0(E);
        if (U.cptr()) {
            Matrix<T> EDE(D.size(),D.size(),T(0));
            EDE.diag(-1) = E;
            EDE.diag(0) = D;
            EDE.diag(1) = E;
            A0 = U * EDE * U.adjoint();
        }
#ifdef LAP
        Vector<RT> D2 = D;
        Vector<RT> E2 = E;
        Matrix<T> U2(D.size(),D.size());
        if (U.cptr()) { 
            U2 = U;
            NonLapEigenFromTridiagonal<T>(U2.view(),D2.view(),E2.view());
        } else {
            NonLapEigenFromTridiagonal<T>(U,D2.view(),E2.view());
        }
#endif // LAP
#endif // XDEBUG

        // Before running the normal algorithms, rescale D,E by the maximum
        // value to help avoid overflow and underflow.
        RT scale = TMV_MAX(D.maxAbs2Element(),E.maxAbs2Element());
        if (TMV_Underflow(scale)) {
            // Hopeless case.  Just zero out D,E and call it done.
            D.setZero();
            E.setZero();
            return;
        }
        D /= scale;
        E /= scale;

        //std::cout<<"After scaling"<<std::endl;

#ifdef LAP
        LapEigenFromTridiagonal(U,D,E);
#else 
        NonLapEigenFromTridiagonal(U,D,E);
#endif

        //std::cout<<"After LAP"<<std::endl;

        // Now A = U * D * Ut
        // Technically, singular values should be positive, but we allow them
        // to be negative, since these are the eigenvalues of A - no sense
        // killing that.  Also, to make them positive, we'd have to break the
        // V = U relationship.  So just keep that in mind later when we use S.

        // Sort output singular values by absolute value:
        AlignedArray<ptrdiff_t> sortp(D.size());
        D.sort(sortp.get(),Descend,AbsComp);
        if (U.cptr()) U.permuteCols(sortp.get());

        // Now undo the scaling
        D *= scale;

        //std::cout<<"After undo scaling"<<std::endl;

#ifdef XDEBUG
        if (U.cptr()) {
            cout<<"Done EigenFromTridiag: Norm(U) = "<<Norm(U)<<endl;
            cout<<"D = "<<D<<endl;
            Matrix<T> UDU = U * DiagMatrixViewOf(D)*U.adjoint();
            cout<<"Norm(UDUt) = "<<Norm(UDU)<<endl;
            cout<<"Norm(UDUt-A0) = "<<Norm(UDU-A0)<<endl;
            cout<<"Norm(U) = "<<Norm(U)<<endl;
            cout<<"Norm(UD) = "<<Norm(U*DiagMatrixViewOf(D))<<endl;
            cout<<"Norm(A0U) = "<<Norm(A0*U)<<endl;
            cout<<"Norm(UD-A0U) = "<<Norm(U*DiagMatrixViewOf(D)-A0*U)<<endl;
            cout<<"Norm(UUt-1) = "<<Norm(U*U.adjoint()-T(1))<<endl;
            cout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
            cout<<"UUt.diag = "<<U*(U.adjoint()).diag()<<endl;
            cout<<"UtU.diag = "<<(U.adjoint()*U).diag()<<endl;
            if (!(Norm(UDU-A0) < THRESH*Norm(A0))) {
                cerr<<"EigenFromTridiagonal:\n";
                cerr<<"D = "<<D0<<endl;
                cerr<<"E = "<<E0<<endl;
                cerr<<"Done: D = "<<D<<endl;
                cerr<<"U = "<<U<<endl;
#ifdef LAP
                cerr<<"U2 = "<<U2<<endl;
                cerr<<"diff = "<<(U-U2)<<endl;
#endif
                cerr<<"Norm(UDU-A0) = "<<Norm(UDU-A0)<<endl;
                cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
                abort();
            }
        }
#endif // XDEBUG
    }

    template <class T> 
    void UnsortedHermEigen(
        MatrixView<T> U, VectorView<RT> SS)
    {
#ifdef XDEBUG
        Matrix<T> A0 = U;
        A0.upperTri() = A0.lowerTri().adjoint();
#endif
        // Decompose Hermitian A (input as lower tri of U) into U S Ut
        // where S is a diagonal real matrix, and U is a unitary matrix.
        // U,S are N x N
        TMVAssert(SS.size() == U.colsize());
        TMVAssert(SS.size() == U.rowsize());

        if (U.isconj()) {
            UnsortedHermEigen(U.conjugate(),SS);
            U.conjugateSelf();
            return;
        }

        TMVAssert(U.ct() == NonConj);

        const ptrdiff_t N = U.colsize();
        if (N == 0) return;

        // First we reduce A to tridiagonal form: A = U * T * Ut
        // using a series of Householder transformations.
        // The diagonal of the Tridiagonal Matrix T is stored in D.
        // The subdiagonal is stored in E.
        Vector<RT> E(N-1);
#ifdef LAP
        Vector<T> Ubeta(N);
#else
        Vector<T> Ubeta(N-1);
#endif
        T signdet(0);
        Tridiagonalize(
            HermMatrixViewOf(U,Lower),Ubeta.view(),SS,E.view(),signdet);

        // Now U stores Householder vectors for U in lower diagonal columns.
        for(ptrdiff_t j=N-1;j>0;--j) U.col(j,j,N) = U.col(j-1,j,N);
        U.col(0).makeBasis(0);
        U.row(0,1,N).setZero();
        GetQFromQR(U.subMatrix(1,N,1,N),Ubeta.subVector(0,N-1));
#ifdef XDEBUG
        Matrix<T> TT = A0;
        TT.setToIdentity();
        TT.diag() = SS;
        TT.diag(1) = E.conjugate();
        TT.diag(-1) = E;
        cout<<"A0 = "<<A0<<std::endl;
        cout<<"UTU = "<<U*TT*U.adjoint()<<std::endl;
        cout<<"Norm(A0-UTU) = "<<Norm(A0-U*TT*U.adjoint())<<std::endl;
#endif

        EigenFromTridiagonal<T>(U,SS,E.view());

#ifdef XDEBUG
        Matrix<T> A2 = U * DiagMatrixViewOf(SS) * U.adjoint();
        cout<<"Norm(A0-USU) = "<<Norm(A0-A2)<<std::endl;
        if (!(Norm(A0-A2) < THRESH * Norm(U) * Norm(SS) * Norm(U))) {
            cerr<<"UnsortedHermEigen:\n";
            cerr<<"A = "<<A0<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"S = "<<SS<<endl;
            cerr<<"USUt = "<<A2<<endl;
            cerr<<"Norm(A-USUt) = "<<Norm(A0-A2)<<std::endl;
            cerr<<"cf "<<THRESH<<
                " * "<<Norm(U)<<" * "<<Norm(SS)<<" * "<<Norm(U)<<
                " = "<<THRESH*Norm(U)*Norm(SS)*Norm(U)<<std::endl;
            abort();
        }
#endif
    }

    // This version does not accumulate U
    template <class T> 
    void UnsortedEigen(SymMatrixView<T> A, VectorView<RT> SS)
    {
        TMVAssert(SS.size() == A.size());

        if (A.isupper()) return UnsortedEigen(A.transpose(),SS);
        if (A.isconj()) return UnsortedEigen(A.conjugate(),SS);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.ct() == NonConj);

        const ptrdiff_t N = A.size();
        if (N == 0) return;

        Vector<RT> E(N-1);
#ifdef LAP
        Vector<T> Ubeta(N);
#else
        Vector<T> Ubeta(N-1);
#endif
        T signdet(0);
        Tridiagonalize(A,Ubeta.view(),SS,E.view(),signdet);
        MatrixView<T> U(0,0,0,1,1,NonConj);
        EigenFromTridiagonal<T>(U,SS,E.view());
    }

    template <class T> 
    void HermSV_Decompose(MatrixView<T> U, DiagMatrixView<RT> SS)
    {
        TMVAssert(U.rowsize() == SS.size());
        TMVAssert(U.colsize() == SS.size());
        TMVAssert(U.ct() == NonConj);
        TMVAssert(SS.diag().ct() == NonConj);

#ifdef XDEBUG
        Matrix<T> A0(U);
        A0.upperTri() = A0.lowerTri().adjoint();
        std::cout<<"Start HermSV_Decompose\n";
        std::cout<<"U = "<<U<<endl;
        std::cout<<"A0 = "<<A0<<endl;
#endif

        UnsortedHermEigen(U,SS.diag());
        AlignedArray<ptrdiff_t> sortp(SS.size());
        SS.diag().sort(sortp.get(),Descend,AbsComp);
        U.permuteCols(sortp.get());

#ifdef XDEBUG
        Matrix<T> A2 = U * SS * U.adjoint();
        std::cout<<"Done HermSV_Decompose\n";
        std::cout<<"U = "<<U<<endl;
        std::cout<<"SS = "<<SS<<endl;
        std::cout<<"A0 = "<<A0<<endl;
        std::cout<<"A2 = "<<A2<<endl;
        std::cout<<"Norm(A0-A2) = "<<Norm(A0-A2)<<std::endl;
        if (!(Norm(A0-A2) < THRESH * TMV_NORM(Norm(U)) * Norm(SS))) {
            cerr<<"HermSV_Decompose:\n";
            cerr<<"A = "<<A0<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"S = "<<SS.diag()<<endl;
            cerr<<"USV = "<<A2<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    void SymSV_Decompose(
        MatrixView<T> U, DiagMatrixView<RT> SS, 
        MatrixView<T> Vt, RT& logdet, T& signdet)
    {
        TMVAssert(U.rowsize() == SS.size());
        TMVAssert(U.colsize() == SS.size());
        TMVAssert(U.ct() == NonConj);
        if (Vt.cptr()) {
            TMVAssert(Vt.rowsize() == SS.size());
            TMVAssert(Vt.colsize() == SS.size());
            TMVAssert(Vt.ct() == NonConj);
        }
        TMVAssert(SS.diag().ct() == NonConj);

#ifdef XDEBUG
        Matrix<T> A0(U);
        A0.upperTri() = A0.lowerTri().transpose();
        std::cout<<"Start SymSV_Decompose\n";
        std::cout<<"U = "<<U<<endl;
        std::cout<<"A0 = "<<A0<<endl;
#endif

        TMVAssert(isComplex(T()));
        // Decompose complex symmetric A (input as lower tri of U) into U S Vt
        // where S is a diagonal real matrix, and U,Vt are unitary matrices.
        // U,S,Vt are N x N
        // If Vt = 0, then U,Vt are not formed.  Only S,det are accurate on return.
        const ptrdiff_t N = U.colsize();
        if (N == 0) return;

        // First we reduce A to tridiagonal form: A = U * T * UT
        // using a series of Householder transformations.
        // The diagonal of the Tridiagonal Matrix T is stored in D.
        // The subdiagonal is stored in E.
        Vector<T> D(N);
        Vector<RT> E(N-1);
#ifdef LAP
        Vector<T> Ubeta(N);
#else
        Vector<T> Ubeta(N-1);
#endif
        Tridiagonalize(
            SymMatrixViewOf(U,Lower),Ubeta.view(),
            D.view(),E.view(),signdet);
        // Now U stores Householder vectors for U in lower diagonal columns.

        BandMatrix<T,ColMajor> B(N,N,1,1);
        B.diag() = D;
        B.diag(-1) = E;
        B.diag(1) = E;

        for(ptrdiff_t j=N-1;j>0;--j) U.col(j,j,N) = U.col(j-1,j,N);
        U.col(0).makeBasis(0);
        U.row(0,1,N).setZero();
        GetQFromQR(U.subMatrix(1,N,1,N),Ubeta.subVector(0,N-1));
        if (Vt.cptr()) Vt = U.transpose();
        Matrix<T,ColMajor> U1(N,N);
        Matrix<T,ColMajor> Vt1(N,N);
        SV_Decompose<T>(B,U1.view(),SS,Vt1.view(),logdet,signdet);
        U = U*U1;
        if (Vt.cptr()) Vt = Vt1*Vt;

#ifdef XDEBUG
        if (Vt.cptr()) {
            Matrix<T> A2 = U * SS * Vt;
            std::cout<<"Done SymSV_Decompose\n";
            std::cout<<"U = "<<U<<endl;
            std::cout<<"SS = "<<SS<<endl;
            std::cout<<"Vt = "<<Vt<<endl;
            std::cout<<"A0 = "<<A0<<endl;
            std::cout<<"A2 = "<<A2<<endl;
            std::cout<<"Norm(A0-A2) = "<<Norm(A0-A2)<<std::endl;
            if (!(Norm(A0-A2) < THRESH * Norm(U) * Norm(SS) * Norm(Vt))) {
                cerr<<"SymSV_Decompose:\n";
                cerr<<"A = "<<A0<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"S = "<<SS.diag()<<endl;
                cerr<<"Vt = "<<Vt<<endl;
                cerr<<"USVt = "<<A2<<endl;
                abort();
            }
        }
#endif
    }

    // This version does not accumulate U or Vt
    template <class T> 
    void SV_Decompose(SymMatrixView<T> A, DiagMatrixView<RT> SS)
    {
        TMVAssert(SS.size() == A.size());

        if (A.isherm()) {
            UnsortedEigen(A,SS.diag());
            for(ptrdiff_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) {
                SS(i) = -SS(i);
            }
            SS.diag().sort(Descend);
        } else {
            TMVAssert(isComplex(T()));
            const ptrdiff_t N = A.size();
            if (N == 0) return;
            Vector<T> D(N);
            Vector<RT> E(N-1);
#ifdef LAP
            Vector<T> Ubeta(N);
#else
            Vector<T> Ubeta(N-1);
#endif
            T signdet(0);
            Tridiagonalize(A,Ubeta.view(),D.view(),E.view(),signdet);

            BandMatrix<T,ColMajor> B(N,N,1,1);
            B.diag() = D;
            B.diag(-1) = E;
            B.diag(1) = E;

            RT logdet(0);
            MatrixView<T> U(0,0,0,1,1,NonConj);
            MatrixView<T> Vt(0,0,0,1,1,NonConj);
            SV_Decompose<T>(B,U,SS,Vt,logdet,signdet);
        }
    }

    template <class T> 
    void Eigen(const GenSymMatrix<T>& A, MatrixView<T> U, VectorView<RT> SS)
    {
        TMVAssert(A.isherm());
        TMVAssert(A.size() == SS.size());
        TMVAssert(A.size() == U.colsize());
        TMVAssert(A.size() == U.rowsize());

        if (U.isconj()) Eigen(A.conjugate(),U.conjugate(),SS);
        else {
            U.lowerTri() = A.lowerTri();
            UnsortedHermEigen(U,SS);
            AlignedArray<ptrdiff_t> sortp(A.size());
            SS.sort(sortp.get(),Ascend);
            U.permuteCols(sortp.get());
        }
    }

    template <class T> 
    void Eigen(const GenSymMatrix<T>& A, VectorView<RT> SS)
    {
        TMVAssert(A.isherm());
        TMVAssert(A.size() == SS.size());

        HermMatrix<T,Lower|ColMajor> A2 = A;
        UnsortedEigen(A2.view(),SS);
        SS.sort(Ascend);
    }

    template <class T> 
    void SV_Decompose(
        const GenSymMatrix<T>& A, MatrixView<T> U,
        DiagMatrixView<RT> SS, MatrixView<T> Vt)
    {
        TMVAssert(A.size() == U.colsize());
        TMVAssert(A.size() == U.rowsize());
        TMVAssert(A.size() == SS.size());
        TMVAssert(A.size() == Vt.colsize());
        TMVAssert(A.size() == Vt.rowsize());

        if (U.isconj()) {
            if (Vt.isconj()) {
                SV_Decompose(A.conjugate(),U.conjugate(),SS,Vt.conjugate());
            } else {
                SV_Decompose(A.conjugate(),U.conjugate(),SS,Vt);
                Vt.conjugateSelf();
            }
        } else {
            if (Vt.isconj()) {
                SV_Decompose(A,U,SS,Vt.conjugate());
                Vt.conjugateSelf();
            } else {
                U.lowerTri() = A.lowerTri();
                if (A.isherm()) {
                    HermSV_Decompose<T>(U,SS);
                    Vt = U.adjoint();
                    for(ptrdiff_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) {
                        SS(i) = -SS(i);
                        Vt.row(i) = -Vt.row(i);
                    }
                } else {
                    RT ld(0);
                    T d(0);
                    SymSV_Decompose<T>(U,SS,Vt,ld,d);
                }
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenSymMatrix<T>& A, MatrixView<T> U, DiagMatrixView<RT> SS)
    {
        TMVAssert(A.size() == U.colsize());
        TMVAssert(A.size() == U.rowsize());
        TMVAssert(A.size() == SS.size());

        if (U.isconj()) SV_Decompose(A.conjugate(),U.conjugate(),SS);
        else {
            U.lowerTri() = A.lowerTri();
            if (A.isherm()) {
                HermSV_Decompose<T>(U,SS);
                for(ptrdiff_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) {
                    SS(i) = -SS(i);
                }
            } else {
                RT ld(0);
                T d(0);
                MatrixView<T> Vt(0,0,0,1,1,NonConj);
                SymSV_Decompose<T>(U,SS,Vt,ld,d);
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenSymMatrix<T>& A, DiagMatrixView<RT> SS, MatrixView<T> Vt)
    {
        TMVAssert(A.size() == SS.size());
        TMVAssert(A.size() == Vt.colsize());
        TMVAssert(A.size() == Vt.rowsize());

        if (A.isherm()) SV_Decompose(A,Vt.adjoint(),SS);
        else SV_Decompose(A,Vt.transpose(),SS);
    }

    template <class T> 
    void PolarDecompose(MatrixView<T> U, SymMatrixView<T> P)
    {
        // Decompose A = UP
        // A is input in the place of U.
        //
        // TODO: This isn't the most efficient way to do this.
        // There is an algorithm from Higham etal (2003) that is supposedly
        // significantly faster. 
        // They iterate the process:
        // A <- A/5 + 8A(5AtA + 7 - 16(5AtA+3)^-1)^-1
        // which leads to A = U.  Then P = UtA.

        // The easier (but slower) algorithm is:
        // A = W S Vt
        // U = W Vt
        // P = V S Vt
#ifdef XDEBUG
        Matrix<T> A0 = U;
        std::cout<<"Start PolarDecompose:\n";
        std::cout<<"U = "<<TMV_Text(U)<<std::endl;
        std::cout<<"P = "<<TMV_Text(P)<<std::endl;
        //std::cout<<"A0 = "<<A0<<std::endl;
#endif
        Matrix<T> Vt(U.rowsize(),U.rowsize());
        DiagMatrix<RT> S(U.rowsize());
        SV_Decompose(U.view(),S.view(),Vt.view(),true);
        //std::cout<<"S = "<<S.diag()<<std::endl;
        RT thresh = TMV_Epsilon<T>()*S.size()*S(0);
        for(ptrdiff_t i=0;i<S.size();i++) if (S(i) < thresh) S(i) = RT(0);
        //std::cout<<"S => "<<S.diag()<<std::endl;
        U *= Vt;
        Matrix<T> VS = Vt.adjoint() * S;
        SymMultMM<false>(T(1),VS,Vt,P);
#ifdef XDEBUG
        Matrix<T> A2 = U*P;
        std::cout<<"Done PolarDecompose:\n";
        std::cout<<"Norm(A0 - USVt) = "<<Norm(A0-U*S*Vt)<<std::endl;
        std::cout<<"Norm(A2-A0) = "<<Norm(A2-A0)<<std::endl;
        if (!(Norm(A2-A0) < THRESH*Norm(A0))) {
            cerr<<"PolarDecompose "<<TMV_Text(U)<<"  "<<A0<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
            cerr<<"P = "<<P<<endl;
            cerr<<"UP = "<<A2<<endl;
            cerr<<"Norm(A2-A0) = "<<Norm(A2-A0)<<endl;
            cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    void PolarDecompose(
        const GenBandMatrix<T>& A, MatrixView<T> U, SymMatrixView<T> P)
    {
        Matrix<T> Vt(A.rowsize(),A.rowsize());
        DiagMatrix<RT> S(A.rowsize());
        SV_Decompose(A,U.view(),S.view(),Vt.view());
        U *= Vt;
        Matrix<T> VS = Vt.adjoint() * S;
        SymMultMM<false>(T(1),VS,Vt,P);
#ifdef XDEBUG
        std::cout<<"Band PolarDecompose "<<TMV_Text(A)<<"  "<<A<<endl;
        //std::cout<<"U = "<<U<<endl;
        std::cout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
        //std::cout<<"P = "<<P<<endl;
        Matrix<T> A2 = U*P;
        //std::cout<<"UP = "<<A2<<endl;
        std::cout<<"Norm(A2-A0) = "<<Norm(A2-A)<<endl;
        std::cout<<"Norm(A0) = "<<Norm(A)<<endl;
        if (!(Norm(A2-A) < THRESH*Norm(A))) {
            cerr<<"Band PolarDecompose "<<TMV_Text(A)<<"  "<<A<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
            cerr<<"P = "<<P<<endl;
            cerr<<"UP = "<<A2<<endl;
            cerr<<"Norm(A2-A0) = "<<Norm(A2-A)<<endl;
            cerr<<"Norm(A0) = "<<Norm(A)<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    void SquareRoot(SymMatrixView<T> A)
    {
        TMVAssert(A.isherm());
        // A -> A^1/2
        //
        // TODO: Again, there are supposedly faster algorithms than this.
        //
        // A = V D Vt
        // A = V D^1/2 Vt
        Matrix<T> V(A.size(),A.size());
        DiagMatrix<RT> D(A.size());
        Eigen(A,V.view(),D.diag());
        for(ptrdiff_t i=0;i<A.size();i++) {
            if (D(i) < RT(0)) {
#ifdef NOTHROW
                std::cerr<<"Non PosDef SymMatrix found in SquareRoot\n"; 
                exit(1); 
#else
                throw NonPosDef("in SymMatrix SquareRoot");
#endif
            }
            D(i) = TMV_SQRT(D(i));
        }
        Matrix<T> DVt = D*V.adjoint();
        SymMultMM<false>(T(1),V,DVt,A);
    }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymSVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


