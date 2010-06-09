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
#define THRESH 1.e-10
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
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
    { 
        EigenFromTridiagonal_DC(U,D,E,false);
        //EigenFromTridiagonal_QR(U,D,E); 
    }

    template <class T> 
    static void NonLapEigenFromTridiagonal(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        if (U) {
            TMVAssert(U->isSquare());
            TMVAssert(U->rowsize() == D.size());
        }

        const int N = D.size();

        // E may have zeros on entry.
        // This routine finds subproblems in which E is fully non-zero.

        // First chop small elements of D,E
        HermTridiagonalChopSmallElements(D,E);

        // Find sub-problems to solve:
        for(int q = N-1; q>0; ) {
            if (E(q-1) == T(0)) --q;
            else {
                int p=q-1;
                while (p > 0 && (E(p-1) != T(0))) --p; 
                // Set p such that E(p-1) = 0 and 
                // all E(i) with p<=i<q are non-zero.
                if (U) {
                    DoEigenFromTridiagonal_NZ<T>(
                        U->colRange(p,q+1),D.subVector(p,q+1),E.subVector(p,q));
                } else {
                    DoEigenFromTridiagonal_NZ<T>(
                        0,D.subVector(p,q+1),E.subVector(p,q));
                }
                q = p;
            }
        }
    }

#ifdef LAP 
    template <class T> 
    static inline void LapEigenFromTridiagonal(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
    { NonLapEigenFromTridiagonal(U,D,E); }

#ifdef INST_DOUBLE
    static void AltLapEigenFromTridiagonal(
        MVP<double> U, const VectorView<double>& D, 
        const VectorView<double>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
        }
        int n = D.size();
        if (U) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n);
            int ldu = U1.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1+(4+n)*n;
            int liwork = 3+5*n;
            auto_array<double> work(new double[lwork]);
            auto_array<int> iwork(new int[liwork]);
#else
            int lwork = -1;
            int liwork = -1;
            auto_array<double> work(new double[1]);
            auto_array<int> iwork(new int[1]);
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.reset(new double[lwork]);
            iwork.reset(new int[liwork]);
#endif
#endif
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results("dstedc");
#else
            LAP_Results(int(work[0]),n,n,lwork,"dstedc");
#endif
            *U *= U1;
        } else {
            char c = 'N';
            int ldu = n;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1;
            int liwork = 1;
            auto_array<double> work(new double[lwork]);
            auto_array<int> iwork(new int[liwork]);
#else
            int lwork = -1;
            int liwork = -1;
            auto_array<double> work(new double[1]);
            auto_array<int> iwork(new int[1]);
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.reset(new double[lwork]);
            iwork.reset(new int[liwork]);
#endif
#endif
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results("dstedc");
#else
            LAP_Results(int(work[0]),n,n,lwork,"dstedc");
#endif
        }
    }
    static void AltLapEigenFromTridiagonal(
        MVP<std::complex<double> > U, const VectorView<double>& D, 
        const VectorView<double>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
        }
        int n = D.size();
        if (U) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n);
            int ldu = U1.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1+(4+n)*n;
            int liwork = 3+5*n;
            auto_array<double> work(new double[lwork]);
            auto_array<int> iwork(new int[liwork]);
#else
            int lwork = -1;
            int liwork = -1;
            auto_array<double> work(new double[1]);
            auto_array<int> iwork(new int[1]);
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.reset(new double[lwork]);
            iwork.reset(new int[liwork]);
#endif
#endif
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results("dstedc");
#else
            LAP_Results(int(work[0]),n,n,lwork,"dstedc");
#endif
            *U *= U1;
        } else {
            char c = 'N';
            int ldu = n;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1;
            int liwork = 1;
            auto_array<double> work(new double[lwork]);
            auto_array<int> iwork(new int[liwork]);
#else
            int lwork = -1;
            int liwork = -1;
            auto_array<double> work(new double[1]);
            auto_array<int> iwork(new int[1]);
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.reset(new double[lwork]);
            iwork.reset(new int[liwork]);
#endif
#endif
            LAPNAME(dstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results("dstedc");
#else
            LAP_Results(int(work[0]),n,n,lwork,"dstedc");
#endif
        }
    }
#endif
#ifdef INST_FLOAT
    static void AltLapEigenFromTridiagonal(
        MVP<float> U, const VectorView<float>& D, 
        const VectorView<float>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
        }
        int n = D.size();
        if (U) {
            char c = 'I';
            Matrix<float,ColMajor> U1(n,n);
            int ldu = U1.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1+(4+n)*n;
            int liwork = 3+5*n;
            auto_array<float> work(new float[lwork]);
            auto_array<int> iwork(new int[liwork]);
#else
            int lwork = -1;
            int liwork = -1;
            auto_array<float> work(new float[1]);
            auto_array<int> iwork(new int[1]);
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.reset(new float[lwork]);
            iwork.reset(new int[liwork]);
#endif
#endif
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results("sstedc");
#else
            LAP_Results(int(work[0]),n,n,lwork,"sstedc");
#endif
            *U *= U1;
        } else {
            char c = 'N';
            int ldu = n;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1;
            int liwork = 1;
            auto_array<float> work(new float[lwork]);
            auto_array<int> iwork(new int[liwork]);
#else
            int lwork = -1;
            int liwork = -1;
            auto_array<float> work(new float[1]);
            auto_array<int> iwork(new int[1]);
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.reset(new float[lwork]);
            iwork.reset(new int[liwork]);
#endif
#endif
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results("sstedc");
#else
            LAP_Results(int(work[0]),n,n,lwork,"sstedc");
#endif
        }
    }
    static void AltLapEigenFromTridiagonal(
        MVP<std::complex<float> > U, const VectorView<float>& D, 
        const VectorView<float>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
        }
        int n = D.size();
        if (U) {
            char c = 'I';
            Matrix<float,ColMajor> U1(n,n);
            int ldu = U1.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1+(4+n)*n;
            int liwork = 3+5*n;
            auto_array<float> work(new float[lwork]);
            auto_array<int> iwork(new int[liwork]);
#else
            int lwork = -1;
            int liwork = -1;
            auto_array<float> work(new float[1]);
            auto_array<int> iwork(new int[1]);
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.reset(new float[lwork]);
            iwork.reset(new int[liwork]);
#endif
#endif
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results("sstedc");
#else
            LAP_Results(int(work[0]),n,n,lwork,"sstedc");
#endif
            *U *= U1;
        } else {
            char c = 'N';
            int ldu = n;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 1;
            int liwork = 1;
            auto_array<float> work(new float[lwork]);
            auto_array<int> iwork(new int[liwork]);
#else
            int lwork = -1;
            int liwork = -1;
            auto_array<float> work(new float[1]);
            auto_array<int> iwork(new int[1]);
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
            lwork = int(work[0]);
            liwork = iwork[0];
            work.reset(new float[lwork]);
            iwork.reset(new int[liwork]);
#endif
#endif
            LAPNAME(sstedc) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
                LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
                LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results("sstedc");
#else
            LAP_Results(int(work[0]),n,n,lwork,"sstedc");
#endif
        }
    }
#endif
#ifdef NOSTEGR
#ifdef INST_DOUBLE
    template <> 
    inline void LapEigenFromTridiagonal(
        MVP<double> U, const VectorView<double>& D, 
        const VectorView<double>& E)
    { AltLapEigenFromTridiagonal(U,D,E); }
    template <> 
    inline void LapEigenFromTridiagonal(
        MVP<std::complex<double> > U, const VectorView<double>& D, 
        const VectorView<double>& E)
    { AltLapEigenFromTridiagonal(U,D,E); }
#endif
#ifdef INST_FLOAT
    template <> 
    inline void LapEigenFromTridiagonal(
        MVP<float> U, const VectorView<float>& D, 
        const VectorView<float>& E)
    { AltLapEigenFromTridiagonal(U,D,E); }
    template <> 
    inline void LapEigenFromTridiagonal(
        MVP<std::complex<float> > U, const VectorView<float>& D, 
        const VectorView<float>& E)
    { AltLapEigenFromTridiagonal(U,D,E); }
#endif
#else // normal stegr implementaion
#ifdef INST_DOUBLE
    template <> 
    void LapEigenFromTridiagonal(
        MVP<double> U, const VectorView<double>& D, 
        const VectorView<double>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
        }
        int n = D.size();
        if (U) {
            char c1 = 'V';
            char c2 = 'A';
            double junk = 0;
            int ijunk = 0;
            double tol = 0;
            int neigen;
            Vector<double> Din = D;
            Vector<double> E1(n);
            E1.subVector(0,n-1) = E;  E1(n-1) = 0.;
            Vector<double> Dout(n);
            Matrix<double,ColMajor> U1(n,n);
            int ldu = U1.stepj();
            auto_array<int> isuppz(new int[2*n]);
#ifndef LAPNOWORK
            int lwork = 18*n;
            auto_array<double> work(new double[lwork]);
            int liwork = 10*n;
            auto_array<int> iwork(new int[liwork]);
#endif
            LAPNAME(dstegr) (
                LAPCM LAPV(c1),LAPV(c2),LAPV(n),
                LAPP(Din.ptr()),LAPP(E1.ptr()),
                LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk),LAPV(tol),
                LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
                LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
            // MJ: The dstegr algorithm has a bug wherein it sometimes fails
            //     to finish, in which case the output values are incorrect.
            //     It informs that this has happened via the info variable.
            //     Specifically, it sets info = 2 in this case.
            //     We test for all info > 0 just to be sure.
            //     When this happens, we call dstedc instead.
            //
            //     In addition to the above problem, the dstegr routine seems to
            //     not be very careful about nan issues.  So we also check for nan's
            //     in U1, and call dstedc if there are any.
            double nantest = U1.linearView().sumElements();
            if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n || isNaN(nantest)) {
                std::ostringstream ss;
                ss << "Error in LAPACK function dstegr: ";
                if (Lap_info > 0) 
                    ss << "Returned info = "<<Lap_info<<".  ";
                else if (E1(n-1) > 0.F) 
                    ss << "Returned non-zero E(n) = "<<E1(n-1)<<".  ";
                else if (neigen < n) 
                    ss << "Only found "<<neigen<<" / "<<n<<" eigenvectors.  ";
                else 
                    ss <<  "NaN found in eigenvector matrix.  ";
                ss << "Calling dstedc instead.";
                TMV_Warning(ss.str());
                return AltLapEigenFromTridiagonal(U,D,E);
            }
            LAP_Results("dstegr");
            D = Dout;
            *U *= U1;
        } else {
            return AltLapEigenFromTridiagonal(U,D,E);
        }
    }
    template <> 
    void LapEigenFromTridiagonal(
        MVP<std::complex<double> > U, const VectorView<double>& D, 
        const VectorView<double>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
        }
        int n = D.size();
        if (U) {
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
            auto_array<int> isuppz(new int[2*n]);
#ifndef LAPNOWORK
            int lwork = 18*n;
            auto_array<double> work(new double[lwork]);
            int liwork = 10*n;
            auto_array<int> iwork(new int[liwork]);
#endif
            LAPNAME(dstegr) (
                LAPCM LAPV(c1),LAPV(c2),LAPV(n),
                LAPP(Din.ptr()),LAPP(E1.ptr()),
                LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk),LAPV(tol),
                LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
                LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
            double nantest = U1.linearView().sumElements();
            if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n || isNaN(nantest)) {
                std::ostringstream ss;
                ss << "Error in LAPACK function dstegr: ";
                if (Lap_info > 0) 
                    ss << "Returned info = "<<Lap_info<<".  ";
                else if (E1(n-1) > 0.F) 
                    ss << "Returned non-zero E(n) = "<<E1(n-1)<<".  ";
                else if (neigen < n) 
                    ss << "Only found "<<neigen<<" / "<<n<<" eigenvectors.  ";
                else 
                    ss <<  "NaN found in eigenvector matrix.  ";
                ss << "Calling dstedc instead.";
                TMV_Warning(ss.str());
                return AltLapEigenFromTridiagonal(U,D,E);
            }
            LAP_Results("dstegr");
            D = Dout;
            *U *= U1;
        } else {
            return AltLapEigenFromTridiagonal(U,D,E);
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapEigenFromTridiagonal(
        MVP<float> U, const VectorView<float>& D, 
        const VectorView<float>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
        }
        int n = D.size();
        if (U) {
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
            auto_array<int> isuppz(new int[2*n]);
#ifndef LAPNOWORK
            int lwork = 18*n;
            auto_array<float> work(new float[lwork]);
            int liwork = 10*n;
            auto_array<int> iwork(new int[liwork]);
#endif
            LAPNAME(sstegr) (
                LAPCM LAPV(c1),LAPV(c2),LAPV(n),
                LAPP(Din.ptr()),LAPP(E1.ptr()),
                LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk), LAPV(tol),
                LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
                LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork)
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
            float nantest = U1.linearView().sumElements();
            if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n || isNaN(nantest)) {
                std::ostringstream ss;
                ss << "Error in LAPACK function sstegr: ";
                if (Lap_info > 0) 
                    ss << "Returned info = "<<Lap_info<<".  ";
                else if (E1(n-1) > 0.F) 
                    ss << "Returned non-zero E(n) = "<<E1(n-1)<<".  ";
                else if (neigen < n) 
                    ss << "Only found "<<neigen<<" / "<<n<<" eigenvectors.  ";
                else 
                    ss <<  "NaN found in eigenvector matrix.  ";
                ss << "Calling sstedc instead.";
                TMV_Warning(ss.str());
                return AltLapEigenFromTridiagonal(U,D,E);
            }
            LAP_Results("sstegr");
            D = Dout;
            *U *= U1;
        } else {
            return AltLapEigenFromTridiagonal(U,D,E);
        }
    }
    template <> 
    void LapEigenFromTridiagonal(
        MVP<std::complex<float> > U, const VectorView<float>& D, 
        const VectorView<float>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        if (U) {
            TMVAssert(U->colsize() >= U->rowsize());
            TMVAssert(U->rowsize() == D.size());
        }
        int n = D.size();
        if (U) {
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
            auto_array<int> isuppz(new int[2*n]);
#ifndef LAPNOWORK
            int lwork = 18*n;
            auto_array<float> work(new float[lwork]);
            int liwork = 10*n;
            auto_array<int> iwork(new int[liwork]);
#endif
            LAPNAME(sstegr) (
                LAPCM LAPV(c1),LAPV(c2),LAPV(n),
                LAPP(Din.ptr()),LAPP(E1.ptr()),
                LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk), LAPV(tol),
                LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
                LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork) 
                LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
            float nantest = U1.linearView().sumElements();
            if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n || isNaN(nantest)) {
                std::ostringstream ss;
                ss << "Error in LAPACK function sstegr: ";
                if (Lap_info > 0) 
                    ss << "Returned info = "<<Lap_info<<".  ";
                else if (E1(n-1) > 0.F) 
                    ss << "Returned non-zero E(n) = "<<E1(n-1)<<".  ";
                else if (neigen < n) 
                    ss << "Only found "<<neigen<<" / "<<n<<" eigenvectors.  ";
                else 
                    ss <<  "NaN found in eigenvector matrix.  ";
                ss << "Calling sstedc instead.";
                TMV_Warning(ss.str());
                return AltLapEigenFromTridiagonal(U,D,E);
            }
            LAP_Results("sstegr");
            D = Dout;
            *U *= U1;
        } else {
            return AltLapEigenFromTridiagonal(U,D,E);
        }
    }
#endif // FLOAT
#endif // NOSTEGR
#endif // LAP

    template <class T> 
    void EigenFromTridiagonal(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        if (U) {
            TMVAssert(U->colsize() == U->rowsize());
            TMVAssert(U->rowsize() == D.size());
            TMVAssert(U->ct()==NonConj);
        }

        if (D.size() == 0) return;

#ifdef XDEBUG
        Matrix<T> A0(D.size(),D.size());
        Vector<RT> D0(D);
        Vector<RT> E0(E);
        if (U) {
            Matrix<T> EDE(D.size(),D.size(),T(0));
            EDE.diag(-1) = E;
            EDE.diag(0) = D;
            EDE.diag(1) = E;
            A0 = *U * EDE * U->adjoint();
        }
#ifdef LAP
        Vector<RT> D2 = D;
        Vector<RT> E2 = E;
        Matrix<T> U2(D.size(),D.size());
        if (U) { 
            U2 = *U;
            NonLapEigenFromTridiagonal<T>(U2.view(),D2.view(),E2.view());
        } else {
            NonLapEigenFromTridiagonal<T>(0,D2.view(),E2.view());
        }
#endif // LAP
#endif // XDEBUG

        // Before running the normal algorithms, rescale D,E by the maximum
        // value to help avoid overflow and underflow.
        RT scale = TMV_MAX(D.maxAbsElement(),E.maxAbsElement());
        if (scale * TMV_Epsilon<T>() == RT(0)) {
            // Hopeless case.  Just zero out D,E and call it done.
            D.setZero();
            E.setZero();
            return;
        }
        D /= scale;
        E /= scale;

#ifdef LAP
        LapEigenFromTridiagonal(U,D,E);
#else 
        NonLapEigenFromTridiagonal(U,D,E);
#endif

        // Now A = U * D * Ut
        // Technically, singular values should be positive, but we allow them
        // to be negative, since these are the eigenvalues of A - no sense
        // killing that.  Also, to make them positive, we'd have to break the
        // V = Ut relationship.  So just keep that in mind later when we use S.

        // Sort output singular values by absolute value:
        auto_array<int> sortp(new int[D.size()]);
        D.sort(sortp.get(),Descend,AbsComp);
        if (U) U->permuteCols(sortp.get());

        // Now undo the scaling
        D *= scale;

#ifdef XDEBUG
        if (U) {
            //cout<<"Done EigenFromTridiag: Norm(U) = "<<Norm(*U)<<endl;
            //cout<<"D = "<<D<<endl;
            Matrix<T> UDU = *U * DiagMatrixViewOf(D)*U->adjoint();
            //cout<<"Norm(UDUt) = "<<Norm(UDU)<<endl;
            //cout<<"Norm(UDUt-A0) = "<<Norm(UDU-A0)<<endl;
            //cout<<"Norm(U) = "<<Norm(*U)<<endl;
            //cout<<"Norm(UD) = "<<Norm(*U*DiagMatrixViewOf(D))<<endl;
            //cout<<"Norm(A0U) = "<<Norm(A0*(*U))<<endl;
            //cout<<"Norm(UD-A0U) = "<<Norm((*U)*DiagMatrixViewOf(D)-A0*(*U))<<endl;
            if (!(Norm(UDU-A0) < THRESH*Norm(A0))) {
                cerr<<"EigenFromTridiagonal:\n";
                cerr<<"D = "<<D0<<endl;
                cerr<<"E = "<<E0<<endl;
                cerr<<"Done: D = "<<D<<endl;
                cerr<<"U = "<<*U<<endl;
#ifdef LAP
                cerr<<"U2 = "<<U2<<endl;
                cerr<<"diff = "<<(*U-U2)<<endl;
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
        const MatrixView<T>& U, const VectorView<RT>& SS)
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

        const int N = U.colsize();
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
        for(int j=N-1;j>0;--j) U.col(j,j,N) = U.col(j-1,j,N);
        U.col(0).makeBasis(0);
        U.row(0,1,N).setZero();
        GetQFromQR(U.subMatrix(1,N,1,N),Ubeta.subVector(0,N-1));

        EigenFromTridiagonal<T>(U,SS,E.view());

#ifdef XDEBUG
        Matrix<T> A2 = U * DiagMatrixViewOf(SS) * U.adjoint();
        if (!(Norm(A0-A2) < THRESH * Norm(U) * Norm(SS) * Norm(U))) {
            cerr<<"UnsortedHermEigen:\n";
            cerr<<"A = "<<A0<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"S = "<<SS<<endl;
            cerr<<"USUt = "<<A2<<endl;
            abort();
        }
#endif
    }

    // This version does not accumulate U
    template <class T> 
    void UnsortedEigen(
        const SymMatrixView<T>& A, const VectorView<RT>& SS)
    {
        TMVAssert(SS.size() == A.size());

        if (A.isupper()) return UnsortedEigen(A.transpose(),SS);
        if (A.isconj()) return UnsortedEigen(A.conjugate(),SS);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.ct() == NonConj);

        const int N = A.size();
        if (N == 0) return;

        Vector<RT> E(N-1);
#ifdef LAP
        Vector<T> Ubeta(N);
#else
        Vector<T> Ubeta(N-1);
#endif
        T signdet(0);
        Tridiagonalize(A,Ubeta.view(),SS,E.view(),signdet);
        EigenFromTridiagonal<T>(0,SS,E.view());
    }

    template <class T> 
    void HermSV_Decompose(
        const MatrixView<T>& U, const DiagMatrixView<RT>& SS)
    {
        TMVAssert(U.rowsize() == SS.size());
        TMVAssert(U.colsize() == SS.size());
        TMVAssert(U.ct() == NonConj);
        TMVAssert(SS.diag().ct() == NonConj);

#ifdef XDEBUG
        Matrix<T> A0(U);
        A0.upperTri() = A0.lowerTri().adjoint();
        std::cout<<"Start SymSV_Decompose\n";
        std::cout<<"U = "<<U<<endl;
        std::cout<<"A0 = "<<A0<<endl;
#endif

        UnsortedHermEigen(U,SS.diag());
        auto_array<int> sortp(new int[SS.size()]);
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
        const MatrixView<T>& U, const DiagMatrixView<RT>& SS, 
        MVP<T> V, RT& logdet, T& signdet)
    {
        TMVAssert(U.rowsize() == SS.size());
        TMVAssert(U.colsize() == SS.size());
        TMVAssert(U.ct() == NonConj);
        if (V) {
            TMVAssert(V->rowsize() == SS.size());
            TMVAssert(V->colsize() == SS.size());
            TMVAssert(V->ct() == NonConj);
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
        // Decompose complex symmetric A (input as lower tri of U) into U S V
        // where S is a diagonal real matrix, and U,V are unitary matrices.
        // U,S,V are N x N
        // If V = 0, then U,V are not formed.  Only S,det are accurate on return.
        const int N = U.colsize();
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

        for(int j=N-1;j>0;--j) U.col(j,j,N) = U.col(j-1,j,N);
        U.col(0).makeBasis(0);
        U.row(0,1,N).setZero();
        GetQFromQR(U.subMatrix(1,N,1,N),Ubeta.subVector(0,N-1));
        if (V) *V = U.transpose();
        Matrix<T,ColMajor> U1(N,N);
        Matrix<T,ColMajor> V1(N,N);
        SV_Decompose<T>(B,U1.view(),SS,V1.view(),logdet,signdet);
        U = U*U1;
        if (V) *V = V1*(*V);

#ifdef XDEBUG
        if (V) {
            Matrix<T> A2 = U * SS * (*V);
            std::cout<<"Done SymSV_Decompose\n";
            std::cout<<"U = "<<U<<endl;
            std::cout<<"SS = "<<SS<<endl;
            std::cout<<"V = "<<*V<<endl;
            std::cout<<"A0 = "<<A0<<endl;
            std::cout<<"A2 = "<<A2<<endl;
            std::cout<<"Norm(A0-A2) = "<<Norm(A0-A2)<<std::endl;
            if (!(Norm(A0-A2) < THRESH * Norm(U) * Norm(SS) * Norm(*V))) {
                cerr<<"SymSV_Decompose:\n";
                cerr<<"A = "<<A0<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"S = "<<SS.diag()<<endl;
                cerr<<"V = "<<*V<<endl;
                cerr<<"USV = "<<A2<<endl;
                abort();
            }
        }
#endif
    }

    // This version does not accumulate U or V
    template <class T> 
    void SV_Decompose(
        const SymMatrixView<T>& A, const DiagMatrixView<RT>& SS)
    {
        TMVAssert(SS.size() == A.size());

        if (A.isherm()) {
            UnsortedEigen(A,SS.diag());
            for(size_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) {
                SS(i) = -SS(i);
            }
            SS.diag().sort(Descend);
        } else {
            TMVAssert(isComplex(T()));
            const int N = A.size();
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
            SV_Decompose<T>(B,0,SS,0,logdet,signdet);
        }
    }

    template <class T> 
    void Eigen(
        const GenSymMatrix<T>& A, const MatrixView<T>& U,
        const VectorView<RT>& SS)
    {
        TMVAssert(A.isherm());
        TMVAssert(A.size() == SS.size());
        TMVAssert(A.size() == U.colsize());
        TMVAssert(A.size() == U.rowsize());

        if (U.isconj()) Eigen(A.conjugate(),U.conjugate(),SS);
        else {
            U.lowerTri() = A.lowerTri();
            UnsortedHermEigen(U,SS);
            auto_array<int> sortp(new int[A.size()]);
            SS.sort(sortp.get(),Ascend);
            U.permuteCols(sortp.get());
        }
    }

    template <class T> 
    void Eigen(
        const SymMatrixView<T>& A, const VectorView<RT>& SS)
    {
        TMVAssert(A.isherm());
        TMVAssert(A.size() == SS.size());

        if (A.isconj())
            UnsortedEigen(A.conjugate(),SS);
        else
            UnsortedEigen(A,SS);
        SS.sort(Ascend);
    }

    template <class T> 
    void SV_Decompose(
        const GenSymMatrix<T>& A, const MatrixView<T>& U,
        const DiagMatrixView<RT>& SS, const MatrixView<T>& V)
    {
        TMVAssert(A.size() == U.colsize());
        TMVAssert(A.size() == U.rowsize());
        TMVAssert(A.size() == SS.size());
        TMVAssert(A.size() == V.colsize());
        TMVAssert(A.size() == V.rowsize());

        if (U.isconj()) {
            if (V.isconj()) {
                SV_Decompose(A.conjugate(),U.conjugate(),SS,V.conjugate());
            } else {
                SV_Decompose(A.conjugate(),U.conjugate(),SS,V);
                V.conjugateSelf();
            }
        } else {
            if (V.isconj()) {
                SV_Decompose(A,U,SS,V.conjugate());
                V.conjugateSelf();
            } else {
                U.lowerTri() = A.lowerTri();
                if (A.isherm()) {
                    HermSV_Decompose<T>(U,SS);
                    V = U.adjoint();
                    for(size_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) {
                        SS(i) = -SS(i);
                        V.row(i) = -V.row(i);
                    }
                } else {
                    RT ld(0);
                    T d(0);
                    SymSV_Decompose<T>(U,SS,V,ld,d);
                }
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenSymMatrix<T>& A,
        const MatrixView<T>& U, const DiagMatrixView<RT>& SS)
    {
        TMVAssert(A.size() == U.colsize());
        TMVAssert(A.size() == U.rowsize());
        TMVAssert(A.size() == SS.size());

        if (U.isconj()) SV_Decompose(A.conjugate(),U.conjugate(),SS);
        else {
            U.lowerTri() = A.lowerTri();
            if (A.isherm()) {
                HermSV_Decompose<T>(U,SS);
                for(size_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) {
                    SS(i) = -SS(i);
                }
            } else {
                RT ld(0);
                T d(0);
                SymSV_Decompose<T>(U,SS,0,ld,d);
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenSymMatrix<T>& A,
        const DiagMatrixView<RT>& SS, const MatrixView<T>& V)
    {
        TMVAssert(A.size() == SS.size());
        TMVAssert(A.size() == V.colsize());
        TMVAssert(A.size() == V.rowsize());

        if (A.isherm()) SV_Decompose(A,V.adjoint(),SS);
        else SV_Decompose(A,V.transpose(),SS);
    }

    template <class T> 
    void PolarDecompose(const MatrixView<T>& U,
              const SymMatrixView<T>& P)
    {
        // Decompose A = UP
        // A is input in the place of U.
        //
        // MJ: This isn't the most efficient way to do this.
        // There is an algorithm from Higham etal (2003) that is supposedly
        // significantly faster. 
        // They iterate the process:
        // A <- A/5 + 8A(5AtA + 7 - 16(5AtA+3)^-1)^-1
        // which leads to A = U.  Then P = UtA.

        // The easier (but slower) algorithm is:
        // A = W S V
        // U = W V
        // P = Vt S V
#ifdef XDEBUG
        Matrix<T> A0 = U;
#endif
        Matrix<T> V(U.rowsize(),U.rowsize());
        DiagMatrix<RT> S(U.rowsize());
        SV_Decompose(U.view(),S.view(),V.view(),true);
        RT thresh = TMV_Epsilon<T>()*S.size()*S(0);
        for(size_t i=0;i<S.size();i++) if (S(i) < thresh) S(i) = RT(0);
        U *= V;
        Matrix<T> VtS = V.adjoint() * S;
        SymMultMM<false>(T(1),VtS,V,P);
#ifdef XDEBUG
        Matrix<T> A2 = U*P;
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
    void PolarDecompose(const GenBandMatrix<T>& A,
              const MatrixView<T>& U, const SymMatrixView<T>& P)
    {
        Matrix<T> V(A.rowsize(),A.rowsize());
        DiagMatrix<RT> S(A.rowsize());
        SV_Decompose(A,U.view(),S.view(),V.view());
        U *= V;
        Matrix<T> VtS = V.adjoint() * S;
        SymMultMM<false>(T(1),VtS,V,P);
#ifdef XDEBUG
        std::cout<<"Band PolarDecompose "<<TMV_Text(A)<<"  "<<A<<endl;
        std::cout<<"U = "<<U<<endl;
        std::cout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
        //std::cout<<"P = "<<P<<endl;
        Matrix<T> A2 = U*P;
        //std::cout<<"UP = "<<A2<<endl;
        std::cout<<"Norm(A2-A0) = "<<Norm(A2-A)<<endl;
        std::cout<<"Norm(A0) = "<<Norm(A)<<endl;
        if (Norm(A2-A) > THRESH*Norm(A)) {
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
    void SquareRoot(const SymMatrixView<T>& A)
    {
        TMVAssert(A.isherm());
        // A -> A^1/2
        //
        // Again, there are supposedly faster algorithms than this.
        //
        // A = V D Vt
        // A = V D^1/2 Vt
        Matrix<T> V(A.size(),A.size());
        DiagMatrix<RT> D(A.size());
        Eigen(A,V.view(),D.diag());
        for(size_t i=0;i<A.size();i++) {
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

#define InstFile "TMV_SymSVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


