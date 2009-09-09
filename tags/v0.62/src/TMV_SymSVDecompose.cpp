///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT RealType(T)

  template <class T> static void DoEigen_From_Tridiagonal_NZ(
      MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
  { 
    Eigen_From_Tridiagonal_DC(U,D,E,false);
    //Eigen_From_Tridiagonal_QR(U,D,E); 
  }

  template <class T> static void NonLapEigen_From_Tridiagonal(
      MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
  {
    TMVAssert(D.size() == E.size()+1);
    if (U) {
      TMVAssert(U->IsSquare());
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
        // Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.
        if (U) {
          DoEigen_From_Tridiagonal_NZ<T>(U->Cols(p,q+1),D.SubVector(p,q+1),
              E.SubVector(p,q));
        } else {
          DoEigen_From_Tridiagonal_NZ<T>(0,D.SubVector(p,q+1),E.SubVector(p,q));
        }
        q = p;
      }
    }
  }

#ifdef LAP 
  template <class T> static inline void LapEigen_From_Tridiagonal(
      MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
  { NonLapEigen_From_Tridiagonal(U,D,E); }

#ifdef INST_DOUBLE
  static void AltLapEigen_From_Tridiagonal(
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
      LAPNAME(dstedc) (LAPCM LAPV(c),LAPV(n),
          LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
          LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
          LAPINFO LAP1);
      lwork = int(work[0]);
      liwork = iwork[0];
      work.reset(new double[lwork]);
      iwork.reset(new int[liwork]);
#endif
#endif
      LAPNAME(dstedc) (LAPCM LAPV(c),LAPV(n),
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
      LAPNAME(dstedc) (LAPCM LAPV(c),LAPV(n),
          LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
          LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
          LAPINFO LAP1);
      lwork = int(work[0]);
      liwork = iwork[0];
      work.reset(new double[lwork]);
      iwork.reset(new int[liwork]);
#endif
#endif
      LAPNAME(dstedc) (LAPCM LAPV(c),LAPV(n),
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
  static void AltLapEigen_From_Tridiagonal(
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
      LAPNAME(dstedc) (LAPCM LAPV(c),LAPV(n),
          LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
          LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
          LAPINFO LAP1);
      lwork = int(work[0]);
      liwork = iwork[0];
      work.reset(new double[lwork]);
      iwork.reset(new int[liwork]);
#endif
#endif
      LAPNAME(dstedc) (LAPCM LAPV(c),LAPV(n),
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
      LAPNAME(dstedc) (LAPCM LAPV(c),LAPV(n),
          LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
          LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
          LAPINFO LAP1);
      lwork = int(work[0]);
      liwork = iwork[0];
      work.reset(new double[lwork]);
      iwork.reset(new int[liwork]);
#endif
#endif
      LAPNAME(dstedc) (LAPCM LAPV(c),LAPV(n),
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
  static void AltLapEigen_From_Tridiagonal(
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
      LAPNAME(sstedc) (LAPCM LAPV(c),LAPV(n),
          LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
          LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
          LAPINFO LAP1);
      lwork = int(work[0]);
      liwork = iwork[0];
      work.reset(new float[lwork]);
      iwork.reset(new int[liwork]);
#endif
#endif
      LAPNAME(sstedc) (LAPCM LAPV(c),LAPV(n),
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
      LAPNAME(sstedc) (LAPCM LAPV(c),LAPV(n),
          LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
          LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
          LAPINFO LAP1);
      lwork = int(work[0]);
      liwork = iwork[0];
      work.reset(new float[lwork]);
      iwork.reset(new int[liwork]);
#endif
#endif
      LAPNAME(sstedc) (LAPCM LAPV(c),LAPV(n),
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
  static void AltLapEigen_From_Tridiagonal(
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
      LAPNAME(sstedc) (LAPCM LAPV(c),LAPV(n),
          LAPP(D.ptr()),LAPP(E.ptr()),LAPP(U1.ptr()),LAPV(ldu)
          LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
          LAPINFO LAP1);
      lwork = int(work[0]);
      liwork = iwork[0];
      work.reset(new float[lwork]);
      iwork.reset(new int[liwork]);
#endif
#endif
      LAPNAME(sstedc) (LAPCM LAPV(c),LAPV(n),
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
      LAPNAME(sstedc) (LAPCM LAPV(c),LAPV(n),
          LAPP(D.ptr()),LAPP(E.ptr()),0,LAPV(ldu)
          LAPWK(work.get()) LAPVWK(lwork) LAPWK(iwork.get()) LAPVWK(liwork) 
          LAPINFO LAP1);
      lwork = int(work[0]);
      liwork = iwork[0];
      work.reset(new float[lwork]);
      iwork.reset(new int[liwork]);
#endif
#endif
      LAPNAME(sstedc) (LAPCM LAPV(c),LAPV(n),
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
#ifndef NOSTEGR
#ifdef INST_DOUBLE
  template <> inline void LapEigen_From_Tridiagonal(
      MVP<double> U, const VectorView<double>& D, 
      const VectorView<double>& E)
  { AltLapEigen_From_Tridiagonal(U,D,E); }
  template <> inline void LapEigen_From_Tridiagonal(
      MVP<std::complex<double> > U, const VectorView<double>& D, 
      const VectorView<double>& E)
  { AltLapEigen_From_Tridiagonal(U,D,E); }
#endif
#ifdef INST_FLOAT
  template <> inline void LapEigen_From_Tridiagonal(
      MVP<float> U, const VectorView<float>& D, 
      const VectorView<float>& E)
  { AltLapEigen_From_Tridiagonal(U,D,E); }
  template <> inline void LapEigen_From_Tridiagonal(
      MVP<std::complex<float> > U, const VectorView<float>& D, 
      const VectorView<float>& E)
  { AltLapEigen_From_Tridiagonal(U,D,E); }
#endif
#else // Normal stegr implementaion
#ifdef INST_DOUBLE
  template <> void LapEigen_From_Tridiagonal(
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
      //cout<<"Lap Eigen real with U\n";
      char c1 = 'V';
      char c2 = 'A';
      double junk = 0;
      int ijunk = 0;
      double tol = 0;
      int neigen;
      Vector<double> Din = D;
      Vector<double> E1(n);
      E1.SubVector(0,n-1) = E;  E1(n-1) = 0.;
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
      //cout<<"c1 = "<<c1<<std::endl;
      //cout<<"c2 = "<<c2<<std::endl;
      //cout<<"n = "<<n<<std::endl;
      //cout<<"tol = "<<tol<<std::endl;
      //cout<<"ldu = "<<ldu<<std::endl;
      //cout<<"Din = "<<Din<<std::endl;
      //cout<<"E1 = "<<E1<<std::endl;
      //cout<<"Before dstegr Norm(U) = "<<Norm(*U)<<endl;
      //cout<<"Before dstegr Norm(U1) = "<<Norm(U1)<<endl;
      LAPNAME(dstegr) (LAPCM LAPV(c1),LAPV(c2),LAPV(n),
          LAPP(Din.ptr()),LAPP(E1.ptr()),
          LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk),
          LAPV(tol),LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
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
      double nantest = U1.LinearView().SumElements();
      if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n || IsNaN(nantest)) {
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
        return AltLapEigen_From_Tridiagonal(U,D,E);
      }
      LAP_Results("dstegr");
      D = Dout;
      *U *= U1;
      //cout<<"After U *= U1 Norm(U) = "<<Norm(*U)<<endl;
    } else {
      return AltLapEigen_From_Tridiagonal(U,D,E);
    }
  }
  template <> void LapEigen_From_Tridiagonal(
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
      //cout<<"Lap Eigen complex with U\n";
      char c1 = 'V';
      char c2 = 'A';
      double junk = 0;
      int ijunk = 0;
      double tol = 0; // This is automatically bumped up to the minimum
      int neigen;
      Vector<double> Din = D;
      Vector<double> E1(n);
      E1.SubVector(0,n-1) = E;  E1(n-1) = 0.;
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
      //cout<<"dstegr (for complex): \n";
      //cout<<"c1 = "<<c1<<std::endl;
      //cout<<"c2 = "<<c2<<std::endl;
      //cout<<"n = "<<n<<std::endl;
      //cout<<"tol = "<<tol<<std::endl;
      //cout<<"ldu = "<<ldu<<std::endl;
      //cout<<"Din = "<<Din<<std::endl;
      //cout<<"E1 = "<<E1<<std::endl;
      //cout<<"Before dstegr Norm(U) = "<<Norm(*U)<<endl;
      LAPNAME(dstegr) (LAPCM LAPV(c1),LAPV(c2),LAPV(n),
          LAPP(Din.ptr()),LAPP(E1.ptr()),
          LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk),
          LAPV(tol),LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
          LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork) 
          LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
      double nantest = U1.LinearView().SumElements();
      if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n || IsNaN(nantest)) {
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
        return AltLapEigen_From_Tridiagonal(U,D,E);
      }
      //cout<<"After dstegr Norm(U) = "<<Norm(*U)<<endl;
      //cout<<"After dstegr Norm(U1) = "<<Norm(U1)<<endl;
      LAP_Results("dstegr");
      D = Dout;
      *U *= U1;
      //cout<<"After *= U1: Norm(U) = "<<Norm(*U)<<endl;
    } else {
      return AltLapEigen_From_Tridiagonal(U,D,E);
    }
  }
#endif
#ifdef INST_FLOAT
  template <> void LapEigen_From_Tridiagonal(
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
      E1.SubVector(0,n-1) = E;  E1(n-1) = 0.F;
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
      LAPNAME(sstegr) (LAPCM LAPV(c1),LAPV(c2),LAPV(n),
          LAPP(Din.ptr()),LAPP(E1.ptr()),
          LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk),
          LAPV(tol),LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
          LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork)
          LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
      float nantest = U1.LinearView().SumElements();
      if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n || IsNaN(nantest)) {
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
        return AltLapEigen_From_Tridiagonal(U,D,E);
      }
      LAP_Results("sstegr");
      D = Dout;
      *U *= U1;
    } else {
      return AltLapEigen_From_Tridiagonal(U,D,E);
    }
  }
  template <> void LapEigen_From_Tridiagonal(
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
      E1.SubVector(0,n-1) = E;  E1(n-1) = 0.F;
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
      LAPNAME(sstegr) (LAPCM LAPV(c1),LAPV(c2),LAPV(n),
          LAPP(Din.ptr()),LAPP(E1.ptr()),
          LAPV(junk),LAPV(junk),LAPV(ijunk),LAPV(ijunk),
          LAPV(tol),LAPP(&neigen),LAPP(Dout.ptr()),LAPP(U1.ptr()),LAPV(ldu),
          LAPP(isuppz.get()) LAPWK(work.get()) LAPVWK(lwork) 
          LAPWK(iwork.get()) LAPVWK(liwork) LAPINFO LAP1 LAP1);
      float nantest = U1.LinearView().SumElements();
      if (Lap_info > 0 || E1(n-1) > 0.F || neigen < n || IsNaN(nantest)) {
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
        return AltLapEigen_From_Tridiagonal(U,D,E);
      }
      LAP_Results("sstegr");
      D = Dout;
      *U *= U1;
    } else {
      return AltLapEigen_From_Tridiagonal(U,D,E);
    }
  }
#endif // FLOAT
#endif // NOSTEGR
#endif // LAP

  template <class T> void Eigen_From_Tridiagonal(
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
      A0 = *U * EDE * U->Adjoint();
    }
#ifdef LAP
    Vector<RT> D2 = D;
    Vector<RT> E2 = E;
    Matrix<T> U2(D.size(),D.size());
    if (U) { 
      U2 = *U;
      NonLapEigen_From_Tridiagonal<T>(U2.View(),D2.View(),
          E2.View());
    } else {
      NonLapEigen_From_Tridiagonal<T>(0,D2.View(),E2.View());
    }
#endif // LAP
#endif // XDEBUG

#ifdef LAP
    //if (U) cout<<"Before LapEigen Norm(U) = "<<Norm(*U)<<endl;
    LapEigen_From_Tridiagonal(U,D,E);
    //if (U) cout<<"After LapEigen Norm(U) = "<<Norm(*U)<<endl;
#else 
    //if (U) cout<<"Before NonLapEigen Norm(U) = "<<Norm(*U)<<endl;
    NonLapEigen_From_Tridiagonal(U,D,E);
    //if (U) cout<<"After NonLapEigen Norm(U) = "<<Norm(*U)<<endl;
#endif

    // Now A = U * D * Ut
    // Technically, singular values should be positive, but we allow them
    // to be negative, since these are the eigenvalues of A - no sense
    // killing that.  Also, to make them positive, we'd have to break the
    // V = Ut relationship.  So just keep that in mind later when we use S.

    // Sort output singular values by absolute value:
    auto_array<int> sortp(new int[D.size()]);
    D.Sort(sortp.get(),DESCEND,ABS_COMP);
    if (U) U->PermuteCols(sortp.get());
    //if (U) cout<<"After Permute Norm(U) = "<<Norm(*U)<<endl;

#ifdef XDEBUG
#if 0
#ifdef LAP
    auto_array<int> sortp2(new int[D.size()]);
    D2.Sort(sortp2.get(),DESCEND,ABS_COMP);
    U2.PermuteCols(sortp2.get());
    if (Norm(D-D2) > 0.001*Norm(D)) {
      cerr<<"Eigen_From_Tridiagonal:\n";
      if (D.size() < 100) {
        cerr<<"D = "<<D0<<endl;
        cerr<<"E = "<<E0<<endl;
        cerr<<"Done: D = "<<D<<endl;
        cerr<<"NonLap version: D = "<<D2<<endl;
        if (U) {
          cerr<<"U = "<<*U<<endl;
          cerr<<"UDUt = "<<(*U*DiagMatrixViewOf(D)*U->Adjoint())<<endl;
          cerr<<"A0 = "<<A0<<endl;
          cerr<<"NonLap U = "<<U2<<endl;
          cerr<<"UDUt = "<<(U2*DiagMatrixViewOf(D2)*U2.Adjoint())<<endl;
        }
      } else {
        cerr<<"D = "<<D<<std::endl;
        cerr<<"D2 = "<<D2<<std::endl;
        cerr<<"D-D2 = "<<D-D2<<std::endl;
        cerr<<"Norm(D-D2) = "<<Norm(D-D2)<<std::endl;
        int imax;
        cerr<<"(D-D2)_max = "<<(D-D2).MaxAbsElement(&imax)<<std::endl;
        cerr<<"imax = "<<imax<<std::endl;
      }
      abort();
    }
#endif
#endif
    if (U) {
      cout<<"Done Eigen_From_Tridiag: Norm(U) = "<<Norm(*U)<<endl;
      cout<<"D = "<<D<<endl;
      Matrix<T> UDU = *U * DiagMatrixViewOf(D)*U->Adjoint();
      cout<<"Norm(UDUt) = "<<Norm(UDU)<<endl;
      cout<<"Norm(UDUt-A0) = "<<Norm(UDU-A0)<<endl;
      cout<<"Norm(U) = "<<Norm(*U)<<endl;
      cout<<"Norm(UD) = "<<Norm(*U*DiagMatrixViewOf(D))<<endl;
      cout<<"Norm(A0U) = "<<Norm(A0*(*U))<<endl;
      cout<<"Norm(UD-A0U) = "<<Norm((*U)*DiagMatrixViewOf(D)-A0*(*U))<<endl;
      if (!(Norm(UDU-A0) < 0.001*Norm(A0))) {
        cerr<<"Eigen_From_Tridiagonal:\n";
        cerr<<"D = "<<D0<<endl;
        cerr<<"E = "<<E0<<endl;
        cerr<<"Done: D = "<<D<<endl;
        cerr<<"U = "<<*U<<endl;
        cerr<<"U2 = "<<U2<<endl;
        cerr<<"diff = "<<(*U-U2)<<endl;
        //cerr<<"U = "<<*U<<endl;
        //cerr<<"UD = "<<*U*DiagMatrixViewOf(D)<<endl;
        //cerr<<"DUt = "<<DiagMatrixViewOf(D)*U->Adjoint()<<endl;
        //cerr<<"UDUt = "<<UDU<<endl;
        //cerr<<"A0 = "<<A0<<endl;
        cerr<<"Norm(UDU-A0) = "<<Norm(UDU-A0)<<endl;
        cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
        throw int(1);
      }
    }
#endif // XDEBUG
  }

  template <class T> void UnsortedHermEigen(
      const MatrixView<T>& U, const VectorView<RT>& SS)
  {
#ifdef XDEBUG
    //cout<<"Start UnsortedHermEigen"<<endl;
    //cout<<"U = "<<TypeText(U)<<"  "<<U<<endl;
    //cout<<"S = "<<TypeText(SS)<<"  "<<SS<<endl;
    Matrix<T> A0 = U;
    A0.UpperTri() = A0.LowerTri().Adjoint();
#endif
    // Decompose Hermitian A (input as lower tri of U) into U S Ut
    // where S is a diagonal real matrix, and U is a unitary matrix.
    // U,S are N x N
    TMVAssert(SS.size() == U.colsize());
    TMVAssert(SS.size() == U.rowsize());

    if (U.isconj()) {
      UnsortedHermEigen(U.Conjugate(),SS);
      U.ConjugateSelf();
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
    Tridiagonalize(HermMatrixViewOf(U,Lower),Ubeta.View(),SS,E.View(),signdet);
    //cout<<"After Tridiag: D,E = \n"<<SS<<endl<<E<<endl;
    //cout<<"U = "<<U.LowerTri()<<endl;
    //cout<<"Ubeta = "<<Ubeta<<endl;

    // Now U stores Householder vectors for U in lower diagonal columns.
    for(int j=N-1;j>0;--j) U.col(j,j,N) = U.col(j-1,j,N);
    //cout<<"U slide columns = "<<U<<endl;
    U.col(0).MakeBasis(0);
    U.row(0,1,N).Zero();
    //cout<<"U.clear first row,col = "<<U<<endl;
    GetQFromQR(U.SubMatrix(1,N,1,N),Ubeta.SubVector(0,N-1));
    //cout<<"U from GetQ = "<<U<<endl;
    //cout<<"U -> "<<U<<endl;
    //BandMatrix<RT> TT = TriDiagMatrix(E,SS,E);
    //cout<<"T = "<<TT<<endl;
    //cout<<"UTUt = "<<U*TT*U.Adjoint()<<endl;
    //cout<<"cf "<<A0<<endl;
    //cout<<"Norm(diff) = "<<Norm(U*TT*U.Adjoint()-A0)<<endl;

    Eigen_From_Tridiagonal<T>(U,SS,E.View());
    //cout<<"After Eigen From Tridiag: S = "<<SS<<std::endl;
    //cout<<"USUt = "<<U*DiagMatrixViewOf(SS)*U.Adjoint()<<endl;
    //cout<<"cf "<<A0<<endl;
    //cout<<"Norm(diff) = "<<Norm(U*DiagMatrixViewOf(SS)*U.Adjoint()-A0)<<endl;

#ifdef XDEBUG
    Matrix<T> A2 = U * DiagMatrixViewOf(SS) * U.Adjoint();
    //cout<<"Done: SS = "<<SS<<std::endl;
    if (!(Norm(A0-A2) < 0.0001 * Norm(U) * Norm(SS) * Norm(U))) {
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
  template <class T> void UnsortedEigen(
      const SymMatrixView<T>& A, const VectorView<RT>& SS)
  {
    TMVAssert(SS.size() == A.size());

    if (A.isupper()) return UnsortedEigen(A.Transpose(),SS);
    if (A.isconj()) return UnsortedEigen(A.Conjugate(),SS);
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
    Tridiagonalize(A,Ubeta.View(),SS,E.View(),signdet);
    //cout<<"After Tridiag: D,E = \n"<<SS<<endl<<E<<endl;
    Eigen_From_Tridiagonal<T>(0,SS,E.View());
    //cout<<"After Eigen From Tridiag: S = \n"<<SS<<endl;
  }

  template <class T> void HermSV_Decompose(
      const MatrixView<T>& U, const DiagMatrixView<RT>& SS)
  {
    TMVAssert(U.rowsize() == SS.size());
    TMVAssert(U.colsize() == SS.size());
    TMVAssert(U.ct() == NonConj);
    TMVAssert(SS.diag().ct() == NonConj);

#ifdef XDEBUG
    //cout<<"Start HermSVDecompose"<<endl;
    Matrix<T> A0(U);
    A0.UpperTri() = A0.LowerTri().Adjoint();
    //cout<<"A0 = "<<A0<<endl;
#endif

    UnsortedHermEigen(U,SS.diag());
    //cout<<"After UnsortedHermEigen: SS = "<<SS.diag()<<endl;
    auto_array<int> sortp(new int[SS.size()]);
    SS.diag().Sort(sortp.get(),DESCEND,ABS_COMP);
    //cout<<"After sort: SS = "<<SS.diag()<<endl;
    U.PermuteCols(sortp.get());

#ifdef XDEBUG
    Matrix<T> A2 = U * SS * U.Adjoint();
    //cout<<"Done S = "<<SS.diag()<<endl;
    if (!(Norm(A0-A2) < 0.0001 * NORM(Norm(U)) * Norm(SS))) {
      cerr<<"HermSV_Decompose:\n";
      cerr<<"A = "<<A0<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"S = "<<SS.diag()<<endl;
      cerr<<"USV = "<<A2<<endl;
      abort();
    }
#endif
  }

  template <class T> void SymSV_Decompose(
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
    //cout<<"Start SymSVDecompose"<<endl;
    Matrix<T> A0(U);
    A0.UpperTri() = A0.LowerTri().Transpose();
    //cout<<"A0 = "<<A0<<endl;
#endif

    TMVAssert(IsComplex(T()));
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
    Tridiagonalize(SymMatrixViewOf(U,Lower),Ubeta.View(),D.View(),E.View(),
        signdet);
    // Now U stores Householder vectors for U in lower diagonal columns.

    BandMatrix<T,ColMajor> B(N,N,1,1);
    B.diag() = D;
    B.diag(-1) = E;
    B.diag(1) = E;

    for(int j=N-1;j>0;--j) U.col(j,j,N) = U.col(j-1,j,N);
    U.col(0).MakeBasis(0);
    U.row(0,1,N).Zero();
    GetQFromQR(U.SubMatrix(1,N,1,N),Ubeta.SubVector(0,N-1));
    if (V) *V = U.Transpose();
    Matrix<T,ColMajor> U1(N,N);
    Matrix<T,ColMajor> V1(N,N);
    SV_Decompose<T>(B,U1.View(),SS,V1.View(),logdet,signdet);
    U = U*U1;
    if (V) *V = V1*(*V);

#ifdef XDEBUG
    //cout<<"Done: S = "<<SS.diag()<<endl;
    if (V) {
      Matrix<T> A2 = U * SS * (*V);
      //std::cout<<"Norm(A2) = "<<Norm(A2)<<std::endl;
      //std::cout<<"Norm(U) = "<<Norm(U)<<std::endl;
      //std::cout<<"Norm(SS) = "<<Norm(SS)<<std::endl;
      //std::cout<<"Norm(*V) = "<<Norm(*V)<<std::endl;
      //std::cout<<"Norm(A0-A2) = "<<Norm(A0-A2)<<std::endl;
      if (!(Norm(A0-A2) < 0.0001 * Norm(U) * Norm(SS) * Norm(*V))) {
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
  template <class T> void SV_Decompose(
      const SymMatrixView<T>& A, const DiagMatrixView<RT>& SS)
  {
    TMVAssert(SS.size() == A.size());

    if (A.isherm()) {
#if 0
      cout<<"SV_Decompose: A = "<<A;
      Matrix<T> U(A.size(),A.size());
      U.LowerTri() = A.LowerTri();
      DiagMatrix<RT> S2(A.size());
      HermSV_Decompose<T>(U.View(),S2.View());
      cout<<"S1 = "<<S2.diag()<<endl;

      U.LowerTri() = A.LowerTri();
      UnsortedHermEigen(U.View(),S2.diag());
      cout<<"S2 = "<<S2.diag()<<endl;

      U = A;
      SV_Decompose(U.View(),S2.View(),false);
      cout<<"S3 = "<<S2.diag()<<endl;
#endif

      UnsortedEigen(A,SS.diag());
      //cout<<"After UnsortedEigen: SS = "<<SS.diag()<<endl;
      for(size_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) {
        SS(i) = -SS(i);
      }
      //cout<<"Positive SS = "<<SS.diag()<<endl;
      SS.diag().Sort(DESCEND);
      //cout<<"Sorted SS = "<<SS.diag()<<endl;
    } else {
      TMVAssert(IsComplex(T()));
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
      Tridiagonalize(A,Ubeta.View(),D.View(),E.View(),signdet);

      BandMatrix<T,ColMajor> B(N,N,1,1);
      B.diag() = D;
      B.diag(-1) = E;
      B.diag(1) = E;

      RT logdet(0);
      SV_Decompose<T>(B,0,SS,0,logdet,signdet);
    }
  }

  template <class T> void Eigen(
      const GenSymMatrix<T>& A, const MatrixView<T>& U,
      const VectorView<RT>& SS)
  {
    TMVAssert(A.isherm());
    TMVAssert(A.size() == SS.size());
    TMVAssert(A.size() == U.colsize());
    TMVAssert(A.size() == U.rowsize());

    if (U.isconj()) Eigen(A.Conjugate(),U.Conjugate(),SS);
    else {
      U.LowerTri() = A.LowerTri();
      UnsortedHermEigen(U,SS);
      auto_array<int> sortp(new int[A.size()]);
      SS.Sort(sortp.get(),ASCEND);
      U.PermuteCols(sortp.get());
    }
  }

  template <class T> void Eigen(
      const SymMatrixView<T>& A, const VectorView<RT>& SS)
  {
    TMVAssert(A.isherm());
    TMVAssert(A.size() == SS.size());

    if (A.isconj())
      UnsortedEigen(A.Conjugate(),SS);
    else
      UnsortedEigen(A,SS);
    SS.Sort(ASCEND);
  }

  template <class T> void SV_Decompose(
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
        SV_Decompose(A.Conjugate(),U.Conjugate(),SS,V.Conjugate());
      } else {
        SV_Decompose(A.Conjugate(),U.Conjugate(),SS,V);
        V.ConjugateSelf();
      }
    } else {
      if (V.isconj()) {
        SV_Decompose(A,U,SS,V.Conjugate());
        V.ConjugateSelf();
      } else {
        U.LowerTri() = A.LowerTri();
        if (A.isherm()) {
          HermSV_Decompose<T>(U,SS);
          V = U.Adjoint();
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

  template <class T> void SV_Decompose(
      const GenSymMatrix<T>& A,
      const MatrixView<T>& U, const DiagMatrixView<RT>& SS)
  {
    TMVAssert(A.size() == U.colsize());
    TMVAssert(A.size() == U.rowsize());
    TMVAssert(A.size() == SS.size());

    if (U.isconj()) SV_Decompose(A.Conjugate(),U.Conjugate(),SS);
    else {
      U.LowerTri() = A.LowerTri();
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

  template <class T> void SV_Decompose(
      const GenSymMatrix<T>& A,
      const DiagMatrixView<RT>& SS, const MatrixView<T>& V)
  {
    //std::cout<<"SV_Decompose: \n";
    //std::cout<<"A = "<<TypeText(A)<<"  "<<A<<std::endl;
    //std::cout<<"SS = "<<TypeText(SS)<<"  "<<SS<<std::endl;
    //std::cout<<"V = "<<TypeText(V)<<"  "<<V<<std::endl;
    TMVAssert(A.size() == SS.size());
    TMVAssert(A.size() == V.colsize());
    TMVAssert(A.size() == V.rowsize());

    if (A.isherm()) SV_Decompose(A,V.Adjoint(),SS);
    else SV_Decompose(A,V.Transpose(),SS);
    //std::cout<<"SS -> "<<SS<<std::endl;
    //std::cout<<"V -> "<<V<<std::endl;
    //std::cout<<"Vt S S V = "<<V.Adjoint()*SS*SS*V<<std::endl;
    //std::cout<<"At A = "<<A.Adjoint()*A<<std::endl;

    //Matrix<T> A2(A);
    //DiagMatrix<RT> S2(SS);
    //Matrix<T> V2(V);
    //SV_Decompose(A2.View(),S2.View(),V2.View());
    //std::cout<<"Regular SVD:\n";
    //std::cout<<"U = "<<A2<<std::endl;
    //std::cout<<"S = "<<S2<<std::endl;
    //std::cout<<"V = "<<V2<<std::endl;
    //std::cout<<"Vt S S V = "<<V2.Adjoint()*S2*S2*V2<<std::endl;
  }

  template <class T> void Polar_Decompose(const MatrixView<T>& U,
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
    //cout<<"Polar_Decompose: A0 = "<<A0<<endl;
#endif
    Matrix<T> V(U.rowsize(),U.rowsize());
    DiagMatrix<RT> S(U.rowsize());
    SV_Decompose(U.View(),S.View(),V.View(),true);
    RT thresh = Epsilon<T>()*S.size()*S(0);
    for(size_t i=0;i<S.size();i++) if (S(i) < thresh) S(i) = RT(0);
    U *= V;
    Matrix<T> VtS = V.Adjoint() * S;
    SymMultMM<false>(T(1),VtS,V,P);
#ifdef XDEBUG
    Matrix<T> A2 = U*P;
    //cout<<"Polar_Decomp: A0 = "<<TypeText(U)<<"  "<<A0<<endl;
    //cout<<"A2 = "<<A2<<endl;
    if (!(Norm(A2-A0) < 0.001*Norm(A0))) {
      cerr<<"Polar_Decompose "<<TypeText(U)<<"  "<<A0<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"Norm(UtU-1) = "<<Norm(U.Adjoint()*U-T(1))<<endl;
      cerr<<"P = "<<P<<endl;
      cerr<<"UP = "<<A2<<endl;
      cerr<<"Norm(A2-A0) = "<<Norm(A2-A0)<<endl;
      cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
      abort();
    }
#endif
  }

  template <class T> void Polar_Decompose(const GenBandMatrix<T>& A,
      const MatrixView<T>& U, const SymMatrixView<T>& P)
  {
    Matrix<T> V(A.rowsize(),A.rowsize());
    DiagMatrix<RT> S(A.rowsize());
    SV_Decompose(A,U.View(),S.View(),V.View());
    U *= V;
    Matrix<T> VtS = V.Adjoint() * S;
    SymMultMM<false>(T(1),VtS,V,P);
#ifdef XDEBUG
    Matrix<T> A2 = U*P;
    //cout<<"Band Polar_Decomp: A0 = "<<TypeText(A)<<"  "<<A<<endl;
    //cout<<"A2 = "<<A2<<endl;
    if (Norm(A2-A) > 0.001*Norm(A)) {
      cerr<<"Polar_Decompose "<<TypeText(A)<<"  "<<A<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"Norm(UtU-1) = "<<Norm(U.Adjoint()*U-T(1))<<endl;
      cerr<<"P = "<<P<<endl;
      cerr<<"UP = "<<A2<<endl;
      cerr<<"Norm(A2-A0) = "<<Norm(A2-A)<<endl;
      cerr<<"Norm(A0) = "<<Norm(A)<<endl;
      abort();
    }
#endif
  }

  template <class T> void SquareRoot(const SymMatrixView<T>& A)
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
    Eigen(A,V.View(),D.diag());
    for(size_t i=0;i<A.size();i++) {
      if (D(i) < RT(0)) 
#ifdef NOTHROW
      { std::cerr<<"Non PosDef SymMatrix found in SquareRoot\n"; exit(1); }
#else
      throw NonPosDef("in SymMatrix SquareRoot");
#endif
      D(i) = SQRT(D(i));
    }
    Matrix<T> DVt = D*V.Adjoint();
    SymMultMM<false>(T(1),V,DVt,A);
  }

#undef RT

#define InstFile "TMV_SymSVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


