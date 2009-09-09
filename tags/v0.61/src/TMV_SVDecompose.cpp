///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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



#include "TMV_Blas.h"
#include "TMV_SVDiv.h"
#include "TMV_SVD.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_QRDiv.h"
#include "TMV_MatrixArith.h"
#include "TMV_VectorArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include "TMV_DiagMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

//#define TIME

#ifdef TIME
#include <sys/time.h>
#include <iostream>
using std::cout;
using std::endl;
#endif

namespace tmv {

#define RT RealType(T)

  template <class T> static inline void DoSV_Decompose_From_Bidiagonal_NZ(
      MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V,
      bool UisI, bool VisI)
  { SV_Decompose_From_Bidiagonal_DC<T>(U,D,E,V,UisI,VisI); }
  //{ SV_Decompose_From_Bidiagonal_QR<T>(U,D,E,V); }

  template <class T> void DoSV_Decompose_From_Bidiagonal(
      MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V,
      bool UisI, bool VisI)
  {
    //cout<<"Start Decompose from Bidiagonal:\n";
    //cout<<"D = "<<D<<std::endl;
    //cout<<"E = "<<E<<std::endl;
    const int N = D.size();

    // D and E may have zeros on entry.  
    // This routine does the trivial deflations to get to subproblems
    // in which D and E are fully non-zero

    // First chop any small elements in D,E
    BidiagonalChopSmallElements(D,E);

    // Find sub-problems to solve:
    for(int q = N-1; q>0; ) {
      if (E(q-1) == T(0)) --q;
      else if (D(q) == T(0)) {
	// We have the end looking like:
	//   ? ?
	//     ? x
	//       0
	// So we need to find a p where all E(i) with p<=i<q are non-zero.
	int p = q-1;
	while (p>0 && E(p-1) != T(0)) --p;
	// Now Zero out the last column:
	if (V) BidiagonalZeroLastCol<T>(D.SubVector(p,q),E.SubVector(p,q),
	    V->Rows(p,q+1));
	else BidiagonalZeroLastCol<T>(D.SubVector(p,q),E.SubVector(p,q),0);
	VisI = false;
	--q;
      } else {
	// Find first p before q with either E(p) = 0 or D(p) = 0
	int p=q-1;
	while (p > 0 && (E(p-1) != T(0)) && D(p) != T(0)) --p; 
	if (p > 0 && D(p) == T(0)) {
	  // We have a block looking like:
	  //   0 x
	  //     x x 
	  //       x x
	  //         x
	  if (U) BidiagonalZeroFirstRow<T>(U->Cols(p,q+1),
	      D.SubVector(p+1,q+1),E.SubVector(p,q));
	  else BidiagonalZeroFirstRow<T>(0,
	      D.SubVector(p+1,q+1),E.SubVector(p,q));
	  UisI = false;
	  ++p;
	}
	if (q > p) {
	  if (U)
	    if (V) 
	      DoSV_Decompose_From_Bidiagonal_NZ<T>(
		  U->Cols(p,q+1),D.SubVector(p,q+1),
		  E.SubVector(p,q),V->Rows(p,q+1),
		  UisI && p==0 && q+1==N, VisI && p==0 && q+1==N);
	    else 
	      DoSV_Decompose_From_Bidiagonal_NZ<T>(
		  U->Cols(p,q+1),D.SubVector(p,q+1),
		  E.SubVector(p,q),0, 
		  UisI && p==0 && q+1==N, false);
	  else
	    if (V) 
	      DoSV_Decompose_From_Bidiagonal_NZ<T>(
		  0,D.SubVector(p,q+1),
		  E.SubVector(p,q),V->Rows(p,q+1),
		  false, VisI && p==0 && q+1==N);
	    else 
	      DoSV_Decompose_From_Bidiagonal_NZ<T>(
		  0,D.SubVector(p,q+1),E.SubVector(p,q),0,
		  false, false);
	}
	q = p;
      }
    }
  }

  template <class T> static void NonLapSV_Decompose_From_Bidiagonal(
      MVP<T> U, const VectorView<RT>& D, 
      const VectorView<RT>& E, MVP<T> V, bool SetUV)
  {
#ifdef XDEBUG
    //cout<<"Start Decompose from Bidiag:\n";
    //if (U) cout<<"U = "<<Type(*U)<<"  "<<*U<<endl;
    //cout<<"D = "<<Type(D)<<"  step "<<D.step()<<"  "<<D<<endl;
    //cout<<"E = "<<Type(E)<<"  step "<<E.step()<<"  "<<E<<endl;
    //if (V) cout<<"V = "<<Type(*V)<<"  "<<*V<<endl;
    //cout<<"SetUV = "<<SetUV<<endl;
    Matrix<RT> B(D.size(),D.size(),RT(0));
    B.diag() = D;
    B.diag(1) = E;
    Matrix<T> A0(U&&V ? U->colsize() : D.size(),D.size());
    if (U && V && !SetUV) A0 = (*U) * B * (*V);
    else A0 = B;
    //cout<<"A0 = "<<A0<<endl;
#endif

    const int N = D.size();

    if (SetUV) {
      TMVAssert(U && V);
      U->SetToIdentity();
      V->SetToIdentity();
    }

    DoSV_Decompose_From_Bidiagonal<T>(U,D,E,V,SetUV,SetUV);

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

    // Now A = U * S * V
    // Sort output singular values 
    auto_array<int> sortp(new int[N]);
    D.Sort(sortp.get(),DESCEND);
    if (U) U->PermuteCols(sortp.get());
    if (V) V->PermuteRows(sortp.get());

#ifdef XDEBUG
    if (U && V) {
      Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
      if (Norm(A0-AA) > 0.001*Norm(A0)) {
	cerr<<"SV_DecomposeFromBidiagonal: \n";
	cerr<<"input B = "<<B<<endl;
	cerr<<"UBV = "<<A0<<endl;
	cerr<<"USV = "<<AA<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"S = "<<D<<endl;
	cerr<<"V = "<<*V<<endl;
	abort();
      }
    }
#endif
  }

#ifdef LAP 
  template <class T> static inline void LapSV_Decompose_From_Bidiagonal(
      MVP<T> U, const VectorView<RT>& D, 
      const VectorView<RT>& E, MVP<T> V, bool SetUV)
  { NonLapSV_Decompose_From_Bidiagonal<T>(U,D,E,V,SetUV); }
#ifdef INST_DOUBLE
  template <> void LapSV_Decompose_From_Bidiagonal(
      MVP<double> U, const VectorView<double>& D, 
      const VectorView<double>& E, MVP<double> V, bool SetUV)
  {
#ifdef TIME
    timeval tp;
    gettimeofday(&tp,0);
    double t0 = tp.tv_sec + tp.tv_usec/1.e6;
    double t1,t2;
    cout<<"LapSV_Decompose_From_Bidiagonal\n";
    cout<<"  U,V,SetUV = "<<bool(U)<<','<<bool(V)<<','<<SetUV<<endl;
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
      if (U && SetUV) TMVAssert(U->stor() == V->stor());
    }
    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      if (U->iscm()) {
	TMVAssert(V->iscm());
	int ldu = U->stepj();
	int ldv = V->stepj();
#ifdef TIME
	gettimeofday(&tp,0);
	t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	    LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(U->ptr()),LAPV(ldu),LAPP(V->ptr()),LAPV(ldv),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      } else {
	u = 'L';
	TMVAssert(U->isrm());
	TMVAssert(V->isrm());
	int ldu = U->stepi();
	int ldv = V->stepi();
#ifdef TIME
	gettimeofday(&tp,0);
	t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	    LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(V->ptr()),LAPV(ldv),LAPP(U->ptr()),LAPV(ldu),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      }
    } else if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
#ifdef TIME
      gettimeofday(&tp,0);
      t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
      gettimeofday(&tp,0);
      t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
#ifdef TIME
      gettimeofday(&tp,0);
      t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
      gettimeofday(&tp,0);
      t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
    }
    LAP_Results("dbdsdc");
#ifdef TIME
    gettimeofday(&tp,0);
    double t3 = tp.tv_sec + tp.tv_usec/1.e6;
    cout<<"  dbdsdc time = "<<t2-t1<<endl;
    cout<<"  other time = "<<t3-t2+t1-t0<<endl;
#endif
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      MVP<std::complex<double> > U, const VectorView<double>& D, 
      const VectorView<double>& E, MVP<std::complex<double> > V, 
      bool SetUV)
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
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      if (SetUV) {
	if (U) *U = U1;
	if (V) *V = V1;
      } else {
	if (U) *U = *U*U1;
	if (V) *V = V1*(*V);
      }
    } else {
      TMVAssert(!SetUV);
      int ldu = n;
      int ldv = n;
      char c = 'N';
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
    }
    LAP_Results("dbdsdc");
  }
#endif
#ifdef INST_FLOAT
  template <> void LapSV_Decompose_From_Bidiagonal(
      MVP<float> U, const VectorView<float>& D, 
      const VectorView<float>& E, MVP<float> V, bool SetUV)
  {
#ifdef TIME
    timeval tp;
    gettimeofday(&tp,0);
    float t0 = tp.tv_sec + tp.tv_usec/1.e6;
    float t1,t2;
    cout<<"LapSV_Decompose_From_Bidiagonal\n";
    cout<<"  U,V,SetUV = "<<bool(U)<<','<<bool(V)<<','<<SetUV<<endl;
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
      if (U && SetUV) TMVAssert(U->stor() == V->stor());
    }
    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      if (U->iscm()) {
	TMVAssert(V->iscm());
	int ldu = U->stepj();
	int ldv = V->stepj();
#ifdef TIME
	gettimeofday(&tp,0);
	t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	    LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(U->ptr()),LAPV(ldu),LAPP(V->ptr()),LAPV(ldv),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      } else {
	u = 'L';
	TMVAssert(U->isrm());
	TMVAssert(V->isrm());
	int ldu = U->stepi();
	int ldv = V->stepi();
#ifdef TIME
	gettimeofday(&tp,0);
	t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(V->ptr()),LAPV(ldv),LAPP(U->ptr()),LAPV(ldu),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      }
    } else if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
#ifdef TIME
      gettimeofday(&tp,0);
      t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
      gettimeofday(&tp,0);
      t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
#ifdef TIME
      gettimeofday(&tp,0);
      t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
      gettimeofday(&tp,0);
      t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
    }
    LAP_Results("sbdsdc");
#ifdef TIME
    gettimeofday(&tp,0);
    float t3 = tp.tv_sec + tp.tv_usec/1.e6;
    cout<<"  sbdsdc time = "<<t2-t1<<endl;
    cout<<"  other time = "<<t3-t2+t1-t0<<endl;
#endif
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      MVP<std::complex<float> > U, const VectorView<float>& D, 
      const VectorView<float>& E, MVP<std::complex<float> > V, 
      bool SetUV)
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
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      if (SetUV) {
	TMVAssert(U && V);
	*U = U1;
	*V = V1;
      } else {
	if (U) *U = *U*U1;
	if (V) *V = V1*(*V);
      }
    } else {
      TMVAssert(!SetUV);
      int ldu = n;
      int ldv = n;
      char c = 'N';
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
    }
    LAP_Results("sbdsdc");
  }
#endif // FLOAT
#endif // LAP

  template <class T> void SV_Decompose_From_Bidiagonal(
      MVP<T> U, const VectorView<RT>& D, 
      const VectorView<RT>& E, MVP<T> V, bool SetUV)
  {
#ifdef XDEBUG
    Matrix<RT> B(D.size(),D.size(),RT(0));
    B.diag() = D;
    B.diag(1) = E;
    Matrix<T> A0(U&&V ? U->colsize() : D.size(),D.size());
    if (U && V && !SetUV) A0 = (*U) * B * (*V);
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
    TMVAssert(!U || !V || !SetUV || U->stor()==V->stor());

    if (D.size() > 0) {
#ifdef LAP
      LapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
      RT Dmax = MaxAbsElement(D);
      D.Clip(Epsilon<T>()*Dmax);
#else 
      NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
#endif
    }
#ifdef XDEBUG
    if (U && V) {
      Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
      //cout<<"SVDecomposeFromBidiag: Norm(A0-AA) = "<<Norm(A0-AA)<<std::endl;
      //cout<<"cf "<<0.001*Norm(A0)<<std::endl;
      if (Norm(A0-AA) > 0.001*Norm(A0)) {
	cerr<<"SV_DecomposeFromBidiagonal: \n";
	cerr<<"input B = "<<B<<endl;
	cerr<<"UBV = "<<A0<<endl;
	cerr<<"USV = "<<AA<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"S = "<<D<<endl;
	cerr<<"V = "<<*V<<endl;
#ifdef LAP
	Matrix<T,ColMajor> U2 = U0;
	Matrix<T,ColMajor> V2 = V0;
	Vector<RT> D2 = D0;
	Vector<RT> E2 = E0;
	NonLapSV_Decompose_From_Bidiagonal<T>(U2.View(),D2.View(),E2.View(),
	    V2.View(),SetUV);
	cerr<<"U = "<<*U<<endl;
	cerr<<"NonLap U: "<<U2<<endl;
	cerr<<"Diff U = "<<Matrix<T>(*U-U2).Clip(1.e-3);
	cerr<<"V = "<<*V<<endl;
	cerr<<"NonLap V: "<<V2<<endl;
	cerr<<"Diff V = "<<Matrix<T>(*V-V2).Clip(1.e-3);
	cerr<<"D = "<<D<<endl;
	cerr<<"NonLap D: = "<<D2<<endl;
	cerr<<"Diff D = "<<Vector<T>(D-D2).Clip(1.e-3);
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
  
  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const DiagMatrixView<RT>& S, 
      MVP<T> V, RT& logdet, T& signdet, bool StoreU)
  {
#ifdef XDEBUG
    Matrix<T> A0(U);
    //cout<<"SVDecompose:\n";
    //cout<<"A0 = "<<A0<<std::endl;
    //cout<<"StoreU = "<<StoreU<<std::endl;
#endif
#ifdef TIME
    cout<<"Start SV_Decompose #1\n";
    timeval tp;
    gettimeofday(&tp,0);
    double t0 = tp.tv_sec + tp.tv_usec/1.e6;
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

    // If M is much larger than N (technically M > 5/3 N), then it is quicker
    // to start by doing a QR decomposition and then do SVD on the square
    // R matrix.  Thus, the final U of the SVD is Q (from the QR decomp)
    // times U from R's SVD.
    if (M > 5*N/3) {
      if (StoreU) {
	Matrix<T,ColMajor> R(N,N);
	R.LowerTri().OffDiag().Zero();
	QR_Decompose(U,R.UpperTri(),signdet);
#ifdef TIME
	gettimeofday(&tp,0);
	double t1 = tp.tv_sec + tp.tv_usec/1.e6;
	cout<<"QR_Decompose: "<<t1-t0<<" seconds\n";
#endif
	SV_Decompose(R.View(),S,V,logdet,signdet,StoreU);
#ifdef TIME
	gettimeofday(&tp,0);
	double t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	// Now R is a Unitary Matrix U'.  Need to multiply U by U'
	U = U*R;
#ifdef TIME
	gettimeofday(&tp,0);
	double t3 = tp.tv_sec + tp.tv_usec/1.e6;
	cout<<"Mult U*R: "<<t3-t2<<" seconds\n";
#endif
      } else {
	Vector<T> Qbeta(N);
	QR_Decompose(U,Qbeta.View(),signdet);
	if (N > 1)
	  LowerTriMatrixViewOf(U.Rows(0,N)).OffDiag().Zero();
	SV_Decompose(U.Rows(0,N),S,V,logdet,signdet,StoreU);
      }
    } else {
      // First we reduce A to bidiagonal form: A = U * B * V
      // using a series of Householder transformations.
      // The diagonal of the Bidiagonal Matrix B is stored in D.
      // The superdiagonal is stored in E.
      Vector<RT> E(N-1);
      Vector<T> Ubeta(N);
      Vector<T> Vbeta(N-1);
      //cout<<"Before Bidiagonalize: \n";
      Bidiagonalize(U,Ubeta.View(),Vbeta.View(),S.diag(),E.View(),
	  signdet);
      //cout<<"After Bidiagonalize: \n";
      //cout<<"U = "<<U<<endl;
      //cout<<"Ubeta = "<<Ubeta<<endl;
      //cout<<"Vbeta = "<<Vbeta<<endl;
      //cout<<"D = "<<S.diag()<<endl;
      //cout<<"E = "<<E<<endl;
#ifdef TIME
      gettimeofday(&tp,0);
      double t1 = tp.tv_sec + tp.tv_usec/1.e6;
      cout<<"Bidiagonalize: "<<t1-t0<<" seconds\n";
#endif
      // The determinant of B is just the product of the diagonal elements:
      if (signdet != T(0)) {
	RT s;
	logdet += S.LogDet(&s);
	signdet *= s;
      }
#ifdef TIME
      gettimeofday(&tp,0);
      double t1b = tp.tv_sec + tp.tv_usec/1.e6;
      cout<<"Calc det: "<<t1b-t1<<" seconds\n";
#endif

      // Now UV stores Householder vectors for U in lower diagonal columns 
      // (HLi) and Householder vectors for V in upper diagonal rows (HRi)
      // The Householder matrices for U are actually the adjoints of the 
      // matrices that bidiagonalize A, and for V are the transposes:
      // U = HLn-1t ... HL1t HL0t A HR0T HR1T ... HRn-2T
      // Using the fact that H Ht = I, we get A = U B V with:
      // U = HL0 ... HLn-1 
      if (V) {
	V->row(0).MakeBasis(0);
	V->Rows(1,N) = U.Rows(0,N-1);
	V->col(0,1,N).Zero();
	GetQFromQR(V->SubMatrix(1,N,1,N).Transpose(),Vbeta);
	//cout<<"V -> "<<V<<endl;
      }
      if (StoreU) {
	GetQFromQR(U,Ubeta);
	//cout<<"U -> "<<U<<endl;
      }
#ifdef TIME
      gettimeofday(&tp,0);
      double t2 = tp.tv_sec + tp.tv_usec/1.e6;
      cout<<"Make U,V: "<<t2-t1b<<" seconds\n";
#endif

      //cout<<"Before SV_Decompose_From_Bidiag\n";
      if (StoreU) SV_Decompose_From_Bidiagonal<T>(U,S.diag(),E.View(),V);
      else SV_Decompose_From_Bidiagonal<T>(0,S.diag(),E.View(),V);
      //cout<<"Done: S = "<<S.diag()<<endl;

#ifdef TIME
      gettimeofday(&tp,0);
      double t3 = tp.tv_sec + tp.tv_usec/1.e6;
      cout<<"Decompose From Bidiagonal: "<<t3-t2<<" seconds\n";
#endif
    }
#ifdef TIME
    gettimeofday(&tp,0);
    double tx = tp.tv_sec + tp.tv_usec/1.e6;
    cout<<"Total time: "<<tx-t0<<" seconds\n";
#endif
#ifdef XDEBUG
    if (StoreU && V && S.size()>0) {
      Matrix<T> A2 = U * S * (*V);
      //cout<<"SVDecompose: Norm(A0-A2) = "<<Norm(A0-A2)<<std::endl;
      //cout<<"cf "<<0.001*Norm(U)*Norm(S)*Norm(*V)<<std::endl;
      if (Norm(A0-A2) > 0.0001 * Norm(U) * Norm(S) * Norm(*V)) {
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
  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const DiagMatrixView<RT>& SS, 
      const MatrixView<T>& V, bool StoreU)
  { 
    if (U.isconj()) {
      if (V.isconj()) {
	SV_Decompose(U.Conjugate(),SS,V.Conjugate(),StoreU);
      } else {
	SV_Decompose(U.Conjugate(),SS,V,StoreU);
	V.ConjugateSelf();
      }
    } else {
      if (V.isconj()) {
	SV_Decompose(U,SS,V.Conjugate(),StoreU);
	V.ConjugateSelf();
      } else {
	RT ld(0);
	T d(0);
	SV_Decompose<T>(U,SS,V,ld,d,StoreU); 
      }
    }
  }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const DiagMatrixView<RT>& SS, bool StoreU)
  { 
    if (U.isconj()) {
      SV_Decompose(U.Conjugate(),SS,StoreU);
    } else {
      RT ld(0);
      T d(0);
      SV_Decompose<T>(U,SS,0,ld,d,StoreU); 
    }
  }

#undef RT

#define InstFile "TMV_SVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


