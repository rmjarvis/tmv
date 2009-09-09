///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_BandMatrix.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_BandMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
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

  template <bool add, class T, class Ta, class Tb> inline void RowMultMM(
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

    size_t j1=0;
    size_t k=A.nlo();
    size_t j2=A.nhi()+1;

    for(size_t i=0;i<C.colsize(); ++i) {
      // C.row(i) (+)= alpha * A.row(i,j1,j2) * B.Rows(j1,j2);
      MultMV<add>(alpha,B.Rows(j1,j2).Transpose(),A.row(i,j1,j2),C.row(i));
      if (k>0) --k; else ++j1;
      if (j2<A.rowsize()) ++j2;
      else if (j1==A.rowsize()) {
	if (!add) C.Rows(i+1,C.colsize()).Zero();
	break;
      }
    }
  }

  template <bool add, class T, class Ta, class Tb> inline void OPMultMM(
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

    size_t i1=0;
    size_t k=A.nhi();
    size_t i2=A.nlo()+1;

    if (!add) C.Zero();
    for(size_t j=0;j<A.rowsize();++j) {
      C.Rows(i1,i2) += alpha * A.col(j,i1,i2) ^ B.row(j);
      if (k>0) --k; else ++i1;
      if (i2<A.colsize()) ++i2;
      else if (i1==A.colsize()) break;
    }
  }

  template <bool add, class T, class Ta, class Tb> inline void ColMultMM(
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

    for(size_t j=0;j<B.rowsize();j++)
      // C.col(j) (+)= alpha * A * B.col(j);
      MultMV<add>(alpha,A,B.col(j),C.col(j));
  }

  template <bool add, class T, class Ta, class Tb> 
    inline void NonLapTriDiagMultMM(
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

      const size_t N = A.diag().size();
      const size_t M = A.rowsize()>A.colsize() ? N : N-1;

      const Ta* di = A.cptr();
      const Ta* dui = A.diag(1).cptr();
      const Ta* dli = A.diag(-1).cptr()-1;

      for(size_t i=0;i<N;++i,++di,++dui,++dli) {
	if (add) C.row(i) += *di*B.row(i);
	else C.row(i) = *di*B.row(i);
	if (i>0) C.row(i) += *dli*B.row(i-1);
	if (i<M) C.row(i) += *dui*B.row(i+1);
      }
      if (A.colsize() > A.rowsize()) {
	if (add) C.row(N) += *dli*B.row(N-1);
	else C.row(N) = *dli*B.row(N-1);
      }
    }

#ifdef ELAP
  template <class T, class Ta, class Tb> inline void LapTriDiagMultMM(
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B, 
      const int beta, const MatrixView<T>& C)
  {
    if (beta == 1) NonLapTriDiagMultMM<true>(A,B,C); 
    else NonLapTriDiagMultMM<false>(A,B,C); 
  }
#ifdef INST_DOUBLE
  template <> inline void LapTriDiagMultMM(
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
    TMVAssert(A.IsSquare());
    TMVAssert(B.iscm());
    TMVAssert(C.iscm());

    int n = A.colsize();
    int nrhs = B.rowsize();
    double a(1);
    int ldB = B.stepj();
    double xbeta(beta);
    int ldC = C.stepj();
    LAPNAMEX(dlagtm) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPV(a),LAPP(A.cptr()+A.stepj()),
	LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
	LAPP(B.cptr()),LAPV(ldB),LAPV(xbeta),LAPP(C.ptr()),LAPV(ldC) LAP1);
  }
  template <> inline void LapTriDiagMultMM(
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
    TMVAssert(A.IsSquare());
    TMVAssert(B.iscm());
    TMVAssert(C.iscm());
    TMVAssert(B.isconj() == C.isconj());

    int n = A.colsize();
    int nrhs = B.rowsize();
    double a(1);
    int ldB = B.stepj();
    double xbeta(beta);
    int ldC = C.stepj();
    LAPNAMEX(zlagtm) (LAPCM B.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPV(a),LAPP(A.cptr()+A.stepj()),
	LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
	LAPP(B.cptr()),LAPV(ldB),LAPV(xbeta),LAPP(C.ptr()),LAPV(ldC) LAP1);
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapTriDiagMultMM(
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
    TMVAssert(A.IsSquare());
    TMVAssert(B.iscm());
    TMVAssert(C.iscm());

    int n = A.colsize();
    int nrhs = B.rowsize();
    float a(1);
    int ldB = B.stepj();
    float xbeta(beta);
    int ldC = C.stepj();
    LAPNAMEX(slagtm) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPV(a),LAPP(A.cptr()+A.stepj()),
	LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
	LAPP(B.cptr()),LAPV(ldB),LAPV(xbeta),LAPP(C.ptr()),LAPV(ldC) LAP1);
  }
  template <> inline void LapTriDiagMultMM(
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
    TMVAssert(A.IsSquare());
    TMVAssert(B.iscm());
    TMVAssert(C.iscm());
    TMVAssert(B.isconj() == C.isconj());

    int n = A.colsize();
    int nrhs = B.rowsize();
    float a(1);
    int ldB = B.stepj();
    float xbeta(beta);
    int ldC = C.stepj();
    LAPNAMEX(clagtm) (LAPCM B.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPV(a),LAPP(A.cptr()+A.stepj()),
	LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
	LAPP(B.cptr()),LAPV(ldB),LAPV(xbeta),LAPP(C.ptr()),LAPV(ldC) LAP1);
  }
#endif // FLOAT
#endif // ELAP

  template <bool add, class T, class Ta, class Tb> inline void DoTriDiagMultMM(
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef ELAP
    if (A.IsSquare() && B.iscm() && C.iscm() && B.isconj() == C.isconj())
      LapTriDiagMultMM(A,B,add?1:0,C);
    else
#endif
      NonLapTriDiagMultMM<add>(A,B,C);
  }

  template <bool add, class T, class Ta, class Tb> inline void TriDiagMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    if (alpha == T(1) && A.isdm()) {
      if (A.isconj()) 
	DoTriDiagMultMM<add>(A.Conjugate(),B.Conjugate(),C.Conjugate());
      else
	DoTriDiagMultMM<add>(A,B,C);
    } else if (IMAG(alpha) == RealType(T)(0)) {
      BandMatrix<Ta,DiagMajor> A1 = REAL(alpha)*A;
      DoTriDiagMultMM<add>(A1,B,C);
    } else {
      BandMatrix<T,DiagMajor> A1 = alpha*A;
      DoTriDiagMultMM<add>(A1,B,C);
    }
  }

  // MJ: Put in a recursive block calculation here.  (Also in Band*Band)
  template <bool add, class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"Start DoMultMM:\n";
    //cerr<<"alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A.cptr()<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B.cptr()<<"  "<<B<<endl;
    //cerr<<"add = "<<add<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C.cptr();
    //if (add) cerr<<"  "<<C;
    //cerr<<endl;
    Matrix<Tb> B0 = B;
    Matrix<Ta> A0 = A;
    Matrix<T> C0 = C;
    Matrix<T> C2 = C;
    MultMM<add>(alpha,A0,B0,C2.View());
    //cerr<<"C2 = "<<C2<<endl;
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
	TriDiagMultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	    C.Conjugate());
      else
	TriDiagMultMM<add>(alpha,A,B,C);
    }
    else if (C.colsize() < C.rowsize()) RowMultMM<add>(alpha,A,B,C);
    else ColMultMM<add>(alpha,A,B,C);
#ifdef XDEBUG
    //cerr<<"C -> "<<C<<endl;
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"DoMultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A.cptr()<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B.cptr()<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C.cptr()<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <bool add, class T, class Ta, class Tb> inline void FullTempMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    } else {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    }
  }

  template <bool add, class T, class Ta, class Tb> inline void BlockTempMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = std::min(C.rowsize(),j+MM_BLOCKSIZE);
      if (IMAG(alpha) == RealType(T)(0)) {
	if (C.isrm()) {
	  Matrix<Tb,RowMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
	} else  {
	  Matrix<Tb,ColMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
	}
      } else {
	if (C.isrm()) {
	  Matrix<T,RowMajor> B2 = alpha * B.Cols(j,j2);
	  DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
	} else  {
	  Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	  DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
	}
      }
      j=j2;
    }
  }

  template <bool add, class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C (+)= alpha * A * B
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

#ifdef XDEBUG
    //cerr<<"Start MultMM:\n";
    //cerr<<"alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A.cptr()<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B.cptr()<<"  "<<B<<endl;
    //cerr<<"add = "<<add<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C.cptr();
    //if (add) cerr<<"  "<<C;
    //cerr<<endl;
    Matrix<Tb> B0 = B;
    Matrix<Ta> A0 = A;
    Matrix<T> C0 = C;
    Matrix<T> C2 = C;
    MultMM<add>(alpha,A0,B0,C2.View());
    //cerr<<"C2 = "<<C2<<endl;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (A.rowsize() == 0 || alpha == T(0)) {
	if (!add) C.Zero();
      } else if (A.rowsize() > A.colsize()+A.nhi()) {
	MultMM<add>(alpha,A.Cols(0,A.colsize()+A.nhi()),
	    B.Rows(0,A.colsize()+A.nhi()),C);
      } else if (A.colsize() > A.rowsize()+A.nlo()) {
	MultMM<add>(alpha,A.Rows(0,A.rowsize()+A.nlo()),
	    B,C.Rows(0,A.rowsize()+A.nlo()));
	if (!add) C.Rows(A.rowsize()+A.nlo(),A.colsize()).Zero();
      } else if (C.isconj()) {
	MultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	    C.Conjugate());
      } else if (SameStorage(A,C)) {
	FullTempMultMM<add>(alpha,A,B,C);
      } else if (SameStorage(B,C)) {
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM<add>(alpha,A,B,C);
	else 
	  FullTempMultMM<add>(alpha,A,B,C);
      } else {
	DoMultMM<add>(alpha, A, B, C);
      }
    }
#ifdef XDEBUG
    //cerr<<"C -> "<<C<<endl;
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A.cptr()<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B.cptr()<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C.cptr()<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_BandMatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv
