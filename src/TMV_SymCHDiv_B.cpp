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
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_SymMatrix.h"
#include "TMV_SymCHDiv.h"
#include "TMV_SymCHDiv_B.h"
#include "TMV_MatrixArith.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_SymMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <class T> inline void NonLapCHInverse(const SymMatrixView<T>& sinv)
  {
    TMVAssert(sinv.isherm());
    // inv = (L Lt)^-1 = Lt^-1 L^-1
    LowerTriMatrixView<T> L = sinv.LowerTri();
    L.InvertSelf();
    sinv = L.Adjoint() * L;
  }

#ifdef ALAP
  template <class T> inline void LapCHInverse(const SymMatrixView<T>& sinv)
  { NonLapCHInverse(sinv); }
#ifdef INST_DOUBLE
  template <> inline void LapCHInverse(const SymMatrixView<double>& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    int n = sinv.size();
    int lda = sinv.stepj();
    LAPNAME(dpotri) (LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
	LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("dpotri");
  }
  template <> inline void LapCHInverse(
      const SymMatrixView<std::complex<double> >& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    int n = sinv.size();
    int lda = sinv.stepj();
    LAPNAME(zpotri) (LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
	LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("zpotri");
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapCHInverse(const SymMatrixView<float>& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    int n = sinv.size();
    int lda = sinv.stepj();
    LAPNAME(spotri) (LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
	LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("spotri");
  }
  template <> inline void LapCHInverse(
      const SymMatrixView<std::complex<float> >& sinv)
  {
    TMVAssert(sinv.isherm());
    TMVAssert(sinv.iscm());

    int n = sinv.size();
    int lda = sinv.stepj();
    LAPNAME(cpotri) (LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
	LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("cpotri");
  }
#endif // FLOAT
#endif // LAP

  template <class T, class T1> void HermCH_Inverse(const GenSymMatrix<T>& LLx,
      const SymMatrixView<T1>& sinv) 
  {
#ifdef XTEST
    TMVAssert(LLx.HermOK());
#endif
    TMVAssert(LLx.size() == sinv.size());

#ifdef XDEBUG
    Matrix<T> A = LLx.LowerTri() * LLx.UpperTri();
#endif

    if (sinv.size() > 0) {
      if (
#ifdef ALAP
	  !SameType(T(),T1()) || 
#endif
	  !(sinv.iscm() || sinv.isrm())) {
	HermMatrix<T,Lower,ColMajor> temp(sinv.size());
	HermCH_Inverse(LLx,temp.View());
	sinv = temp;
#ifdef ALAP
      } else if (sinv.isrm()) {
	HermCH_Inverse(LLx.Transpose(),sinv.Transpose());
#endif
      } else {
	sinv = LLx;
#ifdef ALAP
	LapCHInverse(sinv);
#else
	NonLapCHInverse(sinv);
#endif
      }
    }

#ifdef XDEBUG
    Matrix<T1> eye = A * sinv;
    RealType(T) kappa = Norm(A) * Norm(sinv);
    if (Norm(eye-T1(1)) > 0.0001*kappa*sinv.size()) {
      cerr<<"A = "<<A<<endl;
      cerr<<"sinv = "<<sinv<<endl;
      cerr<<"A*sinv = "<<A*sinv<<endl;
      cerr<<"sinv*A = "<<sinv*A<<endl;
      cerr<<"Norm(A*sinv-1) = "<<Norm(A*sinv-T1(1))<<endl;
      cerr<<"kappa = "<<kappa<<endl;
      abort();
    }
#endif
#ifdef XTEST
    TMVAssert(sinv.HermOK());
#endif
  }

  template <bool herm, class T> void SymATASquare(const MatrixView<T>& A)
  {
    const size_t N = A.colsize();
    if (N == 1) {
      const T A00 = *A.ptr();
      if (herm)
	*A.ptr() = NORM(REAL(A00));
      else 
	*A.ptr() = SQR(A00);
    } else {
      const size_t K = N/2;
      MatrixView<T> A00 = A.SubMatrix(0,K,0,K);
      MatrixView<T> A10 = A.SubMatrix(K,N,0,K);
      MatrixView<T> A01 = A.SubMatrix(0,K,K,N);
      MatrixView<T> A11 = A.SubMatrix(K,N,K,N);
      MatrixView<T> A10t = herm ? A10.Adjoint() : A10.Transpose();

      // [ A00 A10t ] [ A00 A10t ] 
      // [ A10 A11  ] [ A10 A11  ]
      // = [ A00^2 + A10t A10    A00 A10t + A10t A11 ]
      //   [ A10 A00 + A11 A10   A10 A10t + A11^2    ]

      // A10 stores the actual data for A10
      // We can therefore write to A01, sort of as a temp matrix.
      
      A01 = A00 * A10t;
      A01 += A10t * A11;

      SymATASquare<herm>(A00);
      A00 += A10t*A10;
      SymATASquare<herm>(A11);
      A11 += A10*A10t;

      A10t = A01;
    }
  }
  
#define InstFile "TMV_SymCHDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


