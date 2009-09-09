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
#include "TMV_DiagMatrix.h"
#include "TMV_LUDiv.h"
#include "TMV_TriMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <class T> inline void NonLapLUInverse(const MatrixView<T>& minv)
  {
    // m = P L U
    // m^-1 = U^-1 L^-1 Pt
    UpperTriMatrixView<T> U = UpperTriMatrixViewOf(minv);
    LowerTriMatrixView<T> L = LowerTriMatrixViewOf(minv,UnitDiag);
    UpperTriMatrix<T> U0 = U;
    LowerTriMatrix<T> L0 = L;
    U.InvertSelf();
    L.InvertSelf();
    //cerr<<"U*Uinv = "<<U0*U<<endl;
    //cerr<<"L*Linv = "<<L0*L<<endl;
    //cerr<<"LUUinvLinv = "<<L0*U0*U*L<<endl;
    minv = U*L;
    //cerr<<"L*U = "<<L*U<<endl;
    //cerr<<"minv = "<<minv<<endl;
    //cerr<<"LUminv = "<<L*U*minv<<endl;
    // Do Pt back in LU_Inverse
  }

#ifdef ALAP
  template <class T> inline void LapLUInverse(const MatrixView<T>& minv)
  { NonLapLUInverse(minv); }
#ifdef INST_DOUBLE
  template <> inline void LapLUInverse(const MatrixView<double>& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    LAPNAME(dgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("dgetri");
  }
  template <> inline void LapLUInverse(
      const MatrixView<std::complex<double> >& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    std::complex<double>* work = LAP_ZWork(lwork);
#endif
    LAPNAME(zgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("zgetri");
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapLUInverse(const MatrixView<float>& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    LAPNAME(sgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("sgetri");
  }
  template <> inline void LapLUInverse(
      const MatrixView<std::complex<float> >& minv)
  {
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.iscm());

    int n = minv.colsize();
    int lda = minv.stepj();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    std::complex<float>* work = LAP_CWork(lwork);
#endif
    LAPNAME(cgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("cgetri");
  }
#endif // FLOAT
#endif // ALAP

  template <class T, class T1> void LU_Inverse(const GenMatrix<T>& LUx,
      const size_t* P, const MatrixView<T1>& minv)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(minv.IsSquare());
    TMVAssert(minv.colsize() == LUx.colsize());
#ifdef XDEBUG
    Matrix<T> m = LowerTriMatrixViewOf(LUx,UnitDiag) *
      UpperTriMatrixViewOf(LUx);
    m.ReversePermuteRows(P);
#ifdef ALAP
    Matrix<T,ColMajor> minv2(minv.colsize(),minv.colsize());
    minv2 = LUx;
    NonLapLUInverse(minv2.View());
    minv2.ReversePermuteCols(P);
#endif
#endif

    if (minv.colsize() > 0) {
      if ( !(minv.iscm()
#ifndef ALAP
	    || minv.isrm()
#endif
	   )) {
	Matrix<T,ColMajor> temp(minv.colsize(),minv.colsize());
	LU_Inverse(LUx,P,temp.View());
	minv = temp;
      } else {
	minv = LUx;
#ifdef ALAP
	LapLUInverse(minv);
#else
	NonLapLUInverse(minv);
#endif
	//cerr<<"before permute: "<<minv<<endl;
	minv.ReversePermuteCols(P);
	//cerr<<"after permute: "<<minv<<endl;
      }
    }

#ifdef XDEBUG
    RealType(T) normdiff = Norm(m*minv - T(1));
    RealType(T) kappa = Norm(m)*Norm(minv);
    if (normdiff > 0.001*kappa*minv.colsize()) {
      cerr<<"LUInverse:\n";
      cerr<<"m = "<<m<<endl;
      cerr<<"LUx = "<<LUx<<endl;
      cerr<<"P = ";
      for(size_t i=0;i<LUx.colsize();i++) cerr<<P[i]<<" ";
      cerr<<endl;
      cerr<<"minv = "<<minv<<endl;
      cerr<<"m*minv = "<<m*minv<<endl;
      cerr<<"minv*m = "<<minv*m<<endl;
#ifdef ALAP
      cerr<<"Non-lap inverse = "<<minv2<<endl;
      cerr<<"m*minv2 = "<<m*minv2<<endl;
      cerr<<"minv2*m = "<<minv2*m<<endl;
#endif
      cerr<<"Norm(m*minv - 1) = "<<normdiff<<endl;
      cerr<<"kappa = "<<kappa<<endl;

      abort();
    }
#endif
  }


#define InstFile "TMV_LUDiv_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


