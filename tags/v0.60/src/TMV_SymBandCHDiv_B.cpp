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



#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_TriDiv.h"
#include "TMV_DiagMatrix.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_SymBandCHDiv.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_DiagMatrixArith.h"
#include "TMV_SymBandMatrixArith.h"
#include "TMV_SymMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_BandMatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define TRI_DIV_BLOCKSIZE TMV_BLOCKSIZE
#else
#define TRI_DIV_BLOCKSIZE 64
#endif

  // This is a copy of RecursiveInverse in TMV_TriDiv_D.cpp
  // but using the fact that U is banded
  template <bool unit, class T> inline void RecursiveBandInverse(
      const UpperTriMatrixView<T>& U, size_t nhi)
  {
    TMVAssert(U.iscm() || U.isrm());
    TMVAssert(unit == U.isunit());

    const size_t N = U.size();
    const size_t nb = TRI_DIV_BLOCKSIZE;

    if (N == 1) {
      if (!unit) {
	T*const Uptr = U.ptr();
	if (*Uptr == T(0))
	  throw SingularUpperTriMatrix<T>(U);
	*Uptr = RealType(T)(1) / (*Uptr);
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      UpperTriMatrixView<T> U00 = U.SubTriMatrix(0,k);
      MatrixView<T> U01 = U.SubMatrix(0,k,k,N);
      // U01 only has non-zero elements in the lower-left nhixnhi triangle
      MatrixView<T> U01x = (k > nhi && N-k > nhi) ? 
	U01.Cols(0,nhi) : U01;
      UpperTriMatrixView<T> U11 = U.SubTriMatrix(k,N);

      // U00 U01' + U01 U11' = 0
      // U00 U01' = -U01 U11'
      // U01' = -U00' U01 U11'

      RecursiveBandInverse<unit>(U00,nhi);
      RecursiveBandInverse<unit>(U11,nhi);
      U01x = -U00 * U01x;
      U01 *= U11;
    }
  }
								       

  template <class T> inline void DoHermBandCHInverse(
      const SymMatrixView<T>& sinv, int nlo)
  {
    TMVAssert(sinv.isherm());
    // inv = (L Lt)^-1 = Lt^-1 L^-1
    LowerTriMatrixView<T> L = sinv.LowerTri();
    RecursiveBandInverse<false>(L.Transpose(),nlo);
    sinv = L.Adjoint() * L;
  }

  template <class T> inline void SimpleLDLt_AddXtDX(
      const SymMatrixView<T>& sinv, 
      const GenMatrix<T>& X, const GenDiagMatrix<T>& D)
  {
    const size_t N = D.size();
    if (N==1) {
      sinv += D(0) * X.row(0).Conjugate() ^ X.row(0);
    } else {
      size_t No2 = N/2;
      SimpleLDLt_AddXtDX(sinv,X.Rows(0,No2),D.SubDiagMatrix(0,No2));
      SimpleLDLt_AddXtDX(sinv,X.Rows(No2,N),D.SubDiagMatrix(No2,N));
    }
  }
  
  template <class T> inline void SimpleLDLt_CombineInverse(
      const SymMatrixView<T>& sinv)
  {
    // This is basically the same algorithm as the LDLt_CombineInverse
    // in TMV_SymLUDiv_C.cpp for the Bunch-Kauffman inverse.
    // The difference, and why this is the "simple" version, is that
    // D really is a regular DiagMatrix - not the screwy block diagonal
    // version that we had to deal with there.

    const size_t N = sinv.size();
    if (N > 1) {
      size_t No2 = N/2;
      MatrixView<T> X10 = sinv.SubMatrix(No2,N,0,No2);
      LowerTriMatrixView<T> X11 = sinv.LowerTri(UnitDiag).SubTriMatrix(No2,N);
      DiagMatrixView<T> D(sinv.diag());
      
      SimpleLDLt_CombineInverse(sinv.SubSymMatrix(0,No2));
      SimpleLDLt_AddXtDX(sinv.SubSymMatrix(0,No2),X10,D.SubDiagMatrix(No2,N));
      X10 = D.SubDiagMatrix(No2,N) * X10;
      X10 = X11.Adjoint() * X10;
      SimpleLDLt_CombineInverse(sinv.SubSymMatrix(No2,N));
    }
  }
	  
  template <class T> inline void DoHermTriDiagCHInverse(
      const SymMatrixView<T>& sinv)
  {
    LowerTriMatrixView<T> L = sinv.LowerTri(UnitDiag);
    DiagMatrixView<T> D(sinv.diag());
    RecursiveBandInverse<true>(L.Transpose(),1);
    D.Real().InvertSelf();
    SimpleLDLt_CombineInverse(sinv);
  }


  template <class T, class T1> void HermBandCH_Inverse(
      const GenSymBandMatrix<T>& LLx,
      const SymMatrixView<T1>& sinv) 
  {
#ifdef XTEST
    TMVAssert(LLx.HermOK());
#endif
    TMVAssert(LLx.size() == sinv.size());

#ifdef XDEBUG
    Matrix<T> A(LLx.size(),LLx.size());
    if (LLx.nlo() <= 1) {
      BandMatrix<T> L = LLx.LowerBand();
      L.diag().SetAllTo(T(1));
      DiagMatrix<T> D(LLx.diag());
      A = L;
      A *= D;
      A *= L.Adjoint();
    } else {
      BandMatrix<T> L = LLx.LowerBand();
      A = L * L.Adjoint();
    }
#endif

    if (sinv.size() > 0) {
      if (LLx.nlo() == 0) {
	sinv.Zero();
	sinv.diag() = LLx.diag();
	DiagMatrixViewOf(sinv.diag().Real()).InvertSelf();
      } else if (!( sinv.iscm() || sinv.isrm() )) {
	HermMatrix<T,Lower,ColMajor> temp(sinv.size());
	HermBandCH_Inverse(LLx,temp.View());
	sinv = temp;
      } else {
	sinv = LLx;
	if (LLx.nlo() == 1)
	  DoHermTriDiagCHInverse(sinv);
	else
	  DoHermBandCHInverse(sinv,LLx.nlo());
      }
    }

#ifdef XDEBUG
    Matrix<T1> eye = A * sinv;
    RealType(T) kappa = Norm(A) * Norm(sinv);
    if (Norm(eye-T(1)) > 0.0001*kappa*sinv.size()) {
      cerr<<"A = "<<A<<endl;
      cerr<<"sinv = "<<sinv<<endl;
      cerr<<"A*sinv = "<<A*sinv<<endl;
      cerr<<"sinv*A = "<<sinv*A<<endl;
      cerr<<"Norm(A*sinv-1) = "<<Norm(A*sinv-T(1))<<endl;
      cerr<<"kappa = "<<kappa<<endl;
      abort();
    }
#endif
#ifdef XTEST
    TMVAssert(sinv.HermOK());
#endif
  }

  
#define InstFile "TMV_SymBandCHDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


