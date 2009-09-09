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


#include "tmv/TMV_QRD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "TMV_Householder.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define QR_BLOCKSIZE TMV_BLOCKSIZE
#else
#define QR_BLOCKSIZE 64
#endif

  //
  // QR Update
  //

  template <class T> static void NonBlockQR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);
    // Given that A0 = Q0 R0
    // Find R1, so that [ A0 ] = Q1 R1
    //                  [ A  ] 
    // Input R is R0, output is R1

    const int N = A.rowsize();

    T* Rdiag = R.ptr();
    const int ds = R.stepi()+R.stepj();
    T det(0);

    for(int j=0;j<N;++j,Rdiag+=ds) {
      // Apply the Householder Reflection for this column
      const VectorView<T> v = A.col(j);
      T beta = Householder_Reflect(*Rdiag,v,det);
      if (beta != T(0))
        Householder_LMult(v,beta,R.row(j,j+1,N),A.Cols(j+1,N));
    }
  }

  template <class T> static void RecursiveQR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A,
      const UpperTriMatrixView<T>& Z, bool makeZ)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const int N = A.rowsize();
    T det(0);

    TMVAssert(!R.isconj());
    TMVAssert(!Z.isconj());

    if (N==1) {
      T b = Householder_Reflect(*R.ptr(),A.col(0),det);
#ifdef TMVFLDEBUG
      TMVAssert(Z.ptr() >= Z.first);
      TMVAssert(Z.ptr() < Z.last);
#endif
      *Z.ptr() = CONJ(b);
    } else if (N==2) {
      T* R00 = R.ptr(); // = R(0,0)
      T* R01 = R00 + R.stepj(); // = R(0,1)
      T* R11 = R01 + R.stepi(); // = R(1,1)
      T* Z00 = Z.ptr(); // = Z(0,0)
      T* Z01 = Z00 + Z.stepj(); // = Z(0,1)
      T* Z11 = Z01 + Z.stepi(); // = Z(1,1)
      T b0 = Householder_Reflect(*R00,A.col(0),det);
      if (b0 != T(0)) {
        T temp = b0*(A.col(0).Conjugate()*A.col(1) + *R01);
#ifdef TMVFLDEBUG
        TMVAssert(R01 >= R.first);
        TMVAssert(R01 < R.last);
#endif
        *R01 -= temp;
        A.col(1) -= temp * A.col(0);
      }
#ifdef TMVFLDEBUG
      TMVAssert(Z00 >= Z.first);
      TMVAssert(Z00 < Z.last);
      TMVAssert(Z11 >= Z.first);
      TMVAssert(Z11 < Z.last);
#endif
      *Z00 = CONJ(b0);
      T b1 = Householder_Reflect(*R11,A.col(1),det);
      *Z11 = CONJ(b1);

      if (makeZ) {
        T temp = A.col(0).Conjugate()*A.col(1);
#ifdef TMVFLDEBUG
        TMVAssert(Z01 >= Z.first);
        TMVAssert(Z01 < Z.last);
#endif
        *Z01 = -(*Z00 * *Z11)*temp;
      }
    } else {
      int j1 = N/2;

      UpperTriMatrixView<T> R1 = R.SubTriMatrix(0,j1);
      MatrixView<T> Rx = R.SubMatrix(0,j1,j1,N);
      UpperTriMatrixView<T> R2 = R.SubTriMatrix(j1,N);

      MatrixView<T> A1 = A.Cols(0,j1);
      MatrixView<T> A2 = A.Cols(j1,N);

      UpperTriMatrixView<T> Z1 = Z.SubTriMatrix(0,j1);
      MatrixView<T> Zx = Z.SubMatrix(0,j1,j1,N);
      UpperTriMatrixView<T> Z2 = Z.SubTriMatrix(j1,N);

      RecursiveQR_Update(R1,A1,Z1,true);

      // Zx is a temporary here - it happens to be the right shape.
      Zx = A1.Adjoint() * A2; 
      Zx += Rx;
      Zx = Z1.Adjoint()*Zx;
      Rx -= Zx;
      A2 -= A1 * Zx;

      RecursiveQR_Update(R2,A2,Z2,makeZ);

      if (makeZ) {
        Zx = A1.Adjoint() * A2; 
        Zx = -Z1*Zx;
        Zx *= Z2;
      }
    }
  }

  template <class T> static void BlockQR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const int N = A.rowsize();

    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(
        MIN(QR_BLOCKSIZE,N));
    for(int j1=0;j1<N;) {
      int j2 = MIN(N,j1+QR_BLOCKSIZE);
      MatrixView<T> A1 = A.Cols(j1,j2);
      UpperTriMatrixView<T> R1 = R.SubTriMatrix(j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,j2-j1);

      RecursiveQR_Update(R1,A1,Z,j2<N);

      if (j2 < N) {
        Matrix<T,ColMajor> ZtYtm = A.Cols(j1,j2).Adjoint() * A.Cols(j2,N);
        ZtYtm += R.SubMatrix(j1,j2,j2,N);
        ZtYtm = Z.Adjoint() * ZtYtm;
        R.SubMatrix(j1,j2,j2,N) -= ZtYtm;
        A.Cols(j2,N) -= A.Cols(j1,j2) * ZtYtm;
      }
      j1 = j2;
    }
  }

  template <class T> void QR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

#ifdef XDEBUG
    Matrix<T> R0(R);
    Matrix<T> A0(A);
    UpperTriMatrix<T> R2(R);
    Matrix<T> A2(A);
    NonBlockQR_Update(R2.View(),A2.View());
#endif
    if (A.rowsize() > 0) {
      if (A.rowsize() > QR_BLOCKSIZE)
        BlockQR_Update(R,A);
      else {
        UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
        RecursiveQR_Update(R,A,Z.View(),false);
      }
    }
#ifdef XDEBUG
    if (Norm(R2-R) > 1.e-5*Norm(R0)*Norm(A0)) {
      cerr<<"QR_Update\n";
      cerr<<"R0 = "<<TypeText(R)<<"  "<<R0<<endl;
      cerr<<"A0 = "<<TypeText(A)<<"  "<<A0<<endl;
      cerr<<"R -> "<<R<<endl;
      cerr<<"NonBlock R -> "<<R2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_QRUpdate.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


