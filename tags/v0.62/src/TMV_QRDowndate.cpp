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
#include "TMV_Householder.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef NOTHROW
#include <iostream>
#endif
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
  // QR Downdate
  //

#ifndef NOTHROW
  template <class T> class Bad_QR_Downdate :
    public NonPosDef
  {
  public:
    mutable auto_ptr<UpperTriMatrix<T> > R;
    mutable auto_ptr<Matrix<T> > A;

    Bad_QR_Downdate(
        const GenUpperTriMatrix<T>& _R, const GenMatrix<T>& _A) :
      NonPosDef("QR Downdate."), 
      R(new UpperTriMatrix<T>(_R)), A(new Matrix<T>(_A)) {}
    Bad_QR_Downdate(const Bad_QR_Downdate<T>& rhs) :
      R(rhs.R), A(rhs.A) {}
    ~Bad_QR_Downdate() throw() {}

    void Write(std::ostream& os) const throw()
    {
      os<<"TMV NonPosDef: QR Downdate found that the resulting "<<std::endl;
      os<<"down-dated RtR is not positive definite. "<<std::endl;
      os<<"(and hence the down date is impossible)"<<std::endl;
      os<<"The partially downdated matrix is \n"<<*R<<std::endl;
      os<<"The matrix attempting to be down-dated was \n"<<*A<<std::endl;
    }
  };
#endif

  template <class T> static void NonBlockQR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);
    // Given that A0 = Q0 R0
    // Given that [ A0 ] = Q1 R1
    //            [ A  ] 
    // Find R0 so that A0 = Q0 R0
    // Input R is R1, output is R0

    const int N = A.rowsize();

    T* Rdiag = R.ptr();
    const int ds = R.stepi()+R.stepj();

    for(int j=0;j<N;++j,Rdiag+=ds) {
      // Apply the Householder Reflection for this column
      const VectorView<T> v = A.col(j);
      T beta;
      if (!(Householder_UnReflect(*Rdiag,v,beta))) 
#ifdef NOTHROW
      { std::cerr<<"Bad QR Downdate\n"; exit(1); }
#else
      throw Bad_QR_Downdate<T>(R,A);
#endif
      TMVAssert(beta != T(1));
      VectorView<T> m0 = R.row(j,j+1,N);
      MatrixView<T> mx = A.Cols(j+1,N);

      // m0' = m0 - beta m0 - beta vtmx
      // m0 = (m0' + beta btmx)/(1-beta)
      Vector<T> bvtmx = beta*v.Conjugate()*mx;
      m0 += bvtmx;
      m0 /= T(1)-beta;

      // mx' = mx - beta v (m0 + vtmx)
      bvtmx += beta*m0;
      mx -= v ^ bvtmx;
    }
  }

  template <class T> static void RecursiveQR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A,
      const UpperTriMatrixView<T>& Z, bool makeZ)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const int N = A.rowsize();

    if (N==1) {
      T b;
      if (!(Householder_UnReflect(*R.ptr(),A.col(0),b))) 
#ifdef NOTHROW
      { std::cerr<<"Bad QR Downdate\n"; exit(1); }
#else
      throw Bad_QR_Downdate<T>(R,A);
#endif
#ifdef TMVFLDEBUG
      TMVAssert(Z.ptr() >= Z.first);
      TMVAssert(Z.ptr() < Z.last);
#endif
      *Z.ptr() = CONJ(b);
    } else if (N==2) {
      T* R00 = R.ptr();
      T* R01 = R00+R.stepj();
      T* R11 = R01+R.stepi();
      T* Z00 = Z.ptr();
      T* Z01 = Z00+Z.stepj();
      T* Z11 = Z01+Z.stepi();
      T b0;
      if (!(Householder_UnReflect(*R00,A.col(0),b0))) 
#ifdef NOTHROW
      { std::cerr<<"Bad QR Downdate\n"; exit(1); }
#else
      throw Bad_QR_Downdate<T>(R,A);
#endif
#ifdef TMVFLDEBUG
      TMVAssert(Z00 >= Z.first);
      TMVAssert(Z00 < Z.last);
#endif
      *Z00 = CONJ(b0);
      if (b0 != T(0)) {
        TMVAssert(b0 != T(1));
        T vtmx = A.col(0).Conjugate() * A.col(1);
#ifdef TMVFLDEBUG
        TMVAssert(R01 >= R.first);
        TMVAssert(R01 < R.last);
#endif
        *R01 = (*R01 + b0*vtmx)/(T(1)-b0);
        A.col(1) -= b0*(vtmx + *R01) * A.col(0);
      }

      T b1;
      if (!(Householder_UnReflect(*R11,A.col(1),b1))) 
#ifdef NOTHROW
      { std::cerr<<"Bad QR Downdate\n"; exit(1); }
#else
      throw Bad_QR_Downdate<T>(R,A);
#endif
#ifdef TMVFLDEBUG
      TMVAssert(Z11 >= Z.first);
      TMVAssert(Z11 < Z.last);
#endif
      *Z11 = CONJ(b1);

      if (makeZ) {
        T vtmx = A.col(0).Conjugate() * A.col(1);
#ifdef TMVFLDEBUG
        TMVAssert(Z01 >= Z.first);
        TMVAssert(Z01 < Z.last);
#endif
        *Z01 = -CONJ(b0*b1)*vtmx;
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

#ifndef NOTHROW
      try {
#endif
        RecursiveQR_Downdate(R1,A1,Z1,true);
#ifndef NOTHROW
      }
      catch (Bad_QR_Downdate<T>) {
        throw Bad_QR_Downdate<T>(R,A);
      }
#endif

      Zx = A1.Adjoint() * A2;
      Zx = Z1.Adjoint() * Zx;

      Rx += Zx;
      LowerTriMatrix<T> ImZt = T(1)-Z1.Adjoint();
      Rx /= ImZt;

      Zx += Z1.Adjoint() * Rx;
      A2 -= A1 * Zx;

#ifndef NOTHROW
      try {
#endif
        RecursiveQR_Downdate(R2,A2,Z2,makeZ);
#ifndef NOTHROW
      }
      catch (Bad_QR_Downdate<T>) {
        throw Bad_QR_Downdate<T>(R,A);
      }
#endif

      if (makeZ) {
        Zx = A1.Adjoint() * A2; 
        Zx = -Z1*Zx;
        Zx *= Z2;
      }
    }
  }

  template <class T> static void BlockQR_Downdate(
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

#ifndef NOTHROW
      try {
#endif
        RecursiveQR_Downdate(R1,A1,Z,j2<N);
#ifndef NOTHROW
      }
      catch (Bad_QR_Downdate<T>) {
        throw Bad_QR_Downdate<T>(R,A);
      }
#endif

      if (j2 < N) {
        // m0' = m0 - Zt(Ytmx+m0)
        // m0' + ZtYtmx = (I-Zt) m0;
        MatrixView<T> m0 = R.SubMatrix(j1,j2,j2,N);
        MatrixView<T> Y = A.Cols(j1,j2);
        MatrixView<T> mx = A.Cols(j2,N);

        Matrix<T,ColMajor> ZtYtm = Y.Adjoint() * mx;
        ZtYtm = Z.Adjoint() * ZtYtm;

        m0 += ZtYtm;
        LowerTriMatrix<T> ImZt = T(1)-Z.Adjoint();
        m0 /= ImZt;

        ZtYtm += Z.Adjoint() * m0;
        mx -= Y * ZtYtm;
      }
      j1 = j2;
    }
  }

  template <class T> void QR_Downdate(
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
    NonBlockQR_Downdate(R2.View(),A2.View());
#endif

    if (A.rowsize() > 0) {
      if (A.rowsize() > QR_BLOCKSIZE)
        BlockQR_Downdate(R,A);
      else {
        UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
        RecursiveQR_Downdate(R,A,Z.View(),false);
      }
    }

#ifdef XDEBUG
    if (Norm(R2-R) > 1.e-5*Norm(A0)*Norm(R0)) {
      cerr<<"QR_Downdate\n";
      cerr<<"R0 = "<<TypeText(R)<<"  "<<R0<<endl;
      cerr<<"A0 = "<<TypeText(A)<<"  "<<A0<<endl;
      cerr<<"R -> "<<R<<endl;
      cerr<<"NonBlock R -> "<<R2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_QRDowndate.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


