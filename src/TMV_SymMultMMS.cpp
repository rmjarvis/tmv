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


#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include <iostream>
#include "tmv/TMV_SymMatrixArith.h"
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_R2K_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_R2K_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_R2K_BLOCKSIZE 64
#define SYM_R2K_BLOCKSIZE2 1
#endif

  // 
  // SymMultMM
  // A += alpha * x * y
  // where x,y are Matrices, and the product x*y is assumed to be symmetric.
  //

  template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
  static void RecursiveSymMultMM(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.rowsize());
    TMVAssert(x.rowsize() == y.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(ha == A.isherm());
    TMVAssert(a1 == (alpha == T(1)));

    const int nb = SYM_R2K_BLOCKSIZE;
    int N = A.size();

    if (N <= SYM_R2K_BLOCKSIZE2) {
      if (N == 1) {
        T temp = x.row(0) * y.col(0);
        if (!a1) temp *= alpha;
#ifdef TMVFLDEBUG
        TMVAssert(A.ptr() >= A.first);
        TMVAssert(A.ptr() < A.last);
#endif
        if (ha)
          if (add) *(A.ptr()) += REAL(temp);
          else *(A.ptr()) = REAL(temp);
        else
          if (add) *(A.ptr()) += temp;
          else *(A.ptr()) = temp;
      } else {
        if (A.isrm()) {
          for (int i=0;i<N;++i) {
            if (add) 
              A.row(i,0,i+1) += alpha * x.row(i) * y.Cols(0,i+1);
            else 
              A.row(i,0,i+1) = alpha * x.row(i) * y.Cols(0,i+1);
          }
        } else {
          for (int j=0;j<N;++j) {
            if (add) 
              A.col(j,j,N) += alpha * x.Rows(j,N) * y.col(j);
            else 
              A.col(j,j,N) = alpha * x.Rows(j,N) * y.col(j);
          }
        }
        if (ha && IsComplex(T())) {
#ifdef XTEST
          TMVAssert(NormInf(A.diag().Imag()) < A.size()*Epsilon<T>()*(Norm(A)+Norm(x)*Norm(y)));
#endif
          A.diag().Imag().Zero();
        }
      }
    } else { // Not <= BLOCKSIZE2, so do recurse...
      int k = N/2;
      if (k > nb) k = k/nb*nb;

      RecursiveSymMultMM<ha,a1,add>(alpha,x.Rows(0,k),y.Cols(0,k),
          A.SubSymMatrix(0,k));

      if (add) 
        A.SubMatrix(k,N,0,k) += alpha * x.Rows(k,N) * y.Cols(0,k);
      else 
        A.SubMatrix(k,N,0,k) = alpha * x.Rows(k,N) * y.Cols(0,k);

      RecursiveSymMultMM<ha,a1,add>(alpha,x.Rows(k,N),y.Cols(k,N),
          A.SubSymMatrix(k,N));
    }
  }

  template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
  static void RecursiveInPlaceSymMultMM(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
  {
    TMVAssert(SameStorage(x,A) || SameStorage(y,A));
    TMVAssert(A.size() > 0);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.ct() == NonConj);

    int N = A.size();
    if (N == 1) {
      Tx x00 = x(0,0);
      Ty y00 = y(0,0);
      if (ha) {
        RealType(T) temp = a1 ? REAL(x00*y00) : REAL(alpha*x00*y00);
#ifdef TMVFLDEBUG
        TMVAssert(A.ptr() >= A.first);
        TMVAssert(A.ptr() < A.last);
#endif
        if (add) *A.ptr() += temp;
        else *A.ptr() = temp;
      } else {
        T temp = a1 ? (x00 * y00) : (alpha * (x00 * y00));
#ifdef TMVFLDEBUG
        TMVAssert(A.ptr() >= A.first);
        TMVAssert(A.ptr() < A.last);
#endif
        if (add) *A.ptr() += temp;
        else *A.ptr() = temp;
      }
    } else {
      const int k = N/2;
      const ConstMatrixView<Tx> x00 = x.SubMatrix(0,k,0,k);
      const ConstMatrixView<Tx> x10 = x.SubMatrix(k,N,0,k);
      const ConstMatrixView<Tx> x01 = x.SubMatrix(0,k,k,N);
      const ConstMatrixView<Tx> x11 = x.SubMatrix(k,N,k,N);
      const ConstMatrixView<Ty> y00 = y.SubMatrix(0,k,0,k);
      const ConstMatrixView<Ty> y10 = y.SubMatrix(k,N,0,k);
      const ConstMatrixView<Ty> y01 = y.SubMatrix(0,k,k,N);
      const ConstMatrixView<Ty> y11 = y.SubMatrix(k,N,k,N);
      SymMatrixView<T> A00 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A11 = A.SubSymMatrix(k,N);
      MatrixView<T> A10 = A.SubMatrix(k,N,0,k);

      Matrix<T> tempA10 = x10 * y00;
      tempA10 += x11 * y10;

      RecursiveInPlaceSymMultMM<ha,a1,add>(alpha,x11,y11,A11);
      RecursiveSymMultMM<ha,a1,true>(alpha,x10,y01,A11);
      RecursiveInPlaceSymMultMM<ha,a1,add>(alpha,x00,y00,A00);
      RecursiveSymMultMM<ha,a1,true>(alpha,x01,y10,A00);

      if (add) A10 += alpha * tempA10;
      else A10 = alpha * tempA10;
    }
  }

  template <bool add, class T, class Tx, class Ty> 
  static inline void InPlaceSymMultMM(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
  {
    if (A.isherm())
      if (alpha == T(1))
        RecursiveInPlaceSymMultMM<true,true,add>(alpha,x,y,A); 
      else
        RecursiveInPlaceSymMultMM<true,false,add>(alpha,x,y,A); 
    else
      if (alpha == T(1))
        RecursiveInPlaceSymMultMM<false,true,add>(alpha,x,y,A); 
      else
        RecursiveInPlaceSymMultMM<false,false,add>(alpha,x,y,A); 
  }

  template <bool add, class T, class Tx, class Ty> static void DoSymMultMM(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
  { 
    if (SameStorage(x,A) || SameStorage(y,A)) {
      const int N = A.size();
      TMVAssert(int(x.colsize()) == N);
      TMVAssert(int(y.rowsize()) == N);
      const int K = x.rowsize();
      TMVAssert(int(y.colsize()) == K);

      if (K >= N) {
        InPlaceSymMultMM<add>(alpha,x.Cols(0,N),y.Rows(0,N),A);
        if (K > N) DoSymMultMM<add>(alpha,x.Cols(N,K),y.Rows(N,K),A);
      } else { 
        DoSymMultMM<add>(alpha,x.Rows(K,N),y.Cols(K,N),
            A.SubSymMatrix(K,N));
        if ((y.stepi() < y.stepj()) != (A.stepi() < A.stepj())) {
          if ((x.stepi() < x.stepj()) == (A.stepi() < A.stepj())) {
            // Then both x and y overlap with A(K:N,0:K)
            // Need a temporary
            if (A.iscm()) {
              Matrix<T,ColMajor> temp = x.Rows(K,N) * y.Cols(0,K);
              if (add) A.SubMatrix(K,N,0,K) += alpha * temp;
              else A.SubMatrix(K,N,0,K) = alpha * temp;
            } else {
              Matrix<T,RowMajor> temp = x.Rows(K,N) * y.Cols(0,K);
              if (add) A.SubMatrix(K,N,0,K) += alpha * temp;
              else A.SubMatrix(K,N,0,K) = alpha * temp;
            }
          } else {
            MultMM<add>(alpha,x.SubMatrix(K,N,K,N),y.SubMatrix(K,N,0,K),
                A.SubMatrix(K,N,0,K) );
            MultMM<true>(alpha,x.SubMatrix(K,N,0,K),y.SubMatrix(0,K,0,K),
                A.SubMatrix(K,N,0,K) );
          }
        }
        else {
          MultMM<add>(alpha,x.Rows(K,N),y.Cols(0,K),A.SubMatrix(K,N,0,K) );
        }
        InPlaceSymMultMM<add>(alpha,x.Rows(0,K),y.Rows(0,K),
            A.SubSymMatrix(0,K));
      }
    }
    else 
      if (A.isherm())
        if (alpha == T(1))
          RecursiveSymMultMM<true,true,add>(alpha,x,y,A); 
        else
          RecursiveSymMultMM<true,false,add>(alpha,x,y,A); 
      else
        if (alpha == T(1))
          RecursiveSymMultMM<false,true,add>(alpha,x,y,A); 
        else
          RecursiveSymMultMM<false,false,add>(alpha,x,y,A); 
  }

  template <bool add, class T, class Tx, class Ty> void SymMultMM(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
  // A (+)= alpha * x * y
  {
#ifdef XTEST
    TMVAssert(!add || A.HermOK());
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.rowsize());
    TMVAssert(x.rowsize() == y.colsize());

#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<Tx> x0 = x;
    Matrix<Ty> y0 = y;
    Matrix<T> A2 = A;
    if (add) A2 += alpha*x*y;
    else A2 = alpha*x*y;
    //cout<<"Start SymMultMM: alpha = "<<alpha<<endl;
    //cout<<"add = "<<add<<endl;
    //cout<<"x = "<<TypeText(x)<<"  "<<x0<<endl;
    //cout<<"y = "<<TypeText(y)<<"  "<<y0<<endl;
    //cout<<"A = "<<TypeText(A)<<"  "<<A0<<endl;
#endif

    if (alpha != T(0) && A.size() > 0) {
      if (A.uplo() == Upper) {
        if (A.isherm()) SymMultMM<add>(alpha,x,y,A.Adjoint());
        else SymMultMM<add>(alpha,x,y,A.Transpose());
      }
      else if (A.isconj()) 
        SymMultMM<add>(CONJ(alpha),x.Conjugate(),y.Conjugate(),A.Conjugate());
      else DoSymMultMM<add>(alpha,x,y,A);
    }

#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*(ABS(alpha)*Norm(x0)*Norm(y0)+
          (add?Norm(A0):RealType(T)(0)))) {
      cerr<<"SymMultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"x = "<<TypeText(x)<<"  "<<x0<<endl;
      cerr<<"y = "<<TypeText(y)<<"  "<<y0<<endl;
      cerr<<"A = "<<TypeText(A)<<"  "<<A0<<endl;
      cerr<<"x*y = "<<x0*y0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
  }

#define InstFile "TMV_SymMultMMS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


