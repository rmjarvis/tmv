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


#include "TMV_Blas.h"
#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_RK_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_RK_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_RK_BLOCKSIZE 64
#define SYM_RK_BLOCKSIZE2 1
#endif

  // 
  // S += L*Lt
  //

  template <bool ha, bool uu, bool a1, class T, class TL> static void RecursiveRankKUpdate(
      const T alpha, const GenLowerTriMatrix<TL>& L, const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == L.size());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(ha == A.isherm());
    TMVAssert(uu == L.isunit());
    TMVAssert(a1 == (alpha == T(1)));

    int N = A.size();

    if (N == 1) {
#ifdef TMVFLDEBUG
      TMVAssert(A.ptr() >= A.first);
      TMVAssert(A.ptr() < A.last);
#endif
      if (uu)
        *(A.ptr()) += a1 ? T(1) : alpha;
      else 
        if (a1)
          *(A.ptr()) += ha ? NORM(*(L.cptr())) : SQR(*(L.cptr()));
        else
          *(A.ptr()) += alpha * (ha ? NORM(*(L.cptr())) : SQR(*(L.cptr())));
    } else {
      //
      // [ A11 A12 ] += alpha [ L11  0  ] [ L11t L21t ]
      // [ A21 A22 ]          [ L21 L22 ] [  0   L22t ]
      //              = alpha [ L11 L11t        L11 L21t      ]
      //                      [ L21 L11t  L21 L21t + L22 L22t ]
      int k = N/2;
      const int nb = SYM_RK_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;
      SymMatrixView<T> A11 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A22 = A.SubSymMatrix(k,N);
      MatrixView<T> A21 = A.SubMatrix(k,N,0,k);
      ConstLowerTriMatrixView<TL> L11 = L.SubTriMatrix(0,k);
      ConstLowerTriMatrixView<TL> L22 = L.SubTriMatrix(k,N);
      ConstMatrixView<TL> L21 = L.SubMatrix(k,N,0,k);

      RecursiveRankKUpdate<ha,uu,a1>(alpha,L22,A22);

      RankKUpdate<true>(alpha,L21,A22);

      A21 += alpha * L21 * (ha ? L11.Adjoint() : L11.Transpose());

      RecursiveRankKUpdate<ha,uu,a1>(alpha,L11,A11);
    }
  }

  template <bool a1, class T, class TL> static inline void DoRankKUpdate(
      const T alpha, const GenLowerTriMatrix<TL>& L, const SymMatrixView<T>& A)
  {
    if (A.isherm())
      if (L.isunit())
        RecursiveRankKUpdate<true,true,a1>(alpha,L,A);
      else
        RecursiveRankKUpdate<true,false,a1>(alpha,L,A);
    else
      if (L.isunit())
        RecursiveRankKUpdate<false,true,a1>(alpha,L,A);
      else
        RecursiveRankKUpdate<false,false,a1>(alpha,L,A);
  }

  template <bool ha, class T> static void RecursiveSetLLt(const SymMatrixView<T>& A)
  {
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(ha == A.isherm());

    int N = A.size();

    if (N == 1) {
#ifdef TMVFLDEBUG
      TMVAssert(A.ptr() >= A.first);
      TMVAssert(A.ptr() < A.last);
#endif
      *(A.ptr()) = ha ? NORM(*(A.cptr())) : SQR(*(A.cptr()));
    } else {
      // [ A11 A12 ] = [ L11  0  ] [ L11t L21t ]
      // [ A21 A22 ]   [ L21 L22 ] [  0   L22t ]
      //             = [ L11 L11t        L11 L21t      ]
      //               [ L21 L11t  L21 L21t + L22 L22t ]
      int k = N/2;
      const int nb = SYM_RK_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;
      SymMatrixView<T> A11 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A22 = A.SubSymMatrix(k,N);
      MatrixView<T> A21 = A.SubMatrix(k,N,0,k);
      ConstLowerTriMatrixView<T> L11 = A11.LowerTri();

      RecursiveSetLLt<ha>(A22);

      RankKUpdate<true>(T(1),A21,A22);

      A21 *= (ha ? L11.Adjoint() : L11.Transpose());

      RecursiveSetLLt<ha>(A11);
    }
  }

  template <class T> static inline void SetLLt(const SymMatrixView<T>& A)
  {
    if (A.isherm()) RecursiveSetLLt<true>(A);
    else RecursiveSetLLt<false>(A);
  }

  template <bool add, class T, class TL> void RankKUpdate(
      const T alpha, const GenLowerTriMatrix<TL>& L, 
      const SymMatrixView<T>& A)
  // A = A + alpha * L * LT
  {
#ifdef XTEST
    TMVAssert(!add || A.HermOK());
#endif
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<TL> L0 = L;
    Matrix<TL> Lt = A.isherm() ? L.Adjoint() : L.Transpose();
    Matrix<T> A2 = A;
    if (add) A2 += alpha*L*Lt;
    else A2 = alpha*L*Lt;
    //cout<<"Start RankKUpdate: A = "<<TypeText(A)<<"  "<<A<<endl;
    //cout<<"alpha = "<<alpha<<", L = "<<TypeText(L)<<"  "<<L<<endl;
#endif

    TMVAssert(A.size() == L.size());
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());

    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
        RankKUpdate<add>(CONJ(alpha),L.Conjugate(),A.Conjugate());
      else if (A.uplo() == Upper) 
        if (A.isherm()) RankKUpdate<add>(alpha,L,A.Adjoint());
        else RankKUpdate<add>(alpha,L,A.Transpose());
      else if (!add) {
        A.LowerTri() = L;
        SetLLt(A);
        A *= alpha;
      }
      else
        if (alpha == T(1)) DoRankKUpdate<true>(alpha,L,A);
        else DoRankKUpdate<false>(alpha,L,A);
    }

#ifdef XDEBUG
    //cout<<"Done RankK\n";
    if (Norm(A-A2) > 0.001*(ABS(alpha)*SQR(Norm(L0))+
          (add?Norm(A0):RealType(T)(0)))) {
      cerr<<"RankKUpdate: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"L = "<<TypeText(L)<<"  "<<L0<<endl;
      cerr<<"A = "<<TypeText(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
  }

#define InstFile "TMV_RankK_LUS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


