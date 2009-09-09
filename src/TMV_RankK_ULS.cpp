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
  // S += U*Ut
  //

  template <bool ha, bool uu, bool a1, class T, class TU> static void RecursiveRankKUpdate(
      const T alpha, const GenUpperTriMatrix<TU>& U, const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == U.size());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Upper);
    TMVAssert(ha == A.isherm());
    TMVAssert(uu == U.isunit());
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
          *(A.ptr()) += ha ? NORM(*(U.cptr())) : SQR(*(U.cptr()));
        else
          *(A.ptr()) += alpha * (ha ? NORM(*(U.cptr())) : SQR(*(U.cptr())));
    } else {
      //
      // [ A11 A12 ] += alpha [ U11 U12 ] [ U11t  0   ]
      // [ A21 A22 ]          [  0  U22 ] [ U12t U22t ]
      //              = alpha [ U11 U11t + U12 U12t    U12 U22t ]
      //                      [       U22 U12t         U22 U22t ]
      int k = N/2;
      const int nb = SYM_RK_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;
      SymMatrixView<T> A11 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A22 = A.SubSymMatrix(k,N);
      MatrixView<T> A12 = A.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<TU> U11 = U.SubTriMatrix(0,k);
      ConstUpperTriMatrixView<TU> U22 = U.SubTriMatrix(k,N);
      ConstMatrixView<TU> U12 = U.SubMatrix(0,k,k,N);

      RecursiveRankKUpdate<ha,uu,a1>(alpha,U11,A11);

      RankKUpdate<true>(alpha,U12,A11);

      A12 += alpha * U12 * (ha ? U22.Adjoint() : U22.Transpose());

      RecursiveRankKUpdate<ha,uu,a1>(alpha,U22,A22);
    }
  }

  template <bool a1, class T, class TU> static inline void DoRankKUpdate(
      const T alpha, const GenUpperTriMatrix<TU>& U, const SymMatrixView<T>& A)
  {
    if (A.isherm())
      if (U.isunit())
        RecursiveRankKUpdate<true,true,a1>(alpha,U,A);
      else
        RecursiveRankKUpdate<true,false,a1>(alpha,U,A);
    else
      if (U.isunit())
        RecursiveRankKUpdate<false,true,a1>(alpha,U,A);
      else
        RecursiveRankKUpdate<false,false,a1>(alpha,U,A);
  }

  template <bool ha, class T> static void RecursiveSetUUt(const SymMatrixView<T>& A)
  {
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Upper);
    TMVAssert(ha == A.isherm());

    int N = A.size();

    if (N == 1) {
#ifdef TMVFLDEBUG
      TMVAssert(A.ptr() >= A.first);
      TMVAssert(A.ptr() < A.last);
#endif
      *(A.ptr()) = ha ? NORM(*(A.cptr())) : SQR(*(A.cptr()));
    } else {
      //
      // [ A11 A12 ] = alpha [ U11 U12 ] [ U11t  0   ]
      // [ A21 A22 ]         [  0  U22 ] [ U12t U22t ]
      //             = alpha [ U11 U11t + U12 U12t    U12 U22t ]
      //                     [       U22 U12t         U22 U22t ]
      int k = N/2;
      const int nb = SYM_RK_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;
      SymMatrixView<T> A11 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A22 = A.SubSymMatrix(k,N);
      MatrixView<T> A12 = A.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<T> U22 = A22.UpperTri();

      RecursiveSetUUt<ha>(A11);

      // The Transpose's here are because these RankKUpdate routines
      // want A11 to be stored in the Lower triangle.
      if (ha)
        RankKUpdate<true>(T(1),A12.Conjugate(),A11.Transpose());
      else
        RankKUpdate<true>(T(1),A12,A11.Transpose());

      A12 *= (ha ? U22.Adjoint() : U22.Transpose());

      RecursiveSetUUt<ha>(A22);
    }
  }

#ifdef AELAP
  template <class T> static inline void LapSetUUt(const SymMatrixView<T>& A)
  { RecursiveSetUUt<false>(A); }
#ifdef INST_DOUBLE
  template<> void LapSetUUt(const SymMatrixView<double>& A)
  {
    TMVAssert(A.issym());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Upper);

    int n = A.size();
    int lda = A.isrm() ? A.stepi() : A.stepj();
    LAPNAME(dlauum) (LAPCM A.isrm()?LAPCH_LO:LAPCH_UP,LAPV(n),
        LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("dlauum");
  }
  template<> void LapSetUUt(
      const SymMatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.issym());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Upper);

    int n = A.size();
    int lda = A.isrm() ? A.stepi() : A.stepj();
    LAPNAME(zlauum) (LAPCM A.isrm()?LAPCH_LO:LAPCH_UP,LAPV(n),
        LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("zlauum");
  }
#endif
#ifdef INST_FLOAT
  template<> void LapSetUUt(const SymMatrixView<float>& A)
  {
    TMVAssert(A.issym());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Upper);

    int n = A.size();
    int lda = A.isrm() ? A.stepi() : A.stepj();
    LAPNAME(slauum) (LAPCM A.isrm()?LAPCH_LO:LAPCH_UP,LAPV(n),
        LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("slauum");
  }
  template<> void LapSetUUt(
      const SymMatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.issym());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Upper);

    int n = A.size();
    int lda = A.isrm() ? A.stepi() : A.stepj();
    LAPNAME(clauum) (LAPCM A.isrm()?LAPCH_LO:LAPCH_UP,LAPV(n),
        LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("clauum");
  }
#endif // FLOAT
#endif // AELAP

  template <class T> static inline void SetUUt(const SymMatrixView<T>& A)
  {
#ifdef AELAP
    if (A.issym() && ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)))
      LapSetUUt(A);
    else
#endif
      if (A.isherm()) RecursiveSetUUt<true>(A);
      else RecursiveSetUUt<false>(A);
  }

  template <bool add, class T, class TU> void RankKUpdate(
      const T alpha, const GenUpperTriMatrix<TU>& U,
      const SymMatrixView<T>& A)
  // A = A + alpha * U * UT
  {
#ifdef XTEST
    TMVAssert(!add || A.HermOK());
#endif
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<TU> U0 = U;
    Matrix<TU> Ut = A.isherm() ? U.Adjoint() : U.Transpose();
    Matrix<T> A2 = A;
    if (add) A2 += alpha*U*Ut;
    else A2 = alpha*U*Ut;
    //cout<<"Start RankKUpdate: A = "<<TypeText(A)<<"  "<<A<<endl;
    //cout<<"alpha = "<<alpha<<", U = "<<TypeText(U)<<"  "<<U<<endl;
#endif

    TMVAssert(A.size() == U.size());
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());

    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
        RankKUpdate<add>(CONJ(alpha),U.Conjugate(),A.Conjugate());
      else if (A.uplo() == Lower) 
        if (A.issym()) RankKUpdate<add>(alpha,U,A.Transpose());
        else RankKUpdate<add>(CONJ(alpha),U.Conjugate(),A.Transpose());
      else if (!add) {
        A.UpperTri() = U;
        SetUUt(A);
        if (alpha != T(1)) A *= alpha;
      }
      else
        if (alpha == T(1)) DoRankKUpdate<true>(alpha,U,A);
        else DoRankKUpdate<false>(alpha,U,A);
    }

#ifdef XDEBUG
    //cout<<"Done RankK\n";
    if (Norm(A-A2) > 0.001*(ABS(alpha)*SQR(Norm(U0))+
          add?Norm(A0):RealType(T)(0))) {
      cerr<<"RankKUpdate: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"U = "<<TypeText(U)<<"  "<<U0<<endl;
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

#define InstFile "TMV_RankK_ULS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


