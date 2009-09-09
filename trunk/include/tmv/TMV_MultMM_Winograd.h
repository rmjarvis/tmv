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


#ifndef TMV_MultMM_Winograd_H
#define TMV_MultMM_Winograd_H

#include "TMV_Matrix.h"
#include "TMV_AddMM.h"
#include "TMV_MultXM.h"

namespace tmv {

  template <class T, class M1, class M2, class M3>
  static void MultMM_Winograd_Final(
      const bool add, const T x, const M1& m1, const M2& m2, M3& m3)
  {
#ifdef TMV_USE_RECURSIVE_BLOCK
    DoRecursiveBlockMultMM(add,x,m1,m2,m3);
#else
    DoBlockMultMM(add,x,m1,m2,m3);
#endif
  }

#ifdef TMV_CompilingWinograd
#define TMV_STATIC
#else
#define TMV_STATIC static
#endif

  template <class T1, bool CC1, class T2, bool CC2, class T3>
  static void DoMultMM_Winograd(
      const bool add, const T3& x, 
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<T2,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<T3>& m3)
  {
    // This is a recursive algorithm that avoids some of the 
    // matrix multiplication steps by using some clever algebra.
    // This particular implementation of the Winograd algorithm
    // is applicable to non-square matrices and transitions to 
    // a regular algorithm at a particular recursion threshold.
    //
    // In this code, each matrix is divided into 4 parts:
    // m1 == A = ( A0  A1 )
    //           ( A2  A3 )
    // m2 == B = ( B0  B1 )
    //           ( B2  B3 )
    // m3 == C = ( C0  C1 )
    //           ( C2  C3 )
    // 
    // The basic algorithm in pseudo-code, taken directly from the
    // Alberto and Nicolau paper (Algorithm 1) is technically
    // only applicable if M,N,K are all even.  But odd values are
    // allowed if the definition of + for matrices is extended to 
    // apply to different sized matrices, where the extra row or 
    // column is just copied to the destination matrix.
    //
    // Note: There are some typos in both their listing in Algorithm 1
    // and the version found in the appendix, although the 
    // non-psuedo-code listing in the appendix is accurate.
    // The listing below has (I think) corrected all the typos.
    //
    // Notation: (+=) means = if add=false, or += if add=true
    //
    // if (min(m,n,k) < TMV_Q4) { Call direct algorithm }
    // else {
    //   X = A2 + A3
    //   Y = B1 - B0
    //   Z = X * Y
    //   C3 (+=) Z
    //   C1 (+=) Z
    // 
    //   Z = A0 * B0
    //   C0 (+=) Z
    //
    //   C0 += A1 * B2
    //   X = X - A0
    //   Y = B3 - Y
    //   Z += X * Y
    //   C1 += Z
    //   
    //   X = A1 - X
    //   C1 += X * B3
    // 
    //   Y = B2 - Y
    //   C2 (+=) A3 * Y
    //
    //   X = A0 - A2
    //   Y = B3 - B1
    //   Z += X * Y
    //   C3 += Z
    //   C2 += Z
    // }
    //
    // X = A0-A2
    // Y = B3-B1
    // Z = A0*B0 + (A2+A3-A0)*(B0+B3-B1) + (A0-A2)*(B3-B1)
    // If you track through the calculation, the values of C are:
    // C0 (+=) A0*B0 + A1*B2
    // C1 (+=) (A2+A3)*(B1-B0) + A0*B0 + (A2+A3-A0)*(B0+B3-B1)
    //            + (A0+A1-A2-A3)*B3
    // C2 (+=) A3*(B1+B2-B0-B3) + A0*B0 + (A2+A3-A0)*(B0+B3-B1) 
    //            + (A0-A2)*(B3-B1)
    // C3 (+=) (A2+A3)*(B1-B0) + A0*B0 + (A2+A3-A0)*(B0+B3-B1) 
    //            + (A0-A2)*(B3-B1)
    //
    // It is straightforward to check that these equations simplify
    // to the right answers:
    // C0 (+=) A0*B0 + A1*B2
    // C1 (+=) (A2+A3)*(B1-B0) + A0*B0 + (A2+A3-A0)*(B0+B3-B1)
    //            + (A0+A1-A2-A3)*B3
    //    = A2*B1 + A3*B1 - A2*B0 - A3*B0 + A0*B0 + A2*B0 + A2*B3 - A2*B1
    //      + A3*B0 + A3*B3 - A3*B1 - A0*B0 - A0*B3 + A0*B1
    //      + A0*B3 + A1*B3 - A2*B3 - A3*B3
    //    = A0*B1 + A1*B3
    // C2 (+=) A3*(B1+B2-B0-B3) + A0*B0 + (A2+A3-A0)*(B0+B3-B1) 
    //            + (A0-A2)*(B3-B1)
    //    = A3*B1 + A3*B2 - A3*B0 - A3*B3 + A0*B0 + A2*B0 + A2*B3 - A2*B1
    //      + A3*B0 + A3*B3 - A3*B1 - A0*B0 - A0*B3 + A0*B1
    //      + A0*B3 - A0*B1 - A2*B3 + A2*B1
    //    = A2*B0 + A3*B2
    // C3 (+=) (A2+A3)*(B1-B0) + A0*B0 + (A2+A3-A0)*(B0+B3-B1) 
    //            + (A0-A2)*(B3-B1)
    //    = A2*B1 + A3*B1 - A2*B0 - A3*B0 + A0*B0 + A2*B0 + A2*B3 - A2*B1
    //      + A3*B0 + A3*B3 - A3*B1 - A0*B0 - A0*B3 + A0*B1
    //      + A0*B3 - A0*B1 - A2*B3 + A2*B1
    //    = A2*B1 + A3*B3
    //

    const int M = m3.colsize();
    const int N = m3.rowsize();
    const int K = m1.rowsize();
    TMVAssert(m1.colsize() == M);
    TMVAssert(m2.rowsize() == N);
    TMVAssert(m2.colsize() == K);

    const int XX = UNKNOWN;
    typedef ConstMatrixView<T1,XX,XX,CC1> M1x;
    typedef ConstMatrixView<T2,XX,XX,CC2> M2x;
    typedef MatrixView<T3> M3x;
    typedef typename Traits2<T1,T2>::type PT3;
    PT3 one(1);

    const bool recurse = (M >= TMV_Q4 && N >= TMV_Q4 && K >= TMV_Q4);

    const int Mb = M>>1;
    const int Ma = M-Mb;
    const int Nb = N>>1;
    const int Na = N-Nb;
    const int Kb = K>>1;
    const int Ka = K-Kb;

#ifdef PRINTALGO_MM_WIN
    std::cout<<"Winograd recursive branch: \n";
    std::cout<<"M = "<<Ma<<" + "<<Mb<<std::endl;
    std::cout<<"N = "<<Na<<" + "<<Nb<<std::endl;
    std::cout<<"K = "<<Ka<<" + "<<Kb<<std::endl;
#endif

    // X,Y,Z are temporaries.
    // We use the subscripts 1,2,3 to match the size of the 
    // corrsponding A, B or C submatrix.  0 matches the full size.
    typedef Matrix<T1,RowMajor> X_type;
    typedef Matrix<T2,ColMajor> Y_type;
    typedef Matrix<PT3,ColMajor> Z_type;
    X_type X(Ma,Ka,T1(0)); 
    Y_type Y(Ka,Na,T2(0));
    Z_type Z(Ma,Na,PT3(0));

    typename X_type::view_type      X0 = X.View();
    typename X_type::cols_type      X1 = X.Cols(0,Kb);
    typename X_type::rows_type      X2 = X.Rows(0,Mb);
    typename X_type::submatrix_type X3 = X.SubMatrix(0,Mb,0,Kb);

    typename Y_type::view_type      Y0 = Y.View();
    typename Y_type::cols_type      Y1 = Y.Cols(0,Nb);
    typename Y_type::rows_type      Y2 = Y.Rows(0,Kb);
    typename Y_type::submatrix_type Y3 = Y.SubMatrix(0,Kb,0,Nb);

    typename Z_type::view_type      Z0 = Z.View();
    typename Z_type::cols_type      Z1 = Z.Cols(0,Nb);
    typename Z_type::rows_type      Z2 = Z.Rows(0,Mb);
    typename Z_type::submatrix_type Z3 = Z.SubMatrix(0,Mb,0,Nb);

    M1x A0 = m1.SubMatrix(0,Ma,0,Ka);
    M1x A1 = m1.SubMatrix(0,Ma,Ka,K);
    M1x A2 = m1.SubMatrix(Ma,M,0,Ka);
    M1x A3 = m1.SubMatrix(Ma,M,Ka,K);
    M2x B0 = m2.SubMatrix(0,Ka,0,Na);
    M2x B1 = m2.SubMatrix(0,Ka,Na,N);
    M2x B2 = m2.SubMatrix(Ka,K,0,Na);
    M2x B3 = m2.SubMatrix(Ka,K,Na,N);
    M3x C0 = m3.SubMatrix(0,Ma,0,Na);
    M3x C1 = m3.SubMatrix(0,Ma,Na,N);
    M3x C2 = m3.SubMatrix(Ma,M,0,Na);
    M3x C3 = m3.SubMatrix(Ma,M,Na,N);

    // X = A2 + A3
    X2 = A2;
    X3 += A3;

    // Y = B1 - B0
    Y1 = B1;
    Y -= B0;

    // Z = X * Y
    if (recurse) {
      ConstMatrixView<T1> X2x = X2.XView();
      ConstMatrixView<T2> Y0x = Y0.XView();
      MatrixView<PT3> Z2x = Z2.XView();
      MultMM_Winograd(false,one,X2x,Y0x,Z2x);
    }
    else
      MultMM_Winograd_Final(false,one,X2,Y0,Z2);

    // C3 (+=) Z
    if (add) C3 += x * Z3;
    else C3 = x * Z3;

    // C1 (+=) Z
    if (add) C1 += x * Z1;
    else C1 = x * Z1;

    // Z = A0 * B0
    if (recurse) {
      MatrixView<PT3> Z0x = Z0.XView();
      MultMM_Winograd(false,one,A0,B0,Z0x);
    }
    else
      MultMM_Winograd_Final(false,one,A0,B0,Z0);

    // C0 (+=) Z
    if (add) C0 += x * Z;
    else C0 = x * Z;

    // C0 += A1 * B2
    if (recurse) 
      MultMM_Winograd(true,x,A1,B2,C0);
    else
      MultMM_Winograd_Final(true,x,A1,B2,C0);

    // X = X - A0
    X -= A0;

    // Y = B3 - Y
    Y *= -1;
    Y3 += B3;

    // Z += X * Y
    if (recurse) {
      ConstMatrixView<T1> X0x = X0.XView();
      ConstMatrixView<T2> Y0x = Y0.XView();
      MatrixView<PT3> Z0x = Z0.XView();
      MultMM_Winograd(true,one,X0x,Y0x,Z0x);
    }
    else
      MultMM_Winograd_Final(true,one,X0,Y0,Z0);

    // C1 += Z
    C1 += x * Z1;

    // X = A1 - X
    X *= -1;
    X1 += A1;

    // C1 += X * B3
    if (recurse) {
      ConstMatrixView<T1> X1x = X1.XView();
      MultMM_Winograd(true,x,X1x,B3,C1);
    }
    else
      MultMM_Winograd_Final(true,x,X1,B3,C1);

    // Y = B2 - Y
    Y *= -1;
    Y2 += B2;

    // C2 (+=) A3 * Y
    if (recurse) {
      ConstMatrixView<T2> Y2x = Y2.XView();
      MultMM_Winograd(add,x,A3,Y2x,C2);
    }
    else 
      MultMM_Winograd_Final(add,x,A3,Y2,C2);

    // X = A0 - A2
    X = A0;
    X2 -= A2;

    // Y = B3 - B1
    Y1 = -B1;
    Y3 += B3;

    // Z += X * Y
    if (recurse) {
      ConstMatrixView<T1> X0x = X0.XView();
      ConstMatrixView<T2> Y1x = Y1.XView();
      MatrixView<PT3> Z1x = Z1.XView();
      MultMM_Winograd(true,one,X0x,Y1x,Z1x);
    }
    else
      MultMM_Winograd_Final(true,one,X0,Y1,Z1);

    // C3 += Z
    C3 += x * Z3;

    // C2 += Z
    C2 += x * Z2;
  }

  // Generic version for non instantiated types:
  template <class T1, bool CC1, class T2, bool CC2, class T3>
  static void MultMM_Winograd(
      const bool add, const T3& x, 
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<T2,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<T3>& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }


#if !defined(TMV_CompilingWinograd) && defined(TMV_INST_FLOAT)
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const float& x, 
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<float>& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<float>& x, 
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<float> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<float>& x, 
      const ConstMatrixView<std::complex<float>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<float> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<float>& x, 
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<float>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<float> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<float>& x, 
      const ConstMatrixView<std::complex<float>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<float>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<float> >& m3);
#elif \
  (defined(TMV_CompilingWinograd) && defined(TMV_INST_FLOAT)) || \
  (!defined(TMV_CompilingWinograd) && !defined(TMV_INST_FLOAT)) 
  // comp   inst_float
  // !comp  !inst_float
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const float& x, 
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<float>& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<float>& x, 
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<float> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<float>& x, 
      const ConstMatrixView<std::complex<float>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<float> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<float>& x, 
      const ConstMatrixView<float,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<float>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<float> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<float>& x, 
      const ConstMatrixView<std::complex<float>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<float>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<float> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
#endif

#if !defined(TMV_CompilingWinograd) && defined(TMV_INST_DOUBLE)
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const double& x, 
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<double>& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<double>& x, 
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<double> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<double>& x, 
      const ConstMatrixView<std::complex<double>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<double> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<double>& x, 
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<double>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<double> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<double>& x, 
      const ConstMatrixView<std::complex<double>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<double>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<double> >& m3);
#elif \
  (defined(TMV_CompilingWinograd) && defined(TMV_INST_DOUBLE)) || \
  (!defined(TMV_CompilingWinograd) && !defined(TMV_INST_DOUBLE)) 
  // comp   inst_double
  // !comp  !inst_double
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const double& x, 
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<double>& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<double>& x, 
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<double> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<double>& x, 
      const ConstMatrixView<std::complex<double>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<double> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<double>& x, 
      const ConstMatrixView<double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<double>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<double> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<double>& x, 
      const ConstMatrixView<std::complex<double>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<double>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<double> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
#endif
#if !defined(TMV_CompilingWinograd) && defined(TMV_INST_INT)
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const int& x, 
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<int>& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<int>& x, 
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<int> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<int>& x, 
      const ConstMatrixView<std::complex<int>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<int> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<int>& x, 
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<int>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<int> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<int>& x, 
      const ConstMatrixView<std::complex<int>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<int>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<int> >& m3);
#elif \
  (defined(TMV_CompilingWinograd) && defined(TMV_INST_INT)) || \
  (!defined(TMV_CompilingWinograd) && !defined(TMV_INST_INT)) 
  // comp   inst_int
  // !comp  !inst_int
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const int& x, 
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<int>& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<int>& x, 
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<int> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<int>& x, 
      const ConstMatrixView<std::complex<int>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<int> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<int>& x, 
      const ConstMatrixView<int,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<int>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<int> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<int>& x, 
      const ConstMatrixView<std::complex<int>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<int>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<int> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
#endif
#if !defined(TMV_CompilingWinograd) && defined(TMV_INST_LONGDOUBLE)
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const long double& x, 
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<long double>& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<long double>& x, 
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<long double> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<long double>& x, 
      const ConstMatrixView<std::complex<long double>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<long double> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<long double>& x, 
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<long double>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<long double> >& m3);
  template <bool CC1, bool CC2>
  void MultMM_Winograd(
      const bool add, const std::complex<long double>& x, 
      const ConstMatrixView<std::complex<long double>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<long double>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<long double> >& m3);
#elif \
  (defined(TMV_CompilingWinograd) && defined(TMV_INST_LONGDOUBLE)) || \
  (!defined(TMV_CompilingWinograd) && !defined(TMV_INST_LONGDOUBLE)) 
  // comp   inst_long double
  // !comp  !inst_long double
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const long double& x, 
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<long double>& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<long double>& x, 
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<long double> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<long double>& x, 
      const ConstMatrixView<std::complex<long double>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<long double> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<long double>& x, 
      const ConstMatrixView<long double,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<long double>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<long double> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
  template <bool CC1, bool CC2>
  TMV_STATIC void MultMM_Winograd(
      const bool add, const std::complex<long double>& x, 
      const ConstMatrixView<std::complex<long double>,UNKNOWN,UNKNOWN,CC1>& m1,
      const ConstMatrixView<std::complex<long double>,UNKNOWN,UNKNOWN,CC2>& m2,
      MatrixView<std::complex<long double> >& m3)
  { DoMultMM_Winograd(add,x,m1,m2,m3); }
#endif

#undef TMV_STATIC

} // namespace tmv

#endif 
