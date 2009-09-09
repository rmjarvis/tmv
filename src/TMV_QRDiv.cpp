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


#include "TMV_QRDiv.h"
#include "tmv/TMV_QRD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"

namespace tmv {

  //
  // LDiv
  //

  template <class T, class T1, class T2> void QR_LDiv(
      const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const int* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, int N1)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(beta.size() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(x.colsize() == QRx.rowsize());
    TMVAssert(x.rowsize() == m.rowsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    // Solve Q R P x = m
    // where Q and R are stored in QRx, and beta are the beta

    // First Solve Q y = m
    if (QRx.IsSquare()) {
      x = m;
      Q_LDivEq(QRx,beta,x);
    } else {
      if (m.isrm()) {
        Matrix<T,RowMajor> m1 = m;
        // Q is Q1 [ I ]
        //         [ 0 ]
        // where Q1 is the part of Q that is stored in QRx and beta
        // m1 = Q^-1 m1
        Q_LDivEq(QRx,beta,m1.View());
        // y = [ I 0 ] m1
        x = m1.Rows(0,x.colsize()); // x = y here
      } else {
        Matrix<T,ColMajor> m1 = m;
        Q_LDivEq(QRx,beta,m1.View());
        x = m1.Rows(0,x.colsize()); // x = y here
      }
    }

    // Now solve R z = y
    x.Rows(N1,x.colsize()).Zero();
    //x.Rows(0,N1) /= QRx.UpperTri().SubTriMatrix(0,N1);
    QRx.UpperTri().SubTriMatrix(0,N1).LDivEq(x.Rows(0,N1));

    // Finally P x = z
    if (P) x.ReversePermuteRows(P);
  }

  //
  // LDivEq
  //

  template <class T, class T1> void QR_LDivEq(
      const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const int* P,
      const MatrixView<T>& m, int N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(beta.size() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    // Solves Q R P x = m in place (m <- x)
    Q_LDivEq(QRx,beta,m);
    m.Rows(N1,m.colsize()).Zero();
    //m.Rows(0,N1) /= QRx.UpperTri().SubTriMatrix(0,N1);
    QRx.UpperTri().SubTriMatrix(0,N1).LDivEq(m.Rows(0,N1));;
    if (P) m.ReversePermuteRows(P);
  }


  //
  // RDiv
  //

  template <class T, class T1, class T2> void QR_RDiv(
      const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const int* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, int N1)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(beta.size() == QRx.rowsize());
    TMVAssert(x.rowsize() == QRx.colsize());
    TMVAssert(m.rowsize() == QRx.rowsize());
    TMVAssert(x.colsize() == m.colsize());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    // Solve x Q R P = m
    // where Q and R are stored in QRx, and beta are the beta

    // First solve y P = m
    x.Cols(0,m.rowsize()) = m;
    if (P) x.Cols(0,m.rowsize()).PermuteCols(P);

    // Next solve z R = y by forward substitution
    x.Cols(N1,x.rowsize()).Zero();
    //x.Cols(0,N1) %= QRx.UpperTri().SubTriMatrix(0,N1);
    QRx.UpperTri().SubTriMatrix(0,N1).RDivEq(x.Cols(0,N1));;

    // Finally solve x Q = z
    // Q = Q1 [ I ]
    //        [ 0 ]
    // where Q1 is the part of Q that is stored in QRx and beta
    // We've already dealt with the first part by zeroing out the 
    // right columns of x.
    Q_RDivEq(QRx,beta,x);
  }

  //
  // RDivEq
  //

  template <class T, class T1> void QR_RDivEq(
      const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const int* P,
      const MatrixView<T>& m, int N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(beta.size() == QRx.rowsize());
    TMVAssert(m.rowsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    // Solve x Q R P = m in place (m <- x)

    if (P) m.PermuteCols(P);
    m.Cols(N1,m.rowsize()).Zero();
    //m.Cols(0,N1) %= QRx.UpperTri().SubTriMatrix(0,N1);
    QRx.UpperTri().SubTriMatrix(0,N1).RDivEq(m.Cols(0,N1));
    Q_RDivEq(QRx,beta,m);
  }

#define InstFile "TMV_QRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


