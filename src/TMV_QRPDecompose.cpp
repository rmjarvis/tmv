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
#ifdef NOGEQP3
#ifdef LAP
#undef LAP
#endif
#endif

#include "TMV_Blas.h"
#include "TMV_QRPDiv.h"
#include "tmv/TMV_QRPD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "TMV_QRDiv.h"
#include "TMV_Householder.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_PackedQ.h"
#ifdef LAP
#include "TMV_ConvertIndex.h"
#endif

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_VIt.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv 
{

#ifdef TMV_BLOCKSIZE
#define QRP_BLOCKSIZE TMV_BLOCKSIZE
#else
#define QRP_BLOCKSIZE 32
#endif

#define RT RealType(T)

  //
  // QRP Decompose
  //

  template <class T> static void NonBlockQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, int* P, T& det, bool strict)
  {
#ifdef XDEBUG
    Matrix<T> A0(A);
    cout<<"Start NonBlock QRPD with strict = "<<strict<<endl;
    cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    cout<<"beta = "<<Type(beta)<<"  "<<beta<<endl;
#endif
    // Decompose A into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (output of A), 
    // with the beta's for the Householder matrices returned in beta.

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    const int M = A.colsize();
    const int N = A.rowsize();
    const int Astepj = A.stepj();

    // Keep track of the norm of each column
    // When considering column j, these are actually just the norm
    // of each column from j:M, not 0:M.
    RT scale = RT(1) / A.MaxAbsElement(); // for more stable NormSq
    Vector<RT> colnormsq(N);
    for(int j=0;j<N;++j) colnormsq(j) = A.col(j).NormSq(scale);
    //cout<<"colnormsq's are "<<colnormsq<<endl;
    RT anormsq = colnormsq.SumElements();
    RT thresh = RT(N) * SQR(Epsilon<T>()) * anormsq;
    //cout<<"thresh = "<<thresh<<endl;
    // Set to 0 any diag element whose norm is < epsilon * |A|
    RT recalcthresh(0);
    // recalcthresh is the threshold for recalculating the norm to account for
    // rounding errors in the subtractions which keep track of it..
    // The is set to sqrt(Epsilon) * the largest normsq whenever we 
    // recalculate the norms. 

    T* bj = beta.ptr();
    for(int j=0;j<N;++j,++bj) {
      //cout<<"j = "<<j<<" colnormsq = "<<colnormsq(j)<<endl;

      if (strict || j==0 || colnormsq(j) < recalcthresh) {
        // Find the column with the largest norm
        int jpiv;
        RT maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
        //cout<<"maxnormsq = "<<maxnormsq<<" at "<<jpiv+j<<std::endl;
        if (j==0) recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;
        //cout<<"recalcthresh = "<<recalcthresh<<std::endl;
        // Note: jpiv is relative to the SubVector(j,N)

        // If the largest colnormsq is lower than the recalulation threshold,
        // then recalc all colnormsq's, and redetermine max.
        if (maxnormsq < recalcthresh) {
          for(int k=j;k<N;++k) colnormsq(k) = A.col(k,j,M).NormSq(scale);
          //cout<<"colnormsq -> "<<colnormsq<<std::endl;
          maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
          //cout<<"maxnormsq = "<<maxnormsq<<" at "<<jpiv+j<<std::endl;
          recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;
          if (recalcthresh < thresh) recalcthresh = thresh;
          //cout<<"recalcthresh = "<<recalcthresh<<std::endl;
        }

        // If maxnormsq = 0 (technically < thresh to account for rounding)
        // then the rest of the R matrix is 0, and the Householder matrices 
        // are identities (indicated by 0's in the Q part of the matrix).
        if (maxnormsq < thresh) {
          //cout<<"Do Zero"<<std::endl;
          A.SubMatrix(j,M,j,N).Zero();
          // Already essentially zero - make it exact
          beta.SubVector(j,N).Zero();
          // Set the Householder matrices for these to identities
          for(;j<N;j++) P[j] = j;
          break;
        } else {
          // Swap the column with the largest norm into the current column
          if (jpiv != 0) {
            //cout<<"Do Swap"<<std::endl;
            // Add j to get real index
            jpiv += j;
            TMVAssert(jpiv < int(A.rowsize()));
            colnormsq.Swap(j,jpiv);
            A.SwapCols(j,jpiv);
            det = -det;
            P[j] = jpiv;
          } else {
            P[j] = j;
          }
        }
      } else P[j] = j;

      // Apply the Householder Reflection for this column
#ifdef TMVFLDEBUG
      TMVAssert(bj >= beta.first);
      TMVAssert(bj < beta.last);
#endif

      //cout<<"Apply Householder reflection"<<std::endl;
      *bj = Householder_Reflect(A.SubMatrix(j,M,j,N),det);

      // And update the norms for use with the next column
      const T* Ajk = A.row(j,j+1,N).cptr();
      for(int k=j+1;k<N;++k,Ajk+=Astepj) {
        colnormsq(k) -= NORM(*Ajk*scale);
      }
    }
#ifdef XDEBUG
    cout<<"Done NonBlock QRPDecompose"<<std::endl;
    cout<<"A -> "<<A<<std::endl;
    Matrix<T> Q(A);
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*UpperTriMatrixViewOf(A);
    AA.ReversePermuteCols(P);
    cout<<"Norm(AA-A0) = "<<Norm(AA-A0)<<std::endl;
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"StrictBlockQRP_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"P = ";
      for(int i=0;i<int(A.rowsize());i++) cerr<<P[i]<<" ";
      cerr<<endl;
      cerr<<"QRP = "<<AA<<endl;
      abort(); 
    }
#endif
  }

  template <class T> static void StrictBlockQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, int* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in beta.

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

#ifdef XDEBUG
    Matrix<T> A0(A);
#endif

    const int M = A.colsize();
    const int N = A.rowsize();
    const int Astepj = A.stepj();

    RT scale = RT(1) / A.MaxAbsElement(); // for more stable NormSq
    Vector<RT> colnormsq(N);
    for(int j=0;j<N;++j) colnormsq(j) = A.col(j).NormSq(scale);
    RT anormsq = colnormsq.SumElements();
    RT thresh = RT(N) * SQR(Epsilon<T>()) * anormsq;
    RT recalcthresh(0);

    Matrix<T,RowMajor> ZYtA(MIN(QRP_BLOCKSIZE,M),N);
    // We keep track of the product ZYtA(j1:M,j1:N) [ stored as ZYtA ]
    // since this is the product that we need.  We update this one 
    // row at a time.

    //Matrix<T> A2 = A;

    T* bj = beta.ptr();
    for(int j1=0;j1<N;) {
      int j2 = MIN(N,j1+QRP_BLOCKSIZE);
      //cout<<"Start block "<<j1<<","<<j2<<endl;
      //cout<<"P = ";
      //for(int i=0;i<int(A.rowsize());i++) cout<<P[i]<<" ";
      //cout<<endl;
      //cout<<"thresh = "<<thresh<<", recalcthresh = "<<recalcthresh<<std::endl;

      for(int j=j1,jmj1=0; j<j2; ) {
        //cout<<"j = "<<j<<endl;
        //cout<<"beta = "<<beta<<endl;
        //cout<<"A.diag() = "<<A.diag()<<std::endl;
        //cout<<"colnormsq = "<<colnormsq<<std::endl;
        int jpiv;
        RT maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
        //cout<<"jpiv = "<<jpiv<<endl;
        if (recalcthresh == RT(0)) {
          recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;
          //cout<<"recalcthresh -> "<<recalcthresh<<std::endl;
        }

        if (maxnormsq < recalcthresh) {
          // If need to recalculate norms, then do so and end the block here.
          //cout<<"need to recalc"<<std::endl;
          // The columns are only correct if j == j1.  Otherwise,
          // indicate to end the block, which will update the columns.
          // Then next time through, we can do the recalc as needed.
          if (j==j1) {
            for(int k=j;k<N;++k) colnormsq(k) = A.col(k,j,M).NormSq(scale);
            recalcthresh = RT(0);
          }
          j2 = j;
        } else if (maxnormsq < thresh) {
          //cout<<"everything else is below thresh"<<std::endl;
          if (j==j1) {
            // If first in set, then just zero the rest out and 
            // indicate that we are done.
            A.SubMatrix(j,M,j,N).Zero();
            beta.SubVector(j,N).Zero();
            for(;j<N;j++) P[j] = j;
            //cout<<"P("<<j<<":"<<N<<") = j"<<std::endl;
            //cout<<"beta("<<j<<":"<<N<<") = 0"<<std::endl;
            j2 = N; 
          } else {
            // Otherwise, do the block Householder transforms on the 
            // previous columns first.  (The next time through the block
            // loop should still result in maxnormsq < thresh.)
            j2 = j; 
          }
        } else {
          // Normal case with columns still left to do.

          // Pivot
          if (jpiv != 0) {
            jpiv += j;
            TMVAssert(jpiv < int(A.rowsize()));
            colnormsq.Swap(j,jpiv);
            A.SwapCols(j,jpiv);
            //A2.SwapCols(j,jpiv);
            ZYtA.Rows(0,jmj1).SwapCols(j,jpiv);
            det = -det;
            P[j] = jpiv;
          } else {
            P[j] = j;
          }
          //cout<<"P["<<j<<"] = "<<P[j]<<std::endl;

          // Update the pivot column with Block Householder so far:
          // A(j1:M,j) -= Y Z Yt A(j1:M,j)
          // A(j1:j,j) has already been updated, so we only need to do A(j:M,j)
          // A(j:M,j) -= Y(j:M,0:j) (ZYtA)(0:j,j)
          A.col(j,j,M) -= A.SubMatrix(j,M,j1,j) * ZYtA.col(j,0,jmj1);

          // Find Householder matrix for this column
#ifdef TMVFLDEBUG
          TMVAssert(bj >= beta.first);
          TMVAssert(bj < beta.last);
#endif
          //T d2(0);
          //*bj = Householder_Reflect(A2.SubMatrix(j,M,j,N),d2);
          //cout<<"nonblock bj = "<<*bj<<std::endl;
          //cout<<"v = "<<A2.col(j,j,M)<<endl;
          *bj = Householder_Reflect(A.col(j,j,M),det);
          //cout<<"beta["<<j<<"] = "<<*bj<<std::endl;
          //cout<<"v = "<<A.col(j,j,M)<<endl;

          //cout<<"Before update:"<<std::endl;
          //for(int k=j+1;k<N;++k) {
          //cout<<"colnormsq("<<k<<") = "<<colnormsq[k]<<endl;
          //cout<<"Actual colnormsq = "<<A2.col(k,j,M).NormSq(scale)<<endl;
          //}

          // Update ZYtA:
          if (*bj != T(0)) {
            // (I-beta v vt)(I-Y Z Yt) A
            // = I - Y (ZYtA) - v (beta vt A) + v (beta vt Y (ZYtA))
            // The augmented Y now includes v in the j column, 
            // so the augmented ZYtA now has to include in the j row:
            // beta (vt A - vt Y ZYtA)
            VectorView<T> vt = A.col(j,j+1,M).Conjugate();
            // Remember, this doesn't include the implicit 1 at the top of v.
            ZYtA.row(jmj1,j1,j+1).Zero();
            ZYtA.row(jmj1,j+1,N) = vt * A.SubMatrix(j+1,M,j+1,N);
            ZYtA.row(jmj1,j+1,N) += A.row(j,j+1,N);
            Vector<T> vtY = vt * A.SubMatrix(j+1,M,j1,j);
            vtY += A.row(j,j1,j);
            ZYtA.row(jmj1,j1,N) -= vtY * ZYtA.SubMatrix(0,jmj1,j1,N);
            ZYtA.row(jmj1,j1,N) *= *bj;
          } else ZYtA.row(jmj1,j1,N).Zero();

          // Update row j of the rest of the matrix:
          // A(j,j+1:N) -= (Y ZYtA)(j,j+1:N) = Y(j,j1:j+1) ZYtA(j1:j+1,j+1:N)
          VectorView<T> Arowj = A.row(j,j+1,N);
          Arowj -= A.row(j,j1,j)*ZYtA.SubMatrix(0,jmj1,j+1,N);
          Arowj -= ZYtA.row(jmj1,j+1,N);

          //cout<<"A.row("<<j<<") = "<<A.row(j,j+1,N)<<endl;
          //cout<<"A2.row("<<j<<") = "<<A2.row(j,j+1,N)<<endl;

          //cout<<"After Householder reflect"<<std::endl;
          //for(int k=j+1;k<N;++k) {
          //cout<<"Actual colnormsq("<<k<<") = "<<A2.col(k,j,M).NormSq(scale)<<endl;
          //cout<<"New colnormsq("<<k<<") = "<<A2.col(k,j+1,M).NormSq(scale)<<endl;
          //cout<<"NORM(A2("<<j<<","<<k<<")) = "<<NORM(A2(j,k))<<std::endl;
          //cout<<"NORM(A("<<j<<","<<k<<")) = "<<NORM(A(j,k))<<std::endl;
          //}

          // Update the colnormsq values
          const T* Ajk = Arowj.cptr();
          for(int k=j+1;k<N;++k,Ajk+=Astepj) {
            //cout<<"colnormsq("<<k<<") -= "<<NORM(*Ajk);
            colnormsq(k) -= NORM(*Ajk*scale);
            //cout<<" = "<<colnormsq(k)<<endl;
          }
          ++j; ++jmj1; ++bj;
        }
      }
      //cout<<"after j loop "<<std::endl;
      // Do the Block Householder update of the rest of the matrix:
      // A(j2:M,j2:N) -= Y(j2:M,j1:j2) ZYtA(j1:j2,j1:N)
      A.SubMatrix(j2,M,j2,N) -= A.SubMatrix(j2,M,j1,j2) * 
      ZYtA.SubMatrix(0,j2-j1,j2,N);
      j1 = j2;
      //cout<<"after A update "<<std::endl;
    }
    //cout<<"done A.diag = "<<A.diag()<<endl;

#ifdef XDEBUG
    Matrix<T> Q(A);
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*A.UpperTri();
    AA.ReversePermuteCols(P);
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"StrictBlockQRP_Decompose: A = "<<TypeText(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"P = ";
      for(int i=0;i<int(A.rowsize());i++) cerr<<P[i]<<" ";
      cerr<<endl;
      cerr<<"QRP = "<<AA<<endl;
      abort(); 
    }
#endif
  }

  template <class T> static void MoveLowColsToEnd(
      Vector<RT>& colnormsq, RT thresh,
      int j1, int& j2, int& j3, const MatrixView<T>& A, int* P)
  {
    // Move all columns of A (whose norms are in colnormsq) with norms
    // less than thresh to the end.  j1 is the first column we need to 
    // look at.  j2 is one past the last column we need to look at.
    // j3 is one past the last good column overall.
    // On output j2 and j3 are updated lower if necessary.

    --j3; // temporarily, j3 is the last good column.
    while (j3 > j1 && colnormsq(j3) < thresh) --j3;
    if (j3==j1) {
      if (colnormsq(j1) < thresh) j2 = j1;
      else { j2 = j3 = j1+1; P[j1] = j1; }
      return;
    }
    if (j3 < j2) j2 = j3+1;
    for(int i=j1;i<j2;++i) {
      if (colnormsq(i) < thresh) {
        //cerr<<"  "<<i<<" is low ("<<colnormsq(i)<<") - swap with "<<j3<<endl;
        TMVAssert(j3 < int(A.rowsize()));
        A.SwapCols(i,j3);
        colnormsq.Swap(i,j3);
        P[i] = j3;
        while (colnormsq(--j3) < thresh);
        if (j3 < j2) j2 = j3+1;
      } else P[i] = i;
    }
    ++j3; // j3 back to being first bad column.
  }

#ifdef XDEBUG
  static void CheckIndex(const Vector<double>& index, const int* P, int j1)
  {
    const int N = index.size();
    Vector<double> index2(N);
    for(int k=0;k<N;k++) index2(k) = double(k);
    for(int k=0;k<j1;k++) index2.Swap(k,P[k]);
    if (Norm(index-index2) > 0.01) {
      cerr<<"index = "<<index<<endl;
      cerr<<"index2 = "<<index2<<endl;
      cerr<<"norm(diff) = "<<Norm(index-index2)<<endl;
      abort();
    }
  }
#endif

  template <class T> static void LooseBlockQRP_Decompose(
      const MatrixView<T>& A, const VectorView<T>& beta, int* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in beta.
    //
    // This loose version doesn't sort the diagonal of R exactly.
    // It only sorts them enough to make sure the 0's fall at the end.

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

#ifdef XDEBUG
    Matrix<T> A0(A);
#endif

    const int M = A.colsize();
    const int N = A.rowsize();
    const int Astepj = A.stepj();

    RT scale = RT(1) / A.MaxAbsElement(); // for more stable NormSq
    Vector<RT> colnormsq(N);
    for(int j=0;j<N;++j) colnormsq(j) = A.col(j).NormSq(scale);
    //cerr<<"colnormsq = "<<colnormsq<<endl;
    RT anormsq = colnormsq.SumElements();
    //cerr<<"anormsq = "<<anormsq<<endl;
    //cerr<<"eps = "<<Epsilon<T>()<<endl;
    RT thresh = RT(N) * SQR(Epsilon<T>()) * anormsq;
    //cerr<<"thresh = "<<thresh<<endl;

#ifdef XDEBUG
    Vector<double> index(N);
    for(int k=0;k<N;k++) index(k) = double(k);
#endif

    for (int j1 = 0; j1 < N;) {
      //cerr<<"Start outer loop: j1 = "<<j1<<endl;
      // Do as many columns as possible such that none have to have 
      // their norms recalculated.
      int j3=N; // j3 will be how many we have done in this loop
      // Invariant: all columns from j3..N are known to have norms that
      // need to be recalculated.
      // The recalculation is done at the end of the loop.

      int jpiv0;
      RT maxnormsq = colnormsq.SubVector(j1,N).MaxElement(&jpiv0);
      //cerr<<"maxnormsq = "<<maxnormsq<<endl;

      if (maxnormsq < thresh) {
        //cerr<<"Do Zero"<<std::endl;
        // Zero the rest out and we are done.
        A.SubMatrix(j1,M,j1,N).Zero();
        beta.SubVector(j1,N).Zero();
        for(;j1<N;j1++) P[j1] = j1;
        break;
      } 

      // Move max column to the front:
      if (jpiv0 != 0) {
        jpiv0 += j1;
        //cerr<<"swap j1 with jpiv0 = "<<jpiv0<<endl;
        TMVAssert(jpiv0 < int(A.rowsize()));
        A.SwapCols(j1,jpiv0);
        colnormsq.Swap(j1,jpiv0);
        P[j1] = jpiv0;
#ifdef XDEBUG
        index.Swap(j1,jpiv0);
#endif
      } else P[j1] = j1;
#ifdef XDEBUG
      CheckIndex(index,P,j1+1);
#endif

      RT recalcthresh = RT(N)*SqrtEpsilon<T>()*maxnormsq;
      if (recalcthresh < thresh) recalcthresh = thresh;
      //cerr<<"recalcthresh = "<<recalcthresh<<endl;

      TMVAssert(j1<j3);
      int j1x = j1+1; 
      // The first pass through, we don't want to include j1 in the 
      // MoveLowColsToEnd call.

      // Work on this one block at a time:
      while (j1 < j3) {
        //cerr<<"Start inner loop: j1 = "<<j1<<", j3 = "<<j3<<endl;
        int j2 = MIN(j3,j1+QRP_BLOCKSIZE);
        //cerr<<"j2 = "<<j2<<endl;
        TMVAssert(j1 - j2 < 0);
        MoveLowColsToEnd(colnormsq,recalcthresh,j1x,j2,j3,A,P);
        //cerr<<"After MoveLowColsToEnd: "<<std::endl;
        //cerr<<"j2, j3 = "<<j2<<','<<j3<<endl;
#ifdef XDEBUG
        for(int k=j1x;k<j2;k++) index.Swap(k,P[k]);
        CheckIndex(index,P,j2);
#endif

        int origj2 = j2;
        UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(j2-j1);
        //cerr<<"Block "<<j1<<".."<<j2<<endl;

        //cerr<<"colnormsq = "<<colnormsq.SubVector(j1,j2)<<endl;

        T* bj = beta.ptr()+j1*beta.step();
        for(int j=j1; j<j2; ++j, ++bj) {

          if (colnormsq(j) < recalcthresh) {
#ifdef XDEBUG
            CheckIndex(index,P,origj2);
#endif
            --j2;
            //cerr<<" "<<j<<" is low ("<<colnormsq(j)<<") - swap with "<<j2<<endl;
            //cerr<<" Current Ps are "<<P[j]<<"  "<<P[j2]<<endl;
            if (j==j2) break;
            A.SwapCols(j,j2);
            colnormsq.Swap(j,j2);
#ifdef XDEBUG
            index.Swap(j,j2);
#endif
            if (P[j2] > j2) {
              if (P[j] > j2) {
                TMVAssert(P[j] < int(A.rowsize()));
                TMVAssert(P[j2] < int(A.rowsize()));
                __TMV_SWAP(P[j],P[j2]);
                A.SwapCols(P[j],P[j2]);
                colnormsq.Swap(P[j],P[j2]);
#ifdef XDEBUG
                index.Swap(P[j],P[j2]);
#endif
              } else {
                P[j] = P[j2];
              }
            } else {
              if (P[j] > j2) P[j2] = P[j];
              P[j] = j2;
            }
#ifdef XDEBUG
            CheckIndex(index,P,origj2);
#endif
          }

          // Find Householder matrix for this column
          // This multiplies through to the end of the original block.
          // This way, when we are done, the whole block has had the
          // same Householder reflections applied to it.
#ifdef TMVFLDEBUG
          TMVAssert(bj >= beta.first);
          TMVAssert(bj < beta.last);
#endif
          *bj = Householder_Reflect(A.SubMatrix(j,M,j,origj2),det);

          // Update Z:
          BlockHouseholder_Augment(A.SubMatrix(j1,M,j1,j+1),
              Z.SubTriMatrix(0,j-j1+1),CONJ(*bj));

          // Update the colnormsq values within this block
          // (No need to go all the way to origj2, since the j2..origj2
          // columns are those with low norm already - we don't need
          // those values until we recalculate them from scratch anyway.)
          const T* Ajk = A.row(j,j+1,j2).cptr();
          for(int k=j+1;k<j2;++k,Ajk+=Astepj) 
            colnormsq(k) -= tmv::NORM(*Ajk*scale);
        }
        //cerr<<"done: j2 = "<<j2<<endl;

        if (j1 < j2) {

          // Do the Block Householder update of the rest of the matrix:
          BlockHouseholder_LDiv(A.SubMatrix(j1,M,j1,j2),
              Z.SubTriMatrix(0,j2-j1),A.SubMatrix(j1,M,origj2,N));

          // Update the colnormsq values for the rest of the matrix:
          if (M-j2 > j2-j1)
            for(int k=origj2;k<N;++k) colnormsq(k) -= 
              A.col(k,j1,j2).NormSq(scale);
          else 
            for(int k=origj2;k<N;++k) colnormsq(k) = 
              A.col(k,j2,M).NormSq(scale);
        }

        if (j2 < origj2) {
#ifdef XDEBUG
          CheckIndex(index,P,origj2);
#endif
          // Put the bad columns back where they started before this loop:
          for(int j=j2; j<origj2; ++j) if (P[j] > j2) {
            //cerr<<"Sending "<<j<<" back to "<<P[j]<<endl;
            TMVAssert(P[j] < int(A.rowsize()));
            A.SwapCols(j,P[j]);
            colnormsq.Swap(j,P[j]);
#ifdef XDEBUG
            index.Swap(j,P[j]);
#endif
          }
#ifdef XDEBUG
          CheckIndex(index,P,j2);
#endif
        }

        j1 = j1x = j2;
        //cerr<<"For next loop: j1 = "<<j1<<", j2 = "<<j2<<endl;
      }
      //cerr<<"Done main loop pass"<<std::endl;

      if (j3 < N) {
        // Then need to recalculate some of the colnorms:
        for(int k=j3;k<N;++k) colnormsq(k) = A.col(k,j3,M).NormSq(scale);
        //cerr<<"New colnorms are "<<colnormsq.SubVector(j3,N)<<endl;
      }
    }

    if (det != T(0)) {
      for(int i=0;i<N;++i) if (P[i] != i) det = -det;
    }

#ifdef XDEBUG
    CheckIndex(index,P,N);
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*A.UpperTri();
    AA.ReversePermuteCols(P);
    if (Norm(AA-A0) > 0.001*MAX(RT(1),Norm(A0))) {
      cerr<<"LooseBlockQRP_Decompose: "<<std::endl;
      cerr<<"A = "<<TypeText(A)<<endl;
      if (N < 100) {
        cerr<<"  "<<A0<<endl;
        cerr<<"-> "<<A<<endl;
        cerr<<"beta = "<<beta<<endl;
        cerr<<"P = ";
        for(int i=0;i<N;i++) cerr<<P[i]<<" ";
        cerr<<endl;
        cerr<<"QRP = "<<AA<<endl;
        Matrix<T> diff = AA-A0;
        diff.Clip(0.0001);
        cerr<<"diff = "<<diff<<endl;
      }
      cerr<<"Rdiag = "<<A.diag()<<endl;
      cerr<<"Norm(A-QRP) = "<<Norm(AA-A0)<<endl;
      abort(); 
    }
#endif
  }

#undef RT

  template <class T> static inline void NonLapQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, int* P, T& det, bool strict)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in beta.

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

    if (A.rowsize() > QRP_BLOCKSIZE)
      if (strict)
        StrictBlockQRP_Decompose(A,beta,P,det);
      else
        LooseBlockQRP_Decompose(A,beta,P,det);
    else
      NonBlockQRP_Decompose(A,beta,P,det,strict);
  }

#ifdef LAP
  template <class T> static inline void LapQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, int* P, T& det)
  { NonLapQRP_Decompose(A,beta,P,det,true); }
#ifdef INST_DOUBLE
  template <> void LapQRP_Decompose(const MatrixView<double>& A,
      const VectorView<double>& beta, int* P, double& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
    int lda = A.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
    int lwork = 3*n+1;
    auto_array<double> work(new double[lwork]);
#else
    int lwork = -1;
    auto_array<double> work(new double[1]);
    //std::cout<<"Before lwork query"<<std::endl;
    LAPNAME(dgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
        LAPP(lap_p.get()),LAPP(beta.ptr()) 
        LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
    lwork = int(work[0]);
    work.reset(new double[lwork]);
    //std::cout<<"After lwork query: lwork = "<<lwork<<std::endl;
#endif
#endif
    LAPNAME(dgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
        LAPP(lap_p.get()),LAPP(beta.ptr()) 
        LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
    //std::cout<<"After dgeqp3 query: lwork = "<<lwork<<std::endl;
#ifdef LAPNOWORK
    LAP_Results("dgeqp3");
#else
    LAP_Results(int(work[0]),m,n,lwork,"dgeqp3");
#endif
    double thresh = Epsilon<double>()*A.NormF();
    //std::cout<<"thresh = "<<thresh<<std::endl;
    for(int i=0;i<n;++i) {
      if (ABS(A(i,i)) < thresh) {
        TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
        A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
        beta.SubVector(i,n).Zero();
        //std::cout<<"After clear from "<<i<<" to "<<n-1<<std::endl;
        break;
      }
    }
#ifndef CLAP
    for(int i=0;i<n;++i) --((lap_p.get())[i]);
    //std::cout<<"updated lap_p for Fortran BLAS"<<std::endl;
#endif
    ConvertIndexToPermute(n,lap_p.get(),P);
    //std::cout<<"After convert to permute"<<std::endl;
    const double* bi = beta.cptr();
    for(int i=0;i<n;++i,++bi) {
      if (det) {
        if (*bi != 0.) det = -det;
        if (P[i] != i) det = -det;
      }
    }
    //std::cout<<"Update det = "<<det<<std::endl;
  }
  template <> void LapQRP_Decompose(
      const MatrixView<std::complex<double> >& A,
      const VectorView<std::complex<double> >& beta, int* P,
      std::complex<double>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
    int lda = A.stepj();
#ifndef LAPNOWORK
    auto_array<double> rwork(new double[2*n]);
#ifdef NOWORKQUERY
    int lwork = n+1;
    auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#else
    int lwork = -1;
    auto_array<std::complex<double> > work(new std::complex<double>[1]);
    LAPNAME(zgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
        LAPP(lap_p.get()),LAPP(beta.ptr()) 
        LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
    lwork = int(std::real(work[0]));
    work.reset(new std::complex<double>[lwork]);
#endif
#endif
    LAPNAME(zgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
        LAPP(lap_p.get()),LAPP(beta.ptr()) 
        LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
#ifdef LAPNOWORK
    LAP_Results("zgeqp3");
#else
    LAP_Results(int(std::real(work[0])),m,n,lwork,"zgeqp3");
#endif
    beta.ConjugateSelf();
    double thresh = Epsilon<double>()*A.NormF();
    const std::complex<double>* Aii = A.diag().cptr();
    const int Ads = A.diag().step();
    for(int i=0;i<n;++i,Aii+=Ads) {
      TMVAssert(ABS(IMAG(*Aii)) < Epsilon<double>());
      if (ABS(REAL(*Aii)) < thresh) {
        TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
        A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
        beta.SubVector(i,n).Zero();
        break;
      }
    }
#ifndef CLAP
    for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
    ConvertIndexToPermute(n,lap_p.get(),P);
    const std::complex<double>* bi = beta.cptr();
    for(int i=0;i<n;++i,++bi) {
      if (det!=0.) {
        if (IMAG(*bi) != 0.) 
          det *= -CONJ(*bi* *bi)/norm(*bi);
        else if (REAL(*bi) != 0.)
          det = -det;
        if (P[i] != i) det = -det;
      }
    }
  }
#endif
#ifdef INST_FLOAT
  template <> void LapQRP_Decompose(const MatrixView<float>& A,
      const VectorView<float>& beta, int* P, float& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
    int lda = A.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
    int lwork = 3*n+1;
    auto_array<float> work(new float[lwork]);
#else
    int lwork = -1;
    auto_array<float> work(new float[1]);
    LAPNAME(sgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
        LAPP(lap_p.get()),LAPP(beta.ptr()) 
        LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
    lwork = int(work[0]);
    work.reset(new float[lwork]);
#endif
#endif
    LAPNAME(sgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
        LAPP(lap_p.get()),LAPP(beta.ptr()) 
        LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
    LAP_Results("sgeqp3");
#else
    LAP_Results(int(work[0]),m,n,lwork,"sgeqp3");
#endif
    float thresh = Epsilon<float>()*A.NormF();
    for(int i=0;i<n;++i) {
      if (ABS(A(i,i)) < thresh) {
        TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
        A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
        beta.SubVector(i,n).Zero();
        break;
      }
    }
#ifndef CLAP
    for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
    ConvertIndexToPermute(n,lap_p.get(),P);
    const float* bi = beta.cptr();
    for(int i=0;i<n;++i,++bi) {
      if (det) {
        if (*bi != 0.) det = -det;
        if (P[i] != i) det = -det;
      }
    }
  }
  template <> void LapQRP_Decompose(
      const MatrixView<std::complex<float> >& A,
      const VectorView<std::complex<float> >& beta, int* P,
      std::complex<float>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
    int lda = A.stepj();
#ifndef LAPNOWORK
    auto_array<float> rwork(new float[2*n]);
#ifdef NOWORKQUERY
    int lwork = n+1;
    auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#else
    int lwork = -1;
    auto_array<std::complex<float> > work(new std::complex<float>[1]);
    LAPNAME(cgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
        LAPP(lap_p.get()),LAPP(beta.ptr()) 
        LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
    lwork = int(std::real(work[0]));
    work.reset(new std::complex<float>[lwork]);
#endif
#endif
    LAPNAME(cgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
        LAPP(lap_p.get()),LAPP(beta.ptr()) 
        LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
#ifdef LAPNOWORK
    LAP_Results("cgeqp3");
#else
    LAP_Results(int(std::real(work[0])),m,n,lwork,"cgeqp3");
#endif
    float thresh = Epsilon<float>()*A.NormF();
    const std::complex<float>* Aii = A.diag().cptr();
    const int Ads = A.diag().step();
    for(int i=0;i<n;++i,Aii+=Ads) {
      TMVAssert(ABS(IMAG(*Aii)) < Epsilon<float>());
      if (ABS(REAL(*Aii)) < thresh) {
        TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
        A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
        beta.SubVector(i,n).Zero();
        break;
      }
    }
    beta.ConjugateSelf();
#ifndef CLAP
    for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
    ConvertIndexToPermute(n,lap_p.get(),P);
    const std::complex<float>* bi = beta.cptr();
    for(int i=0;i<n;++i,++bi) {
      if (det!=0.F) {
        if (IMAG(*bi) != 0.) 
          det *= -CONJ(*bi* *bi)/norm(*bi);
        else if (REAL(*bi) != 0.)
          det = -det;
        if (P[i] != i) det = -det;
      }
    }
  }
#endif // FLOAT
#endif // LAP

  template <class T> void QRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, int* P, T& det, bool strict)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

#ifdef XDEBUG
    Matrix<T> A0(A);
    std::cout<<"Start QRP_Decompose:"<<std::endl;
    std::cout<<"A = "<<A0<<std::endl;
    std::cout<<"strict = "<<strict<<std::endl;
#ifdef LAP
    Matrix<T> A2(A);
    Vector<T> beta2(beta);
    auto_array<int> P2(new int[beta.size()]);
    T det2=det;
    NonLapQRP_Decompose(A2.View(),beta2.View(),P2.get(),det2,strict);
    std::cout<<"NonLap QRP = "<<A2<<std::endl;
#endif
#endif

    if (A.rowsize() > 0) {
#ifdef LAP
      if (strict && A.iscm()) {
        LapQRP_Decompose(A,beta,P,det);
      } else {
        NonLapQRP_Decompose(A,beta,P,det,strict);
      }
#else
      NonLapQRP_Decompose(A,beta,P,det,strict);
#endif
    }
#ifdef XDEBUG
    Matrix<T> Q(A);
    std::cout<<"Q = "<<Q<<std::endl;
    GetQFromQR(Q.View(),beta);
    std::cout<<"Q => "<<Q<<std::endl;
    Matrix<T> R(UpperTriMatrixViewOf(A));
    std::cout<<"R = "<<R<<std::endl;
    Matrix<T> AA = Q*R;
    std::cout<<"AA = "<<AA<<std::endl;
    AA.ReversePermuteCols(P);
    std::cout<<"Norm(AA -A0) = "<<Norm(AA-A0)<<std::endl;
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"QRP_Decompose: A = "<<TypeText(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"P = ";
      for(int i=0;i<int(A.rowsize());i++) cerr<<P[i]<<" ";
      cerr<<endl;
      cerr<<"QRP = "<<AA<<endl;
#ifdef LAP
      cerr<<"NonLap: A = "<<A2<<endl;
      cerr<<"beta2 = "<<beta2<<endl;
      cerr<<"det2 = "<<det2<<endl;
#endif
      cerr<<"Norm(AA-A0) = "<<Norm(A0-AA)<<endl;
      abort(); 
    }
#endif
  }

  //
  // QRP Decompose - Unpacked
  //

  template <class T> void QRP_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R, int* P,
      bool strict)
  {
    // Decompose A (input as Q) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation

    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(R.colsize() == Q.rowsize());
    TMVAssert(R.rowsize() == Q.rowsize());

    if (Q.isconj()) {
      if (R.isconj()) {
        QRP_Decompose(Q.Conjugate(),R.Conjugate(),P,strict);
      } else {
        QRP_Decompose(Q.Conjugate(),R,P,strict);
        R.ConjugateSelf();
      }
    } else {
      if (R.isconj()) {
        QRP_Decompose(Q,R.Conjugate(),P,strict);
        R.ConjugateSelf();
      } else {
        Vector<T> beta(Q.rowsize());
        T d(0);
        QRP_Decompose(Q,beta.View(),P,d,strict);
        R = Q.UpperTri();
        GetQFromQR(Q,beta.View());
      }
    }
  }

  template <class T> void QRP_Decompose(const MatrixView<T>& A, bool strict)
  {
    // Decompose A into Q R P, but don't keep Q or P.
    // R is returned as A.UpperTri().

    TMVAssert(A.colsize() >= A.rowsize());

    Vector<T> beta(A.rowsize());
    T d(0);
    auto_array<int> P(new int[A.rowsize()]);
    if (A.isconj())
      QRP_Decompose(A.Conjugate(),beta.View(),P.get(),d,strict);
    else
      QRP_Decompose(A,beta.View(),P.get(),d,strict);
  }


#define InstFile "TMV_QRPDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


