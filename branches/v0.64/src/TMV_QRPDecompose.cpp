///////////////////////////////////////////////////////////////////////////////
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

#define RT TMV_RealType(T)

    //
    // QRP Decompose
    //

    template <class T> 
    static void NonBlockQRPDecompose(
        const MatrixView<T>& A, const VectorView<T>& beta, int* P, T& det,
        bool strict)
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
        RT scale = RT(1) / A.maxAbsElement(); // for more stable normSq
        Vector<RT> colnormsq(N);
        for(int j=0;j<N;++j) colnormsq(j) = A.col(j).normSq(scale);
        RT anormsq = colnormsq.sumElements();
        RT thresh = RT(N) * TMV_SQR(TMV_Epsilon<T>()) * anormsq;
        // Set to 0 any diag element whose norm is < epsilon * |A|
        RT recalcthresh(0);
        // recalcthresh is the threshold for recalculating the norm to account 
        // for rounding errors in the subtractions which keep track of it..
        // The is set to sqrt(TMV_Epsilon) * the largest normsq whenever we 
        // recalculate the norms. 

        T* bj = beta.ptr();
        for(int j=0;j<N;++j,++bj) {
            if (strict || j==0 || colnormsq(j) < recalcthresh) {
                // Find the column with the largest norm
                int jpiv;
                RT maxnormsq = colnormsq.subVector(j,N).maxElement(&jpiv);
                if (j==0) recalcthresh = 4*TMV_SqrtEpsilon<T>() * maxnormsq;
                // Note: jpiv is relative to the subVector(j,N)

                // If the largest colnormsq is lower than the recalulation 
                // threshold, then recalc all colnormsq's, and redetermine max.
                if (maxnormsq < recalcthresh) {
                    for(int k=j;k<N;++k) 
                        colnormsq(k) = A.col(k,j,M).normSq(scale);
                    maxnormsq = colnormsq.subVector(j,N).maxElement(&jpiv);
                    recalcthresh = 4*TMV_SqrtEpsilon<T>() * maxnormsq;
                    if (recalcthresh < thresh) recalcthresh = thresh;
                }

                // If maxnormsq = 0 (technically < thresh to account for
                // rounding) then the rest of the R matrix is 0, and the 
                // Householder matrices are identities (indicated by 0's 
                // in the Q part of the matrix).
                if (maxnormsq < thresh) {
                    A.subMatrix(j,M,j,N).setZero();
                    // Already essentially zero - make it exact
                    beta.subVector(j,N).setZero();
                    // Set the Householder matrices for these to identities
                    for(;j<N;j++) P[j] = j;
                    break;
                } else {
                    // Swap the column with the largest norm into the current 
                    // column
                    if (jpiv != 0) {
                        // Add j to get real index
                        jpiv += j;
                        TMVAssert(jpiv < int(A.rowsize()));
                        colnormsq.swap(j,jpiv);
                        A.swapCols(j,jpiv);
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

            *bj = HouseholderReflect(A.subMatrix(j,M,j,N),det);

            // And update the norms for use with the next column
            const T* Ajk = A.row(j,j+1,N).cptr();
            for(int k=j+1;k<N;++k,Ajk+=Astepj) {
                colnormsq(k) -= TMV_NORM(*Ajk*scale);
            }
        }
#ifdef XDEBUG
        cout<<"Done NonBlock QRPDecompose"<<std::endl;
        cout<<"A -> "<<A<<std::endl;
        Matrix<T> Q(A);
        GetQFromQR(Q.view(),beta);
        Matrix<T> AA = Q*UpperTriMatrixViewOf(A);
        AA.reversePermuteCols(P);
        cout<<"Norm(AA-A0) = "<<Norm(AA-A0)<<std::endl;
        if (Norm(AA-A0) > 0.001*Norm(A0)) {
            cerr<<"StrictBlockQRPDecompose: A = "<<Type(A)<<"  "<<A0<<endl;
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

    template <class T> 
    static void StrictBlockQRPDecompose(
        const MatrixView<T>& A, const VectorView<T>& beta, int* P, T& det)
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

        RT scale = RT(1) / A.maxAbsElement(); // for more stable normSq
        Vector<RT> colnormsq(N);
        for(int j=0;j<N;++j) colnormsq(j) = A.col(j).normSq(scale);
        RT anormsq = colnormsq.sumElements();
        RT thresh = RT(N) * TMV_SQR(TMV_Epsilon<T>()) * anormsq;
        RT recalcthresh(0);

        Matrix<T,RowMajor> ZYtA(TMV_MIN(QRP_BLOCKSIZE,M),N);
        // We keep track of the product ZYtA(j1:M,j1:N) [ stored as ZYtA ]
        // since this is the product that we need.  We update this one 
        // row at a time.

        T* bj = beta.ptr();
        for(int j1=0;j1<N;) {
            int j2 = TMV_MIN(N,j1+QRP_BLOCKSIZE);
            for(int j=j1,jmj1=0; j<j2; ) {
                int jpiv;
                RT maxnormsq = colnormsq.subVector(j,N).maxElement(&jpiv);
                if (recalcthresh == RT(0)) {
                    recalcthresh = 4*TMV_SqrtEpsilon<T>() * maxnormsq;
                }

                if (maxnormsq < recalcthresh) {
                    // If need to recalculate norms, then do so and end the 
                    // block here.
                    // The columns are only correct if j == j1.  Otherwise,
                    // indicate to end the block, which will update the columns.
                    // Then next time through, we can do the recalc as needed.
                    if (j==j1) {
                        for(int k=j;k<N;++k) 
                            colnormsq(k) = A.col(k,j,M).normSq(scale);
                        recalcthresh = RT(0);
                    }
                    j2 = j;
                } else if (maxnormsq < thresh) {
                    if (j==j1) {
                        // If first in set, then just zero the rest out and 
                        // indicate that we are done.
                        A.subMatrix(j,M,j,N).setZero();
                        beta.subVector(j,N).setZero();
                        for(;j<N;j++) P[j] = j;
                        j2 = N; 
                    } else {
                        // Otherwise, do the block Householder transforms 
                        // on the previous columns first.  
                        // (The next time through the block loop should 
                        // still result in maxnormsq < thresh.)
                        j2 = j; 
                    }
                } else {
                    // normal case with columns still left to do.

                    // Pivot
                    if (jpiv != 0) {
                        jpiv += j;
                        TMVAssert(jpiv < int(A.rowsize()));
                        colnormsq.swap(j,jpiv);
                        A.swapCols(j,jpiv);
                        //A2.swapCols(j,jpiv);
                        ZYtA.rowRange(0,jmj1).swapCols(j,jpiv);
                        det = -det;
                        P[j] = jpiv;
                    } else {
                        P[j] = j;
                    }

                    // Update the pivot column with Block Householder so far:
                    // A(j1:M,j) -= Y Z Yt A(j1:M,j)
                    // A(j1:j,j) has already been updated, so we only 
                    // need to do A(j:M,j)
                    // A(j:M,j) -= Y(j:M,0:j) (ZYtA)(0:j,j)
                    A.col(j,j,M) -= A.subMatrix(j,M,j1,j) * ZYtA.col(j,0,jmj1);

                    // Find Householder matrix for this column
#ifdef TMVFLDEBUG
                    TMVAssert(bj >= beta.first);
                    TMVAssert(bj < beta.last);
#endif
                    *bj = HouseholderReflect(A.col(j,j,M),det);

                    // Update ZYtA:
                    if (*bj != T(0)) {
                        // (I-beta v vt)(I-Y Z Yt) A
                        // = I - Y (ZYtA) - v (beta vt A) + v (beta vt Y (ZYtA))
                        // The augmented Y now includes v in the j column, 
                        // so the augmented ZYtA now has to include in the 
                        // j row:
                        // beta (vt A - vt Y ZYtA)
                        VectorView<T> vt = A.col(j,j+1,M).conjugate();
                        // Remember, this doesn't include the implicit 1 
                        // at the top of v.
                        ZYtA.row(jmj1,j1,j+1).setZero();
                        ZYtA.row(jmj1,j+1,N) = vt * A.subMatrix(j+1,M,j+1,N);
                        ZYtA.row(jmj1,j+1,N) += A.row(j,j+1,N);
                        Vector<T> vtY = vt * A.subMatrix(j+1,M,j1,j);
                        vtY += A.row(j,j1,j);
                        ZYtA.row(jmj1,j1,N) -= 
                            vtY * ZYtA.subMatrix(0,jmj1,j1,N);
                        ZYtA.row(jmj1,j1,N) *= *bj;
                    } else ZYtA.row(jmj1,j1,N).setZero();

                    // Update row j of the rest of the matrix:
                    // A(j,j+1:N) -= (Y ZYtA)(j,j+1:N) 
                    //             = Y(j,j1:j+1) ZYtA(j1:j+1,j+1:N)
                    VectorView<T> Arowj = A.row(j,j+1,N);
                    Arowj -= A.row(j,j1,j)*ZYtA.subMatrix(0,jmj1,j+1,N);
                    Arowj -= ZYtA.row(jmj1,j+1,N);

                    // Update the colnormsq values
                    const T* Ajk = Arowj.cptr();
                    for(int k=j+1;k<N;++k,Ajk+=Astepj) {
                        colnormsq(k) -= TMV_NORM(*Ajk*scale);
                    }
                    ++j; ++jmj1; ++bj;
                }
            }
            // Do the Block Householder update of the rest of the matrix:
            // A(j2:M,j2:N) -= Y(j2:M,j1:j2) ZYtA(j1:j2,j1:N)
            A.subMatrix(j2,M,j2,N) -= A.subMatrix(j2,M,j1,j2) * 
                ZYtA.subMatrix(0,j2-j1,j2,N);
            j1 = j2;
        }

#ifdef XDEBUG
        Matrix<T> Q(A);
        GetQFromQR(Q.view(),beta);
        Matrix<T> AA = Q*A.upperTri();
        AA.reversePermuteCols(P);
        if (Norm(AA-A0) > 0.001*Norm(A0)) {
            cerr<<"StrictBlockQRPDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
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

    template <class T> 
    static void moveLowColsToEnd(
        Vector<RT>& colnormsq, RT thresh, int j1, int& j2, int& j3,
        const MatrixView<T>& A, int* P)
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
                TMVAssert(j3 < int(A.rowsize()));
                A.swapCols(i,j3);
                colnormsq.swap(i,j3);
                P[i] = j3;
                while (colnormsq(--j3) < thresh);
                if (j3 < j2) j2 = j3+1;
            } else P[i] = i;
        }
        ++j3; // j3 back to being first bad column.
    }

#ifdef XDEBUG
    static void checkIndex(const Vector<double>& index, const int* P, int j1)
    {
        const int N = index.size();
        Vector<double> index2(N);
        for(int k=0;k<N;k++) index2(k) = double(k);
        for(int k=0;k<j1;k++) index2.swap(k,P[k]);
        if (Norm(index-index2) > 0.01) {
            cerr<<"index = "<<index<<endl;
            cerr<<"index2 = "<<index2<<endl;
            cerr<<"Norm(diff) = "<<Norm(index-index2)<<endl;
            abort();
        }
    }
#endif

    template <class T> 
    static void LooseBlockQRPDecompose(
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

        RT scale = RT(1) / A.maxAbsElement(); // for more stable normSq
        Vector<RT> colnormsq(N);
        for(int j=0;j<N;++j) colnormsq(j) = A.col(j).normSq(scale);
        RT anormsq = colnormsq.sumElements();
        RT thresh = RT(N) * TMV_SQR(TMV_Epsilon<T>()) * anormsq;

#ifdef XDEBUG
        Vector<double> index(N);
        for(int k=0;k<N;k++) index(k) = double(k);
#endif

        for (int j1 = 0; j1 < N;) {
            // Do as many columns as possible such that none have to have 
            // their norms recalculated.
            int j3=N; // j3 will be how many we have done in this loop
            // Invariant: all columns from j3..N are known to have norms that
            // need to be recalculated.
            // The recalculation is done at the end of the loop.

            int jpiv0;
            RT maxnormsq = colnormsq.subVector(j1,N).maxElement(&jpiv0);

            if (maxnormsq < thresh) {
                // Zero the rest out and we are done.
                A.subMatrix(j1,M,j1,N).setZero();
                beta.subVector(j1,N).setZero();
                for(;j1<N;j1++) P[j1] = j1;
                break;
            } 

            // Move max column to the front:
            if (jpiv0 != 0) {
                jpiv0 += j1;
                TMVAssert(jpiv0 < int(A.rowsize()));
                A.swapCols(j1,jpiv0);
                colnormsq.swap(j1,jpiv0);
                P[j1] = jpiv0;
#ifdef XDEBUG
                index.swap(j1,jpiv0);
#endif
            } else P[j1] = j1;
#ifdef XDEBUG
            checkIndex(index,P,j1+1);
#endif

            RT recalcthresh = RT(N)*TMV_SqrtEpsilon<T>()*maxnormsq;
            if (recalcthresh < thresh) recalcthresh = thresh;

            TMVAssert(j1<j3);
            int j1x = j1+1; 
            // The first pass through, we don't want to include j1 in the 
            // moveLowColsToEnd call.

            // Work on this one block at a time:
            while (j1 < j3) {
                int j2 = TMV_MIN(j3,j1+QRP_BLOCKSIZE);
                TMVAssert(j1 - j2 < 0);
                moveLowColsToEnd(colnormsq,recalcthresh,j1x,j2,j3,A,P);
#ifdef XDEBUG
                for(int k=j1x;k<j2;k++) index.swap(k,P[k]);
                checkIndex(index,P,j2);
#endif

                int origj2 = j2;
                UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(j2-j1);


                T* bj = beta.ptr()+j1*beta.step();
                for(int j=j1; j<j2; ++j, ++bj) {

                    if (colnormsq(j) < recalcthresh) {
#ifdef XDEBUG
                        checkIndex(index,P,origj2);
#endif
                        --j2;
                        if (j==j2) break;
                        A.swapCols(j,j2);
                        colnormsq.swap(j,j2);
#ifdef XDEBUG
                        index.swap(j,j2);
#endif
                        if (P[j2] > j2) {
                            if (P[j] > j2) {
                                TMVAssert(P[j] < int(A.rowsize()));
                                TMVAssert(P[j2] < int(A.rowsize()));
                                TMV_SWAP(P[j],P[j2]);
                                A.swapCols(P[j],P[j2]);
                                colnormsq.swap(P[j],P[j2]);
#ifdef XDEBUG
                                index.swap(P[j],P[j2]);
#endif
                            } else {
                                P[j] = P[j2];
                            }
                        } else {
                            if (P[j] > j2) P[j2] = P[j];
                            P[j] = j2;
                        }
#ifdef XDEBUG
                        checkIndex(index,P,origj2);
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
                    *bj = HouseholderReflect(A.subMatrix(j,M,j,origj2),det);

                    // Update Z:
                    BlockHouseholderAugment(
                        A.subMatrix(j1,M,j1,j+1),
                        Z.subTriMatrix(0,j-j1+1),TMV_CONJ(*bj));

                    // Update the colnormsq values within this block
                    // (No need to go all the way to origj2, since the 
                    // j2..origj2 columns are those with low norm already -
                    // we don't need those values until we recalculate them 
                    // from scratch anyway.)
                    const T* Ajk = A.row(j,j+1,j2).cptr();
                    for(int k=j+1;k<j2;++k,Ajk+=Astepj) 
                        colnormsq(k) -= tmv::TMV_NORM(*Ajk*scale);
                }

                if (j1 < j2) {

                    // Do the Block Householder update of the rest of the 
                    // matrix:
                    BlockHouseholderLDiv(
                        A.subMatrix(j1,M,j1,j2),
                        Z.subTriMatrix(0,j2-j1),A.subMatrix(j1,M,origj2,N));

                    // Update the colnormsq values for the rest of the matrix:
                    if (M-j2 > j2-j1)
                        for(int k=origj2;k<N;++k) colnormsq(k) -= 
                            A.col(k,j1,j2).normSq(scale);
                    else 
                        for(int k=origj2;k<N;++k) colnormsq(k) = 
                            A.col(k,j2,M).normSq(scale);
                }

                if (j2 < origj2) {
#ifdef XDEBUG
                    checkIndex(index,P,origj2);
#endif
                    // Put the bad columns back where they started before this 
                    // loop:
                    for(int j=j2; j<origj2; ++j) if (P[j] > j2) {
                        TMVAssert(P[j] < int(A.rowsize()));
                        A.swapCols(j,P[j]);
                        colnormsq.swap(j,P[j]);
#ifdef XDEBUG
                        index.swap(j,P[j]);
#endif
                    }
#ifdef XDEBUG
                    checkIndex(index,P,j2);
#endif
                }

                j1 = j1x = j2;
            }

            if (j3 < N) {
                // Then need to recalculate some of the colnorms:
                for(int k=j3;k<N;++k) 
                    colnormsq(k) = A.col(k,j3,M).normSq(scale);
            }
        }

        if (det != T(0)) {
            for(int i=0;i<N;++i) if (P[i] != i) det = -det;
        }

#ifdef XDEBUG
        checkIndex(index,P,N);
        Matrix<T> Q = A;
        GetQFromQR(Q.view(),beta);
        Matrix<T> AA = Q*A.upperTri();
        AA.reversePermuteCols(P);
        if (Norm(AA-A0) > 0.001*TMV_MAX(RT(1),Norm(A0))) {
            cerr<<"LooseBlockQRPDecompose: "<<std::endl;
            cerr<<"A = "<<TMV_Text(A)<<endl;
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

    template <class T> 
    static inline void NonLapQRPDecompose(
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
                StrictBlockQRPDecompose(A,beta,P,det);
            else
                LooseBlockQRPDecompose(A,beta,P,det);
        else
            NonBlockQRPDecompose(A,beta,P,det,strict);
    }

#ifdef LAP
    template <class T> 
    static inline void LapQRPDecompose(
        const MatrixView<T>& A,
        const VectorView<T>& beta, int* P, T& det)
    { NonLapQRPDecompose(A,beta,P,det,true); }
#ifdef INST_DOUBLE
    template <> 
    void LapQRPDecompose(
        const MatrixView<double>& A,
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
        LAPNAME(dgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new double[lwork]);
#endif
#endif
        LAPNAME(dgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("dgeqp3");
#else
        LAP_Results(int(work[0]),m,n,lwork,"dgeqp3");
#endif
        double thresh = TMV_Epsilon<double>()*A.normF();
        for(int i=0;i<n;++i) {
            if (TMV_ABS(A(i,i)) < thresh) {
                TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
                A.subMatrix(i,A.colsize(),i,A.rowsize()).setZero();
                beta.subVector(i,n).setZero();
                break;
            }
        }
#ifndef CLAP
        for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
        ConvertIndexToPermute(n,lap_p.get(),P);
        const double* bi = beta.cptr();
        for(int i=0;i<n;++i,++bi) {
            if (det) {
                if (*bi != 0.) det = -det;
                if (P[i] != i) det = -det;
            }
        }
    }
    template <> 
    void LapQRPDecompose(
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
        LAPNAME(zgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<double>[lwork]);
#endif
#endif
        LAPNAME(zgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("zgeqp3");
#else
        LAP_Results(int(std::real(work[0])),m,n,lwork,"zgeqp3");
#endif
        beta.conjugateSelf();
        double thresh = TMV_Epsilon<double>()*A.normF();
        const std::complex<double>* Aii = A.diag().cptr();
        const int Ads = A.diag().step();
        for(int i=0;i<n;++i,Aii+=Ads) {
            TMVAssert(TMV_ABS(TMV_IMAG(*Aii)) < TMV_Epsilon<double>());
            if (TMV_ABS(TMV_REAL(*Aii)) < thresh) {
                TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
                A.subMatrix(i,A.colsize(),i,A.rowsize()).setZero();
                beta.subVector(i,n).setZero();
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
                if (TMV_IMAG(*bi) != 0.) 
                    det *= -TMV_CONJ(*bi* *bi)/norm(*bi);
                else if (TMV_REAL(*bi) != 0.)
                    det = -det;
                if (P[i] != i) det = -det;
            }
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapQRPDecompose(
        const MatrixView<float>& A,
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
        LAPNAME(sgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new float[lwork]);
#endif
#endif
        LAPNAME(sgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("sgeqp3");
#else
        LAP_Results(int(work[0]),m,n,lwork,"sgeqp3");
#endif
        float thresh = TMV_Epsilon<float>()*A.normF();
        for(int i=0;i<n;++i) {
            if (TMV_ABS(A(i,i)) < thresh) {
                TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
                A.subMatrix(i,A.colsize(),i,A.rowsize()).setZero();
                beta.subVector(i,n).setZero();
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
    template <> 
    void LapQRPDecompose(
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
        LAPNAME(cgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<float>[lwork]);
#endif
#endif
        LAPNAME(cgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("cgeqp3");
#else
        LAP_Results(int(std::real(work[0])),m,n,lwork,"cgeqp3");
#endif
        float thresh = TMV_Epsilon<float>()*A.normF();
        const std::complex<float>* Aii = A.diag().cptr();
        const int Ads = A.diag().step();
        for(int i=0;i<n;++i,Aii+=Ads) {
            TMVAssert(TMV_ABS(TMV_IMAG(*Aii)) < TMV_Epsilon<float>());
            if (TMV_ABS(TMV_REAL(*Aii)) < thresh) {
                TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
                A.subMatrix(i,A.colsize(),i,A.rowsize()).setZero();
                beta.subVector(i,n).setZero();
                break;
            }
        }
        beta.conjugateSelf();
#ifndef CLAP
        for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
        ConvertIndexToPermute(n,lap_p.get(),P);
        const std::complex<float>* bi = beta.cptr();
        for(int i=0;i<n;++i,++bi) {
            if (det!=0.F) {
                if (TMV_IMAG(*bi) != 0.) 
                    det *= -TMV_CONJ(*bi* *bi)/norm(*bi);
                else if (TMV_REAL(*bi) != 0.)
                    det = -det;
                if (P[i] != i) det = -det;
            }
        }
    }
#endif // FLOAT
#endif // LAP

    template <class T> 
    void QRP_Decompose(
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
        std::cout<<"Start QRPDecompose:"<<std::endl;
        std::cout<<"A = "<<A0<<std::endl;
        std::cout<<"strict = "<<strict<<std::endl;
#ifdef LAP
        Matrix<T> A2(A);
        Vector<T> beta2(beta);
        auto_array<int> P2(new int[beta.size()]);
        T det2=det;
        NonLapQRPDecompose(A2.view(),beta2.view(),P2.get(),det2,strict);
        std::cout<<"NonLap QRP = "<<A2<<std::endl;
#endif
#endif

        if (A.rowsize() > 0) {
#ifdef LAP
            if (strict && A.iscm()) {
                LapQRPDecompose(A,beta,P,det);
            } else {
                NonLapQRPDecompose(A,beta,P,det,strict);
            }
#else
            NonLapQRPDecompose(A,beta,P,det,strict);
#endif
        }
#ifdef XDEBUG
        Matrix<T> Q(A);
        std::cout<<"Q = "<<Q<<std::endl;
        GetQFromQR(Q.view(),beta);
        std::cout<<"Q => "<<Q<<std::endl;
        Matrix<T> R(UpperTriMatrixViewOf(A));
        std::cout<<"R = "<<R<<std::endl;
        Matrix<T> AA = Q*R;
        std::cout<<"AA = "<<AA<<std::endl;
        AA.reversePermuteCols(P);
        std::cout<<"Norm(AA -A0) = "<<Norm(AA-A0)<<std::endl;
        if (Norm(AA-A0) > 0.001*Norm(A0)) {
            cerr<<"QRPDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
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

    template <class T> 
    void QRP_Decompose(
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
                QRP_Decompose(Q.conjugate(),R.conjugate(),P,strict);
            } else {
                QRP_Decompose(Q.conjugate(),R,P,strict);
                R.conjugateSelf();
            }
        } else {
            if (R.isconj()) {
                QRP_Decompose(Q,R.conjugate(),P,strict);
                R.conjugateSelf();
            } else {
                Vector<T> beta(Q.rowsize());
                T d(0);
                QRP_Decompose(Q,beta.view(),P,d,strict);
                R = Q.upperTri();
                GetQFromQR(Q,beta.view());
            }
        }
    }

    template <class T> 
    void QRP_Decompose(const MatrixView<T>& A, bool strict)
    {
        // Decompose A into Q R P, but don't keep Q or P.
        // R is returned as A.upperTri().

        TMVAssert(A.colsize() >= A.rowsize());

        Vector<T> beta(A.rowsize());
        T d(0);
        auto_array<int> P(new int[A.rowsize()]);
        if (A.isconj())
            QRP_Decompose(A.conjugate(),beta.view(),P.get(),d,strict);
        else
            QRP_Decompose(A,beta.view(),P.get(),d,strict);
    }


#define InstFile "TMV_QRPDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


