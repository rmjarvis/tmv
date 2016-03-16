///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//#define XDEBUG

#ifdef NOGEQP3
#ifndef NOLAP
#define NOLAP
#endif
#endif

#include "TMV_Blas.h"
#include "tmv/TMV_QRPD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "TMV_QRDiv.h"
#include "tmv/TMV_Householder.h"
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
        MatrixView<T> A, VectorView<T> beta, ptrdiff_t* P, T& det,
        bool strict)
    {
#ifdef XDEBUG
        Matrix<T> A0(A);
        cout<<"Start NonBlock QRPD with strict = "<<strict<<endl;
        cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        cout<<"beta = "<<TMV_Text(beta)<<"  "<<beta<<endl;
        cout<<"Norm(A) = "<<Norm(A)<<std::endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
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
        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();
        const ptrdiff_t Astepj = A.stepj();
        const RT sqrteps = TMV_SQRT(TMV_Epsilon<T>());

        // Keep track of the norm of each column
        // When considering column j, these are actually just the norm
        // of each column from j:M, not 0:M.
        RT scale = RT(1) / A.maxAbs2Element(); // for more stable normSq
        //cout<<"scale = "<<scale<<std::endl;
        //cout<<"maxAbs = "<<A.maxAbsElement()<<std::endl;
        //cout<<"maxAbs2 = "<<A.maxAbs2Element()<<std::endl;
        Vector<RT> colnormsq(N);
        for(ptrdiff_t j=0;j<N;++j) colnormsq(j) = A.col(j).normSq(scale);
        //cout<<"colnormsq = "<<colnormsq<<std::endl;
        RT anormsq = colnormsq.sumElements();
        RT thresh = RT(N) * TMV_SQR(TMV_Epsilon<T>()) * anormsq;
        //cout<<"anormsq = "<<anormsq<<std::endl;
        //cout<<"thresh = "<<thresh<<std::endl;
        // Set to 0 any diag element whose norm is < epsilon * |A|
        RT recalcthresh(0);
        // recalcthresh is the threshold for recalculating the norm to account 
        // for rounding errors in the subtractions which keep track of it..
        // The is set to sqrt(TMV_Epsilon) * the largest normsq whenever we 
        // recalculate the norms. 

        T* bj = beta.ptr();
        for(ptrdiff_t j=0;j<N;++j,++bj) {
            //cout<<"j = "<<j<<std::endl;
            //cout<<"colnormsq = "<<colnormsq<<std::endl;
            //cout<<"recalcthresh = "<<recalcthresh<<std::endl;
            if (strict || j==0 || colnormsq(j) < recalcthresh) {
                // Find the column with the largest norm
                ptrdiff_t jpiv;
                RT maxnormsq = colnormsq.subVector(j,N).maxElement(&jpiv);
                //cout<<"jpiv = "<<jpiv<<std::endl;
                if (j==0) recalcthresh = 4*sqrteps * maxnormsq;
                //cout<<"maxnormsq = "<<maxnormsq<<std::endl;
                //cout<<"recalcthresh = "<<recalcthresh<<std::endl;
                // Note: jpiv is relative to the subVector(j,N)

                // If the largest colnormsq is lower than the recalulation 
                // threshold, then recalc all colnormsq's, and redetermine max.
                if (maxnormsq < recalcthresh) {
                    //cout<<"do recalc\n";
                    for(ptrdiff_t k=j;k<N;++k) 
                        colnormsq(k) = A.col(k,j,M).normSq(scale);
                    //cout<<"colnormsq => "<<colnormsq<<std::endl;
                    maxnormsq = colnormsq.subVector(j,N).maxElement(&jpiv);
                    recalcthresh = 4*sqrteps* maxnormsq;
                    if (recalcthresh < thresh) recalcthresh = thresh;
                    //cout<<"maxnormsq => "<<maxnormsq<<std::endl;
                    //cout<<"recalcthresh => "<<recalcthresh<<std::endl;
                }

                // If maxnormsq = 0 (technically < thresh to account for
                // rounding) then the rest of the R matrix is 0, and the 
                // Householder matrices are identities (indicated by 0's 
                // in the Q part of the matrix).
                if (maxnormsq < thresh) {
                    //cout<<"maxnormsq < thresh\n";
                    A.subMatrix(j,M,j,N).setZero();
                    // Already essentially zero - make it exact
                    beta.subVector(j,N).setZero();
                    // Set the Householder matrices for these to identities
                    for(;j<N;j++) P[j] = j;
                    break;
                } else {
                    //cout<<"apply jpiv = "<<jpiv<<std::endl;
                    // Swap the column with the largest norm into the current 
                    // column
                    if (jpiv != 0) {
                        // Add j to get real index
                        jpiv += j;
                        TMVAssert(jpiv < A.rowsize());
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
            TMVAssert(bj >= beta._first);
            TMVAssert(bj < beta._last);
#endif

            //cout<<"Before HouseholderReflect: A.col(j) = "<<A.col(j,j,M)<<std::endl;
            *bj = HouseholderReflect(A.subMatrix(j,M,j,N),det);
            //cout<<"bj = "<<*bj<<std::endl;

            // And update the norms for use with the next column
            const T* Ajk = A.row(j,j+1,N).cptr();
            for(ptrdiff_t k=j+1;k<N;++k,Ajk+=Astepj) {
                colnormsq(k) -= TMV_NORM(*Ajk*scale);
            }
        }
#ifdef XDEBUG
        cout<<"Done NonBlock QRPDecompose"<<std::endl;
        cout<<"A -> "<<A<<std::endl;
        cout<<"Norm(A) = "<<Norm(A)<<std::endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
        Matrix<T> Q(A);
        GetQFromQR(Q.view(),beta);
        cout<<"Q = "<<Q<<std::endl;
        Matrix<T> AA = Q*A.upperTri();
        cout<<"R = "<<A.upperTri()<<std::endl;
        cout<<"QR = "<<AA<<std::endl;
        AA.reversePermuteCols(P);
        cout<<"QRP = "<<AA<<std::endl;
        cout<<"Norm(AA-A0) = "<<Norm(AA-A0)<<std::endl;
        if (!(Norm(AA-A0) < 0.001*Norm(A0))) {
            cerr<<"NonBlockQRPDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> "<<A<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"P = ";
            for(ptrdiff_t i=0;i<A.rowsize();i++) cerr<<P[i]<<" ";
            cerr<<endl;
            cerr<<"QRP = "<<AA<<endl;
            cerr<<"A0 = "<<A0<<endl;
            cerr<<"diff = "<<Matrix<T>(AA-A0).clip(1.e-5)<<endl;
            cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
            abort(); 
        }
#endif
    }

    template <class T> 
    static void StrictBlockQRPDecompose(
        MatrixView<T> A, VectorView<T> beta, ptrdiff_t* P, T& det)
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
        cout<<"Start StrictBlockQRPDecompose: "<<std::endl;
        cout<<"Norm(A) = "<<Norm(A)<<std::endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
#endif

        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();
        const ptrdiff_t Astepj = A.stepj();
        const RT sqrteps = TMV_SQRT(TMV_Epsilon<T>());

        RT scale = RT(1) / A.maxAbs2Element(); // for more stable normSq
        //cout<<"scale = "<<scale<<std::endl;
        //cout<<"maxAbs = "<<A.maxAbsElement()<<std::endl;
        //cout<<"maxAbs2 = "<<A.maxAbs2Element()<<std::endl;
        Vector<RT> colnormsq(N);
        for(ptrdiff_t j=0;j<N;++j) colnormsq(j) = A.col(j).normSq(scale);
        RT anormsq = colnormsq.sumElements();
        RT thresh = RT(N) * TMV_SQR(TMV_Epsilon<T>()) * anormsq;
        RT recalcthresh(0);

        Matrix<T,RowMajor> ZYtA(TMV_MIN(QRP_BLOCKSIZE,int(M)),N);
        // We keep track of the product ZYtA(j1:M,j1:N) [ stored as ZYtA ]
        // since this is the product that we need.  We update this one 
        // row at a time.

        T* bj = beta.ptr();
        for(ptrdiff_t j1=0;j1<N;) {
            ptrdiff_t j2 = TMV_MIN(N,j1+QRP_BLOCKSIZE);
            for(ptrdiff_t j=j1,jmj1=0; j<j2; ) {
                ptrdiff_t jpiv;
                RT maxnormsq = colnormsq.subVector(j,N).maxElement(&jpiv);
                if (recalcthresh == RT(0)) {
                    recalcthresh = 4*sqrteps * maxnormsq;
                }

                if (maxnormsq < recalcthresh) {
                    // If need to recalculate norms, then do so and end the 
                    // block here.
                    // The columns are only correct if j == j1.  Otherwise,
                    // indicate to end the block, which will update the columns.
                    // Then next time through, we can do the recalc as needed.
                    if (j==j1) {
                        for(ptrdiff_t k=j;k<N;++k) 
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
                        TMVAssert(jpiv < A.rowsize());
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
                    TMVAssert(bj >= beta._first);
                    TMVAssert(bj < beta._last);
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
                    for(ptrdiff_t k=j+1;k<N;++k,Ajk+=Astepj) {
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
        cout<<"Done StrictBlockQRPDecompose\n";
        cout<<"Norm(A) = "<<Norm(A)<<std::endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
        if (!(Norm(AA-A0) < 0.001*Norm(A0))) {
            cerr<<"StrictBlockQRPDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> "<<A<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"P = ";
            for(ptrdiff_t i=0;i<A.rowsize();i++) cerr<<P[i]<<" ";
            cerr<<endl;
            cerr<<"QRP = "<<AA<<endl;
            cerr<<"A0 = "<<A0<<endl;
            cerr<<"diff = "<<Matrix<T>(AA-A0).clip(1.e-5)<<endl;
            cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
            abort(); 
        }
#endif
    }

    template <class T> 
    static void moveLowColsToEnd(
        Vector<RT>& colnormsq, RT thresh, ptrdiff_t j1, ptrdiff_t& j2, ptrdiff_t& j3,
        MatrixView<T> A, ptrdiff_t* P)
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
        for(ptrdiff_t i=j1;i<j2;++i) {
            if (colnormsq(i) < thresh) {
                TMVAssert(j3 < A.rowsize());
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
    static void checkIndex(const Vector<double>& index, const ptrdiff_t* P, ptrdiff_t j1)
    {
        const ptrdiff_t N = index.size();
        Vector<double> index2(N);
        for(ptrdiff_t k=0;k<N;k++) index2(k) = double(k);
        for(ptrdiff_t k=0;k<j1;k++) index2.swap(k,P[k]);
        if (Norm(index-index2) > 0.01) {
            cout<<"index = "<<index<<endl;
            cout<<"index2 = "<<index2<<endl;
            cout<<"Norm(diff) = "<<Norm(index-index2)<<endl;
            abort();
        }
    }
#endif

    template <class T> 
    static void LooseBlockQRPDecompose(
        MatrixView<T> A, VectorView<T> beta, ptrdiff_t* P, T& det)
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
        cout<<"Start LooseBlockQRPDecompose: "<<std::endl;
        cout<<"Norm(A) = "<<Norm(A)<<std::endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
#endif

        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();
        const ptrdiff_t Astepj = A.stepj();
        const RT sqrteps = TMV_SQRT(TMV_Epsilon<T>());

        RT scale = RT(1) / A.maxAbs2Element(); // for more stable normSq
        //cout<<"scale = "<<scale<<std::endl;
        //cout<<"maxAbs = "<<A.maxAbsElement()<<std::endl;
        //cout<<"maxAbs2 = "<<A.maxAbs2Element()<<std::endl;
        Vector<RT> colnormsq(N);
        for(ptrdiff_t j=0;j<N;++j) colnormsq(j) = A.col(j).normSq(scale);
        //cout<<"colnormsq = "<<colnormsq<<std::endl;
        RT anormsq = colnormsq.sumElements();
        RT thresh = RT(N) * TMV_SQR(TMV_Epsilon<T>()) * anormsq;
        //cout<<"anormsq = "<<anormsq<<std::endl;
        //cout<<"thresh = "<<thresh<<std::endl;
        //cout<<"Norm(A) = "<<Norm(A)<<std::endl;

#ifdef XDEBUG
        Vector<double> index(N);
        for(ptrdiff_t k=0;k<N;k++) index(k) = double(k);
#endif

        for (ptrdiff_t j1 = 0; j1 < N;) {
            // Do as many columns as possible such that none have to have 
            // their norms recalculated.
            ptrdiff_t j3=N; // j3 will be how many we have done in this loop
            // Invariant: all columns from j3..N are known to have norms that
            // need to be recalculated.
            // The recalculation is done at the end of the loop.

            ptrdiff_t jpiv0;
            RT maxnormsq = colnormsq.subVector(j1,N).maxElement(&jpiv0);
            //cout<<"j1 = "<<j1<<", j3 = "<<j3<<", jpiv = "<<jpiv0<<std::endl;
            //cout<<"colnormsq = "<<colnormsq<<std::endl;
            //cout<<"maxnormsq = "<<maxnormsq<<std::endl;
            //cout<<"Norm(A) = "<<Norm(A)<<std::endl;

            if (maxnormsq < thresh) {
                //cout<<"OK, zero the rest and we're done.\n";
                // Zero the rest out and we are done.
                A.subMatrix(j1,M,j1,N).setZero();
                beta.subVector(j1,N).setZero();
                for(;j1<N;j1++) P[j1] = j1;
                break;
            } 

            // Move max column to the front:
            if (jpiv0 != 0) {
                //cout<<"pivot\n";
                jpiv0 += j1;
                TMVAssert(jpiv0 < A.rowsize());
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

            RT recalcthresh = RT(N)*sqrteps*maxnormsq;
            if (recalcthresh < thresh) recalcthresh = thresh;

            TMVAssert(j1<j3);
            ptrdiff_t j1x = j1+1; 
            // The first pass through, we don't want to include j1 in the 
            // moveLowColsToEnd call.

            // Work on this one block at a time:
            while (j1 < j3) {
                ptrdiff_t j2 = TMV_MIN(j3,j1+QRP_BLOCKSIZE);
                //cout<<"j1,j2,j3 = "<<j1<<','<<j2<<','<<j3<<std::endl;
                //cout<<"Norm(A) = "<<Norm(A)<<std::endl;
                TMVAssert(j1 - j2 < 0);
                moveLowColsToEnd(colnormsq,recalcthresh,j1x,j2,j3,A,P);
#ifdef XDEBUG
                for(ptrdiff_t k=j1x;k<j2;k++) index.swap(k,P[k]);
                checkIndex(index,P,j2);
#endif

                ptrdiff_t origj2 = j2;
                UpperTriMatrix<T,NonUnitDiag|ColMajor> Z(j2-j1);

                T* bj = beta.ptr()+j1*beta.step();
                for(ptrdiff_t j=j1; j<j2; ++j, ++bj) {
                    //cout<<"j = "<<j<<std::endl;
                    //cout<<"Norm(A) = "<<Norm(A)<<std::endl;

                    if (colnormsq(j) < recalcthresh) {
                        ////cout<<"swap this column out\n";
#ifdef XDEBUG
                        checkIndex(index,P,origj2);
#endif
                        --j2;
#ifdef TMV_INITIALIZE_NAN
                        // Need this to avoid gratuitously populating
                        // uninitialized memory with nan's.
                        // Sometimes LAPACK or BLAS functions will allocate
                        // memory internally, and screw up if they get
                        // memory with nan's.  The INITIALIZE_NAN stuff
                        // is only designed to check that the TMV code
                        // is able to handle such values.
                        Z.col(j2-j1,0,j2-j1+1).setZero();
#endif
                        if (j==j2) break;
                        A.swapCols(j,j2);
                        colnormsq.swap(j,j2);
#ifdef XDEBUG
                        index.swap(j,j2);
#endif
                        if (P[j2] > j2) {
                            if (P[j] > j2) {
                                TMVAssert(P[j] < A.rowsize());
                                TMVAssert(P[j2] < A.rowsize());
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
                    //cout<<"Norm(A) = "<<Norm(A)<<std::endl;

                    // Find Householder matrix for this column
                    // This multiplies through to the end of the original block.
                    // This way, when we are done, the whole block has had the
                    // same Householder reflections applied to it.
#ifdef TMVFLDEBUG
                    TMVAssert(bj >= beta._first);
                    TMVAssert(bj < beta._last);
#endif
                    *bj = HouseholderReflect(A.subMatrix(j,M,j,origj2),det);
                    //cout<<"After reflect Norm(A) = "<<Norm(A)<<std::endl;

                    // Update Z:
                    BlockHouseholderAugment(
                        A.subMatrix(j1,M,j1,j+1),
                        Z.subTriMatrix(0,j-j1+1),TMV_CONJ(*bj));
                    //cout<<"After Householder Augment Norm(A) = "<<Norm(A)<<std::endl;

                    // Update the colnormsq values within this block
                    // (No need to go all the way to origj2, since the 
                    // j2..origj2 columns are those with low norm already -
                    // we don't need those values until we recalculate them 
                    // from scratch anyway.)
                    const T* Ajk = A.row(j,j+1,j2).cptr();
                    for(ptrdiff_t k=j+1;k<j2;++k,Ajk+=Astepj) 
                        colnormsq(k) -= TMV_NORM(*Ajk*scale);
                }

                if (j1 < j2) {

                    // Do the Block Householder update of the rest of the 
                    // matrix:
                    //cout<<"Before BlockLDiv Norm(A) = "<<Norm(A)<<std::endl;
                    //cout<<"j1,j2 = "<<j1<<" "<<j2<<std::endl;
                    //cout<<"origj2 = "<<origj2<<std::endl;
                    //cout<<"Z = "<<Z.subTriMatrix(0,j2-j1)<<std::endl;
                    //cout<<"Norm(Z) = "<<Norm(Z.subTriMatrix(0,j2-j1))<<std::endl;
                    //cout<<"Norm(A1) = "<<Norm(A.subMatrix(j1,M,j1,j2))<<std::endl;
                    //cout<<"Norm(A2) = "<<Norm(A.subMatrix(j1,M,origj2,N))<<std::endl;
                    BlockHouseholderLDiv(
                        A.subMatrix(j1,M,j1,j2),
                        Z.subTriMatrix(0,j2-j1),A.subMatrix(j1,M,origj2,N));
                    //cout<<"After BlockLDiv Norm(A) = "<<Norm(A)<<std::endl;

                    // Update the colnormsq values for the rest of the matrix:
                    if (M-j2 > j2-j1)
                        for(ptrdiff_t k=origj2;k<N;++k) colnormsq(k) -= 
                            A.col(k,j1,j2).normSq(scale);
                    else 
                        for(ptrdiff_t k=origj2;k<N;++k) colnormsq(k) = 
                            A.col(k,j2,M).normSq(scale);
                }

                if (j2 < origj2) {
#ifdef XDEBUG
                    checkIndex(index,P,origj2);
#endif
                    // Put the bad columns back where they started before this 
                    // loop:
                    for(ptrdiff_t j=j2; j<origj2; ++j) if (P[j] > j2) {
                        TMVAssert(P[j] < A.rowsize());
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
                //cout<<"Done block: Norm(A) = "<<Norm(A)<<std::endl;

                j1 = j1x = j2;
            }

            if (j3 < N) {
                //cout<<"recalculate colnorms\n";
                // Then need to recalculate some of the colnorms:
                for(ptrdiff_t k=j3;k<N;++k) 
                    colnormsq(k) = A.col(k,j3,M).normSq(scale);
                //cout<<"colnormsq = "<<colnormsq<<std::endl;
            }
        }

        if (det != T(0)) {
            for(ptrdiff_t i=0;i<N;++i) if (P[i] != i) det = -det;
        }

#ifdef XDEBUG
        cout<<"Done LooseBlockQRPDecompose\n";
        cout<<"Norm(A) = "<<Norm(A)<<std::endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
        checkIndex(index,P,N);
        Matrix<T> Q = A;
        GetQFromQR(Q.view(),beta);
        Matrix<T> AA = Q*A.upperTri();
        AA.reversePermuteCols(P);
        if (!(Norm(AA-A0) < 0.001*TMV_MAX(RT(1),Norm(A0)))) {
            cerr<<"LooseBlockQRPDecompose: "<<std::endl;
            cerr<<"A = "<<TMV_Text(A)<<endl;
            if (N < 100) {
                cerr<<"  "<<A0<<endl;
                cerr<<"-> "<<A<<endl;
                cerr<<"beta = "<<beta<<endl;
                cerr<<"P = ";
                for(ptrdiff_t i=0;i<N;i++) cerr<<P[i]<<" ";
                cerr<<endl;
                cerr<<"QRP = "<<AA<<endl;
                cerr<<"A0 = "<<A0<<endl;
                cerr<<"diff = "<<Matrix<T>(AA-A0).clip(1.e-5)<<endl;
            }
            cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
            cerr<<"Rdiag = "<<A.diag()<<endl;
            abort(); 
        }
#endif
    }

#undef RT

    template <class T> 
    static inline void NonLapQRPDecompose(
        MatrixView<T> A,
        VectorView<T> beta, ptrdiff_t* P, T& det, bool strict)
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
        MatrixView<T> A, VectorView<T> beta, ptrdiff_t* P, T& det)
    { NonLapQRPDecompose(A,beta,P,det,true); }
#ifdef INST_DOUBLE
    template <> 
    void LapQRPDecompose(
        MatrixView<double> A, VectorView<double> beta, ptrdiff_t* P, double& det)
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
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
        beta.setZero();
        int lda = A.stepj();
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = 3*n+1;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<double> work(1);
        work.get()[0] = 0.;
        LAPNAME(dgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(dgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"dgeqp3");
#else
        LAP_Results(Lap_info,int(work[0]),m,n,lwork,"dgeqp3");
#endif
        double thresh = TMV_Epsilon<double>()*A.normF();
        for(int i=0;i<n;++i) {
            if (TMV_ABS(A(i,i)) < thresh) {
                A.subMatrix(i,A.colsize(),i,A.rowsize()).setZero();
                beta.subVector(i,n).setZero();
                break;
            }
        }
#ifndef CLAP
        for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
        ConvertIndexToPermute(n,lap_p.get(),P);
        double* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi) {
            if (TMV_ABS(*bi-1.) > 1.01) *bi = 0.;
            if (det) {
                if (*bi != 0.) det = -det;
                if (P[i] != i) det = -det;
            }
        }
    }
    template <> 
    void LapQRPDecompose(
        MatrixView<std::complex<double> > A,
        VectorView<std::complex<double> > beta, ptrdiff_t* P,
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
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
        beta.setZero();
        int lda = A.stepj();
        int Lap_info=0;
#ifndef LAPNOWORK
        AlignedArray<double> rwork(2*n);
        VectorViewOf(rwork.get(),2*n).setZero();
#ifdef NOWORKQUERY
        int lwork = n+1;
        AlignedArray<std::complex<double> > work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<std::complex<double> > work(1);
        work.get()[0] = 0.;
        LAPNAME(zgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
        lwork = int(TMV_REAL(work[0]));
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(zgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"zgeqp3");
#else
        LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"zgeqp3");
#endif
        beta.conjugateSelf();
        double thresh = TMV_Epsilon<double>()*A.normF();
        const std::complex<double>* Aii = A.diag().cptr();
        const int Ads = A.diag().step();
        for(int i=0;i<n;++i,Aii+=Ads) {
            TMVAssert(TMV_ABS(TMV_IMAG(*Aii)) < TMV_Epsilon<double>());
            if (TMV_ABS(TMV_REAL(*Aii)) < thresh) {
                A.subMatrix(i,A.colsize(),i,A.rowsize()).setZero();
                beta.subVector(i,n).setZero();
                break;
            }
        }
#ifndef CLAP
        for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
        ConvertIndexToPermute(n,lap_p.get(),P);
        std::complex<double>* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi) {
            if (TMV_ABS(*bi-1.) > 1.01) *bi = 0.;
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
        MatrixView<float> A, VectorView<float> beta, ptrdiff_t* P, float& det)
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
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
        beta.setZero();
        int lda = A.stepj();
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = 3*n+1;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<float> work(1);
        work.get()[0] = 0.;
        LAPNAME(sgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(sgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"sgeqp3");
#else
        LAP_Results(Lap_info,int(work[0]),m,n,lwork,"sgeqp3");
#endif
        float thresh = TMV_Epsilon<float>()*A.normF();
        for(int i=0;i<n;++i) {
            if (TMV_ABS(A(i,i)) < thresh) {
                A.subMatrix(i,A.colsize(),i,A.rowsize()).setZero();
                beta.subVector(i,n).setZero();
                break;
            }
        }
#ifndef CLAP
        for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
        ConvertIndexToPermute(n,lap_p.get(),P);
        float* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi) {
            if (TMV_ABS(*bi-1.F) > 1.01F) *bi = 0.F;
            if (det) {
                if (*bi != 0.) det = -det;
                if (P[i] != i) det = -det;
            }
        }
    }
    template <> 
    void LapQRPDecompose(
        MatrixView<std::complex<float> > A,
        VectorView<std::complex<float> > beta, ptrdiff_t* P,
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
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
        beta.setZero();
        int lda = A.stepj();
        int Lap_info=0;
#ifndef LAPNOWORK
        AlignedArray<float> rwork(2*n);
        VectorViewOf(rwork.get(),2*n).setZero();
#ifdef NOWORKQUERY
        int lwork = n+1;
        AlignedArray<std::complex<float> > work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<std::complex<float> > work(1);
        work.get()[0] = 0.;
        LAPNAME(cgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
        lwork = int(TMV_REAL(work[0]));
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(cgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPWK(rwork.get()) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"cgeqp3");
#else
        LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"cgeqp3");
#endif
        float thresh = TMV_Epsilon<float>()*A.normF();
        const std::complex<float>* Aii = A.diag().cptr();
        const int Ads = A.diag().step();
        for(int i=0;i<n;++i,Aii+=Ads) {
            TMVAssert(TMV_ABS(TMV_IMAG(*Aii)) < TMV_Epsilon<float>());
            if (TMV_ABS(TMV_REAL(*Aii)) < thresh) {
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
        std::complex<float>* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi) {
            if (TMV_ABS(*bi-1.F) > 1.01F) *bi = 0.F;
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
        MatrixView<T> A, VectorView<T> beta, ptrdiff_t* P, T& det, bool strict)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == beta.size());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step()==1);

#ifdef XDEBUG
        Matrix<T> A0(A);
        cout<<"Start QRPDecompose:"<<std::endl;
        cout<<"A = "<<A0<<std::endl;
        cout<<"strict = "<<strict<<std::endl;
        cout<<"Norm(A) = "<<Norm(A)<<std::endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
#ifdef LAP
        Matrix<T> A2(A);
        Vector<T> beta2(beta);
        AlignedArray<ptrdiff_t> P2(beta.size());
        T det2=det;
        NonLapQRPDecompose(A2.view(),beta2.view(),P2.get(),det2,strict);
        cout<<"NonLap QRP = "<<A2<<std::endl;
#endif
#endif

        if (A.rowsize() > 0) {
#ifdef LAP
            if (strict && BlasIsCM(A)) {
                LapQRPDecompose(A,beta,P,det);
            } else {
                NonLapQRPDecompose(A,beta,P,det,strict);
            }
#else
            NonLapQRPDecompose(A,beta,P,det,strict);
#endif
        }
#ifdef XDEBUG
        cout<<"Done QRPDecompose\n";
        cout<<"Norm(A) = "<<Norm(A)<<std::endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
        Matrix<T> Q(A);
        cout<<"Q = "<<Q<<std::endl;
        GetQFromQR(Q.view(),beta);
        cout<<"Q => "<<Q<<std::endl;
        Matrix<T> R = A.upperTri();
        cout<<"R = "<<R<<std::endl;
        Matrix<T> AA = Q*R;
        cout<<"AA = "<<AA<<std::endl;
        AA.reversePermuteCols(P);
        cout<<"Norm(AA -A0) = "<<Norm(AA-A0)<<std::endl;
        if (!(Norm(AA-A0) < 0.001*Norm(A0))) {
            cerr<<"QRPDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> "<<A<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"P = ";
            for(ptrdiff_t i=0;i<A.rowsize();i++) cerr<<P[i]<<" ";
            cerr<<endl;
            cerr<<"QRP = "<<AA<<endl;
            cerr<<"A0 = "<<A0<<endl;
            cerr<<"diff = "<<Matrix<T>(AA-A0).clip(1.e-5)<<endl;
            cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
#ifdef LAP
            cerr<<"NonLap: A = "<<A2<<endl;
            cerr<<"diff = "<<Matrix<T>(A-A2).clip(1.e-5)<<endl;
            cerr<<"beta2 = "<<beta2<<endl;
            cerr<<"det2 = "<<det2<<endl;
#endif
            abort(); 
        }
#endif
    }

    //
    // QRP Decompose - Unpacked
    //

    template <class T> 
    void QRP_Decompose(
        MatrixView<T> Q, UpperTriMatrixView<T> R, ptrdiff_t* P, bool strict)
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
    void QRP_Decompose(MatrixView<T> A, bool strict)
    {
        // Decompose A into Q R P, but don't keep Q or P.
        // R is returned as A.upperTri().

        TMVAssert(A.colsize() >= A.rowsize());

        Vector<T> beta(A.rowsize());
        T d(0);
        AlignedArray<ptrdiff_t> P(A.rowsize());
        if (A.isconj())
            QRP_Decompose(A.conjugate(),beta.view(),P.get(),d,strict);
        else
            QRP_Decompose(A,beta.view(),P.get(),d,strict);
    }


#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_QRPDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


