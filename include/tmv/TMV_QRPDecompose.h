

#ifndef TMV_QRPDecompose_H
#define TMV_QRPDecompose_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Householder.h"
#include "TMV_Permutation.h"
//#include "TMV_MultMM.h"
//#include "TMV_DivVU.h"
//#include "TMV_DivMU.h"

#ifdef PRINTALGO_QR
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_TriMatrixIO.h"
#include "TMV_ProdMM.h"
#include "TMV_MultUL.h"
#include "TMV_PermuteM.h"
#endif

#ifdef XDEBUG_QR
#include "tmv/TMV_UnpackQ.h"
#include "tmv/TMV_AddVV.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_SumVV.h"
#include "tmv/TMV_Vector.h"
#endif

// BLOCKSIZE is the block size to use in algo 21, etc.
#define TMV_QR_BLOCKSIZE 48


namespace tmv {

    // Defined in TMV_QRPDecompose.cpp
    template <class T, class RT>
    void InstQRP_Decompose(
        MatrixView<T> m, VectorView<RT> beta, int* P, bool strict);

    template <int algo, int cs, int rs, class M, class V>
    struct QRPDecompose_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or 1, or N == 0)
    template <int cs, int rs, class M, class V>
    struct QRPDecompose_Helper<0,cs,rs,M,V>
    { static TMV_INLINE void call(M& , V& , int* , bool ) {} };

    // algo 11: Non-block algorithm, loop over n
    // This is basically the same algorithm for strict and not.
    template <int cs, int rs, class M1, class V>
    struct QRPDecompose_Helper<11,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta, int* P, bool strict)
        {
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const int M = cs==TMV_UNKNOWN ? int(A.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(A.rowsize()) : rs;
#ifdef PRINTALGO_QR
            std::cout<<"QRPDecompose algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif

#ifdef XDEBUG_QR
            Matrix<T> A0(A);
            std::cout<<"Start NonBlock QRPD with strict = "<<strict<<std::endl;
            std::cout<<"A = "<<TMV_Text(A)<<"  "<<A<<std::endl;
            std::cout<<"beta = "<<TMV_Text(beta)<<"  "<<beta<<std::endl;
            std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;
            std::cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
#endif
            // Decompose A into A = Q R P
            // where Q is unitary, R is upper triangular, and P is a permutation
            // Q and R are stored in the same matrix (output of A), 
            // with the beta's for the Householder matrices returned in beta.

            const RT sqrteps = TMV_SQRT(TMV_Epsilon<T>());

            typedef typename V::iterator IT;
            typedef typename M1::col_sub_type M1c;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;
            typedef typename VCopyHelper<T,rs>::type V2;
            typedef typename V2::subvector_type V2s;
            V2 tempBase = VectorSizer<T>(N);

            // Keep track of the norm of each column
            // When considering column j, these are actually just the norm
            // of each column from j:M, not 0:M.
            RT scale = RT(1) / A.maxAbs2Element(); // for more stable normSq
            //std::cout<<"scale = "<<scale<<std::endl;
            //std::cout<<"maxAbs = "<<A.maxAbsElement()<<std::endl;
            //std::cout<<"maxAbs2 = "<<A.maxAbs2Element()<<std::endl;
            Vector<RT> colnormsq(N);
            for(int j=0;j<N;++j) colnormsq(j) = A.col(j).normSq(scale);
            //std::cout<<"colnormsq = "<<colnormsq<<std::endl;
            RT anormsq = colnormsq.sumElements();
            RT thresh = RT(N) * TMV_SQR(TMV_Epsilon<T>()) * anormsq;
            //std::cout<<"anormsq = "<<anormsq<<std::endl;
            //std::cout<<"thresh = "<<thresh<<std::endl;
            // Set to 0 any diag element whose norm is < epsilon * |A|
            RT recalcthresh(0);
            // recalcthresh is the threshold for recalculating the norm to account 
            // for rounding errors in the subtractions which keep track of it..
            // The is set to sqrt(TMV_Epsilon) * the largest normsq whenever we 
            // recalculate the norms. 

            IT bj = beta.begin();
            for(int j=0;j<N;++j,++bj) {
                //std::cout<<"j = "<<j<<std::endl;
                //std::cout<<"A = "<<A<<std::endl;
                if (!(Norm(A) >= 0.)) abort();
                //std::cout<<"colnormsq = "<<colnormsq<<std::endl;
                //std::cout<<"recalcthresh = "<<recalcthresh<<std::endl;
                if (strict || j==0 || colnormsq(j) < recalcthresh) {
                    // Find the column with the largest norm
                    int jpiv;
                    RT maxnormsq = colnormsq.subVector(j,N).maxElement(&jpiv);
                    //std::cout<<"maxnormsq = "<<maxnormsq<<std::endl;
                    //std::cout<<"jpiv = "<<jpiv<<std::endl;
                    if (j==0) recalcthresh = 4*sqrteps * maxnormsq;
                    //std::cout<<"recalcthresh = "<<recalcthresh<<std::endl;
                    // Note: jpiv is relative to the subVector(j,N)

                    // If the largest colnormsq is lower than the recalulation 
                    // threshold, then recalc all colnormsq's, and redetermine max.
                    if (maxnormsq < recalcthresh) {
                        //std::cout<<"do recalc\n";
                        for(int k=j;k<N;++k) 
                            colnormsq(k) = A.col(k,j,M).normSq(scale);
                        //std::cout<<"colnormsq => "<<colnormsq<<std::endl;
                        maxnormsq = colnormsq.subVector(j,N).maxElement(&jpiv);
                        recalcthresh = 4*sqrteps* maxnormsq;
                        if (recalcthresh < thresh) recalcthresh = thresh;
                        //std::cout<<"maxnormsq => "<<maxnormsq<<std::endl;
                        //std::cout<<"recalcthresh => "<<recalcthresh<<std::endl;
                    }

                    // If maxnormsq = 0 (technically < thresh to account for
                    // rounding) then the rest of the R matrix is 0, and the 
                    // Householder matrices are identities (indicated by 0's 
                    // in the Q part of the matrix).
                    if (maxnormsq < thresh) {
                        //std::cout<<"maxnormsq < thresh\n";
                        A.subMatrix(j,M,j,N).setZero();
                        // Already essentially zero - make it exact
                        beta.subVector(j,N).setZero();
                        // Set the Householder matrices for these to identities
                        for(;j<N;j++) P[j] = j;
                        break;
                    } else {
                        //std::cout<<"apply jpiv = "<<jpiv<<std::endl;
                        // Swap the column with the largest norm into the current 
                        // column
                        if (jpiv != 0) {
                            // Add j to get real index
                            jpiv += j;
                            TMVAssert(jpiv < int(A.rowsize()));
                            colnormsq.swap(j,jpiv);
                            A.swapCols(j,jpiv);
                            P[j] = jpiv;
                        } else {
                            P[j] = j;
                        }
                    }
                } else P[j] = j;

                // Apply the Householder Reflection for this column
                //std::cout<<"Before HouseholderReflect: A.col(j) = "<<
                    //A.col(j,j,M)<<std::endl;
                //std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;
                M1c u = A.col(j,j+1,M);
                HouseholderReflect(A.ref(j,j),u,*bj);
                //std::cout<<"After HouseholderReflect\n";
                //std::cout<<"bj = "<<*bj<<std::endl;
                //std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;
                M1r A2a = A.row(j,j+1,N);
                M1s A2b = A.subMatrix(j+1,M,j+1,N);
                V2s temp = tempBase.subVector(0,N-j-1);
                HouseholderMultEq(u,*bj,A2a,A2b,temp);
                //std::cout<<"After HouseholderMultEq\n";
                //std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;

                // And update the norms for use with the next column
                for(int k=j+1;k<N;++k) {
                    colnormsq(k) -= TMV_NORM(A.cref(j,k)*scale);
                }
            }
#ifdef XDEBUG_QR
            std::cout<<"Done NonBlock QRPDecompose"<<std::endl;
            std::cout<<"A -> "<<A<<std::endl;
            if (!(Norm(A) >= 0.)) abort();
            std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;
            std::cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
            Matrix<T> Q(A);
            UnpackQ(Q,beta);
            std::cout<<"Q = "<<Q<<std::endl;
            Matrix<T> AA = Q*A.upperTri();
            std::cout<<"R = "<<A.upperTri()<<std::endl;
            std::cout<<"QR = "<<AA<<std::endl;
            AA.reversePermuteCols(P);
            std::cout<<"QRP = "<<AA<<std::endl;
            std::cout<<"Norm(AA-A0) = "<<Norm(AA-A0)<<std::endl;
            if (!(Norm(AA-A0) <= 0.001*Norm(A0))) {
                std::cerr<<"NonBlockQRPDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<std::endl;
                std::cerr<<"-> "<<A<<std::endl;
                std::cerr<<"beta = "<<beta<<std::endl;
                std::cerr<<"P = ";
                for(int i=0;i<int(A.rowsize());i++) std::cerr<<P[i]<<" ";
                std::cerr<<std::endl;
                std::cerr<<"QRP = "<<AA<<std::endl;
                std::cerr<<"A0 = "<<A0<<std::endl;
                std::cerr<<"diff = "<<Matrix<T>(AA-A0).clip(1.e-5)<<std::endl;
                std::cerr<<"Norm(diff) = "<<Norm(AA-A0)<<std::endl;
                abort(); 
            }
#endif
        }
    };

    // algo 21: Block strict algorithm
    template <int cs, int rs, class M1, class V>
    struct QRPDecompose_Helper<21,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta, int* P, bool)
        {
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const int M = cs==TMV_UNKNOWN ? int(A.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(A.rowsize()) : rs;
#ifdef PRINTALGO_QR
            std::cout<<"QRPDecompose algo 21: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            //typedef typename M1::copy_type M1x;
            //typedef typename V::copy_type Vx;
            //M1x Ac = A;
            //Vx bc = beta;
#endif
            // Decompose A (input as A) into A = Q R P
            // where Q is unitary, R is upper triangular, and P is a permutation
            // Q and R are stored in the same matrix (A), with the beta's for
            // the Householder matrices returned in beta.

#ifdef XDEBUG_QR
            Matrix<T> A0(A);
            //std::cout<<"Start StrictBlockQRPDecompose: "<<std::endl;
            //std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;
            //std::cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
#endif
            const RT sqrteps = TMV_SQRT(TMV_Epsilon<T>());

            RT scale = RT(1) / A.maxAbs2Element(); // for more stable normSq
            //std::cout<<"scale = "<<scale<<std::endl;
            //std::cout<<"maxAbs = "<<A.maxAbsElement()<<std::endl;
            //std::cout<<"maxAbs2 = "<<A.maxAbs2Element()<<std::endl;
            Vector<RT> colnormsq(N);
            for(int j=0;j<N;++j) colnormsq(j) = A.col(j).normSq(scale);
            RT anormsq = colnormsq.sumElements();
            RT thresh = RT(N) * TMV_SQR(TMV_Epsilon<T>()) * anormsq;
            RT recalcthresh(0);

            typedef typename V::iterator IT;
            typedef typename M1::col_sub_type M1c;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;

            const int Nx = TMV_QR_BLOCKSIZE;
            const int rs1 = IntTraits2<rs,Nx>::min;
            typedef typename VCopyHelper<T,rs1>::type V2;
            typedef typename V2::subvector_type V2s;
            V2 tempBase = VectorSizer<T>(TMV_MIN(N,Nx));

            // We keep track of the product ZYtA(j1:M,j1:N) [ stored as ZYtA ]
            // since this is the product that we need.  We update this one 
            // row at a time.
            const int cs1 = IntTraits2<cs,Nx>::min;
            // TODO: Double check this.  Is RowMajor really better?
            typedef typename MCopyHelper<T,Rec,cs1,rs,true>::type M3;
            M3 ZYtA = MatrixSizer<T>(TMV_MIN(M,Nx),N);
            // TODO: Would this be faster for the regular QRDecomposition too?
            // i.e. Rather than using the normal BlockHouseholder calls.

            IT bj = beta.begin();
            for(int j1=0;j1<N;) {
                int j2 = TMV_MIN(N,j1+Nx);
                //std::cout<<"Start block "<<j1<<" .. "<<j2<<std::endl;
                //std::cout<<"A = "<<A<<std::endl;
                for(int j=j1,jmj1=0; j<j2; ) {
                    //std::cout<<"Start j = "<<j<<std::endl;
                    //std::cout<<"A = "<<A<<std::endl;
                    //std::cout<<"colnormsq = "<<colnormsq<<std::endl;
                    //std::cout<<"recalcthresh = "<<recalcthresh<<std::endl;
                    int jpiv;
                    RT maxnormsq = colnormsq.subVector(j,N).maxElement(&jpiv);
                    //std::cout<<"maxnormsq = "<<maxnormsq<<std::endl;
                    //std::cout<<"jpiv = "<<jpiv<<std::endl;
                    if (recalcthresh == RT(0)) {
                        recalcthresh = 4*sqrteps * maxnormsq;
                    }
                    //std::cout<<"recalcthresh = "<<recalcthresh<<std::endl;

                    if (maxnormsq < recalcthresh) {
                        //std::cout<<"Need to recalc norms\n";
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
                        //std::cout<<"Done with this block.\n";
                        if (j==j1) {
                            // If first in set, then just zero the rest out and 
                            // indicate that we are done.
                            //std::cout<<"j == j1\n";
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
                        //std::cout<<"Normal case\n";
                        // normal case with columns still left to do.

                        // Pivot
                        if (jpiv != 0) {
                            jpiv += j;
                            TMVAssert(jpiv < int(A.rowsize()));
                            //std::cout<<"jpiv => "<<jpiv<<std::endl;
                            colnormsq.swap(j,jpiv);
                            A.swapCols(j,jpiv);
                            ZYtA.rowRange(0,jmj1).swapCols(j,jpiv);
                            P[j] = jpiv;
                        } else {
                            P[j] = j;
                        }

                        // Update the pivot column with Block Householder so far:
                        // A(j1:M,j) -= Y Z Yt A(j1:M,j)
                        // A(j1:j,j) has already been updated, so we only 
                        // need to do A(j:M,j)
                        // A(j:M,j) -= Y(j:M,0:j) (ZYtA)(0:j,j)
                        //std::cout<<"Update pivot column\n";
                        A.col(j,j,M) -= A.subMatrix(j,M,j1,j) * ZYtA.col(j,0,jmj1);

                        // Find Householder matrix for this column
                        M1c u = A.col(j,j+1,M);
                        HouseholderReflect(A.ref(j,j),u,*bj);
                        //std::cout<<"bj = "<<*bj<<std::endl;

                        // Update ZYtA:
                        if (*bj != T(0)) {
                            // (I-beta u ut)(I-Y Z Yt) A
                            // = I - Y (ZYtA) - u (beta ut A) + u (beta ut Y (ZYtA))
                            // The augmented Y now includes u in the j column, 
                            // so the augmented ZYtA now has to include in the 
                            // j row:
                            // beta (ut A - ut Y ZYtA)
                            typename M1c::conjugate_type ut = u.conjugate();
                            // Remember, this doesn't include the implicit 1 
                            // at the top of u.
                            ZYtA.row(jmj1,j1,j+1).setZero();
                            ZYtA.row(jmj1,j+1,N) = ut * A.subMatrix(j+1,M,j+1,N);
                            ZYtA.row(jmj1,j+1,N) += A.row(j,j+1,N);
                            V2s temp2 = tempBase.subVector(0,jmj1);
                            temp2 = ut * A.subMatrix(j+1,M,j1,j);
                            temp2 += A.row(j,j1,j);
                            ZYtA.row(jmj1,j1,N) -= 
                                temp2 * ZYtA.subMatrix(0,jmj1,j1,N);
                            ZYtA.row(jmj1,j1,N) *= *bj;
                        } else ZYtA.row(jmj1,j1,N).setZero();

                        // Update row j of the rest of the matrix:
                        // A(j,j+1:N) -= (Y ZYtA)(j,j+1:N) 
                        //             = Y(j,j1:j+1) ZYtA(j1:j+1,j+1:N)
                        M1r Arowj = A.row(j,j+1,N);
                        Arowj -= A.row(j,j1,j)*ZYtA.subMatrix(0,jmj1,j+1,N);
                        Arowj -= ZYtA.row(jmj1,j+1,N);

                        // Update the colnormsq values
                        for(int k=j+1;k<N;++k) {
                            colnormsq(k) -= TMV_NORM(A.cref(j,k)*scale);
                        }
                        ++j; ++jmj1; ++bj;
                    }
                }

                if (j2 < N) {
                    // Do the Block Householder update of the rest of the matrix:
                    // A(j2:M,j2:N) -= Y(j2:M,j1:j2) ZYtA(j1:j2,j1:N)
                    A.subMatrix(j2,M,j2,N) -= A.subMatrix(j2,M,j1,j2) * 
                        ZYtA.subMatrix(0,j2-j1,j2,N);
                }
                j1 = j2;
                //std::cout<<"Done block: j1 => "<<j1<<std::endl;
            }

#ifdef PRINTALGO_QR
            //std::cout<<"A -> "<<A<<std::endl;
            //std::cout<<"beta -> "<<beta<<std::endl;
            //QRPDecompose_Helper<11,cs,rs,M1x,Vx>::call(Ac,bc,P,true);
            //std::cout<<"correct A = "<<Ac<<std::endl;
            //std::cout<<"beta = "<<bc<<std::endl;
            //std::cout<<"diff = "<<(A-Ac).copy().clip(1.e-5)<<std::endl;
            //std::cout<<(beta-bc).copy().clip(1.e-5)<<std::endl;
#endif

#ifdef XDEBUG_QR
            Matrix<T> Q(A);
            UnpackQ(Q,beta);
            Matrix<T> AA = Q*A.upperTri();
            AA.reversePermuteCols(P);
            //std::cout<<"Done StrictBlockQRPDecompose\n";
            //std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;
            //std::cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
            if (!(Norm(AA-A0) <= 0.001*Norm(A0))) {
                std::cerr<<"StrictBlockQRPDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<std::endl;
                std::cerr<<"-> "<<A<<std::endl;
                std::cerr<<"beta = "<<beta<<std::endl;
                std::cerr<<"P = ";
                for(int i=0;i<int(A.rowsize());i++) std::cerr<<P[i]<<" ";
                std::cerr<<std::endl;
                std::cerr<<"QRP = "<<AA<<std::endl;
                std::cerr<<"A0 = "<<A0<<std::endl;
                std::cerr<<"diff = "<<Matrix<T>(AA-A0).clip(1.e-5)<<std::endl;
                std::cerr<<"Norm(diff) = "<<Norm(AA-A0)<<std::endl;
                typedef typename M1::copy_type M1c;
                typedef typename V::copy_type Vc;
                M1c A2 = A0;
                Vc beta2 = beta;
                int* P2 = new int[beta.size()];
                QRPDecompose_Helper<11,cs,rs,M1c,Vc>::call(A2,beta2,P2,true);
                std::cout<<"Algo 11:\n";
                std::cout<<"A = "<<A<<std::endl;
                std::cout<<"A2 = "<<A2<<std::endl;
                std::cout<<"diff = "<<A2-A<<std::endl;
                std::cout<<"Norm(A2-A) = "<<Norm(A2-A)<<std::endl;
                std::cout<<"beta = "<<beta<<std::endl;
                std::cout<<"beta2 = "<<beta2<<std::endl;
                std::cout<<"diff = "<<beta2-beta<<std::endl;
                std::cout<<"Norm(b2-b) = "<<Norm(beta2-beta)<<std::endl;
                std::cerr<<"P = ";
                for(int i=0;i<int(A.rowsize());i++) std::cerr<<P[i]<<" ";
                std::cerr<<"\nP2 = ";
                for(int i=0;i<int(A.rowsize());i++) std::cerr<<P2[i]<<" ";
                std::cerr<<"\n";
                abort(); 
            }
#endif
        }
    };

    // algo 22: Block non-strict algorithm
    template <int cs, int rs, class M1, class V>
    struct QRPDecompose_Helper<22,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        typedef typename M1::real_type RT;

        template <class RT, class M> 
        static void moveLowColsToEnd(
            Vector<RT>& colnormsq, RT thresh, int j1, int& j2, int& j3,
            M& A, int* P)
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

#ifdef XDEBUG_QR
        static void checkIndex(
            const Vector<double>& index, const int* P, int j1)
        {
            const int N = index.size();
            Vector<double> index2(N);
            for(int k=0;k<N;k++) index2(k) = double(k);
            for(int k=0;k<j1;k++) index2.swap(k,P[k]);
            if (Norm(index-index2) > 0.01) {
                std::cout<<"index = "<<index<<std::endl;
                std::cout<<"index2 = "<<index2<<std::endl;
                std::cout<<"Norm(diff) = "<<Norm(index-index2)<<std::endl;
                abort();
            }
        }
#endif

        static void call(M1& A, V& beta, int* P, bool)
        {
            const int M = cs==TMV_UNKNOWN ? int(A.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(A.rowsize()) : rs;
#ifdef PRINTALGO_QR
            std::cout<<"QRPDecompose algo 22: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            // Decompose A (input as A) into A = Q R P
            // where Q is unitary, R is upper triangular, and P is a permutation
            // Q and R are stored in the same matrix (A), with the beta's for
            // the Householder matrices returned in beta.
            //
            // This loose version doesn't sort the diagonal of R exactly.
            // It only sorts them enough to make sure the 0's fall at the end.

#ifdef XDEBUG_QR
            Matrix<T> A0(A);
            std::cout<<"Start LooseBlockQRPDecompose: "<<std::endl;
            std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;
            std::cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
#endif

            const RT sqrteps = TMV_SQRT(TMV_Epsilon<T>());

            typedef typename M1::col_sub_type M1c;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;
            typedef typename V::iterator IT;

            const int Nx = TMV_QR_BLOCKSIZE;
            const int s1 = IntTraits2<Nx,rs>::min;
            const int N1 = TMV_MIN(Nx,N);
            typedef typename MCopyHelper<T,UpperTri,s1,s1>::type Ztype;
            typedef typename Ztype::subtrimatrix_type Zs;
            Ztype BaseZ = MatrixSizer<T>(N1,N1);

            typedef typename MCopyHelper<T,Rec,s1,rs>::type M3;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3::submatrix_type M3s;
            M3 tempBase = MatrixSizer<T>(N1,N);

            RT scale = RT(1) / A.maxAbs2Element(); // for more stable normSq
            //std::cout<<"scale = "<<scale<<std::endl;
            //std::cout<<"maxAbs = "<<A.maxAbsElement()<<std::endl;
            //std::cout<<"maxAbs2 = "<<A.maxAbs2Element()<<std::endl;
            Vector<RT> colnormsq(N);
            for(int j=0;j<N;++j) colnormsq(j) = A.col(j).normSq(scale);
            //std::cout<<"colnormsq = "<<colnormsq<<std::endl;
            RT anormsq = colnormsq.sumElements();
            RT thresh = RT(N) * TMV_SQR(TMV_Epsilon<T>()) * anormsq;
            //std::cout<<"anormsq = "<<anormsq<<std::endl;
            //std::cout<<"thresh = "<<thresh<<std::endl;
            //std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;

#ifdef XDEBUG_QR
            Vector<double> index(N);
            for(int k=0;k<N;k++) index(k) = double(k);
#endif

            IT bj = beta.begin();
            for (int j1 = 0; j1 < N;) {
                // Do as many columns as possible such that none have to have 
                // their norms recalculated.
                int j3=N; // j3 will be how many we have done in this loop
                // Invariant: all columns from j3..N are known to have norms that
                // need to be recalculated.
                // The recalculation is done at the end of the loop.

                int jpiv0;
                RT maxnormsq = colnormsq.subVector(j1,N).maxElement(&jpiv0);
                //std::cout<<"j1 = "<<j1<<", j3 = "<<j3<<", jpiv = "<<jpiv0<<std::endl;
                //std::cout<<"colnormsq = "<<colnormsq<<std::endl;
                //std::cout<<"maxnormsq = "<<maxnormsq<<std::endl;
                //std::cout<<"thresh = "<<thresh<<std::endl;

                if (maxnormsq < thresh) {
                    //std::cout<<"OK, zero the rest and we're done.\n";
                    // Zero the rest out and we are done.
                    A.subMatrix(j1,M,j1,N).setZero();
                    beta.subVector(j1,N).setZero();
                    for(;j1<N;j1++) P[j1] = j1;
                    break;
                } 

                // Move max column to the front:
                if (jpiv0 != 0) {
                    //std::cout<<"pivot\n";
                    jpiv0 += j1;
                    TMVAssert(jpiv0 < int(A.rowsize()));
                    A.swapCols(j1,jpiv0);
                    colnormsq.swap(j1,jpiv0);
                    P[j1] = jpiv0;
#ifdef XDEBUG_QR
                    index.swap(j1,jpiv0);
#endif
                } else P[j1] = j1;
#ifdef XDEBUG_QR
                checkIndex(index,P,j1+1);
#endif

                RT recalcthresh = RT(N)*sqrteps*maxnormsq;
                if (recalcthresh < thresh) recalcthresh = thresh;

                TMVAssert(j1<j3);
                int j1x = j1+1; 
                // The first pass through, we don't want to include j1 in the 
                // moveLowColsToEnd call.

                // Work on this one block at a time:
                while (j1 < j3) {
                    int j2 = TMV_MIN(j3,j1+Nx);
                    //std::cout<<"j1,j2,j3 = "<<j1<<','<<j2<<','<<j3<<std::endl;
                    TMVAssert(j1 - j2 < 0);
                    moveLowColsToEnd(colnormsq,recalcthresh,j1x,j2,j3,A,P);
#ifdef XDEBUG_QR
                    for(int k=j1x;k<j2;k++) index.swap(k,P[k]);
                    checkIndex(index,P,j2);
#endif

                    int origj2 = j2;
                    typename Ztype::subtrimatrix_type Z = 
                        BaseZ.subTriMatrix(0,j2-j1);

                    for(int j=j1; j<j2; ++j, ++bj) {
                        //std::cout<<"j = "<<j<<std::endl;
                        //std::cout<<"j1,j2,j3,origj2 = "<<j1<<','<<j2<<','<<j3<<','<<origj2<<std::endl;

                        if (colnormsq(j) < recalcthresh) {
                            //std::cout<<"swap this column out\n";
#ifdef XDEBUG_QR
                            checkIndex(index,P,origj2);
#endif
                            --j2;
                            //std::cout<<"j2 => "<<j2<<std::endl;
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
#ifdef XDEBUG_QR
                            index.swap(j,j2);
#endif
                            if (P[j2] > j2) {
                                if (P[j] > j2) {
                                    TMVAssert(P[j] < int(A.rowsize()));
                                    TMVAssert(P[j2] < int(A.rowsize()));
                                    TMV_SWAP(P[j],P[j2]);
                                    A.swapCols(P[j],P[j2]);
                                    colnormsq.swap(P[j],P[j2]);
#ifdef XDEBUG_QR
                                    index.swap(P[j],P[j2]);
#endif
                                } else {
                                    P[j] = P[j2];
                                }
                            } else {
                                if (P[j] > j2) P[j2] = P[j];
                                P[j] = j2;
                            }
#ifdef XDEBUG_QR
                            checkIndex(index,P,origj2);
#endif
                        }

                        // Find Householder matrix for this column
                        M1c u = A.col(j,j+1,M);
                        HouseholderReflect(A.ref(j,j),u,*bj);
                        //std::cout<<"bj = "<<*bj<<std::endl;

                        // This multiplies through to the end of the original block.
                        // This way, when we are done, the whole block has had the
                        // same Householder reflections applied to it.
                        M1r A2a = A.row(j,j+1,origj2);
                        M1s A2b = A.subMatrix(j+1,M,j+1,origj2);
                        M3c temp = tempBase.col(0,0,origj2-j-1);
                        HouseholderMultEq(u,*bj,A2a,A2b,temp);

                        // Update Z:
                        M1s A3 = A.subMatrix(j1,M,j1,j+1);
                        Zs Z1 = Z.subTriMatrix(0,j-j1+1);
                        BlockHouseholderAugment(A3,Z1,*bj);

                        // Update the colnormsq values within this block
                        // (No need to go all the way to origj2, since the 
                        // j2..origj2 columns are those with low norm already -
                        // we don't need those values until we recalculate them 
                        // from scratch anyway.)
                        for(int k=j+1;k<j2;++k) 
                            colnormsq(k) -= TMV_NORM(A.cref(j,k)*scale);
                        //std::cout<<"colnormsq => "<<colnormsq<<std::endl;
                    }

                    if (j1 < j2) {
                        // Do the Block Householder update of the rest of the 
                        // matrix:
                        //std::cout<<"j1,j2 = "<<j1<<" "<<j2<<std::endl;
                        //std::cout<<"origj2 = "<<origj2<<std::endl;
                        M1s A1 = A.subMatrix(j1,M,j1,j2);
                        //std::cout<<"A1 = "<<A1<<std::endl;
                        Zs Z1 = Z.subTriMatrix(0,j2-j1);
                        //std::cout<<"Z1 = "<<Z1<<std::endl;
                        M1s A4 = A.subMatrix(j1,M,origj2,N);
                        //std::cout<<"A4 = "<<A4<<std::endl;
                        M3s temp = tempBase.subMatrix(0,j2-j1,0,N-origj2);
                        //std::cout<<"temp = "<<temp<<std::endl;
                        BlockHouseholderLDiv(A1,Z1,A4,temp);

                        // Update the colnormsq values for the rest of the matrix:
                        if (M-j2 > j2-j1)
                            for(int k=origj2;k<N;++k) colnormsq(k) -= 
                                A.col(k,j1,j2).normSq(scale);
                        else 
                            for(int k=origj2;k<N;++k) colnormsq(k) = 
                                A.col(k,j2,M).normSq(scale);
                    }

#ifdef XDEBUG_QR
                    checkIndex(index,P,origj2);
#endif
                    // Put the bad columns back where they started before this 
                    // loop:
                    for(int j=j2; j<origj2; ++j) if (P[j] > j2) {
                        TMVAssert(P[j] < int(A.rowsize()));
                        A.swapCols(j,P[j]);
                        colnormsq.swap(j,P[j]);
#ifdef XDEBUG_QR
                        index.swap(j,P[j]);
#endif
                    }
#ifdef XDEBUG_QR
                    checkIndex(index,P,j2);
#endif

                    j1 = j1x = j2;
                }

                if (j3 < N) {
                    //std::cout<<"recalculate colnorms\n";
                    // Then need to recalculate some of the colnorms:
                    for(int k=j3;k<N;++k) 
                        colnormsq(k) = A.col(k,j3,M).normSq(scale);
                    //std::cout<<"colnormsq = "<<colnormsq<<std::endl;
                }
            }

#ifdef PRINTALGO_QR
            //std::cout<<"A -> "<<A<<std::endl;
            //std::cout<<"beta -> "<<beta<<std::endl;
            //QRPDecompose_Helper<11,cs,rs,M1x,Vx>::call(Ac,bc,P,false);
            //std::cout<<"correct A = "<<Ac<<std::endl;
            //std::cout<<"beta = "<<bc<<std::endl;
            //std::cout<<"diff = "<<(A-Ac).copy().clip(1.e-5)<<std::endl;
            //std::cout<<(beta-bc).copy().clip(1.e-5)<<std::endl;
#endif

#ifdef XDEBUG_QR
            std::cout<<"Done LooseBlockQRPDecompose\n";
            std::cout<<"Norm(A) = "<<Norm(A)<<std::endl;
            std::cout<<"Norm(beta) = "<<Norm(beta)<<std::endl;
            checkIndex(index,P,N);
            Matrix<T> Q = A;
            UnpackQ(Q,beta);
            Matrix<T> AA = Q*A.upperTri();
            AA.reversePermuteCols(P);
            if (!(Norm(AA-A0) <= 0.001*TMV_MAX(RT(1),Norm(A0)))) {
                std::cerr<<"LooseBlockQRPDecompose: "<<std::endl;
                std::cerr<<"A = "<<TMV_Text(A)<<std::endl;
                if (N < 100) {
                    std::cerr<<"  "<<A0<<std::endl;
                    std::cerr<<"-> "<<A<<std::endl;
                    std::cerr<<"beta = "<<beta<<std::endl;
                    std::cerr<<"P = ";
                    for(int i=0;i<N;i++) std::cerr<<P[i]<<" ";
                    std::cerr<<std::endl;
                    std::cerr<<"QRP = "<<AA<<std::endl;
                    std::cerr<<"A0 = "<<A0<<std::endl;
                    std::cerr<<"diff = "<<Matrix<T>(AA-A0).clip(1.e-5)<<std::endl;
                }
                std::cerr<<"Norm(diff) = "<<Norm(AA-A0)<<std::endl;
                std::cerr<<"Rdiag = "<<A.diag()<<std::endl;
                abort(); 
            }
#endif
        }
    };

    // algo 31: Decide which algorithm to use from runtime size
    template <int cs, int rs, class M1, class V>
    struct QRPDecompose_Helper<31,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta, int* P, bool strict)
        {
            typedef typename M1::value_type T;

            const int N = rs==TMV_UNKNOWN ? int(A.rowsize()) : rs;
            const int M = cs==TMV_UNKNOWN ? int(A.colsize()) : cs;
#ifdef PRINTALGO_QR
            std::cout<<"QRPDecompose algo 31: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const int l2cache = TMV_L2_CACHE*1024/sizeof(T);

#ifdef PRINTALGO_QR
            std::cout<<"M*N = "<<M*N<<std::endl;
            std::cout<<"l2cache = "<<l2cache<<std::endl;
            std::cout<<"strict = "<<strict<<std::endl;
            std::cout<<"algo = "<<
                (M*N <= l2cache ? 11 : strict ? 21 : 22)<<std::endl;
#endif
            if (M*N <= l2cache)
                QRPDecompose_Helper<11,cs,rs,M1,V>::call(A,beta,P,strict);
            else if (strict)
                QRPDecompose_Helper<21,cs,rs,M1,V>::call(A,beta,P,strict);
            else
                QRPDecompose_Helper<22,cs,rs,M1,V>::call(A,beta,P,strict);
        }
    };

    // algo 32: Call strict or non-strict algorithm according to strict param.
    template <int cs, int rs, class M1, class V>
    struct QRPDecompose_Helper<32,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta, int* P, bool strict)
        {
#ifdef PRINTALGO_QR
            const int M = cs==TMV_UNKNOWN ? int(A.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(A.rowsize()) : rs;
            std::cout<<"QRPDecompose algo 32: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            if (strict)
                QRPDecompose_Helper<21,cs,rs,M1,V>::call(A,beta,P,strict);
            else
                QRPDecompose_Helper<22,cs,rs,M1,V>::call(A,beta,P,strict);
        }
    };

    // algo 81: Copy to colmajor
    template <int cs, int rs, class M, class V>
    struct QRPDecompose_Helper<81,cs,rs,M,V>
    {
        static inline void call(M& m, V& beta, int* P, bool strict)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRPDecompose algo 81: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M::value_type T;
            typedef typename MCopyHelper<T,Rec,cs,rs>::type Mcm;
            Mcm mcm = m;
            QRPDecompose_Helper<-2,cs,rs,Mcm,V>::call(mcm,beta,P,strict);
            NoAliasCopy(mcm,m);
        }
    };

    // algo 90: call InstQRP_Decompose
    template <int cs, int rs, class M, class V>
    struct QRPDecompose_Helper<90,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta, int* P, bool strict)
        { InstQRP_Decompose(m.xView(),beta.xView(),P,strict); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M, class V>
    struct QRPDecompose_Helper<97,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta, int* P, bool strict)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            QRPDecompose_Helper<-2,cs,rs,Mc,V>::call(mc,beta,P,strict);
        }
    };

    // algo -4: No copies or branches
    template <int cs, int rs, class M, class V>
    struct QRPDecompose_Helper<-4,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta, int* P, bool strict)
        {
#if 0
            const int algo = 21;
#else
            typedef typename M::value_type T;
            const int csrs = IntTraits2<cs,rs>::prod;
            const int l2cache = TMV_L2_CACHE*1024/sizeof(T);
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                rs == TMV_UNKNOWN ? 31 :
                cs == TMV_UNKNOWN ? 31 : 
                csrs <= l2cache ? 11 : 32;
#endif
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRPDecompose: \n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<m.colsize()<<"  "<<m.rowsize()<<std::endl;
            //std::cout<<"csrs = "<<csrs<<std::endl;
            //std::cout<<"l2cache = "<<l2cache<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
#endif
            QRPDecompose_Helper<algo,cs,rs,M,V>::call(m,beta,P,strict);
#ifdef PRINTALGO_QR
            //std::cout<<"m => "<<m<<std::endl;
            //std::cout<<"beta => "<<beta<<std::endl;
#endif
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class V>
    struct QRPDecompose_Helper<-3,cs,rs,M1,V>
    {
        static TMV_INLINE void call(M1& m, V& beta, int* P, bool strict)
        {
            const int algo = (
                ( cs != TMV_UNKNOWN && rs != TMV_UNKNOWN &&
                  cs <= 16 && rs <= 16 ) ? -4 :
                ( TMV_OPT >= 2 && !M1::_colmajor ) ? 81 :
                -4 );
#ifdef PRINTALGO_QR
            const int M = cs==TMV_UNKNOWN ? int(m.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m.rowsize()) : rs;
            std::cout<<"QRPDecompose algo -3: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"strict = "<<strict<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            QRPDecompose_Helper<algo,cs,rs,M1,V>::call(m,beta,P,strict);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M, class V>
    struct QRPDecompose_Helper<-2,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta, int* P, bool strict)
        {
            typedef typename M::value_type T;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                M::_conj ? 97 :
                inst ? 90 :
                -3;
            QRPDecompose_Helper<algo,cs,rs,M,V>::call(m,beta,P,strict);
        }
    };

    template <int cs, int rs, class M, class V>
    struct QRPDecompose_Helper<-1,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta, int* P, bool strict)
        { QRPDecompose_Helper<-2,cs,rs,M,V>::call(m,beta,P,strict); }
    };

    template <class M, class V>
    static inline void InlineQRP_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta,
        int* P, bool strict=false)
    {
        const int cs = M::_colsize;
        const int rs = Sizes<M::_rowsize,V::_size>::size;
        typedef typename M::cview_type Mv;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TMV_MAYBE_REF(V,Vv) betav = beta.cView();
        QRPDecompose_Helper<-3,cs,rs,Mv,Vv>::call(mv,betav,P,strict);
    }

    // This is the basic functionality
    template <class M, class V>
    static inline void QRP_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta,
        int* P, bool strict=false)
    {
        TMVStaticAssert(V::isreal);
        TMVStaticAssert((Traits2<
                         typename M::value_type,
                         typename V::value_type>::samebase));
        TMVAssert(m.colsize() >= m.rowsize());
        TMVAssert(beta.size() == m.rowsize());
        TMVStaticAssert((Sizes<M::_rowsize,V::_size>::same));
        TMVAssert(m.rowsize() == beta.size());

        const int cs = M::_colsize;
        const int rs = Sizes<M::_rowsize,V::_size>::size;
        typedef typename M::cview_type Mv;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TMV_MAYBE_REF(V,Vv) betav = beta.cView();
        QRPDecompose_Helper<-2,cs,rs,Mv,Vv>::call(mv,betav,P,strict);
    }

    // This function is a friend of the Permutation class.
    // Note: The strict=false default is in TMV_Permutation.h
    template <class M, class V>
    static inline void QRP_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta,
        Permutation& P, bool strict)
    {
        TMVAssert(P.size() == m.rowsize());
        P.allocateMem();
        QRP_Decompose(m,beta,P.getMem(),strict);
        P.isinv = false;
        P.calcDet();
    }


    // The rest of these below are basically convenience functions
    // to allow the user to provide fewer or different arguments.
    template <class M1, class M2>
    static inline void QRP_Decompose(
        BaseMatrix_Rec_Mutable<M1>& Q, BaseMatrix_Tri_Mutable<M2>& R,
        Permutation& P, bool strict=false)
    {
        TMVStaticAssert((Traits2<
                         typename M1::value_type,
                         typename M2::value_type>::sametype));
        TMVStaticAssert(M2::_upper);
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(R.colsize() == Q.rowsize());

        typedef typename M1::real_type RT;
        Vector<RT> beta(Q.rowsize());
        QRP_Decompose(Q,beta,P,strict);
        NoAliasCopy(Q.upperTri(),R);
        UnpackQ(Q,beta);
    }

    template <class M>
    static inline void QRP_Decompose(BaseMatrix_Rec_Mutable<M>& A,
        bool strict=false)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        typedef typename M::real_type RT;
        Vector<RT> beta(A.rowsize());
        Permutation P(A.rowsize());
        QRP_Decompose(A,beta,P,strict);
    }


    // Allow views as an argument by value (for convenience)
    template <class T, int A, int A2>
    static inline void QRP_Decompose(
        MatrixView<T,A> Q, UpperTriMatrixView<T,A2> R,
        Permutation& P, bool strict=false)
    {
        typedef MatrixView<T,A> M1;
        typedef UpperTriMatrixView<T,A2> M2;
        QRP_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R),P,strict); 
    }

    template <class T, int M, int N, int Si, int Sj, int A, int Si2, int Sj2, int A2>
    static inline void QRP_Decompose(
        SmallMatrixView<T,M,N,Si,Sj,A> Q,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> R,
        Permutation& P, bool strict=false)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        typedef SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> M2;
        QRP_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R),P,strict);
    }

    template <class T, int N, int A, int Si2, int Sj2, int A2>
    static inline void QRP_Decompose(
        MatrixView<T,A> Q,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> R,
        Permutation& P, bool strict=false)
    {
        typedef MatrixView<T,A> M1;
        typedef SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> M2;
        QRP_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R),P,strict);
    }

    template <class T, int M, int N, int Si, int Sj, int A, int A2>
    static inline void QRP_Decompose(
        SmallMatrixView<T,M,N,Si,Sj,A> Q,
        UpperTriMatrixView<T,A2> R,
        Permutation& P, bool strict=false)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        typedef UpperTriMatrixView<T,A2> M2;
        QRP_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R),P,strict);
    }

    template <class T, int A>
    static inline void QRP_Decompose(MatrixView<T,A> m, bool strict=false)
    {
        typedef MatrixView<T,A> M1;
        QRP_Decompose(static_cast<BaseMatrix_Rec_Mutable<M1>&>(m),strict);
    }

    template <class T, int M, int N, int Si, int Sj, int A>
    static inline void QRP_Decompose(
        SmallMatrixView<T,M,N,Si,Sj,A> m, bool strict=false)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        QRP_Decompose(static_cast<BaseMatrix_Rec_Mutable<M1>&>(m),strict);
    }

    // Don't forget the ones that mix *MatrixView with BaseMatrix_*_Mutable
    template <class T, int A, class M2>
    static inline void QRP_Decompose(
        MatrixView<T,A> Q, BaseMatrix_Tri_Mutable<M2>& R,
        Permutation& P, bool strict=false)
    {
        typedef MatrixView<T,A> M1;
        QRP_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),R,P,strict);
    }

    template <class M1, class T, int A2>
    static inline void QRP_Decompose(
        BaseMatrix_Rec_Mutable<M1>& Q, UpperTriMatrixView<T,A2> R,
        Permutation& P, bool strict=false)
    {
        typedef UpperTriMatrixView<T,A2> M2;
        QRP_Decompose(
            Q,static_cast<BaseMatrix_Tri_Mutable<M2>&>(R),P,strict); 
    }

    template <class T, int M, int N, int Si, int Sj, int A, class M2>
    static inline void QRP_Decompose(
        SmallMatrixView<T,M,N,Si,Sj,A> Q, BaseMatrix_Tri_Mutable<M2>& R,
        Permutation& P, bool strict=false)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        QRP_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),R,P,strict);
    }

    template <class M1, class T, int N, int Si2, int Sj2, int A2>
    static inline void QRP_Decompose(
        BaseMatrix_Rec_Mutable<M1>& Q,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> R,
        Permutation& P, bool strict=false)
    {
        typedef SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> M2;
        QRP_Decompose(
            Q,static_cast<BaseMatrix_Tri_Mutable<M2>&>(R),P,strict);
    }

} // namespace tmv

#undef TMV_QR_BLOCKSIZE

#endif

