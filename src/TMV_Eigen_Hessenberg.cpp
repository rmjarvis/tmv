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



#include "TMV_Blas.h"
#include "TMV_Eigen.h"
#include "TMV_Matrix.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define HESS_BLOCKSIZE TMV_BLOCKSIZE
#else
#define HESS_BLOCKSIZE 16
#endif

    //
    // Reduce Matrix to Hessenberg Form (upper tri with one lower sub-diag)
    //

    template <class T> 
    static void NonBlockHessenberg(
        MatrixView<T> A, VectorView<T> Ubeta)
    {
#ifdef XDEBUG
        cout<<"Start NonBlock Hessenberg Reduction: A = "<<A<<endl;
        Matrix<T> A0(A);
#endif
        // Decompose A into U H Ut
        // H is a Hessenberg Matrix
        // U is a Unitary Matrix
        // On output, H is stored in the upper-Hessenberg part of A
        // U is stored in compact form in the rest of A along with 
        // the vector Ubeta.
        const ptrdiff_t N = A.rowsize();

        TMVAssert(A.colsize() == A.rowsize());
        TMVAssert(N > 0);
        TMVAssert(Ubeta.size() == N-1);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(!Ubeta.isconj());
        TMVAssert(Ubeta.step()==1);

        // We use Householder reflections to reduce A to the Hessenberg form:
        T* Uj = Ubeta.ptr();
        T det = 0; // Ignore Householder det calculations
        for(ptrdiff_t j=0;j<N-1;++j,++Uj) {
#ifdef TMVFLDEBUG
            TMVAssert(Uj >= Ubeta._first);
            TMVAssert(Uj < Ubeta._last);
#endif
            *Uj = Householder_Reflect(A.subMatrix(j+1,N,j,N),det);
            if (*Uj != T(0))
                Householder_LMult(A.col(j+2,N),*Uj,A.subMatrix(0,N,j+1,N).adjoint());
        }

#ifdef XDEBUG
        Matrix<T> U(N,N,T(0));
        U.subMatrix(1,N,1,N) = A.subMatrix(1,N,0,N-1);
        U.upperTri().setZero();
        Vector<T> Ubeta2(N);
        Ubeta2.subVector(1,N) = Ubeta;
        Ubeta2(0) = T(0);
        GetQFromQR(U.view(),Ubeta2);
        Matrix<T> H = A;
        if (N>2) LowerTriMatrixViewOf(H).offDiag(2).setZero();
        Matrix<T> AA = U*H*U.adjoint();
        if (Norm(A0-AA) > 0.001*Norm(A0)) {
            cerr<<"NonBlock Hessenberg: A = "<<Type(A)<<"  "<<A0<<endl;
            cerr<<"A = "<<A<<endl;
            cerr<<"Ubeta = "<<Ubeta<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"UHUt = "<<AA<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    static void BlockHessenberg(
        MatrixView<T> A, VectorView<T> Ubeta)
    {
        // Much like the block version of Bidiagonalize, we try to maintain
        // the operation of several successive Householder matrices in
        // a block form, where the net Block Householder is I - YZYt.
        //
        // But as with the bidiagonlization algorithm (and unlike a simple
        // block QR decomposition), we update the matrix from both the left 
        // and the right, so we also need to keep track of the product
        // ZYtm in addition.
        //
        // The block update at the end of the block loop is
        // m' = (I-YZYt) m (I-YZtYt)
        //
        // The Y matrix is stored in the first K columns of m,
        // and the Hessenberg portion of these columns is updated as we go.
        // For the right-hand-side update, m -= mYZtYt, the m on the right
        // needs to be the full original matrix m, including the original
        // versions of these K columns.  Therefore, we can't wait until 
        // the end for this calculation.  
        //
        // Instead, we keep track of mYZt as we progress, so the final update
        // is:
        //
        // m' = (I-YZYt) (m - mYZt Y)
        //
        // We also need to do this same calculation for each column as we
        // progress through the block.
        //
        const ptrdiff_t N = A.rowsize();

#ifdef XDEBUG
        Matrix<T> A0(A);
#endif

        TMVAssert(A.rowsize() == A.colsize());
        TMVAssert(N > 0);
        TMVAssert(Ubeta.size() == N-1);
        TMVAssert(!Ubeta.isconj());
        TMVAssert(Ubeta.step()==1);

        ptrdiff_t ncolmax = MIN(HESS_BLOCKSIZE,N-1);
        Matrix<T,RowMajor> mYZt_full(N,ncolmax);
        UpperTriMatrix<T,NonUnitDiag|ColMajor> Z_full(ncolmax);

        T det(0); // Ignore Householder Determinant calculations
        T* Uj = Ubeta.ptr();
        for(ptrdiff_t j1=0;j1<N-1;) {
            ptrdiff_t j2 = MIN(N-1,j1+HESS_BLOCKSIZE);
            ptrdiff_t ncols = j2-j1;
            MatrixView<T> mYZt = mYZt_full.subMatrix(0,N-j1,0,ncols);
            UpperTriMatrixView<T> Z = Z_full.subTriMatrix(0,ncols);

            for(ptrdiff_t j=j1,jj=0;j<j2;++j,++jj,++Uj) { // jj = j-j1

                // Update current column of A
                //
                // m' = (I - YZYt) (m - mYZt Yt)
                // A(0:N,j)' = A(0:N,j) - mYZt(0:N,0:j) Y(j,0:j)t
                A.col(j,j1+1,N) -= mYZt.Cols(0,j) * A.row(j,0,j).Conjugate();
                //
                // A(0:N,j)'' = A(0:N,j) - Y Z Yt A(0:N,j)'
                // 
                // Let Y = (L)     where L is unit-diagonal, lower-triangular,
                //         (M)     and M is rectangular
                //
                LowerTriMatrixView<T> L = 
                    LowerTriMatrixViewOf(A.subMatrix(j1+1,j+1,j1,j),UnitDiag);
                MatrixView<T> M = A.subMatrix(j+1,N,j1,j);
                // Use the last column of Z as temporary storage for Yt A(0:N,j)'
                VectorView<T> YtAj = Z.col(jj,0,jj);
                YtAj = L.adjoint() * A.col(j,j1+1,j+1);
                YtAj += M.adjoint() * A.col(j,j+1,N);
                YtAj = Z.subTriMatrix(0,jj) * YtAj;
                A.col(j,j1+1,j+1) -= L * YtAj;
                A.col(j,j+1,N) -= M * YtAj;

                // Do the Householder reflection 
                VectorView<T> u = A.col(j,j+1,N);
                T bu = Householder_Reflect(u,det);
#ifdef TMVFLDEBUG
                TMVAssert(Uj >= Ubeta._first);
                TMVAssert(Uj < Ubeta._last);
#endif
                *Uj = bu;

                // Save the top of the u vector, which isn't actually part of u
                T& Atemp = *u.cptr();
                TMVAssert(IMAG(Atemp) == RealType(T)(0));
                RealType(T) Aorig = REAL(Atemp);
                Atemp = RealType(T)(1);

                // Update Z
                VectorView<T> Zj = Z.col(jj,0,jj);
                Zj = -bu * M.adjoint() * u;
                Zj = Z * Zj;
                Z(jj,jj) = -bu;

                // Update mYtZt:
                //
                // mYZt(0:N,j) = m(0:N,0:N) Y(0:N,0:j) Zt(0:j,j)
                //             = m(0:N,j+1:N) Y(j+1:N,j) Zt(j,j)
                //             = bu* m(0:N,j+1:N) u 
                //
                mYZt.col(jj) = CONJ(bu) * A.subMatrix(j1,N,j+1,N) * u;

                // Restore Aorig, which is actually part of the Hessenberg matrix.
                Atemp = Aorig;
            }

            // Update the rest of the matrix:
            // A(j2,j2-1) needs to be temporarily changed to 1 for use in Y
            T& Atemp = *(A.ptr() + j2*A.stepi() + (j2-1)*A.stepj());
            TMVAssert(IMAG(Atemp) == RealType(T)(0));
            RealType(T) Aorig = Atemp;
            Atemp = RealType(T)(1);

            // m' = (I-YZYt) (m - mYZt Y)
            MatrixView<T> m = A.subMatrix(j1,N,j2,N);
            ConstMatrixView<T> Y = A.subMatrix(j2+1,N,j1,j2);
            m -= mYZt * Y.adjoint();
            BlockHouseholder_LMult(Y,Z,m);

            // Restore A(j2,j2-1)
            Atemp = Aorig;
            j1 = j2;
        }

#ifdef XDEBUG
        Matrix<T> U(N,N,T(0));
        U.subMatrix(1,N,1,N) = A.subMatrix(1,N,0,N-1);
        U.upperTri().setZero();
        U(0,0) = T(1);
        Vector<T> Ubeta2(N);
        Ubeta2.subVector(1,N) = Ubeta;
        Ubeta2(0) = T(0);
        GetQFromQR(U.view(),Ubeta2);
        Matrix<T> H = A;
        if (N>2) LowerTriMatrixViewOf(H).offDiag(2).setZero();
        Matrix<T> AA = U*H*U.adjoint();
        if (Norm(A0-AA) > 0.001*Norm(A0)) {
            cerr<<"NonBlock Hessenberg: A = "<<Type(A)<<"  "<<A0<<endl;
            cerr<<"A = "<<A<<endl;
            cerr<<"Ubeta = "<<Ubeta<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"UHUt = "<<AA<<endl;
            Matrix<T,ColMajor> A2 = A0;
            Vector<T> Ub2(Ubeta.size());
            NonBlockHessenberg(A2.view(),Ub2.view());
            cerr<<"cf NonBlock: A -> "<<A2<<endl;
            cerr<<"Ubeta = "<<Ub2<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    static inline void NonLapHessenberg(
        MatrixView<T> A, VectorView<T> Ubeta)
    {
        TMVAssert(A.rowsize() == A.colsize());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(Ubeta.size() == A.rowsize()-1);

#if 0
        if (A.rowsize() > HESS_BLOCKSIZE)
            BlockHessenberg(A,Ubeta,Vbeta,D,E,det);
        else
#endif
            NonBlockHessenberg(A,Ubeta);
    }

#ifdef LAP
    template <class T> 
    static inline void LapHessenberg(
        MatrixView<T> A, VectorView<T> Ubeta)
    { NonLapHessenberg(A,Ubeta,Vbeta,D,E,det); }
#ifdef INST_DOUBLE
    template <> void LapHessenberg(
        MatrixView<double> A, VectorView<double> Ubeta)
    {
        TMVAssert(A.iscm());
        TMVAssert(A.colsize() == A.rowsize());
        TMVAssert(Ubeta.size() == A.rowsize()-1);
        TMVAssert(A.ct()==NonConj);

        int n = A.rowsize();
        int ilo = 1;
        int ihi = n;
        int lda = A.stepj();
        int Lap_info=0;
#ifndef LAPNOWORK
        int lwork = n*LAP_BLOCKSIZE;
        double* work = LAP_DWork(lwork);
#endif
        LAPNAME(dgehrd) (
            LAPCM LAPV(n),LAPV(ilo),LAPV(ihi),
            LAPP(A.ptr()),LAPV(lda),LAPP(Ubeta.ptr())
            LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"dgehrd");
#else
        LAP_Results(Lap_info,int(work[0]),m,n,lwork,"dgehrd");
#endif
    }
    template <> void LapHessenberg(
        MatrixView<std::complex<double> > A, 
        VectorView<std::complex<double> > Ubeta)
    {
        TMVAssert(A.iscm());
        TMVAssert(A.colsize() == A.rowsize());
        TMVAssert(Ubeta.size() == A.rowsize()-1);
        TMVAssert(A.ct()==NonConj);

        int n = A.rowsize();
        int ilo = 1;
        int ihi = n;
        int lda = A.stepj();
        int Lap_info=0;
#ifndef LAPNOWORK
        int lwork = n*LAP_BLOCKSIZE;
        std::complex<double>* work = LAP_ZWork(lwork);
#endif
        LAPNAME(zgehrd) (
            LAPCM LAPV(n),LAPV(ilo),LAPV(ihi),
            LAPP(A.ptr()),LAPV(lda),LAPP(Ubeta.ptr())
            LAPWK(work) LAPVWK(lwork) LAPINFO);
        Ubeta.ConjugateSelf();
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"zgehrd");
#else
        LAP_Results(Lap_info,int(REAL(work[0])),m,n,lwork,"zgehrd");
#endif
    }
#endif
#ifdef INST_FLOAT
    template <> void LapHessenberg(
        MatrixView<float> A, VectorView<float> Ubeta)
    {
        TMVAssert(A.iscm());
        TMVAssert(A.colsize() == A.rowsize());
        TMVAssert(Ubeta.size() == A.rowsize()-1);
        TMVAssert(A.ct()==NonConj);

        int n = A.rowsize();
        int ilo = 1;
        int ihi = n;
        int lda = A.stepj();
        int Lap_info=0;
#ifndef LAPNOWORK
        int lwork = n*LAP_BLOCKSIZE;
        float* work = LAP_SWork(lwork);
#endif
        LAPNAME(sgebrd) (
            LAPCM LAPV(n),LAPV(ilo),LAPV(ihi),
            LAPP(A.ptr()),LAPV(lda),LAPP(Ubeta.ptr())
            LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"sgehrd");
#else
        LAP_Results(Lap_info,int(work[0]),m,n,lwork,"sgehrd");
#endif
    }
    template <> void LapHessenberg(
        MatrixView<std::complex<float> > A, 
        VectorView<std::complex<float> > Ubeta) 
    {
        TMVAssert(A.iscm());
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(Ubeta.size() == A.rowsize());
        TMVAssert(A.ct()==NonConj);

        int n = A.rowsize();
        int ilo = 1;
        int ihi = n;
        int lda = A.stepj();
        int Lap_info=0;
#ifndef LAPNOWORK
        int lwork = n*LAP_BLOCKSIZE;
        std::complex<float>* work = LAP_CWork(lwork);
#endif
        LAPNAME(cgehrd) (
            LAPCM LAPV(n),LAPV(ilo),LAPV(ihi),
            LAPP(A.ptr()),LAPV(lda),LAPP(Ubeta.ptr())
            LAPWK(work) LAPVWK(lwork) LAPINFO);
        Ubeta.ConjugateSelf();
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"cgehrd");
#else
        LAP_Results(Lap_info,int(REAL(work[0])),m,n,lwork,"cgehrd");
#endif
    }
#endif 
#endif // LAP

    template <class T> 
    static inline void Hessenberg(
        MatrixView<T> A, VectorView<T> Ubeta)
    {
        TMVAssert(A.colsize() == A.rowsize());
        TMVAssert(Ubeta.size() == A.rowsize()-1);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(Ubeta.step() == 1);

        if (A.rowsize() > 0) {
#ifdef LAP
            if (A.iscm()) 
                LapHessenberg(A,Ubeta);
            else 
#endif
                NonLapHessenberg(A,Ubeta);
        }
    }

    //
    // Schur_From_Hessenberg: QR method
    //

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_Eigen.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


