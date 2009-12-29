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


#include "TMV_Blas.h"
#include "TMV_SVDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "TMV_Householder.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "TMV_QRDiv.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif


namespace tmv {

#define RT TMV_RealType(T)

#ifdef TMV_BLOCKSIZE
#define BIDIAG_BLOCKSIZE TMV_BLOCKSIZE/4
#else
#define BIDIAG_BLOCKSIZE 16
#endif

    //
    // Bidiagonalize
    //

    template <class T> 
    static void NonBlockBidiagonalize(
        const MatrixView<T>& A, const VectorView<T>& Ubeta,
        const VectorView<T>& Vbeta, const VectorView<RT>& D,
        const VectorView<RT>& E, T& signdet)
    {
#ifdef XDEBUG
        //cout<<"Start NonBlockBidiag: A = "<<A<<endl;
        Matrix<T> A0(A);
#endif
        // Decompose A into U B V
        // The Bidiagonal Matrix B is stored as two vectors: D, E
        // D is the diagonal, E is the super-diagonal
        // A along with Ubeta and Vbeta hold the U and V matrices.
        const int M = A.colsize();
        const int N = A.rowsize();

        TMVAssert(N <= M);
        TMVAssert(N > 0);
        TMVAssert(int(Ubeta.size()) == N);
        TMVAssert(int(Vbeta.size()) == N-1);
        TMVAssert(int(D.size()) == N);
        TMVAssert(int(E.size()) == N-1);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(!Ubeta.isconj());
        TMVAssert(!Vbeta.isconj());

        // We use Householder reflections to reduce A to the bidiagonal form:
        T* Uj = Ubeta.ptr();
        T* Vj = Vbeta.ptr();
        //cout<<"Start Bidiag\n";
        //cout<<"A = "<<A<<endl;
        for(int j=0;j<N-1;++j,++Uj,++Vj) {
#ifdef TMVFLDEBUG
            TMVAssert(Uj >= Ubeta.first);
            TMVAssert(Uj < Ubeta.last);
            TMVAssert(Vj >= Vbeta.first);
            TMVAssert(Vj < Vbeta.last);
#endif
            //cout<<"j = "<<j<<endl;
            *Uj = HouseholderReflect(A.subMatrix(j,M,j,N),signdet);
            //cout<<"U reflect: A->"<<A<<endl;
            *Vj = HouseholderReflect(A.transpose().subMatrix(j+1,N,j,M),signdet);
            //cout<<"V reflect: A->"<<A<<endl;
        }
#ifdef TMVFLDEBUG
        TMVAssert(Uj >= Ubeta.first);
        TMVAssert(Uj < Ubeta.last);
#endif
        //cout<<"j = "<<N-1<<endl;
        *Uj = HouseholderReflect(A.col(N-1,N-1,M),signdet);
        //cout<<"U reflect: A->"<<A<<endl;

        // The bidiagonal of A is the bidiagonal we want, so copy it to D,E
#ifdef XTEST
        if (isComplex(T())) {
            TMVAssert(normInf(A.diag().imag()) == RT(0));
            TMVAssert(normInf(A.diag(1).imag()) == RT(0));
        }
#endif
        D = A.diag().real();
        E = A.diag(1).real();

#ifdef XDEBUG
        Matrix<T> U(A);
        GetQFromQR(U.view(),Ubeta);
        Matrix<T> V(N,N);
        V.setToIdentity();
        V.subMatrix(1,N,1,N) = A.subMatrix(0,N-1,1,N);
        GetQFromQR(V.subMatrix(1,N,1,N).transpose(),Vbeta);
        Matrix<RT> B(N,N,RT(0));
        B.diag() = D;
        B.diag(1) = E;
        Matrix<T> AA = U*B*V;
        if (Norm(A0-AA) > 0.001*Norm(A0)) {
            cerr<<"Bidiagonalize: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"A = "<<A<<endl;
            cerr<<"Ubeta = "<<Ubeta<<endl;
            cerr<<"Vbeta = "<<Vbeta<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"B = "<<B<<endl;
            cerr<<"V = "<<V<<endl;
            cerr<<"UBV = "<<AA<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    static void BlockBidiagonalize(
        const MatrixView<T>& A, const VectorView<T>& Ubeta,
        const VectorView<T>& Vbeta, const VectorView<RT>& D,
        const VectorView<RT>& E, T& signdet)
    {
        // Normally we keep the Z matrix for block Householder matrices where 
        // the block Householder is I - YZYt (and Z is upper triangular).
        //
        // However, since the bidiagonalizing process proceeds along both 
        // the rows and columns, we have a problem.  If we just kept the Z
        // matrix for each set, then we would not be able to update the 
        // resulting submatrix at the end of the block.
        // m' = (I-YZYt) m (I-XtWX)
        //
        // For the left multiply, the m needs to have full height for the 
        // Yt m product.  Likewise, for the right multiply, it needs full width.
        // Since we update the first K rows and columns with the Y matrix, 
        // this doesn't work.  So instead of keeping Z,W we are forced to use
        // a bit more temporary storage and store the products ZYtm and mXtW.
        //
        // Furthermore, the m in these products is maintained such that the
        // it already has the appropriate multiplies from the other side.
        // Then, when we are done with the block, the update becomes just:
        //
        // m' = m' - Y (ZYtm) - (mXtW) X
        //
        const int M = A.colsize();
        const int N = A.rowsize();

#ifdef XDEBUG
        Matrix<T> A0(A);
#endif

        TMVAssert(N <= M);
        TMVAssert(N > 0);
        TMVAssert(int(Ubeta.size()) == N);
        TMVAssert(int(Vbeta.size()) == N-1);
        TMVAssert(int(D.size()) == N);
        TMVAssert(int(E.size()) == N-1);
        TMVAssert(!Ubeta.isconj());
        TMVAssert(!Vbeta.isconj());
        TMVAssert(Ubeta.step()==1);
        TMVAssert(Vbeta.step()==1);
        TMVAssert(E.step()==1);

        Matrix<T,RowMajor> ZYtm(TMV_MIN(BIDIAG_BLOCKSIZE,N-1),N);
        Matrix<T,ColMajor> mXtW(M,TMV_MIN(BIDIAG_BLOCKSIZE,N-1));

        T* Uj = Ubeta.ptr();
        T* Vj = Vbeta.ptr();
        RT* Dj = D.ptr();
        const int Ds = D.step();
        RT* Ej = E.ptr();
        for(int j1=0;j1<N-1;) {
            int j2 = TMV_MIN(N-1,j1+BIDIAG_BLOCKSIZE);
            for(int j=j1,jj=0;j<j2;++j,++jj,++Uj,++Vj,Dj+=Ds,++Ej) { // jj = j-j1

                // Update current column:
                // A(j:M,j) -= Y(j:M,0:j) ZYtm(0:j,j) + mXtW(j:M,0:j) X(0:j,j)
                //
                VectorView<T> u = A.col(j,j,M);
                MatrixView<T> Y0 = A.subMatrix(j,M,j1,j);
                MatrixView<T> ZYtm0 = ZYtm.subMatrix(0,jj,j,N);
                MatrixView<T> mXtW0 = mXtW.subMatrix(j,M,0,jj);
                MatrixView<T> X0 = A.subMatrix(j1,j,j,N);
                if (jj > 0) {
                    u -= Y0 * ZYtm0.col(0);
                    u -= mXtW0 * X0.col(0);
                }

                // Do the Householder reflection for U
                // Copy the reflection into D(j), and set the top of the 
                // Householder vector to be explicitly 1.  (It makes life easier
                // if it's actually 1 rather than dealing with it implicitly.)
                //
                T bu = HouseholderReflect(u,signdet);
#ifdef TMVFLDEBUG
                TMVAssert(Uj >= Ubeta.first);
                TMVAssert(Uj < Ubeta.last);
                TMVAssert(Dj >= D.first);
                TMVAssert(Dj < D.last);
                TMVAssert(u.ptr() >= A.first);
                TMVAssert(u.ptr() < A.last);
#endif
                *Uj = bu;
                TMVAssert(TMV_IMAG(*u.cptr()) == RT(0)); 
                *Dj = TMV_REAL(*u.cptr());
                *u.ptr() = T(1);

                // Update ZYtm:
                //
                // ZYtm(j,j+1:N) = Z(j,0:j+1) Yt(0:j+1,0:M) m'(0:M,j+1:N)
                // m' is the matrix after taking into account the Householder 
                // multiplies that we have already done:
                //
                // m' = m' - Y ZYtm - mXtW X
                //
                // The new ZYtm(j,j+1:N) = (Z Yt m')(j,j+1:N)
                // = Z(j,0:j+1) Yt(0:j+1,0:M) m'(0:M,j+1:N)
                //
                // Z is Upper Triangular, so Z(j,0:j+1) = 0
                // Also Z(j,j) = bu, so:
                // ZYtm(j,j+1:N) = bu Yt(j,0:M) m'(0:M,j+1:N)
                //
                // Y is Lower Unit Trapezoidal, so Yt(j,0:j) = 0
                // The rest, Yt(j,j:M) is just ut, so:
                // ZYtm(j,j+1:N) = bu ut m'(j:M,j+1:N)
                //
                // Finally, expand out the m':
                // m'(j:M,j+1:N) = A(j:M,j+1:N) 
                //                 - Y(j:M,0:j) ZYtm(0:j,j+1:N)
                //                 - mXtW(j:M,0:j) X(0:j,j+1:N)
                //
                VectorView<T> ZYtmj = ZYtm.row(jj,j+1,N);
                VectorView<T> temp = ZYtm.row(jj,j1,j);
                VectorView<T> ut = u.conjugate();
                MatrixView<T> ZYtm0a = ZYtm0.colRange(1,ZYtm0.rowsize());
                MatrixView<T> X0a = X0.colRange(1,X0.rowsize());
                ZYtmj = ut * A.subMatrix(j,M,j+1,N);
                ZYtmj -= (temp = ut*Y0) * ZYtm0a;
                ZYtmj -= (temp = ut*mXtW0) * X0a;
                ZYtmj *= bu;

                // Update the current row:
                // A(j,j+1:N) -= Y(j,0:j+1) ZYtm(0:j+1,j+1:N) + mXtW(j,0:j) X(0:j:j+1,N)
                //
                MatrixView<T> Y1 = A.subMatrix(j,M,j1,j+1);
                MatrixView<T> ZYtm1 = ZYtm.subMatrix(0,jj+1,j+1,N);
                VectorView<T> v = A.row(j,j+1,N);
                v -= Y1.row(0) * ZYtm1;
                v -= mXtW0.row(0) * X0a;

                // Do the Householder reflection for V
                //
                T bv = HouseholderReflect(v,signdet);
#ifdef TMVFLDEBUG
                TMVAssert(Vj >= Vbeta.first);
                TMVAssert(Vj < Vbeta.last);
                TMVAssert(Ej >= E.first);
                TMVAssert(Ej < E.last);
                TMVAssert(v.ptr() >= A.first);
                TMVAssert(v.ptr() < A.last);
#endif
                *Vj = bv;
                TMVAssert(TMV_IMAG(*v.cptr()) == RT(0));
                *Ej = TMV_REAL(*v.cptr());
                *v.ptr() = T(1);

                // Update mXtW:
                //
                // mXtW(j+1:M,j) = m'(j+1:M,0:N) Xt(0:N,0:j+1) W(0:j+1,j)
                // = bv m'(j+1:M,j+1:N) vt
                //
                // And m' is:
                //
                // m'(j+1:M,j+1:N) = A(j+1:M,j+1:N) 
                //                   - Y(j+1:M,0:j+1) ZYtm(0:j+1,j+1:N) 
                //                   - mXtW(j+1:M,0:j) X(0:j,j+1:N)
                //
                VectorView<T> mXtWj = mXtW.col(jj,j+1,M);
                VectorView<T> temp1 = mXtW.col(jj,j1,j+1);
                VectorView<T> temp2 = mXtW.col(jj,j1,j);
                VectorView<T> vt = v.conjugate();
                MatrixView<T> Y1a = Y1.rowRange(1,Y1.colsize());
                MatrixView<T> mXtW0a = mXtW0.rowRange(1,mXtW0.colsize());
                mXtWj = A.subMatrix(j+1,M,j+1,N)*vt;
                mXtWj -= Y1a * (temp1 = ZYtm1*vt);
                mXtWj -= mXtW0a * (temp2 = X0a*vt);
                mXtWj *= bv;
            }

            // Update the rest of the matrix:
            A.subMatrix(j2,M,j2,N) -= A.subMatrix(j2,M,j1,j2) *
                ZYtm.subMatrix(0,j2-j1,j2,N);
            A.subMatrix(j2,M,j2,N) -= mXtW.subMatrix(j2,M,0,j2-j1) *
                A.subMatrix(j1,j2,j2,N);
            j1 = j2;
        }

        // Do the last U Householder vector:
#ifdef TMVFLDEBUG
        TMVAssert(Uj >= Ubeta.first);
        TMVAssert(Uj < Ubeta.last);
        TMVAssert(Dj >= D.first);
        TMVAssert(Dj < D.last);
#endif
        *Uj = HouseholderReflect(A.col(N-1,N-1,M),signdet);
        const T Aend = *(A.cptr()+(N-1)*(A.stepi()+A.stepj()));
        TMVAssert(TMV_IMAG(Aend) == RT(0));
        *Dj = TMV_REAL(Aend);

#ifdef XDEBUG
        Matrix<T> U(A);
        GetQFromQR(U.view(),Ubeta);
        Matrix<T> V(N,N);
        V.setToIdentity();
        V.subMatrix(1,N,1,N) = A.subMatrix(0,N-1,1,N);
        GetQFromQR(V.subMatrix(1,N,1,N).transpose(),Vbeta);
        Matrix<RT> B(N,N,RT(0));
        B.diag() = D;
        B.diag(1) = E;
        Matrix<T> AA = U*B*V;
        if (Norm(A0-AA) > 0.001*Norm(A0)) {
            cerr<<"Bidiagonalize: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"A = "<<A<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"B = "<<B<<endl;
            cerr<<"V = "<<V<<endl;
            cerr<<"UBV = "<<AA<<endl;
            Matrix<T,ColMajor> A2 = A0;
            Vector<T> Ub2(Ubeta.size());
            Vector<T> Vb2(Vbeta.size());
            Vector<RT> D2(D.size());
            Vector<RT> E2(E.size());
            T signdet2(0);
            NonBlockBidiagonalize(
                A2.view(),Ub2.view(),Vb2.view(),D2.view(),E2.view(),signdet2);
            cerr<<"NonBlock: "<<A2<<endl;
            cerr<<"Ubeta = "<<Ubeta<<endl;
            cerr<<"Nonblock: "<<Ub2<<endl;
            cerr<<"Vbeta = "<<Vbeta<<endl;
            cerr<<"Nonblock: "<<Vb2<<endl;
            cerr<<"D = "<<D<<endl;
            cerr<<"D2 = "<<D2<<endl;
            cerr<<"E = "<<E<<endl;
            cerr<<"E2 = "<<E2<<endl;

            abort();
        }
#endif
    }

    template <class T> 
    static inline void NonLapBidiagonalize(
        const MatrixView<T>& A, const VectorView<T>& Ubeta,
        const VectorView<T>& Vbeta, const VectorView<RT>& D,
        const VectorView<RT>& E, T& signdet)
    {
        TMVAssert(A.rowsize() <= A.colsize());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(Ubeta.size() == A.rowsize());
        TMVAssert(Vbeta.size() == A.rowsize()-1);
        TMVAssert(D.size() == A.rowsize());
        TMVAssert(E.size() == A.rowsize()-1);

        if (A.rowsize() > BIDIAG_BLOCKSIZE)
            BlockBidiagonalize(A,Ubeta,Vbeta,D,E,signdet);
        else
            NonBlockBidiagonalize(A,Ubeta,Vbeta,D,E,signdet);
    }

#ifdef LAP
    template <class T> 
    static inline void LapBidiagonalize(
        const MatrixView<T>& A, const VectorView<T>& Ubeta,
        const VectorView<T>& Vbeta, const VectorView<RT>& D,
        const VectorView<RT>& E, T& signdet)
    { NonLapBidiagonalize(A,Ubeta,Vbeta,D,E,signdet); }
#ifdef INST_DOUBLE
    template <> 
    void LapBidiagonalize(
        const MatrixView<double>& A, const VectorView<double>& Ubeta,
        const VectorView<double>& Vbeta, const VectorView<double>& D,
        const VectorView<double>& E, double& signdet)
    {
        TMVAssert(A.iscm());
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(Ubeta.size() == A.rowsize());
        TMVAssert(Vbeta.size() == A.rowsize()-1);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int ldu = A.stepj();
        Vector<double> Vbeta2(n);  
        // Stupid LAPACK requires an extra element in the Vbeta vector
        // which it sets to 0 (!!!) rather than ignores.
        // So we need to create a temporary Vector which is size n.
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = (m+n)*LAP_BLOCKSIZE;
        auto_array<double> work(new double[lwork]);
#else
        int lwork = -1;
        auto_array<double> work(new double[1]);
        LAPNAME(dgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new double[lwork]);
#endif
#endif
        LAPNAME(dgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        Vbeta = Vbeta2.subVector(0,n-1);
#ifdef LAPNOWORK
        LAP_Results("dgebrd");
#else
        LAP_Results(int(work[0]),m,n,lwork,"dgebrd");
#endif
        if (signdet) {
            const double* Ubi = Ubeta.cptr();
            for(int i=0;i<n;++i,++Ubi) if (*Ubi != 0.) 
                signdet = -signdet;
            const double* Vbi = Vbeta.cptr();
            for(int i=0;i<n-1;++i,++Vbi) if (*Vbi != 0.) 
                signdet = -signdet;
        }
    }
    template <> 
    void LapBidiagonalize(
        const MatrixView<std::complex<double> >& A, 
        const VectorView<std::complex<double> >& Ubeta, 
        const VectorView<std::complex<double> >& Vbeta,
        const VectorView<double>& D, const VectorView<double>& E,
        std::complex<double>& signdet)
    {
        TMVAssert(A.iscm());
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(Ubeta.size() == A.rowsize());
        TMVAssert(Vbeta.size() == A.rowsize()-1);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        int m = A.colsize();
        int n = A.rowsize();
        int ldu = A.stepj();
        Vector<std::complex<double> > Vbeta2(n);
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = (m+n)*LAP_BLOCKSIZE;
        auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#else
        int lwork = -1;
        auto_array<std::complex<double> > work(new std::complex<double>[1]);
        LAPNAME(zgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<double>[lwork]);
#endif
#endif
        LAPNAME(zgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        Ubeta.conjugateSelf();
        Vbeta = Vbeta2.subVector(0,n-1);
#ifdef LAPNOWORK
        LAP_Results("zgebrd");
#else
        LAP_Results(int(std::real(work[0])),m,n,lwork,"zgebrd");
#endif
        if (signdet!=0.) {
            const std::complex<double>* Ubi = Ubeta.cptr();
            for(int i=0;i<n;++i,++Ubi) if (*Ubi != 0.) {
                signdet *= TMV_CONJ((*Ubi)*(*Ubi))/norm(*Ubi);
            }
            const std::complex<double>* Vbi = Vbeta.cptr();
            for(int i=0;i<n-1;++i,++Vbi) if (*Vbi != 0.) {
                signdet *= -TMV_CONJ((*Vbi)*(*Vbi))/norm(*Vbi);
            }
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapBidiagonalize(
        const MatrixView<float>& A, const VectorView<float>& Ubeta,
        const VectorView<float>& Vbeta, const VectorView<float>& D,
        const VectorView<float>& E, float& signdet)
    {
        TMVAssert(A.iscm());
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(Ubeta.size() == A.rowsize());
        TMVAssert(Vbeta.size() == A.rowsize()-1);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int ldu = A.stepj();
        Vector<float> Vbeta2(n);  
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = (m+n)*LAP_BLOCKSIZE;
        auto_array<float> work(new float[lwork]);
#else
        int lwork = -1;
        auto_array<float> work(new float[1]);
        LAPNAME(sgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new float[lwork]);
#endif
#endif
        LAPNAME(sgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        Vbeta = Vbeta2.subVector(0,n-1);
#ifdef LAPNOWORK
        LAP_Results("sgebrd");
#else
        LAP_Results(int(work[0]),m,n,lwork,"sgebrd");
#endif
        if (signdet) {
            const float* Ubi = Ubeta.cptr();
            for(int i=0;i<n;++i,++Ubi) if (*Ubi != 0.F) 
                signdet = -signdet;
            const float* Vbi = Vbeta.cptr();
            for(int i=0;i<n-1;++i,++Vbi) if (*Vbi != 0.F) 
                signdet = -signdet;
        }
    }
    template <> 
    void LapBidiagonalize(
        const MatrixView<std::complex<float> >& A, 
        const VectorView<std::complex<float> >& Ubeta, 
        const VectorView<std::complex<float> >& Vbeta,
        const VectorView<float>& D, const VectorView<float>& E,
        std::complex<float>& signdet)
    {
        TMVAssert(A.iscm());
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(Ubeta.size() == A.rowsize());
        TMVAssert(Vbeta.size() == A.rowsize()-1);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        int m = A.colsize();
        int n = A.rowsize();
        int ldu = A.stepj();
        Vector<std::complex<float> > Vbeta2(n);
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = (m+n)*LAP_BLOCKSIZE;
        auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#else
        int lwork = -1;
        auto_array<std::complex<float> > work(new std::complex<float>[1]);
        LAPNAME(cgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<float>[lwork]);
#endif
#endif
        LAPNAME(cgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        Ubeta.conjugateSelf();
        Vbeta = Vbeta2.subVector(0,n-1);
#ifdef LAPNOWORK
        LAP_Results("cgebrd");
#else
        LAP_Results(int(std::real(work[0])),m,n,lwork,"cgebrd");
#endif
        if (signdet!=0.F) {
            const std::complex<float>* Ubi = Ubeta.cptr();
            for(int i=0;i<n;++i,++Ubi) if (*Ubi != 0.F) {
                signdet *= TMV_CONJ((*Ubi)*(*Ubi))/norm(*Ubi);
            }
            const std::complex<float>* Vbi = Vbeta.cptr();
            for(int i=0;i<n-1;++i,++Vbi) if (*Vbi != 0.F) {
                signdet *= -TMV_CONJ((*Vbi)*(*Vbi))/norm(*Vbi);
            }
        }
    }
#endif 
#endif // LAP

    template <class T> 
    void Bidiagonalize(
        const MatrixView<T>& A, const VectorView<T>& Ubeta,
        const VectorView<T>& Vbeta, const VectorView<RT>& D,
        const VectorView<RT>& E, T& signdet)
    {
#ifdef XDEBUG
        Matrix<T> A0(A);
#endif
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == D.size());
        TMVAssert(Ubeta.size() == A.rowsize());
        TMVAssert(Vbeta.size() == A.rowsize()-1);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        TMVAssert(Ubeta.step() == 1);
        TMVAssert(Vbeta.step() == 1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);

        if (A.rowsize() > 0) {
#ifdef LAP
            if (A.iscm()) 
                LapBidiagonalize(A,Ubeta,Vbeta,D,E,signdet);
            else 
#endif
                NonLapBidiagonalize(A,Ubeta,Vbeta,D,E,signdet);
        }
#ifdef XDEBUG
        int N = D.size();
        Matrix<T> U(A);
        GetQFromQR(U.view(),Ubeta);
        Matrix<T> V(N,N);
        V.setToIdentity();
        V.subMatrix(1,N,1,N) = A.subMatrix(0,N-1,1,N);
        GetQFromQR(V.subMatrix(1,N,1,N).transpose(),Vbeta);
        Matrix<RT> B(N,N,RT(0));
        B.diag() = D;
        B.diag(1) = E;
        Matrix<T> AA = U*B*V;
        //cout<<"SVBidiag: Norm(A0-AA) = "<<Norm(A0-AA)<<std::endl;
        //cout<<"cf "<<0.001*Norm(A0)<<std::endl;
        if (Norm(A0-AA) > 0.001*Norm(A0)) {
            cerr<<"Bidiagonalize: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"A = "<<A<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"B = "<<B<<endl;
            cerr<<"V = "<<V<<endl;
            cerr<<"UBV = "<<AA<<endl;
#ifdef LAP
            Matrix<T,ColMajor> A2 = A0;
            Vector<T> Ub2(Ubeta.size());
            Vector<T> Vb2(Vbeta.size());
            Vector<RT> D2(D.size());
            Vector<RT> E2(E.size());
            T signdet2(0);
            NonLapBidiagonalize(
                A2.view(),Ub2.view(),Vb2.view(),D2.view(),E2.view(),signdet2);
            cerr<<"NonLap: "<<A2<<endl;
            cerr<<"Ubeta = "<<Ubeta<<endl;
            cerr<<"NonLap: "<<Ub2<<endl;
            cerr<<"Vbeta = "<<Vbeta<<endl;
            cerr<<"NonLap: "<<Vb2<<endl;
            cerr<<"D = "<<D<<endl;
            cerr<<"D2 = "<<D2<<endl;
            cerr<<"E = "<<E<<endl;
            cerr<<"E2 = "<<E2<<endl;
#endif
            abort();
        }
#endif
    }

#undef RT

#define InstFile "TMV_SVDecompose_Bidiag.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


