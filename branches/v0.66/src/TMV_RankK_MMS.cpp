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
#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#ifdef BLAS
#include "tmv/TMV_SymMatrixArith.h"
#endif

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_SymMatrixArith.h"
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
    // RankKUpdate
    //

    template <bool ha, bool a1, bool add, class T, class Tx> 
    static void RecursiveRankKUpdate(
        const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(alpha != T(0));
        TMVAssert(TMV_IMAG(alpha)==TMV_RealType(T)(0) || !A.isherm());
        TMVAssert(x.colsize() > 0);
        TMVAssert(x.rowsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(ha == A.isherm());
        TMVAssert(a1 == (alpha == T(1)));

        int N = A.size();

        if (N <= SYM_RK_BLOCKSIZE2) {
            if (N == 1) {
                if (ha) {
                    TMV_RealType(T) temp = x.row(0).normSq();
                    TMVAssert(TMV_IMAG(alpha) == TMV_RealType(T)(0));
                    if (!a1) temp *= TMV_REAL(alpha);
#ifdef TMVFLDEBUG
                    TMVAssert(A.ptr() >= A.first);
                    TMVAssert(A.ptr() < A.last);
#endif
                    if (add) *(A.ptr()) += temp;
                    else *(A.ptr()) = temp;
                } else {
                    T temp = x.row(0) * x.row(0);
                    if (!a1) temp *= alpha;
#ifdef TMVFLDEBUG
                    TMVAssert(A.ptr() >= A.first);
                    TMVAssert(A.ptr() < A.last);
#endif
                    if (add) *(A.ptr()) += temp;
                    else *(A.ptr()) = temp;
                }
            } else {
                if (x.isrm()) {
                    if (A.isrm()) {
                        for (int i=0;i<N;++i) {
                            if (add)
                                A.row(i,0,i+1) += alpha * x.row(i) * 
                                    (ha ?
                                     x.rowRange(0,i+1).adjoint() :
                                     x.rowRange(0,i+1).transpose());
                            else
                                A.row(i,0,i+1) = alpha * x.row(i) * 
                                    (ha ?
                                     x.rowRange(0,i+1).adjoint() :
                                     x.rowRange(0,i+1).transpose());
                        }
                    } else {
                        for (int j=0;j<N;++j) {
                            if (add)
                                A.col(j,j,N) += alpha * x.rowRange(j,N) * 
                                    (ha ? x.row(j).conjugate() : x.row(j));
                            else
                                A.col(j,j,N) = alpha * x.rowRange(j,N) * 
                                    (ha ? x.row(j).conjugate() : x.row(j));
                        }
                    }
                } else { // x not row major
                    const int K = x.rowsize();
                    for (int i=0;i<K;i++)
                        Rank1Update<add>(alpha,x.col(i),A);
                }
            }
        } else {
            int k = N/2;
            const int nb = SYM_RK_BLOCKSIZE;
            if (k > nb) k = k/nb*nb;
            RecursiveRankKUpdate<ha,a1,add>(
                alpha,x.rowRange(0,k),A.subSymMatrix(0,k));
            MultMM<add>(
                alpha,x.rowRange(k,N),
                (ha ? x.rowRange(0,k).adjoint() : x.rowRange(0,k).transpose()),
                A.subMatrix(k,N,0,k));
            RecursiveRankKUpdate<ha,a1,add>(
                alpha,x.rowRange(k,N),A.subSymMatrix(k,N));
        }
    }

    template <bool cm, bool ha, bool a1, bool add, class T, class Tx> 
    static void RecursiveInPlaceRankKUpdate(
        const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
    {
        TMVAssert(SameStorage(x,A));
        TMVAssert(A.size() > 0);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(cm == A.iscm());
        TMVAssert(ha == A.isherm());
        TMVAssert(a1 == (alpha == T(1)));

        int N = A.size();
        if (N == 1) {
            Tx x00 = x(0,0);
            if (ha) {
                TMV_RealType(T) temp = TMV_NORM(x00);
                TMVAssert(TMV_IMAG(alpha) == TMV_RealType(T)(0));
                if (!a1) temp *= TMV_REAL(alpha);
#ifdef TMVFLDEBUG
                TMVAssert(A.ptr() >= A.first);
                TMVAssert(A.ptr() < A.last);
#endif
                if (add) *A.ptr() += temp;
                else *A.ptr() = temp;
            } else {
                T temp = x00 * x00;
                if (!a1) temp *= alpha;
#ifdef TMVFLDEBUG
                TMVAssert(A.ptr() >= A.first);
                TMVAssert(A.ptr() < A.last);
#endif
                if (add) *A.ptr() += temp;
                else *A.ptr() = temp;
            }
        } else {
            // [ A00  A01 ] += alpha * [ x00  x01 ] [ x00t  x10t ]
            // [ A10  A11 ]            [ x10  x11 ] [ x01t  x11t ]
            //               = alpha * [ x00 x00t+x01 x01t  x00 x10t+x01 x11t ]
            //                         [ x10 x00t+x11 x01t  x10 x10t+x11 x11t ]
            // Note that there is no order to do these in which overwriting A??
            // won't screw up an x?? which is needed later.
            // Need a temporary.  I choose A10, but it doesn't much matter.
            const int k = N/2;
            const ConstMatrixView<Tx> x00 = x.subMatrix(0,k,0,k);
            const ConstMatrixView<Tx> x10 = x.subMatrix(k,N,0,k);
            const ConstMatrixView<Tx> x01 = x.subMatrix(0,k,k,N);
            const ConstMatrixView<Tx> x11 = x.subMatrix(k,N,k,N);
            SymMatrixView<T> A00 = A.subSymMatrix(0,k);
            SymMatrixView<T> A11 = A.subSymMatrix(k,N);
            MatrixView<T> A10 = A.subMatrix(k,N,0,k);

            Matrix<T,cm?ColMajor:RowMajor> tempA10 = 
                x10 * (ha ? x00.adjoint() : x00.transpose());
            tempA10 += x11 * (ha ? x01.adjoint() : x01.transpose());

            RecursiveInPlaceRankKUpdate<cm,ha,a1,add>(alpha,x11,A11);
            RecursiveRankKUpdate<ha,a1,true>(alpha,x10,A11);
            RecursiveInPlaceRankKUpdate<cm,ha,a1,add>(alpha,x00,A00);
            RecursiveRankKUpdate<ha,a1,true>(alpha,x01,A00);

            if (add) A10 += alpha * tempA10;
            else A10 = alpha * tempA10;
        }
    }

    template <bool add, class T, class Tx> 
    static void InPlaceRankKUpdate(
        const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
    {
        TMVAssert(A.uplo() == Lower);
        if (A.iscm()) {
            if (A.isherm()) {
                if (alpha == T(1))
                    RecursiveInPlaceRankKUpdate<true,true,true,add>(
                        alpha,x,A); 
                else
                    RecursiveInPlaceRankKUpdate<true,true,false,add>(
                        alpha,x,A); 
            } else {
                if (alpha == T(1))
                    RecursiveInPlaceRankKUpdate<true,false,true,add>(
                        alpha,x,A); 
                else
                    RecursiveInPlaceRankKUpdate<true,false,false,add>(
                        alpha,x,A); 
            }
        } else {
            if (A.isherm()) {
                if (alpha == T(1))
                    RecursiveInPlaceRankKUpdate<false,true,true,add>(
                        alpha,x,A); 
                else
                    RecursiveInPlaceRankKUpdate<false,true,false,add>(
                        alpha,x,A); 
            } else {
                if (alpha == T(1))
                    RecursiveInPlaceRankKUpdate<false,false,true,add>(
                        alpha,x,A); 
                else
                    RecursiveInPlaceRankKUpdate<false,false,false,add>(
                        alpha,x,A); 
            }
        }
    }

    template <bool add, class T, class Tx> 
    static void NonBlasRankKUpdate(
        const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
    {
        if (A.uplo() == Upper) {
            return NonBlasRankKUpdate<add>(
                alpha,x,A.issym()?A.transpose():A.adjoint());
        } else if (A.isconj()) {
            return NonBlasRankKUpdate<add>(
                TMV_CONJ(alpha),x.conjugate(),A.conjugate());
        } else if (SameStorage(x,A)) {
            const int N = A.size();
            TMVAssert(int(x.colsize()) == N);
            const int K = x.rowsize();
            if (K >= N) {
                InPlaceRankKUpdate<add>(alpha,x.colRange(0,N),A);
                if (K > N) NonBlasRankKUpdate<add>(alpha,x.colRange(N,K),A);
            } else {
                NonBlasRankKUpdate<add>(
                    alpha,x.rowRange(K,N),A.subSymMatrix(K,N));
                MultMM<add>(
                    alpha,x.rowRange(K,N),
                    A.isherm() ?
                    x.rowRange(0,K).adjoint() :
                    x.rowRange(0,K).transpose(),
                    A.subMatrix(K,N,0,K) );
                InPlaceRankKUpdate<add>(
                    alpha,x.rowRange(0,K),A.subSymMatrix(0,K));
            }
        } else {
            if (A.isherm())
                if (alpha == T(1))
                    RecursiveRankKUpdate<true,true,add>(alpha,x,A); 
                else
                    RecursiveRankKUpdate<true,false,add>(alpha,x,A); 
            else
                if (alpha == T(1))
                    RecursiveRankKUpdate<false,true,add>(alpha,x,A); 
                else
                    RecursiveRankKUpdate<false,false,add>(alpha,x,A); 
        }
    }

#ifdef BLAS
    template <class T, class Tx> 
    static inline void BlasRankKUpdate(
        const T alpha, const GenMatrix<Tx>& x,
        const SymMatrixView<T>& A)
    { NonBlasRankKUpdate<true>(alpha,x,A); }
#ifdef INST_DOUBLE
    template <> 
    void BlasRankKUpdate(
        const double alpha, const GenMatrix<double>& x,
        const SymMatrixView<double>& A)
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(alpha != 0.);
        TMVAssert(x.colsize() > 0);
        TMVAssert(x.rowsize() > 1);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.iscm());
        TMVAssert(x.isrm() || x.iscm());

        int n = A.size();
        int k = x.rowsize();
        int ldx = x.isrm()?x.stepi():x.stepj();
        double beta(1);
        int lda = A.stepj();
        BLASNAME(dsyrk) (
            BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
            x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(x.cptr()),BLASV(ldx),BLASV(beta),BLASP(A.ptr()),
            BLASV(lda) BLAS1 BLAS1);
    }
    template <> 
    void BlasRankKUpdate(
        const std::complex<double> alpha,
        const GenMatrix<std::complex<double> >& x, 
        const SymMatrixView<std::complex<double> >& A)
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(alpha != 0.);
        TMVAssert(TMV_IMAG(alpha)==0. || !A.isherm());
        TMVAssert(x.colsize() > 0);
        TMVAssert(x.rowsize() > 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.iscm());
        TMVAssert(x.isrm() || x.iscm());
        TMVAssert(A.isherm() || !x.isconj());
        TMVAssert(A.issym() || x.isrm() == x.isconj());

        int n = A.size();
        int k = x.rowsize();
        int ldx = x.isrm()?x.stepi():x.stepj();
        int lda = A.stepj();
        if (A.isherm()) {
            double a(TMV_REAL(alpha));
            double beta(1);
            BLASNAME(zherk) (
                BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                x.isrm()?BLASCH_CT:BLASCH_NT,BLASV(n),BLASV(k),BLASV(a),
                BLASP(x.cptr()),BLASV(ldx),BLASV(beta),BLASP(A.ptr()),
                BLASV(lda) BLAS1 BLAS1);
        } else {
            std::complex<double> beta(1);
            BLASNAME(zsyrk) (
                BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(x.cptr()),BLASV(ldx),BLASP(&beta),BLASP(A.ptr()),
                BLASV(lda) BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasRankKUpdate(
        const std::complex<double> alpha,
        const GenMatrix<double>& x, 
        const SymMatrixView<std::complex<double> >& A)
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(alpha != 0.);
        TMVAssert(TMV_IMAG(alpha)==0. || !A.isherm());
        TMVAssert(x.colsize() > 0);
        TMVAssert(x.rowsize() > 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.iscm());
        TMVAssert(x.isrm() || x.iscm());
        TMVAssert(A.isherm() || !x.isconj());
        TMVAssert(A.issym() || x.isrm() == x.isconj());

        SymMatrix<double,Lower,ColMajor> A1(A.size(),0.);
        BlasRankKUpdate(1.,x,A1.view());
        A += alpha*A1;
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasRankKUpdate(
        const float alpha, const GenMatrix<float>& x,
        const SymMatrixView<float>& A)
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(alpha != 0.F);
        TMVAssert(x.colsize() > 0);
        TMVAssert(x.rowsize() > 1);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.iscm());
        TMVAssert(x.isrm() || x.iscm());

        int n = A.size();
        int k = x.rowsize();
        int ldx = x.isrm()?x.stepi():x.stepj();
        float beta(1);
        int lda = A.stepj();
        BLASNAME(ssyrk) (
            BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
            x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(x.cptr()),BLASV(ldx),BLASV(beta),BLASP(A.ptr()),
            BLASV(lda) BLAS1 BLAS1);
    }
    template <> 
    void BlasRankKUpdate(
        const std::complex<float> alpha,
        const GenMatrix<std::complex<float> >& x, 
        const SymMatrixView<std::complex<float> >& A)
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(alpha != 0.F);
        TMVAssert(TMV_IMAG(alpha)==0.F || !A.isherm());
        TMVAssert(x.colsize() > 0);
        TMVAssert(x.rowsize() > 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.iscm());
        TMVAssert(x.isrm() || x.iscm());
        TMVAssert(A.isherm() || !x.isconj());
        TMVAssert(A.issym() || x.isrm() == x.isconj());

        int n = A.size();
        int k = x.rowsize();
        int ldx = x.isrm()?x.stepi():x.stepj();
        int lda = A.stepj();
        if (A.isherm()) {
            float a(TMV_REAL(alpha));
            float beta(1);
            BLASNAME(cherk) (
                BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                x.isrm()?BLASCH_CT:BLASCH_NT,BLASV(n),BLASV(k),BLASV(a),
                BLASP(x.cptr()),BLASV(ldx),BLASV(beta),BLASP(A.ptr()),
                BLASV(lda) BLAS1 BLAS1);
        } else {
            std::complex<float> beta(1);
            BLASNAME(csyrk) (
                BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(x.cptr()),BLASV(ldx),BLASP(&beta),BLASP(A.ptr()),
                BLASV(lda) BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasRankKUpdate(
        const std::complex<float> alpha,
        const GenMatrix<float>& x, 
        const SymMatrixView<std::complex<float> >& A)
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(alpha != 0.F);
        TMVAssert(TMV_IMAG(alpha)==0.F || !A.isherm());
        TMVAssert(x.colsize() > 0);
        TMVAssert(x.rowsize() > 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.iscm());
        TMVAssert(x.isrm() || x.iscm());
        TMVAssert(A.isherm() || !x.isconj());
        TMVAssert(A.issym() || x.isrm() == x.isconj());

        SymMatrix<float,Lower,ColMajor> A1(A.size(),0.F);
        BlasRankKUpdate(1.F,x,A1.view());
        A += alpha*A1;
    }
#endif 
#endif // BLAS

    template <bool add, class T, class Tx> 
    void RankKUpdate(
        const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
    // A (+)= alpha * x * xT
    {
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<Tx> x0 = x;
        Matrix<Tx> xt = A.isherm() ? x.adjoint() : x.transpose();
        Matrix<T> A2 = A;
        Matrix<Tx> xxt = x0 * xt;
        if (add) A2 += alpha*xxt;
        else A2 = alpha * xxt;
        //cout<<"Start RankKUpdate: A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"alpha = "<<alpha<<", x = "<<TMV_Text(x)<<"  "<<x<<endl;
#endif

        TMVAssert(A.size() == x.colsize());
        TMVAssert(TMV_IMAG(alpha)==TMV_RealType(T)(0) || !A.isherm());
        if (alpha != T(0) && x.colsize() > 0 && x.rowsize() > 0) {
            if (x.rowsize() == 1)
                return Rank1Update<add>(alpha,x.col(0),A);
#ifdef BLAS
            else if (A.isrm())
                return RankKUpdate<add>(
                    alpha,x,A.issym()?A.transpose():A.adjoint());
            else if (A.isconj()) 
                return RankKUpdate<add>(
                    TMV_CONJ(alpha),x.conjugate(),A.conjugate());
            else if (A.iscm() && A.stepj()>0) {
                if (!add) A.setZero();
                if (!((x.iscm() && x.stepj()>0) || (x.isrm() && x.stepi()>0)) || 
                    (!A.issym() && x.iscm() == x.isconj()) ||
                    (!A.isherm() && x.isconj()) || SameStorage(x,A)) {
                    Matrix<T,ColMajor> xx = x;
                    BlasRankKUpdate(alpha,xx,A);
                } else {
                    BlasRankKUpdate(alpha,x,A);
                }
            } else {
                if (A.isherm()) {
                    HermMatrix<T,Lower,ColMajor> AA(A.size());
                    RankKUpdate<false>(alpha,x,AA.view());
                    if (add) A += AA;
                    else A = AA;
                } else {
                    SymMatrix<T,Lower,ColMajor> AA(A.size());
                    RankKUpdate<false>(alpha,x,AA.view());
                    if (add) A += AA;
                    else A = AA;
                }
            } 
#else
            NonBlasRankKUpdate<add>(alpha,x,A);
#endif
        }

#ifdef XDEBUG
        TMVAssert(A.isHermOK());
        if (!(Norm(A-A2) < 0.001*(TMV_ABS(alpha)*TMV_SQR(Norm(x0))+
                                  add?Norm(A0):TMV_RealType(T)(0)))) {
            cerr<<"RankKUpdate: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"x = "<<TMV_Text(x)<<"  "<<x0<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> A = "<<A<<endl;
            cerr<<"A2 = "<<A2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_RankK_MMS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

