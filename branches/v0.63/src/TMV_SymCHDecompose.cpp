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
#include "TMV_SymCHDiv.h"
#include "tmv/TMV_SymCHD.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_TriMatrixArith.h"

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include "tmv/TMV_TriMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

#ifdef TMV_BLOCKSIZE
#define CH_BLOCKSIZE TMV_BLOCKSIZE
#define CH_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define CH_BLOCKSIZE 64
#define CH_BLOCKSIZE2 2
#endif

    //
    // Decompose
    //

    template <bool cm, class T> 
    static void NonBlockCH_Decompose(const SymMatrixView<T>& A)
    {
        // Cholesky Decompostion is basically a simpler version of LU Decomp,
        // since U is just Lt for Hermitian matrices.  Cholesky cannot
        // be used for complex symmetric matrices.
        //
        // Also, the matrix must be positive definite.
        // Positive definite by definition means that:
        // - All eignevalues are positive.
        // Necessary, but not sufficient, conditions of positive definite
        // matrices are:
        // - All diagonal elements are positive.
        // - Aii + Ajj > 2*Re[Aij]
        // - det(A) > 0
        // - maxElement(A) is on main diagonal
        // A necessary and sufficient condition is:
        // - det(A(0:i,0:i)) > 0 for all 0<i<=N
        //
        // We want to decompose the matrix (input as A) into L * Lt.
        //
        // The equations for the rows of A are:
        //
        // A(i,0:i)   = L(i,0:i) * Lt(0:i,0:i)
        // A(i,i)     = L(i,0:i) * Lt(0:i,i) + |L(i,i)|^2
        // A(i,i+1:N) = L(i,0:i) * Lt(0:i,i+1:N) + L(i,i) * Lt(i,i+1:N)
        // 
        // These are solved as: (Taking A to be stored in its LowerTriangle part)
        //
        // L(i,0:i) = A(i,0:i) % L(0:i,0:i)t
        // L(i,i) = sqrt( A(i,i) - normSq(L(i,0:i)) )
        // [ ie. row by row from the top. ]
        //
        // or 
        // 
        // L(i,i) = sqrt( A(i,i) - normSq(L(i,0:i)) )
        // L(i+1:N,i) = (A(i+1:N,i) - L(i,0:i) * L(i+1:N,0:i)t) / L(i,i)
        // [ ie. col by col from the left ]
        //
        // The second version is not so good as written, since it accesses L
        // by both rows and columns, which means there will be non-unit strided
        // vectors in the algorithm.  However, we can save it by
        // noting that the second equation can be calculated piecemeal:
        //
        // L(i+1:N,i) = A(i+1:N,i)
        // for(j=0..i) L(i+1:N,i) -= L(i,j) * L(i+1:N,j)*
        // L(i+1:N,i) /= L(i,i)
        //
        // (The row access in the first equation can proceed similarly.)
        // 
        // Then note that for a given j, all the i's with i>j can have 
        // L(i+1:N,i) calculated at once.  So the final column algorithm
        // is:
        //
        // L = A.lowerTri();
        // for(j=0..N) 
        //   L(j,j) = sqrt(L(j,j))
        //   L(j+1:N,j) /= L(j,j)
        //   L(j+1:N,j+1:N) -= L(j+1:N,j) ^ L(j+1:N,j)*
        //
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(isReal(T()) || A.isherm());
        TMVAssert(cm == A.iscm());
        const int N = A.size();
#ifdef XDEBUG
        Matrix<T> A0(A);
#ifdef XTEST
        for(int i=1;i<=N;i++) {
            T d = Matrix<T>(A.subSymMatrix(0,i)).det();
            if (TMV_REAL(d) < 0) {
                cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
                cout<<"det(0.."<<i<<") = "<<d<<endl;
            }
        }
#endif
#endif

        const VectorView<RT> Adiag = A.real().diag();
        if (cm) {
            RT* Ajj= Adiag.ptr();
            const int ds = Adiag.step();
            for(int j=0;j<N-1;++j,Ajj+=ds) {
                if (*Ajj <= RT(0)) {
#ifdef NOTHROW
                    std::cerr<<"Non Posdef HermMatrix found \n"; 
                    exit(1); 
#else
                    throw NonPosDefHermMatrix<T>(A);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(Ajj >= A.real().first);
                TMVAssert(Ajj < A.real().last);
#endif
                *Ajj = TMV_SQRT(*Ajj);
                A.col(j,j+1,N) /= *Ajj;
                A.subSymMatrix(j+1,N) -= 
                    A.col(j,j+1,N) ^ A.col(j,j+1,N).conjugate();
            }
            if (*Ajj <= RT(0)) {
#ifdef NOTHROW
                std::cerr<<"Non Posdef HermMatrix found \n"; 
                exit(1); 
#else
                throw NonPosDefHermMatrix<T>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(Ajj >= A.real().first);
            TMVAssert(Ajj < A.real().last);
#endif
            *Ajj = TMV_SQRT(*Ajj);
        } else {
            RT* Aii = Adiag.ptr();
            const int ds = Adiag.step();
            if (*Aii <= RT(0)) {
#ifdef NOTHROW
                std::cerr<<"Non Posdef HermMatrix found \n"; 
                exit(1); 
#else
                throw NonPosDefHermMatrix<T>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(Aii >= A.real().first);
            TMVAssert(Aii < A.real().last);
#endif
            *Aii = TMV_SQRT(*Aii);
            for(int i=1; i<N; ++i) {
                Aii+=ds;
                A.row(i,0,i) %= A.lowerTri().subTriMatrix(0,i).adjoint();
#ifdef TMVFLDEBUG
                TMVAssert(Aii >= A.real().first);
                TMVAssert(Aii < A.real().last);
#endif
                *Aii -= normSq(A.row(i,0,i));
                if (*Aii <= RT(0)) {
#ifdef NOTHROW
                    std::cerr<<"Non Posdef HermMatrix found \n"; 
                    exit(1); 
#else
                    throw NonPosDefHermMatrix<T>(A);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(Aii >= A.real().first);
                TMVAssert(Aii < A.real().last);
#endif
                *Aii = TMV_SQRT(*Aii);
            }
        }
#ifdef XDEBUG
        RT norml = Norm(A.lowerTri());
        Matrix<T> A2 = A.lowerTri()*A.lowerTri().adjoint();
        if (Norm(A2-A0) > 0.001*TMV_SQR(norml)) {
            cerr<<"CHDecomp: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"Done: A = "<<A<<endl;
            cerr<<"L = "<<A.lowerTri()<<endl;
            cerr<<"Lt = "<<A.lowerTri().adjoint()<<endl;
            cerr<<"L*Lt = "<<A2<<endl;
            cerr<<"Norm(diff) = "<<Norm(A2-A0)<<"  Norm(L) = "<<norml<<endl;
            abort();
        }
#endif
    }

    template <bool cm, class T> 
    static void BlockCH_Decompose(const SymMatrixView<T>& A)
    {
        // If A is large, we can take advantage of Blas Level 3 speed
        // by partitioning the matrix into blocks.
        //
        // A = [ A00  A01 ]
        //     [ A10  A11 ]
        // A00, A11 are Hermitian,
        // A01 = A10t
        // A00 is k by k, A11 is (N-k) by (N-k)
        //
        // Then, the L Lt decomposition of A is:
        //
        // A = [ L00   0  ] [ L00t L10t ]
        //     [ L10  L11 ] [  0   L11t ]
        //   = [ L00 L00t         L00 L10t      ]
        //     [ L10 L00t   L10 L10t + L11 L11t ]
        //
        // So A00 = L00 L00t
        //    A10 = L10 L00t
        //    A11 = L10 L10t + L11 L11t
        //
        // The block routine is:
        //
        // 1) First Cholesky decompose A00 into L00 L00t.  
        // 2) Solve for L10 = A10 L00t^-1.
        // 3) Find A11' = A11 - L10 L10t.
        // 4) Repeat with A11'.
        //
#ifdef XDEBUG
        //cout<<"Start Block CH_Decomp: A = "<<TMV_Text(A)<<"  "<<A<<endl;
        Matrix<T> A0(A);
        HermMatrix<T,Lower,ColMajor> A2 = A;
        NonBlockCH_Decompose(A2.view());
#endif
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(isReal(T()) || A.isherm());
        TMVAssert(cm == A.iscm());

        const int N = A.size();

        for (int jk=0; jk<N; jk+=CH_BLOCKSIZE)
        {
            int jkpk = TMV_MIN(jk+CH_BLOCKSIZE,N);

            SymMatrixView<T> A00 = A.subSymMatrix(jk,jkpk);
            SymMatrixView<T> A11 = A.subSymMatrix(jkpk,N);
            MatrixView<T> A10 = A.subMatrix(jkpk,N,jk,jkpk);

            NonBlockCH_Decompose<cm>(A00);

            UpperTriMatrixView<T> L00t = A00.upperTri();
            A10 %= L00t;

            A11 -= A10 * A10.adjoint();
        }

#ifdef XDEBUG
        ConstLowerTriMatrixView<T> L = A.lowerTri();
        Matrix<T> AA = L * L.adjoint();
        if (Norm(AA-A0) > 0.001*Norm(A0)) {
            cerr<<"Done Block CH: \n";
            cerr<<"A (block) = "<<A<<endl;
            cerr<<"A (nonblock) = "<<A2<<endl;
            cerr<<"Orig A = "<<A0<<endl;
            cerr<<"LLt (block) = "<<AA<<endl;
            cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
            cerr<<"Norm(L-L2) = "<<Norm(L-A2.lowerTri())<<endl;
            cerr<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
            abort();
        }
#endif
    }

    template <bool cm, class T> 
    static void RecursiveCH_Decompose(const SymMatrixView<T>& A)
    {
#ifdef XDEBUG
        //cout<<"Start Recursive CH_Decomp: A = "<<TMV_Text(A)<<"  "<<A<<endl;
        Matrix<T> A0(A);
        HermMatrix<T,Lower,ColMajor> A2(A);
        NonBlockCH_Decompose<true>(A2.view());
#endif
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isherm());
        TMVAssert(cm == A.iscm());
        const int N = A.size();

        if (N > CH_BLOCKSIZE2) {
            int k = N/2;
            if (k > CH_BLOCKSIZE) k = k/CH_BLOCKSIZE*CH_BLOCKSIZE;

            SymMatrixView<T> A00 = A.subSymMatrix(0,k);
            SymMatrixView<T> A11 = A.subSymMatrix(k,N);
            MatrixView<T> A10 = A.subMatrix(k,N,0,k);

            RecursiveCH_Decompose<cm>(A00);

            UpperTriMatrixView<T> L00t = A00.upperTri();
            A10 %= L00t;
            A11 -= A10 * A10.adjoint();

            RecursiveCH_Decompose<cm>(A11);

        } else if (CH_BLOCKSIZE2 > 2 && N > 2) {
            NonBlockCH_Decompose<cm>(A);
        } else if (N > 0) {
            T* Aptr = A.ptr();
            RT A00 = TMV_REAL(*Aptr);
            if (A00 <= RT(0)) {
#ifdef NOTHROW
                std::cerr<<"Non Posdef HermMatrix found \n"; 
                exit(1); 
#else
                throw NonPosDefHermMatrix<T>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(Aptr >= A.first);
            TMVAssert(Aptr < A.last);
#endif
            *Aptr = TMV_SQRT(A00);

            if (N == 2) {
                A00 = TMV_REAL(*Aptr);
                if (cm) {
                    T A10 = (*(Aptr+1) /= A00);
                    T* A11ptr = Aptr+1+A.stepj();
                    RT A11 = TMV_REAL(*A11ptr);
                    A11 -= TMV_NORM(A10);
                    if (A11 <= RT(0)) {
#ifdef NOTHROW
                        std::cerr<<"Non Posdef HermMatrix found \n"; 
                        exit(1); 
#else
                        throw NonPosDefHermMatrix<T>(A);
#endif
                    }
#ifdef TMVFLDEBUG
                    TMVAssert(A11ptr >= A.first);
                    TMVAssert(A11ptr < A.last);
#endif
                    *(A11ptr) = TMV_SQRT(A11);
                } else {
                    const int si = A.stepi();
                    T A10 = (*(Aptr+si) /= A00);
                    T* A11ptr = Aptr+1+si;
                    RT A11 = TMV_REAL(*A11ptr);
                    A11 -= TMV_NORM(A10);
                    if (A11 <= RT(0)) {
#ifdef NOTHROW
                        std::cerr<<"Non Posdef HermMatrix found \n"; 
                        exit(1); 
#else
                        throw NonPosDefHermMatrix<T>(A);
#endif
                    }
#ifdef TMVFLDEBUG
                    TMVAssert(A11ptr >= A.first);
                    TMVAssert(A11ptr < A.last);
#endif
                    *(A11ptr) = TMV_SQRT(A11);
                }
            }
        }
#ifdef XDEBUG
        ConstLowerTriMatrixView<T> L = A.lowerTri();
        Matrix<T> AA = L * L.adjoint();
        if (Norm(AA-A0) > 0.001*Norm(A0)) {
            cerr<<"Done Recursive CH: \n";
            cerr<<"A (block) = "<<A<<endl;
            cerr<<"A (nonblock) = "<<A2<<endl;
            cerr<<"Orig A = "<<A0<<endl;
            cerr<<"LLt (block) = "<<AA<<endl;
            cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
            cerr<<"Norm(L-L2) = "<<Norm(L-A2.lowerTri())<<endl;
            cerr<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    static void NonLapCH_Decompose(const SymMatrixView<T>& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(isReal(T()) || A.isherm());

        if (A.iscm())
            RecursiveCH_Decompose<true>(A);
        else
            RecursiveCH_Decompose<false>(A);
    }

#ifdef ALAP
    template <class T> 
    static inline void LapCH_Decompose(const SymMatrixView<T>& A)
    { NonLapCH_Decompose(A); }
#ifdef INST_DOUBLE
    template <> 
    void LapCH_Decompose(const SymMatrixView<double>& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(A.ct()==NonConj);

        int n = A.size();
        int lda = A.iscm() ? A.stepj() : A.stepi();

        LAPNAME(dpotrf) (
            LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n),
            LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
        if (Lap_info > 0) {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermMatrix found \n"; 
            exit(1); 
#else
            throw NonPosDefHermMatrix<double>(A);
#endif
        }
        LAP_Results("dpotrf");
    }
    template <> 
    void LapCH_Decompose(
        const SymMatrixView<std::complex<double> >& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.isherm());

        int n = A.size();
        int lda = A.iscm() ? A.stepj() : A.stepi();
        LAPNAME(zpotrf) (
            LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n),
            LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
        if (Lap_info > 0) {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermMatrix found \n"; 
            exit(1); 
#else
            throw NonPosDefHermMatrix<std::complex<double> >(A);
#endif
        }
        LAP_Results("zpotrf");
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapCH_Decompose(const SymMatrixView<float>& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(A.ct()==NonConj);

        int n = A.size();
        int lda = A.iscm() ? A.stepj() : A.stepi();
        LAPNAME(spotrf) (
            LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n),
            LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
        if (Lap_info > 0) {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermMatrix found \n"; 
            exit(1); 
#else
            throw NonPosDefHermMatrix<float>(A);
#endif
        }
        LAP_Results("spotrf");
    }
    template <> 
    void LapCH_Decompose(
        const SymMatrixView<std::complex<float> >& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.isherm());

        int n = A.size();
        int lda = A.iscm() ? A.stepj() : A.stepi();
        LAPNAME(cpotrf) (
            LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n),
            LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
        if (Lap_info > 0) {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermMatrix found \n"; 
            exit(1); 
#else
            throw NonPosDefHermMatrix<std::complex<float> >(A);
#endif
        }
        LAP_Results("cpotrf");
    }
#endif 
#endif // ALAP

    template <class T> 
    void CH_Decompose(const SymMatrixView<T>& A)
    {
#ifdef XTEST
        TMVAssert(A.isHermOK());
#endif
        TMVAssert(isReal(T()) || A.isherm());
        TMVAssert(A.isrm() || A.iscm());

        if (A.uplo() == Upper) CH_Decompose(A.adjoint());
        else if (A.isconj()) CH_Decompose(A.conjugate());
        else if (A.size() > 0) {
#ifdef ALAP 
            LapCH_Decompose(A);
#else
            NonLapCH_Decompose(A);
#endif
        }
#ifdef XTEST
        TMVAssert(A.isHermOK());
#endif
    }

#undef RT

#define InstFile "TMV_SymCHDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


