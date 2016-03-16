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


#include "TMV_Blas.h"
#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
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
    // S += L*Lt
    //

    template <bool ha, bool uu, bool a1, class T, class TL> 
    static void RecursiveRankKUpdate(
        const T alpha, const GenLowerTriMatrix<TL>& L,
        SymMatrixView<T> A)
    {
        TMVAssert(A.size() == L.size());
        TMVAssert(alpha != T(0));
        TMVAssert(TMV_IMAG(alpha)==TMV_RealType(T)(0) || !A.isherm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(ha == A.isherm());
        TMVAssert(uu == L.isunit());
        TMVAssert(a1 == (alpha == T(1)));

        ptrdiff_t N = A.size();

        if (N == 1) {
#ifdef TMVFLDEBUG
            TMVAssert(A.ptr() >= A._first);
            TMVAssert(A.ptr() < A._last);
#endif
            if (uu)
                *(A.ptr()) += a1 ? T(1) : alpha;
            else 
                if (a1)
                    *(A.ptr()) += 
                        ha ? TMV_NORM(*(L.cptr())) : TMV_SQR(*(L.cptr()));
                else
                    *(A.ptr()) += alpha * 
                        (ha ? TMV_NORM(*(L.cptr())) : TMV_SQR(*(L.cptr())));
        } else {
            //
            // [ A11 A12 ] += alpha [ L11  0  ] [ L11t L21t ]
            // [ A21 A22 ]          [ L21 L22 ] [  0   L22t ]
            //              = alpha [ L11 L11t        L11 L21t      ]
            //                      [ L21 L11t  L21 L21t + L22 L22t ]
            ptrdiff_t k = N/2;
            const ptrdiff_t nb = SYM_RK_BLOCKSIZE;
            if (k > nb) k = k/nb*nb;
            SymMatrixView<T> A11 = A.subSymMatrix(0,k);
            SymMatrixView<T> A22 = A.subSymMatrix(k,N);
            MatrixView<T> A21 = A.subMatrix(k,N,0,k);
            ConstLowerTriMatrixView<TL> L11 = L.subTriMatrix(0,k);
            ConstLowerTriMatrixView<TL> L22 = L.subTriMatrix(k,N);
            ConstMatrixView<TL> L21 = L.subMatrix(k,N,0,k);

            RecursiveRankKUpdate<ha,uu,a1>(alpha,L22,A22);

            RankKUpdate<true>(alpha,L21,A22);

            A21 += alpha * L21 * (ha ? L11.adjoint() : L11.transpose());

            RecursiveRankKUpdate<ha,uu,a1>(alpha,L11,A11);
        }
    }

    template <bool a1, class T, class TL> 
    static inline void DoRankKUpdate(
        const T alpha, const GenLowerTriMatrix<TL>& L,
        SymMatrixView<T> A)
    {
        if (A.isherm())
            if (L.isunit())
                RecursiveRankKUpdate<true,true,a1>(alpha,L,A);
            else
                RecursiveRankKUpdate<true,false,a1>(alpha,L,A);
        else
            if (L.isunit())
                RecursiveRankKUpdate<false,true,a1>(alpha,L,A);
            else
                RecursiveRankKUpdate<false,false,a1>(alpha,L,A);
    }

    template <bool ha, class T> 
    static void RecursiveSetLLt(SymMatrixView<T> A)
    {
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(ha == A.isherm());

        ptrdiff_t N = A.size();

        if (N == 1) {
#ifdef TMVFLDEBUG
            TMVAssert(A.ptr() >= A._first);
            TMVAssert(A.ptr() < A._last);
#endif
            *(A.ptr()) = ha ? TMV_NORM(*(A.cptr())) : TMV_SQR(*(A.cptr()));
        } else {
            // [ A11 A12 ] = [ L11  0  ] [ L11t L21t ]
            // [ A21 A22 ]   [ L21 L22 ] [  0   L22t ]
            //             = [ L11 L11t        L11 L21t      ]
            //               [ L21 L11t  L21 L21t + L22 L22t ]
            ptrdiff_t k = N/2;
            const ptrdiff_t nb = SYM_RK_BLOCKSIZE;
            if (k > nb) k = k/nb*nb;
            SymMatrixView<T> A11 = A.subSymMatrix(0,k);
            SymMatrixView<T> A22 = A.subSymMatrix(k,N);
            MatrixView<T> A21 = A.subMatrix(k,N,0,k);
            ConstLowerTriMatrixView<T> L11 = A11.lowerTri();

            RecursiveSetLLt<ha>(A22);

            RankKUpdate<true>(T(1),A21,A22);

            A21 *= (ha ? L11.adjoint() : L11.transpose());

            RecursiveSetLLt<ha>(A11);
        }
    }

    template <class T> 
    static inline void setLLt(SymMatrixView<T> A)
    {
        if (A.isherm()) RecursiveSetLLt<true>(A);
        else RecursiveSetLLt<false>(A);
    }

    template <bool add, class T, class TL> 
    void RankKUpdate(
        const T alpha, const GenLowerTriMatrix<TL>& L, SymMatrixView<T> A)
        // A = A + alpha * L * LT
    {
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<TL> L0 = L;
        Matrix<TL> Lt = A.isherm() ? L.adjoint() : L.transpose();
        Matrix<T> A2 = A;
        if (add) A2 += alpha*L*Lt;
        else A2 = alpha*L*Lt;
        //cout<<"Start RankKUpdate: A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"alpha = "<<alpha<<", L = "<<TMV_Text(L)<<"  "<<L<<endl;
#endif

        TMVAssert(A.size() == L.size());
        TMVAssert(TMV_IMAG(alpha)==TMV_RealType(T)(0) || !A.isherm());

        if (alpha != T(0) && A.size() > 0) {
            if (A.isconj()) 
                RankKUpdate<add>(TMV_CONJ(alpha),L.conjugate(),A.conjugate());
            else if (A.uplo() == Upper) 
                if (A.isherm()) RankKUpdate<add>(alpha,L,A.adjoint());
                else RankKUpdate<add>(alpha,L,A.transpose());
            else if (!add) {
                A.lowerTri() = L;
                setLLt(A);
                A *= alpha;
            }
            else
                if (alpha == T(1)) DoRankKUpdate<true>(alpha,L,A);
                else DoRankKUpdate<false>(alpha,L,A);
        }

#ifdef XDEBUG
        TMVAssert(A.isHermOK());
        if (!(Norm(A-A2) < 0.001*(TMV_ABS(alpha)*TMV_SQR(Norm(L0))+
                                  (add?Norm(A0):TMV_RealType(T)(0))))) {
            cerr<<"RankKUpdate: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"L = "<<TMV_Text(L)<<"  "<<L0<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> A = "<<A<<endl;
            cerr<<"A2 = "<<A2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_RankK_LUS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


