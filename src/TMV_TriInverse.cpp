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
#include "tmv/TMV_TriDiv.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define TRI_DIV_BLOCKSIZE TMV_BLOCKSIZE
#else
#define TRI_DIV_BLOCKSIZE 64
#endif

    template <bool unit, class T> 
    static void RecursiveInverse(UpperTriMatrixView<T> U)
    {
        TMVAssert(U.iscm() || U.isrm());
        TMVAssert(unit == U.isunit());

        const ptrdiff_t N = U.size();
        const ptrdiff_t nb = TRI_DIV_BLOCKSIZE;

        if (N == 1) {
            if (!unit) {
                T*const Uptr = U.ptr();
                if (*Uptr == T(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular UpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularUpperTriMatrix<T>(U);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(Uptr >= U._first);
                TMVAssert(Uptr < U._last);
#endif
                *Uptr = TMV_InverseOf(*Uptr);
            }
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            UpperTriMatrixView<T> U00 = U.subTriMatrix(0,k);
            MatrixView<T> U01 = U.subMatrix(0,k,k,N);
            UpperTriMatrixView<T> U11 = U.subTriMatrix(k,N);

            // U00 U01' + U01 U11' = 0
            // U00 U01' = -U01 U11'
            // U01' = -U00' U01 U11'

            RecursiveInverse<unit>(U00);
            RecursiveInverse<unit>(U11);
            U01 = -U00 * U01;
            U01 *= U11;
        }
    }

    template <class T> 
    static inline void NonLapInverse(UpperTriMatrixView<T> U)
    {
#ifndef NOTHROW
        try {
#endif
            if (U.isunit()) RecursiveInverse<true>(U);
            else RecursiveInverse<false>(U);
#ifndef NOTHROW
        } catch (Singular) {
            throw SingularUpperTriMatrix<T>(U);
        }
#endif
    }

#ifdef ALAP
    template <class T> 
    static inline void LapInverse(UpperTriMatrixView<T> m)
    { NonLapInverse(m); }
#ifdef INST_DOUBLE
    template <> 
    void LapInverse(UpperTriMatrixView<double> m)
    {
        int n = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        int Lap_info=0;
        LAPNAME(dtrtri) (
            LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
            m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
            LAP1 LAP1);
        LAP_Results(Lap_info,"dtrtri");
    }
    template <> 
    void LapInverse(
        UpperTriMatrixView<std::complex<double> > m)
    {
        int n = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        int Lap_info=0;
        LAPNAME(ztrtri) (
            LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
            m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
            LAP1 LAP1);
        LAP_Results(Lap_info,"ztrtri");
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapInverse(UpperTriMatrixView<float> m)
    {
        int n = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        int Lap_info=0;
        LAPNAME(strtri) (
            LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
            m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
            LAP1 LAP1);
        LAP_Results(Lap_info,"strtri");
    }
    template <> 
    void LapInverse(
        UpperTriMatrixView<std::complex<float> > m)
    {
        int n = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        int Lap_info=0;
        LAPNAME(ctrtri) (
            LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
            m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
            LAP1 LAP1);
        LAP_Results(Lap_info,"ctrtri");
    }
#endif // FLOAT
#endif // ALaP

    template <class T> 
    void TriInverse(UpperTriMatrixView<T> U)
    {
#ifdef XDEBUG
        Matrix<T> U0(U);
#endif

        if (U.size() > 0) {
            if (!(U.iscm() || U.isrm())) {
                UpperTriMatrix<T> temp = U;
                TriInverse(temp.view());
                U = temp;
            } else {
#ifdef ALAP
                LapInverse(U);
#else
                NonLapInverse(U);
#endif
            }
        }
#ifdef XDEBUG
        Matrix<T> eye = U*U0;
        if (Norm(eye-T(1)) > 0.0001*(Norm(U0)+Norm(U))) {
            cerr<<"UpperTriMatrix Inverse:\n";
            cerr<<"U = "<<TMV_Text(U)<<"  "<<U0<<endl;
            cerr<<"Uinv = "<<U<<endl;
            cerr<<"Uinv*U = "<<U*U0<<endl;
            cerr<<"U*Uinv = "<<U0*U<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_TriInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


