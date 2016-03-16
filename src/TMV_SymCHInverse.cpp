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
#include "TMV_SymCHDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_TriMatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <class T> 
    static void NonLapCHInverse(SymMatrixView<T> sinv)
    {
        TMVAssert(sinv.isherm());
        // inv = (L Lt)^-1 = Lt^-1 L^-1
        LowerTriMatrixView<T> L = sinv.lowerTri();
        L.invertSelf();
        sinv = L.adjoint() * L;
    }

#ifdef ALAP
    template <class T> 
    static inline void LapCHInverse(SymMatrixView<T> sinv)
    { NonLapCHInverse(sinv); }
#ifdef INST_DOUBLE
    template <> 
    void LapCHInverse(SymMatrixView<double> sinv)
    {
        TMVAssert(sinv.isherm());
        TMVAssert(sinv.iscm());

        int n = sinv.size();
        int lda = sinv.stepj();
        int Lap_info=0;
        LAPNAME(dpotri) (
            LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
            LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
        LAP_Results(Lap_info,"dpotri");
    }
    template <> 
    void LapCHInverse(SymMatrixView<std::complex<double> > sinv)
    {
        TMVAssert(sinv.isherm());
        TMVAssert(sinv.iscm());

        int n = sinv.size();
        int lda = sinv.stepj();
        int Lap_info=0;
        LAPNAME(zpotri) (
            LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
            LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
        LAP_Results(Lap_info,"zpotri");
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapCHInverse(SymMatrixView<float> sinv)
    {
        TMVAssert(sinv.isherm());
        TMVAssert(sinv.iscm());

        int n = sinv.size();
        int lda = sinv.stepj();
        int Lap_info=0;
        LAPNAME(spotri) (
            LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
            LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
        LAP_Results(Lap_info,"spotri");
    }
    template <> 
    void LapCHInverse(SymMatrixView<std::complex<float> > sinv)
    {
        TMVAssert(sinv.isherm());
        TMVAssert(sinv.iscm());

        int n = sinv.size();
        int lda = sinv.stepj();
        int Lap_info=0;
        LAPNAME(cpotri) (
            LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
            LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
        LAP_Results(Lap_info,"cpotri");
    }
#endif // FLOAT
#endif // LAP

    template <class T, class T1> 
    void CH_Inverse(const GenSymMatrix<T1>& LLx, SymMatrixView<T> sinv) 
    {
        TMVAssert(LLx.size() == sinv.size());

#ifdef XDEBUG
        Matrix<T> A = LLx.lowerTri() * LLx.upperTri();
#endif

        if (sinv.size() > 0) {
            if (
#ifdef ALAP
                !isSameType(T(),T1()) || 
#endif
                !(sinv.iscm() || sinv.isrm())) {
                HermMatrix<T,Lower|ColMajor> temp(sinv.size());
                CH_Inverse(LLx,temp.view());
                sinv = temp;
#ifdef ALAP
            } else if (sinv.isrm()) {
                CH_Inverse(LLx.transpose(),sinv.transpose());
#endif
            } else {
                sinv = LLx;
#ifdef ALAP
                LapCHInverse(sinv);
#else
                NonLapCHInverse(sinv);
#endif
            }
        }

#ifdef XDEBUG
        Matrix<T> eye = A * sinv;
        TMV_RealType(T) kappa = Norm(A) * Norm(sinv);
        if (Norm(eye-T(1)) > 0.0001*kappa*sinv.size()) {
            cerr<<"A = "<<A<<endl;
            cerr<<"sinv = "<<sinv<<endl;
            cerr<<"A*sinv = "<<A*sinv<<endl;
            cerr<<"sinv*A = "<<sinv*A<<endl;
            cerr<<"Norm(A*sinv-1) = "<<Norm(A*sinv-T(1))<<endl;
            cerr<<"kappa = "<<kappa<<endl;
            abort();
        }
#endif
        TMVAssert(sinv.isHermOK());
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymCHInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


